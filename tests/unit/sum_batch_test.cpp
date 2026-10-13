// SPDX-License-Identifier: Apache-2.0
// The grid reduction's sum kernels (P-6, D-404; core/reducers/sum_batch.hpp): every kernel this CPU and build have
// (scalar, SSE2, AVX2), selected in turn on one machine, bit for bit against the scalar code:
//   - weighted_sums against weighted_sum on random and adversarial surfaces (huge, tiny, subnormal, signed zeros,
//     cancelling terms, +-inf, NaN), with zero, integer and fractional weights, at every run length 1-40 and offsets
//     that leave remainders;
//   - tile_bounds against the bounded argmax's scalar loop;
//   - calibrate, and replicate_estimates with the full-grid and the bounded search, under every kernel against none
//     (the scalar loops), on a real panel's surface;
// and the selection by VCAL_REDUCE_KERNEL (a kernel this machine lacks skips visibly, exit code 77).
#include "core/reducers/sum_batch.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "core/reducers/weighted_sum.hpp"
#include "engine/calibrate.hpp"
#include "resample/bootstrap.hpp"
#include "resample/weights.hpp"
#include "tests/harness/vcal_test.hpp"

namespace {

namespace k = vcal::reducers::sum_batch;
namespace e = vcal::engine;
namespace rs = vcal::resample;
using Backend = vcal::backends::CpuBackend;
using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
using Adaptive = vcal::quadrature::GaussHermiteAdaptive<vcal::PrecisionF64>;
using Obs = Objective::Obs;

const k::Kernel kInitial = k::selected();
constexpr k::Kernel kAll[] = {k::Kernel::Scalar, k::Kernel::Sse2, k::Kernel::Avx2};

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
bool same(double a, double b) { return bits_equal(a, b) || (std::isnan(a) && std::isnan(b)); }

struct Restore {
    k::Kernel saved = k::selected();
    ~Restore() { k::select(saved); }
};

template <class F>
void for_each_kernel(F fn) {
    Restore restore;
    for (const auto kern : kAll) {
        if (!k::select(kern)) {
            vcal::test::note(std::string("reduce kernel ") + k::name(kern) + " not available here");
            continue;
        }
        fn(kern);
    }
}

struct Lcg {  // deterministic inputs
    std::uint64_t s;
    double uniform() {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        return static_cast<double>(s >> 11) * 0x1p-53;
    }
};

struct Surface {
    std::int64_t T, K;
    std::vector<double> L, w;
};

// Kinds: 0 log-likelihood-like, 1 adversarial magnitudes and cancellation, 2 with non-finite values.
Surface make_surface(Lcg& g, std::int64_t T, std::int64_t K, int kind) {
    Surface s{T, K, std::vector<double>(static_cast<std::size_t>(T * K)), std::vector<double>(static_cast<std::size_t>(T))};
    const double specials[] = {1e300, -1e300, 1e-300, -1e-300, std::numeric_limits<double>::denorm_min(), 0.0, -0.0,
                               1e16, -1e16, 1.0, -1.0, 3.0e-17};
    for (double& x : s.L) {
        const double u = g.uniform();
        x = -50.0 + 40.0 * g.uniform();
        if (kind >= 1 && u < 0.3) x = specials[static_cast<int>(g.uniform() * 12.0) % 12];
        if (kind == 2 && u > 0.97) {
            const double v = g.uniform();
            x = v < 0.4 ? HUGE_VAL : (v < 0.8 ? -HUGE_VAL : std::numeric_limits<double>::quiet_NaN());
        }
    }
    for (double& x : s.w) {
        const double u = g.uniform();
        x = u < 0.25 ? 0.0 : (u < 0.7 ? static_cast<double>(static_cast<int>(1.0 + 4.0 * g.uniform())) : 3.0 * g.uniform());
    }
    return s;
}

}  // namespace

// VCAL_REDUCE_KERNEL selects the kernel at start-up; a kernel this CPU or build lacks makes the test skip (exit 77,
// CTest SKIP_RETURN_CODE), never pass without having run that path.
VCAL_TEST(reduce_kernel_selection_follows_vcal_reduce_kernel) {
    const std::string env = vcal::test::env_or_empty("VCAL_REDUCE_KERNEL");
    for (const auto kern : kAll) {
        if (env == k::name(kern) && !k::available(kern)) {
            std::printf("SKIPPED: VCAL_REDUCE_KERNEL=%s, but this CPU or build has no %s kernel\n", env.c_str(),
                        env.c_str());
            std::fflush(stdout);
            std::exit(77);
        }
    }
    k::Kernel want = k::Kernel::Scalar;
    if (env.empty()) {
        if (k::available(k::Kernel::Sse2)) want = k::Kernel::Sse2;
        if (k::available(k::Kernel::Avx2)) want = k::Kernel::Avx2;
    } else {
        want = k::Kernel::Unavailable;
        for (const auto kern : kAll) {
            if (env == k::name(kern) && k::available(kern)) want = kern;
        }
    }
    vcal::test::note("VCAL_REDUCE_KERNEL '" + env + "': selected at start-up " + k::name(kInitial));
    VCAL_CHECK(kInitial == want);
#if defined(__x86_64__) || defined(_M_X64)
    VCAL_CHECK(k::available(k::Kernel::Sse2));
#endif
    Restore restore;
    VCAL_REQUIRE(k::select(k::Kernel::Unavailable));
    double L = 1.0, w = 1.0, out = 0.0;
    VCAL_CHECK(!k::weighted_sums(&L, 1, 1, &w, 0, 1, &out));
    VCAL_CHECK(!k::tile_bounds(&L, 1, 1, &w, 1e-9, &out));
    VCAL_CHECK(!k::select(static_cast<k::Kernel>(7)));
}

// weighted_sums is weighted_sum at every point, every run length and offset; tile_bounds is the scalar loop.
VCAL_TEST(reduce_kernels_equal_the_scalar_sums_bit_for_bit) {
    for_each_kernel([&](k::Kernel kern) {
        Lcg g{0x5u};
        std::size_t bad = 0, sums = 0, bounds = 0;
        for (const std::int64_t T : {1, 2, 3, 13, 40, 100}) {
            for (int kind = 0; kind < 3; ++kind) {
                const std::int64_t K = 37 + 13 * kind;
                const Surface s = make_surface(g, T, K, kind);
                for (std::int64_t n = 1; n <= 40 && n <= K; ++n) {
                    for (std::int64_t k0 = 0; k0 + n <= K; k0 += 7) {
                        double out[40];
                        VCAL_REQUIRE(k::weighted_sums(s.L.data(), T, K, s.w.data(), k0, n, out));
                        for (std::int64_t i = 0; i < n; ++i) {
                            bad += same(out[i], vcal::reducers::weighted_sum(s.L.data(), T, K, s.w.data(), k0 + i)) ? 0 : 1;
                        }
                        sums += static_cast<std::size_t>(n);
                    }
                }
                if (kind < 2) {  // the bounded search's guard keeps non-finite values away from its bounds
                    std::vector<double> bound(static_cast<std::size_t>(K));
                    VCAL_REQUIRE(k::tile_bounds(s.L.data(), T, K, s.w.data(), 1e-9, bound.data()));
                    for (std::int64_t i = 0; i < K; ++i) {
                        double sm = 0.0, mag = 0.0;
                        for (std::int64_t j = 0; j < T; ++j) {
                            if (s.w[static_cast<std::size_t>(j)] == 0.0) continue;
                            const double x = s.w[static_cast<std::size_t>(j)] * s.L[static_cast<std::size_t>(j * K + i)];
                            sm += x;
                            mag += std::fabs(x);
                        }
                        bad += same(bound[static_cast<std::size_t>(i)], sm + 1e-9 * mag) ? 0 : 1;
                        ++bounds;
                    }
                }
            }
        }
        vcal::test::note(std::string(k::name(kern)) + ": " + std::to_string(sums) + " sums, " + std::to_string(bounds) +
                         " tile bounds, " + std::to_string(bad) + " differ");
        VCAL_CHECK(bad == 0);
    });
}

// Calibrate and the bootstrap's replicate estimates (full grid and bounded) are the same under every kernel as with
// none selected, on a real panel's surface over a grid whose last axis leaves a one-point tile column.
VCAL_TEST(calibrate_and_replicate_estimates_are_the_same_under_every_reduce_kernel) {
    const std::int64_t d[] = {2, 5, 1, 9, 3, 0, 4, 12, 6, 2, 1, 3, 7, 15, 2, 4, 0, 5, 8, 3};
    std::vector<Obs> p;
    for (const auto x : d) p.push_back({1000, x});
    const vcal::Grid<2> grid{{{2e-4, 0.1, 41, vcal::AxisScale::Logit}, {5e-3, vcal::kDefaultRhoUpper, 33, vcal::AxisScale::Logit}}};
    const Adaptive primary{vcal::quadrature::gauss_hermite_rule(64)}, check{vcal::quadrature::gauss_hermite_rule(128)};
    const std::int64_t T = static_cast<std::int64_t>(p.size());
    constexpr std::uint32_t B = 64;
    const auto idx = rs::bootstrap_indices(0x5641534943454B31ull, rs::Scheme::IidBootstrap, B, T);
    const auto W = rs::weights_from_indices(idx.data(), B, T, T);
    struct Run {
        e::Estimate2 est;
        std::vector<rs::Replicate2> full, bounded;
    };
    auto run = [&]() {
        Run r{{}, std::vector<rs::Replicate2>(B), std::vector<rs::Replicate2>(B)};
        std::vector<double> L;
        VCAL_REQUIRE(e::calibrate(Backend{}, Objective{}, primary, check, p.data(), T, grid, L, r.est) == e::Status::Ok);
        rs::replicate_estimates(Backend{}, grid, L.data(), T, W.data(), B, r.full.data(), rs::ArgmaxSearch::Full);
        rs::replicate_estimates(Backend{}, grid, L.data(), T, W.data(), B, r.bounded.data(), rs::ArgmaxSearch::Bounded);
        return r;
    };
    auto same_rep = [](const rs::Replicate2& a, const rs::Replicate2& b) {
        return bits_equal(a.value[0], b.value[0]) && bits_equal(a.value[1], b.value[1]) &&
               bits_equal(a.surface_max, b.surface_max) && a.grid_index == b.grid_index && a.flags == b.flags;
    };
    Restore restore;
    VCAL_REQUIRE(k::select(k::Kernel::Unavailable));
    const Run want = run();
    for_each_kernel([&](k::Kernel kern) {
        const Run got = run();
        std::size_t bad = 0;
        bad += (bits_equal(got.est.value[0], want.est.value[0]) && bits_equal(got.est.value[1], want.est.value[1]) &&
                bits_equal(got.est.loglik, want.est.loglik) && bits_equal(got.est.se[0], want.est.se[0]) &&
                bits_equal(got.est.se[1], want.est.se[1]) && got.est.flags == want.est.flags)
                   ? 0
                   : 1;
        for (std::uint32_t b = 0; b < B; ++b) {
            bad += same_rep(got.full[b], want.full[b]) ? 0 : 1;
            bad += same_rep(got.bounded[b], want.bounded[b]) ? 0 : 1;
            bad += same_rep(got.full[b], got.bounded[b]) ? 0 : 1;  // the bounded search finds the full grid's argmax
        }
        vcal::test::note(std::string(k::name(kern)) + ": calibrate and " + std::to_string(B) +
                         " replicates (full and bounded), " + std::to_string(bad) + " differ");
        VCAL_CHECK(bad == 0);
    });
}
