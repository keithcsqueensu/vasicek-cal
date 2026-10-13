// SPDX-License-Identifier: Apache-2.0
// The derivative node-function kernels (P-22, D-402; core/model/derivs_batch.hpp): every kernel this CPU and build
// have (scalar, SSE2, AVX2), selected in turn on one machine, bit for bit against the derivative functor's own
// per-node path (BinomialLogIntegrandDerivs without a node cache). On random periods (0 < d < n, d = 0, d = n) and
// node counts 1 to 256 (whole blocks and remainders), with nodes whose |x| overflows (no cached e: the functor's own
// path), under every erfcx kernel and none (P-7, which decides what the cache holds). Then the functor read through
// the cache, with skipped nodes as the accumulator skips underflowed ones; the per-period log-likelihood with its
// score and Hessian (log_contrib_derivs) under every kernel against none, on the parity rule and its doubled rule;
// and the selection by VCAL_DERIVS_KERNEL (a kernel this machine lacks skips visibly, exit code 77).
#include "core/model/derivs_batch.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "core/model/binomial_mixture.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "core/special/erfcx_batch.hpp"
#include "tests/harness/vcal_test.hpp"

namespace {

namespace k = vcal::model::derivs_batch;
namespace ek = vcal::special::erfcx_batch;
namespace m = vcal::model;

const k::Kernel kInitial = k::selected();
constexpr k::Kernel kAll[] = {k::Kernel::Scalar, k::Kernel::Sse2, k::Kernel::Avx2};

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
bool same(double a, double b) { return bits_equal(a, b) || (std::isnan(a) && std::isnan(b)); }

struct Restore {
    k::Kernel saved = k::selected();
    ek::Kernel saved_e = ek::selected();
    ~Restore() {
        k::select(saved);
        ek::select(saved_e);
    }
};

template <class F>
void for_each_kernel(F fn) {
    Restore restore;
    for (const auto kern : kAll) {
        if (!k::select(kern)) {
            vcal::test::note(std::string("derivs kernel ") + k::name(kern) + " not available here");
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

struct Period {
    double pd, rho;
    std::int64_t n, d;
    int count;
    std::vector<double> z;
};

std::vector<Period> periods() {
    Lcg g{0xD22u};
    std::vector<Period> out;
    for (int rep = 0; rep < 600; ++rep) {
        Period p;
        p.pd = std::exp(-9.0 + 8.9 * g.uniform());
        p.rho = 0.001 + 0.99 * g.uniform();
        p.n = 1 + static_cast<std::int64_t>(std::exp(14.0 * g.uniform()));
        p.d = static_cast<std::int64_t>(g.uniform() * static_cast<double>(p.n + 1));
        if (rep % 5 == 0) p.d = 0;
        if (rep % 5 == 1) p.d = p.n;
        if (p.d > p.n) p.d = p.n;
        p.count = 1 + rep % 256;
        p.z.resize(static_cast<std::size_t>(p.count));
        for (double& zi : p.z) zi = -12.0 + 24.0 * g.uniform();
        if (rep % 7 == 0) p.z[0] = 1e308;  // |x| overflows: no cached e
        if (rep % 11 == 0) p.z[static_cast<std::size_t>(p.count - 1)] = -1e308;
        if (rep % 13 == 0) p.z[static_cast<std::size_t>(p.count / 2)] = 0.0;  // x = c / b, and z = 0 in x_rho
        out.push_back(std::move(p));
    }
    return out;
}

}  // namespace

// VCAL_DERIVS_KERNEL selects the kernel at start-up; a kernel this CPU or build lacks makes the test skip (exit 77,
// CTest SKIP_RETURN_CODE), never pass without having run that path.
VCAL_TEST(derivs_kernel_selection_follows_vcal_derivs_kernel) {
    const std::string env = vcal::test::env_or_empty("VCAL_DERIVS_KERNEL");
    for (const auto kern : kAll) {
        if (env == k::name(kern) && !k::available(kern)) {
            std::printf("SKIPPED: VCAL_DERIVS_KERNEL=%s, but this CPU or build has no %s kernel\n", env.c_str(),
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
    vcal::test::note("VCAL_DERIVS_KERNEL '" + env + "': selected at start-up " + k::name(kInitial));
    VCAL_CHECK(kInitial == want);
#if defined(__x86_64__) || defined(_M_X64)
    VCAL_CHECK(k::available(k::Kernel::Sse2));
#endif
    Restore restore;
    VCAL_REQUIRE(k::select(k::Kernel::Unavailable));
    const k::Params p{0.0, 0.5, 0.5, 0.25, 1.0, 1.0, 1.0, 1.0};
    double z = 0.0, e = 1.0, E = 1.0, f[5];
    VCAL_CHECK(!k::node_functions(p, &z, &e, &E, f, 1));
    VCAL_CHECK(!k::select(static_cast<k::Kernel>(7)));
}

// Every kernel's node functions, filled into the cache after the value's batch, are the functor's own per-node
// values bit for bit at every cached node; and the functor read through the cache (skipping nodes as the accumulator
// does) is the same. Under every erfcx kernel and none, so the cache holds e for all nodes, some, or none.
VCAL_TEST(derivs_kernels_equal_the_functor_bit_for_bit) {
    const std::vector<Period> ps = periods();
    for_each_kernel([&](k::Kernel kern) {
        for (const auto ekern : {ek::Kernel::Scalar, ek::Kernel::Sse2, ek::Kernel::Avx2, ek::Kernel::Unavailable}) {
            if (!ek::select(ekern)) continue;
            std::size_t bad = 0, filled = 0, through = 0;
            for (const Period& p : ps) {
                m::NodeCache cache;
                auto f = m::make_binomial_log_integrand(m::make_vasicek1f(p.pd, p.rho), p.n, p.d);
                f.cache = &cache;
                auto df = m::make_binomial_log_integrand_derivs(p.pd, p.rho, p.n, p.d);
                const auto plain = df;  // no cache: the per-node path
                df.cache = &cache;
                double out[256];
                f.batch(p.z.data(), out, p.count);
                // The kernel's values directly.
                df.fill_node_functions(cache);
                VCAL_REQUIRE(cache.f_state == 1);
                for (int i = 0; i < p.count; ++i) {
                    if (std::isnan(cache.e[i])) continue;
                    double want[5];
                    plain(p.z[static_cast<std::size_t>(i)], want);
                    for (int j = 0; j < 5; ++j) bad += same(cache.f[5 * i + j], want[j]) ? 0 : 1;
                    ++filled;
                }
                // Through operator(), from a fresh batch (f_state reset), skipping some nodes.
                f.batch(p.z.data(), out, p.count);
                VCAL_REQUIRE(cache.f_state == 0);
                const bool skip = (p.count % 3) == 0;
                for (int i = 0; i < p.count; ++i) {
                    if (skip && i % 4 == 1) continue;
                    double got[5], want[5];
                    df(p.z[static_cast<std::size_t>(i)], got);
                    plain(p.z[static_cast<std::size_t>(i)], want);
                    for (int j = 0; j < 5; ++j) bad += same(got[j], want[j]) ? 0 : 1;
                    ++through;
                }
            }
            vcal::test::note(std::string(k::name(kern)) + " (erfcx " + ek::name(ekern) + "): " +
                             std::to_string(filled) + " filled nodes, " + std::to_string(through) +
                             " nodes through the functor, " + std::to_string(bad) + " values differ");
            VCAL_CHECK(bad == 0);
        }
    });
}

// The per-period log-likelihood, score and Hessian through the parity rule (and its doubled check rule) are the same
// under every kernel as with none selected, where every node takes the functor's own path.
VCAL_TEST(log_contrib_derivs_is_the_same_under_every_derivs_kernel) {
    using Obj = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
    const Obj obj{};
    struct Row {
        double l, grad[2], hess[3];
    };
    std::vector<Obj::Obs> obs;
    std::vector<Obj::Theta> th;
    Lcg g{0xC0DE22u};
    for (int i = 0; i < 400; ++i) {
        const std::int64_t n = 1 + static_cast<std::int64_t>(std::exp(13.0 * g.uniform()));
        std::int64_t d = static_cast<std::int64_t>(g.uniform() * static_cast<double>(n + 1));
        if (i % 6 == 0) d = 0;
        if (d > n) d = n;
        obs.push_back({n, d});
        th.push_back({std::exp(-8.0 + 7.8 * g.uniform()), 0.002 + 0.6 * g.uniform()});
    }
    auto run = [&](bool doubled) {
        const auto rule = vcal::quadrature::parity_rule(doubled);
        std::vector<Row> rows(obs.size());
        for (std::size_t i = 0; i < obs.size(); ++i) {
            rows[i].l = obj.log_contrib_derivs(obs[i], th[i], rule, rows[i].grad, rows[i].hess);
        }
        return rows;
    };
    Restore restore;
    for (const bool doubled : {false, true}) {
        VCAL_REQUIRE(k::select(k::Kernel::Unavailable));
        const std::vector<Row> want = run(doubled);
        for_each_kernel([&](k::Kernel kern) {
            const std::vector<Row> got = run(doubled);
            std::size_t bad = 0;
            for (std::size_t i = 0; i < got.size(); ++i) {
                bad += same(got[i].l, want[i].l) ? 0 : 1;
                for (int j = 0; j < 2; ++j) bad += same(got[i].grad[j], want[i].grad[j]) ? 0 : 1;
                for (int j = 0; j < 3; ++j) bad += same(got[i].hess[j], want[i].hess[j]) ? 0 : 1;
            }
            vcal::test::note(std::string(k::name(kern)) + (doubled ? " (doubled rule)" : "") + ": " +
                             std::to_string(got.size()) + " periods, " + std::to_string(bad) + " values differ");
            VCAL_CHECK(bad == 0);
        });
    }
}
