// SPDX-License-Identifier: Apache-2.0
// The erfcx lane kernels (P-7, D-210; core/special/erfcx_batch.hpp): every kernel this CPU and build have (scalar,
// SSE2, AVX2), selected in turn on one machine, bit for bit against the scalar erfcx: a fine sweep over all four
// pieces, every piece boundary and its neighbouring doubles, the signed zeros, the subnormals and huge arguments,
// and random arguments, at lengths that leave remainders. Then the quadrature's node batches: the binomial
// log-integrand's batch against its operator() node by node, and the per-period log-likelihood with its score and
// Hessian (log_contrib, log_contrib_derivs) under every kernel against the path with no kernel selected, which is
// the scalar erfcx node by node. And the selection by VCAL_ERFCX_KERNEL (a kernel this machine lacks skips
// visibly, exit code 77).
#include "core/special/erfcx_batch.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "core/model/binomial_mixture.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "core/special/erfcx.hpp"
#include "core/special/generated/erfcx_chebyshev.hpp"
#include "tests/harness/vcal_test.hpp"

namespace {

namespace k = vcal::special::erfcx_batch;

const k::Kernel kInitial = k::selected();
constexpr k::Kernel kAll[] = {k::Kernel::Scalar, k::Kernel::Sse2, k::Kernel::Avx2};

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

struct Restore {
    k::Kernel saved = k::selected();
    ~Restore() { k::select(saved); }
};

template <class F>
void for_each_kernel(F fn) {
    Restore restore;
    for (const auto kern : kAll) {
        if (!k::select(kern)) {
            vcal::test::note(std::string("erfcx kernel ") + k::name(kern) + " not available here");
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

std::vector<double> erfcx_arguments() {
    std::vector<double> x;
    for (int i = 0; i <= 400000; ++i) x.push_back(40.0 * i / 400000.0);  // all four pieces, finely
    for (const double edge : {1.0, 3.0, 8.0}) {                       // each boundary and its neighbours
        double lo = edge, hi = edge;
        x.push_back(edge);
        for (int i = 0; i < 64; ++i) {
            lo = std::nextafter(lo, 0.0);
            hi = std::nextafter(hi, HUGE_VAL);
            x.push_back(lo);
            x.push_back(hi);
        }
    }
    for (const double v : {0.0, -0.0, std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::min(),
                           1e-300, 1e-17, 1e-8, 0.5, 26.6, 1e8, 1e150, 1e300, std::numeric_limits<double>::max()}) {
        x.push_back(v);
    }
    Lcg g{0x5EED0007u};
    for (int i = 0; i < 400000; ++i) x.push_back(std::exp(-40.0 + 120.0 * g.uniform()));  // log-uniform, 4e-18 to 1e34
    return x;
}

}  // namespace

// VCAL_ERFCX_KERNEL selects the kernel at start-up; a kernel this CPU or build lacks makes the test skip (exit 77,
// CTest SKIP_RETURN_CODE), never pass without having run that path.
VCAL_TEST(erfcx_kernel_selection_follows_vcal_erfcx_kernel) {
    const std::string env = vcal::test::env_or_empty("VCAL_ERFCX_KERNEL");
    for (const auto kern : kAll) {
        if (env == k::name(kern) && !k::available(kern)) {
            std::printf("SKIPPED: VCAL_ERFCX_KERNEL=%s, but this CPU or build has no %s kernel\n", env.c_str(),
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
    vcal::test::note("VCAL_ERFCX_KERNEL '" + env + "': selected at start-up " + k::name(kInitial));
    VCAL_CHECK(kInitial == want);
#if defined(__x86_64__) || defined(_M_X64)
    VCAL_CHECK(k::available(k::Kernel::Sse2));
#endif
    Restore restore;
    VCAL_REQUIRE(k::select(k::Kernel::Unavailable));
    double x = 1.0, y = 0.0;
    VCAL_CHECK(!k::nonnegative(&x, &y, 1));
    VCAL_CHECK(!k::select(static_cast<k::Kernel>(7)));
}

// Every kernel gives generated::erfcx_nonnegative, and so erfcx, bit for bit, at every length 1 to 9 (whole blocks
// and remainders) and over the whole argument set at once.
VCAL_TEST(erfcx_kernels_equal_the_scalar_erfcx_bit_for_bit) {
    const std::vector<double> x = erfcx_arguments();
    std::vector<double> want(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        want[i] = vcal::special::generated::erfcx_nonnegative(x[i]);
        VCAL_REQUIRE(bits_equal(want[i], vcal::special::erfcx(x[i])) || x[i] == HUGE_VAL);
    }
    for_each_kernel([&](k::Kernel kern) {
        std::vector<double> out(x.size(), -1.0);
        VCAL_REQUIRE(k::nonnegative(x.data(), out.data(), x.size()));
        std::size_t bad = 0;
        for (std::size_t i = 0; i < x.size(); ++i) bad += bits_equal(out[i], want[i]) ? 0 : 1;
        for (std::size_t len = 1; len <= 9; ++len) {
            for (std::size_t start = 0; start + len <= 4096; start += 517) {
                double o[9];
                VCAL_REQUIRE(k::nonnegative(x.data() + start, o, len));
                for (std::size_t j = 0; j < len; ++j) bad += bits_equal(o[j], want[start + j]) ? 0 : 1;
            }
        }
        vcal::test::note(std::string(k::name(kern)) + ": " + std::to_string(x.size()) + " arguments, " +
                         std::to_string(bad) + " differ");
        VCAL_CHECK(bad == 0);
    });
}

// The binomial log-integrand's batch is its operator() node by node, bit for bit, for 0 < d < n, d = 0 and d = n,
// at node counts 1 to the largest rule, with nodes far enough out that |x| overflows to inf (the scalar fallback).
VCAL_TEST(binomial_integrand_batch_equals_operator_bit_for_bit) {
    namespace m = vcal::model;
    Lcg g{0xBA7C4u};
    std::vector<double> z(256);
    for_each_kernel([&](k::Kernel kern) {
        Lcg gk = g;
        std::size_t bad = 0, nodes = 0;
        for (int rep = 0; rep < 600; ++rep) {
            const double pd = std::exp(-9.0 + 8.9 * gk.uniform());
            const double rho = 0.001 + 0.99 * gk.uniform();
            const std::int64_t n = 1 + static_cast<std::int64_t>(std::exp(14.0 * gk.uniform()));
            std::int64_t d = static_cast<std::int64_t>(gk.uniform() * static_cast<double>(n + 1));
            if (rep % 5 == 0) d = 0;
            if (rep % 5 == 1) d = n;
            if (d > n) d = n;
            const auto f = m::make_binomial_log_integrand(m::make_vasicek1f(pd, rho), n, d);
            const int count = 1 + rep % 256;
            for (int i = 0; i < count; ++i) z[static_cast<std::size_t>(i)] = -12.0 + 24.0 * gk.uniform();
            if (rep % 7 == 0) z[0] = 1e308;  // x = -inf: the scalar path
            if (rep % 11 == 0) z[static_cast<std::size_t>(count - 1)] = -1e308;
            double out[256];
            f.batch(z.data(), out, count);
            for (int i = 0; i < count; ++i) {
                const double want = f(z[static_cast<std::size_t>(i)]);
                bad += (bits_equal(out[i], want) || (std::isnan(out[i]) && std::isnan(want))) ? 0 : 1;
            }
            nodes += static_cast<std::size_t>(count);
        }
        vcal::test::note(std::string(k::name(kern)) + ": " + std::to_string(nodes) + " nodes, " +
                         std::to_string(bad) + " differ");
        VCAL_CHECK(bad == 0);
    });
}

// The per-period log-likelihood and its score and Hessian through the parity rule (and its doubled check rule) are
// the same under every kernel as with none selected, where every node takes the scalar erfcx.
VCAL_TEST(log_contrib_is_the_same_under_every_kernel) {
    using Obj = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
    const Obj obj{};
    struct Row {
        double l, grad[2], hess[3];
    };
    std::vector<Obj::Obs> obs;
    std::vector<Obj::Theta> th;
    Lcg g{0xC0DE7u};
    for (int i = 0; i < 400; ++i) {
        const std::int64_t n = 1 + static_cast<std::int64_t>(std::exp(13.0 * g.uniform()));
        std::int64_t d = static_cast<std::int64_t>(g.uniform() * static_cast<double>(n + 1));
        if (d > n) d = n;
        obs.push_back({n, d});
        th.push_back({std::exp(-8.0 + 7.8 * g.uniform()), 0.002 + 0.6 * g.uniform()});
    }
    auto run = [&](bool doubled) {
        const auto rule = vcal::quadrature::parity_rule(doubled);
        std::vector<Row> rows(obs.size());
        for (std::size_t i = 0; i < obs.size(); ++i) {
            rows[i].l = obj.log_contrib(obs[i], th[i], rule);
            const double l2 = obj.log_contrib_derivs(obs[i], th[i], rule, rows[i].grad, rows[i].hess);
            VCAL_CHECK(bits_equal(l2, rows[i].l) || (std::isnan(l2) && std::isnan(rows[i].l)));
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
                bad += bits_equal(got[i].l, want[i].l) ? 0 : 1;
                for (int j = 0; j < 2; ++j) bad += bits_equal(got[i].grad[j], want[i].grad[j]) ? 0 : 1;
                for (int j = 0; j < 3; ++j) bad += bits_equal(got[i].hess[j], want[i].hess[j]) ? 0 : 1;
            }
            vcal::test::note(std::string(k::name(kern)) + (doubled ? " (doubled rule)" : "") + ": " +
                             std::to_string(got.size()) + " periods, " + std::to_string(bad) + " values differ");
            VCAL_CHECK(bad == 0);
        });
    }
}

