// SPDX-License-Identifier: Apache-2.0
// The AR(1) filter's transition kernels (P-18, D-199; benchmark_models/ar1/transition.hpp): every kernel this CPU
// and build have (scalar, SSE2, AVX2) selected in turn on one machine, each against a plain reference written from
// the contract, bit for bit: the vector det_exp over its whole range (the subnormal band and the edges included),
// the forward and backward log-sum-exps and the shifted exps at sizes that leave remainders, with -inf, NaN and
// +inf terms; the selection by VCAL_AR1_KERNEL (a kernel this machine lacks skips visibly, exit code 77); and the
// threshold below which exp is exactly +0.
#include "benchmark_models/ar1/transition.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "common/detmath/det_exp.hpp"
#include "tests/harness/vcal_test.hpp"

namespace {

namespace k = vcal::benchmark_models::ar1::kernel;
using vcal::test::describe;

const k::Kernel kInitial = k::selected();
constexpr k::Kernel kAll[] = {k::Kernel::Scalar, k::Kernel::Sse2, k::Kernel::Avx2};

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

template <class F>
void for_each_kernel(F fn) {
    struct Restore {
        k::Kernel saved = k::selected();
        ~Restore() { k::select(saved); }
    } restore;
    for (const auto kern : kAll) {
        if (!k::select(kern)) {
            vcal::test::note(std::string("AR(1) kernel ") + k::name(kern) + " not available here");
            continue;
        }
        fn(kern);
    }
}

struct Lcg {  // deterministic inputs; no dependence on the DGP
    std::uint64_t s;
    double uniform() {
        s = s * 6364136223846793005ull + 1442695040888963407ull;
        return static_cast<double>(s >> 11) * 0x1p-53;
    }
};

// The contract, written plainly (transition.hpp): core's online log-sum-exp with det_exp.
struct RefLse {
    double mx = -HUGE_VAL, sum = 0.0;
    unsigned char flags = 0;
    void add(double x) {
        if (std::isnan(x)) {
            flags |= k::kFlagNan;
        } else if (x == HUGE_VAL) {
            flags |= k::kFlagInf;
        } else if (x == -HUGE_VAL) {
        } else if (x > mx) {
            sum = sum * vcal::detmath::det_exp(mx - x) + 1.0;
            mx = x;
        } else {
            sum += vcal::detmath::det_exp(x - mx);
        }
    }
};

double ratio(double z, double a, double pu, double c, double hlv) {
    const double d = z - pu;
    return a - ((d * d) * c + hlv);
}

}  // namespace

// VCAL_AR1_KERNEL selects the kernel at start-up; a kernel this CPU or build lacks makes the test skip (exit 77,
// CTest SKIP_RETURN_CODE), never pass without having run that path.
VCAL_TEST(ar1_kernel_selection_follows_vcal_ar1_kernel) {
    const std::string env = vcal::test::env_or_empty("VCAL_AR1_KERNEL");
    for (const auto kern : kAll) {
        if (env == k::name(kern) && !k::available(kern)) {
            std::printf("SKIPPED: VCAL_AR1_KERNEL=%s, but this CPU or build has no %s kernel\n", env.c_str(), env.c_str());
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
    vcal::test::note("VCAL_AR1_KERNEL '" + env + "': selected at start-up " + k::name(kInitial));
    VCAL_CHECK(kInitial == want);
#if defined(__x86_64__) || defined(_M_X64)
    VCAL_CHECK(k::available(k::Kernel::Sse2));
#endif
    struct Restore {
        k::Kernel saved = k::selected();
        ~Restore() { k::select(saved); }
    } restore;
    VCAL_REQUIRE(k::select(k::Kernel::Unavailable));
    double x = 0.0, y = 0.0;
    VCAL_CHECK(!k::exp_lanes(&x, 1, &y));
    VCAL_CHECK(!k::select(static_cast<k::Kernel>(7)));
}

// The vector det_exp (2^k applied as two exact-then-rounded products, transition_lanes.hpp) is det_exp bit for bit
// at every kernel: a fine sweep of [-810, 810], the subnormal band [-745.2, -708] densely, the cut-offs at +-800
// and +-801, the signed zeros, NaN and the infinities, and random arguments.
VCAL_TEST(ar1_kernel_exp_is_det_exp_bit_for_bit) {
    std::vector<double> xs;
    for (int i = -81000; i <= 81000; ++i) xs.push_back(i * 0.01);
    for (int i = 0; i <= 200000; ++i) xs.push_back(-745.2 + i * (37.2 / 200000.0));
    for (const double e : {800.0, -800.0, 801.0, -801.0, 709.782712893384, 709.79, -708.396418532264, -745.1332191019411,
                           -745.14, 0.0, -0.0, 1e-300, -1e-300, 1e300, -1e300}) {
        xs.push_back(e);
        xs.push_back(std::nextafter(e, HUGE_VAL));
        xs.push_back(std::nextafter(e, -HUGE_VAL));
    }
    xs.push_back(std::numeric_limits<double>::quiet_NaN());
    xs.push_back(HUGE_VAL);
    xs.push_back(-HUGE_VAL);
    Lcg g{12345};
    for (int i = 0; i < 100000; ++i) xs.push_back(-900.0 + 1800.0 * g.uniform());
    std::vector<double> want(xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i) want[i] = vcal::detmath::det_exp(xs[i]);
    for_each_kernel([&](k::Kernel kern) {
        std::vector<double> got(xs.size());
        VCAL_REQUIRE(k::exp_lanes(xs.data(), xs.size(), got.data()));
        std::size_t bad = 0;
        for (std::size_t i = 0; i < xs.size(); ++i) {
            if (!bits_equal(got[i], want[i]) && !(std::isnan(got[i]) && std::isnan(want[i]))) {
                if (++bad <= 5) {
                    vcal::test::report_failure(__FILE__, __LINE__, std::string(k::name(kern)) + ": exp(" + describe(xs[i]) +
                                                                       ") = " + describe(got[i]) + ", det_exp " +
                                                                       describe(want[i]));
                }
            }
        }
        VCAL_CHECK_EQ(bad, std::size_t{0});
        vcal::test::note(std::string(k::name(kern)) + ": " + std::to_string(xs.size()) + " arguments bit for bit");
    });
    // The non-positive path (lanes_exp_nonpos: the LSE terms and the pair weights), reached through exp_shifted with
    // shift 0 on blocks of arguments <= 0: across k = -1021, where the single exact product gives way to the
    // two-product form, through the subnormal band, to the cut-offs and -inf.
    std::vector<double> nonpos;
    for (int i = 0; i <= 81000; ++i) nonpos.push_back(-i * 0.01);
    for (int i = 0; i <= 200000; ++i) nonpos.push_back(-708.0 - i * (37.2 / 200000.0));
    for (const double e : {-707.7, -707.75, -708.05, -708.4, -0.0, 0.0, -800.0, -801.0, -1e300}) {
        nonpos.push_back(e);
        nonpos.push_back(std::nextafter(e, -HUGE_VAL));
    }
    nonpos.push_back(-HUGE_VAL);
    for (int i = 0; i < 100000; ++i) nonpos.push_back(-850.0 * g.uniform());
    for_each_kernel([&](k::Kernel kern) {
        std::vector<double> got(nonpos.size());
        VCAL_REQUIRE(k::exp_shifted(nonpos.data(), nonpos.size(), 0.0, got.data()));
        std::size_t bad = 0;
        for (std::size_t i = 0; i < nonpos.size(); ++i) {
            const double want_i = vcal::detmath::det_exp(nonpos[i] - 0.0);
            if (!bits_equal(got[i], want_i) && ++bad <= 5) {
                vcal::test::report_failure(__FILE__, __LINE__, std::string(k::name(kern)) + ": nonpos exp(" +
                                                                   describe(nonpos[i]) + ") = " + describe(got[i]) +
                                                                   ", det_exp " + describe(want_i));
            }
        }
        VCAL_CHECK_EQ(bad, std::size_t{0});
    });
}

// The skip below kExpZeroBelow is exact only if exp there is exactly +0, for det_exp and for the platform's
// std::exp (core's OnlineLogSumExp) alike.
VCAL_TEST(ar1_exp_is_exactly_zero_below_the_skip_threshold) {
    for (const double x : {k::kExpZeroBelow, std::nextafter(k::kExpZeroBelow, -HUGE_VAL), k::kExpZeroBelow - 1.0,
                           -800.0, -1e300}) {
        const double a = vcal::detmath::det_exp(x), b = std::exp(x);
        VCAL_CHECK(bits_equal(a, 0.0));
        VCAL_CHECK(bits_equal(b, 0.0));
    }
}

// Forward, backward (with its terms written out) and the shifted exps against the plain reference, at every
// kernel and at sizes that leave remainders after blocks of 2 and 4, with -inf weights among the sources.
VCAL_TEST(ar1_kernels_equal_the_reference_bit_for_bit) {
    const double phi = 0.83, var = 1.0 - phi * phi, c = 1.0 / (2.0 * var), hlv = 0.5 * std::log(var);
    int compared = 0;
    for (const std::size_t n : {std::size_t{1}, std::size_t{2}, std::size_t{3}, std::size_t{5}, std::size_t{7},
                                std::size_t{128}, std::size_t{131}, std::size_t{256}}) {
        Lcg g{n * 7919u + 1u};
        const std::size_t m = n + 3;  // the other side's size
        std::vector<double> pu(m), lhat(m), z(n), a(n), v(n), b(n);
        for (std::size_t j = 0; j < m; ++j) {
            pu[j] = phi * (-6.0 + 12.0 * g.uniform());
            lhat[j] = j % 11 == 3 ? -HUGE_VAL : -30.0 * g.uniform();
        }
        for (std::size_t i = 0; i < n; ++i) {
            z[i] = -7.0 + 14.0 * g.uniform();
            a[i] = 0.5 * (z[i] * z[i]) - 0.2;
            v[i] = -40.0 * g.uniform();
            b[i] = i % 13 == 5 ? -HUGE_VAL : -900.0 * g.uniform();  // far terms exercise the skip
        }
        // Forward over m sources to n destinations; backward over n... sources (pu's first n) to m destinations.
        std::vector<double> fm(n), fs(n), bm(n), bs(n), e(n * n), sh(n * n), shx(n * n);
        std::vector<unsigned char> ff(n), bf(n);
        for (std::size_t kk = 0; kk < n; ++kk) {
            RefLse r;
            for (std::size_t j = 0; j < m; ++j) r.add(lhat[j] + ratio(z[kk], a[kk], pu[j], c, hlv));
            fm[kk] = r.mx;
            fs[kk] = r.sum;
            ff[kk] = r.flags;
        }
        for (std::size_t j = 0; j < n; ++j) {
            RefLse r;
            for (std::size_t kk = 0; kk < n; ++kk) {
                const double x = (ratio(z[kk], a[kk], pu[j], c, hlv) + v[kk]) + b[kk];
                e[j * n + kk] = x;
                r.add(x);
            }
            bm[j] = r.mx;
            bs[j] = r.sum;
            bf[j] = r.flags;
        }
        const double shift = 3.25;
        for (std::size_t i = 0; i < n * n; ++i) {
            shx[i] = e[i] == -HUGE_VAL ? -1000.0 : e[i];
            sh[i] = vcal::detmath::det_exp(shx[i] - shift);
        }
        for_each_kernel([&](k::Kernel kern) {
            std::vector<double> gm(n), gs(n), ge(n * n), gsh(n * n);
            std::vector<unsigned char> gf(n);
            VCAL_REQUIRE(k::forward(pu.data(), lhat.data(), m, z.data(), a.data(), n, c, hlv, gm.data(), gs.data(), gf.data()));
            for (std::size_t i = 0; i < n; ++i) {
                VCAL_CHECK(bits_equal(gm[i], fm[i]) && bits_equal(gs[i], fs[i]) && gf[i] == ff[i]);
            }
            VCAL_REQUIRE(k::backward(pu.data(), n, z.data(), a.data(), v.data(), b.data(), n, c, hlv, gm.data(), gs.data(),
                                     gf.data(), ge.data()));
            for (std::size_t i = 0; i < n; ++i) {
                VCAL_CHECK(bits_equal(gm[i], bm[i]) && bits_equal(gs[i], bs[i]) && gf[i] == bf[i]);
            }
            for (std::size_t i = 0; i < n * n; ++i) VCAL_CHECK(bits_equal(ge[i], e[i]));
            VCAL_REQUIRE(k::exp_shifted(shx.data(), n * n, shift, gsh.data()));
            for (std::size_t i = 0; i < n * n; ++i) VCAL_CHECK(bits_equal(gsh[i], sh[i]));
            // The score's pair terms, against the filter's scalar expression.
            std::vector<double> gt(n * n), u(n);
            for (std::size_t j = 0; j < n; ++j) u[j] = pu[j] / phi;
            VCAL_REQUIRE(k::score_terms(sh.data(), n, n, u.data(), z.data(), phi, var, gt.data()));
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t kk = 0; kk < n; ++kk) {
                    const double d = z[kk] - phi * u[j];
                    const double want = sh[j * n + kk] * (d * u[j] / var - phi * (d * d) / (var * var) + phi / var);
                    VCAL_CHECK(bits_equal(gt[j * n + kk], want));
                }
            }
            ++compared;
            (void)kern;
        });
    }
    // Special terms: a NaN and a +inf among the sources set the flags, as OnlineLogSumExp's result reads them.
    for_each_kernel([&](k::Kernel) {
        const double pu[3] = {0.0, 0.1, 0.2}, z[4] = {0.0, 1.0, -1.0, 2.0}, a[4] = {0.0, 0.5, 0.5, 2.0};
        double mx[4], sum[4];
        unsigned char f[4];
        const double nanw[3] = {0.0, std::numeric_limits<double>::quiet_NaN(), 0.0};
        VCAL_REQUIRE(k::forward(pu, nanw, 3, z, a, 4, 1.0, 0.0, mx, sum, f));
        for (int i = 0; i < 4; ++i) VCAL_CHECK(std::isnan(k::lse_result(mx[i], sum[i], f[i])));
        const double infw[3] = {0.0, HUGE_VAL, 0.0};
        VCAL_REQUIRE(k::forward(pu, infw, 3, z, a, 4, 1.0, 0.0, mx, sum, f));
        for (int i = 0; i < 4; ++i) VCAL_CHECK(k::lse_result(mx[i], sum[i], f[i]) == HUGE_VAL);
        const double ninf[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
        VCAL_REQUIRE(k::forward(pu, ninf, 3, z, a, 4, 1.0, 0.0, mx, sum, f));
        for (int i = 0; i < 4; ++i) VCAL_CHECK(k::lse_result(mx[i], sum[i], f[i]) == -HUGE_VAL);
    });
    vcal::test::note(std::to_string(compared) + " (size, kernel) cases bit for bit");
}
