// SPDX-License-Identifier: Apache-2.0
// engine/hac.hpp: the HAC (Newey-West) sandwich SE, a `native` option (S-48, D-209,
// docs/methodology/hac_se.md).
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/calibrate.hpp"
#include "engine/hac.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

#define VCAL_TEST_STR2(x) #x
#define VCAL_TEST_STR(x) VCAL_TEST_STR2(x)
#if defined(_MSC_VER)
#define VCAL_TEST_COMPILER "MSVC " VCAL_TEST_STR(_MSC_VER)
#elif defined(__clang__)
#define VCAL_TEST_COMPILER "Clang " __clang_version__
#elif defined(__GNUC__)
#define VCAL_TEST_COMPILER "GCC " __VERSION__
#else
#define VCAL_TEST_COMPILER "unknown compiler"
#endif

namespace {

namespace e = vcal::engine;
namespace tol = vcal::tol;
using vcal::test::describe;
using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
using Obs = Objective::Obs;

// The recovery harness's grid and rules (tests/recovery/recovery.hpp), which S-48 fits with.
vcal::Grid<2> grid() {
    return {{{1e-4, 0.2, 61, vcal::AxisScale::Logit}, {1e-3, vcal::kDefaultRhoUpper, 41, vcal::AxisScale::Logit}}};
}

// The M1.4 test panel: 20 periods of 1000 obligors.
std::vector<Obs> panel() {
    const std::int64_t d[] = {2, 5, 1, 9, 3, 0, 4, 12, 6, 2, 1, 3, 7, 15, 2, 4, 0, 5, 8, 3};
    std::vector<Obs> p;
    for (const auto x : d) p.push_back({1000, x});
    return p;
}

e::Estimate2 fit(const std::vector<Obs>& p) {
    static const auto primary = vcal::quadrature::parity_rule();
    static const auto check = vcal::quadrature::parity_rule(true);
    std::vector<double> L;
    e::Estimate2 est{};
    VCAL_REQUIRE(e::calibrate(vcal::backends::CpuBackend{1}, Objective{}, primary, check, p.data(),
                              static_cast<std::int64_t>(p.size()), grid(), L, est) == e::Status::Ok);
    return est;
}

e::HacSe2 hac(const std::vector<Obs>& p, const e::Estimate2& est, std::int32_t lag) {
    static const auto primary = vcal::quadrature::parity_rule();
    return e::hac_se(Objective{}, primary, p.data(), static_cast<std::int64_t>(p.size()), grid(), est, lag);
}

}  // namespace

// The rule floor(4 (T / 100)^(2/9)): 2, 3, 4 at S-48's T, and the step from 4 to 5 between T = 272 and 273
// (4 (T / 100)^(2/9) = 5 at T = 100 * 1.25^4.5 = 272.96).
VCAL_TEST(hac_lag_rule) {
    VCAL_CHECK_EQ(e::newey_west_lag(20), 2);
    VCAL_CHECK_EQ(e::newey_west_lag(40), 3);
    VCAL_CHECK_EQ(e::newey_west_lag(100), 4);
    VCAL_CHECK_EQ(e::newey_west_lag(272), 4);
    VCAL_CHECK_EQ(e::newey_west_lag(273), 5);
}

// S by hand on four score vectors s = (1, 0), (-1, 2), (2, -1), (-2, -1), whose mean is 0:
//   lag 0:  sum s s' = (1 + 1 + 4 + 4, 0 + 4 + 1 + 1, 0 - 2 - 2 + 2) = (10, 6, -2);
//   lag 1:  sum_{t>=1} (s_t s_{t-1}' + s_{t-1} s_t'):
//           t = 1: (-2, 0, 2); t = 2: (-4, -4, 5); t = 3: (-8, 2, 0); total (-14, -2, 7);
//   L = 1:  w_1 = 1/2, S = (10 - 7, 6 - 1, -2 + 3.5) = (3, 5, 1.5).
// Longer lags, by hand:
//   lag 2:  t = 2: s_2 = (2, -1), s_0 = (1, 0): (4, 0, -1); t = 3: s_3 = (-2, -1), s_1 = (-1, 2): (4, -4, -3);
//           total (8, -4, -4);
//   lag 3:  t = 3: s_3 = (-2, -1), s_0 = (1, 0): (-4, 0, -1);
//   L = 3:  w = 3/4, 1/2, 1/4: S = (10 - 10.5 + 4 - 1, 6 - 1.5 - 2 + 0, -2 + 5.25 - 2 - 0.25) = (2.5, 2.5, 1).
// Every value is exact in binary, so the checks are exact. Adding a constant to every score changes
// nothing (the scores are centred); lags of T or more contribute nothing.
VCAL_TEST(newey_west_by_hand) {
    const double s[] = {1, 0, -1, 2, 2, -1, -2, -1};
    double out[3];
    e::newey_west(s, 4, 0, out);
    VCAL_CHECK_EQ(out[0], 10.0);
    VCAL_CHECK_EQ(out[1], 6.0);
    VCAL_CHECK_EQ(out[2], -2.0);
    e::newey_west(s, 4, 1, out);
    VCAL_CHECK_EQ(out[0], 3.0);
    VCAL_CHECK_EQ(out[1], 5.0);
    VCAL_CHECK_EQ(out[2], 1.5);
    e::newey_west(s, 4, 3, out);
    VCAL_CHECK_EQ(out[0], 2.5);
    VCAL_CHECK_EQ(out[1], 2.5);
    VCAL_CHECK_EQ(out[2], 1.0);

    const double shifted[] = {2, 1, 0, 3, 3, 0, -1, 0};
    e::newey_west(shifted, 4, 1, out);
    VCAL_CHECK_EQ(out[0], 3.0);
    VCAL_CHECK_EQ(out[1], 5.0);
    VCAL_CHECK_EQ(out[2], 1.5);

    // L = 7: weights 7/8, 6/8, 5/8 on lags 1-3, and no lag 4 or beyond in four periods.
    double far[3];
    e::newey_west(s, 4, 7, far);
    VCAL_CHECK_EQ(far[0], 10.0 + (7.0 / 8.0) * -14.0 + (6.0 / 8.0) * 8.0 + (5.0 / 8.0) * -4.0);
    VCAL_CHECK_EQ(far[1], 6.0 + (7.0 / 8.0) * -2.0 + (6.0 / 8.0) * -4.0 + (5.0 / 8.0) * 0.0);
    VCAL_CHECK_EQ(far[2], -2.0 + (7.0 / 8.0) * 7.0 + (6.0 / 8.0) * -4.0 + (5.0 / 8.0) * -1.0);
}

// On the M1.4 panel, hac_se against an assembly written out separately: per-period scores in logit units
// by central differences of each period's log-likelihood (log_contrib, not the analytic derivatives the
// engine uses), H by central differences of their sum, and S by the textbook double sum
// sum_{t,r} k(|t - r|) c_t c_r' with the Bartlett weights k(m) = 1 - m / (L + 1) for m <= L. The difference
// steps (1e-4 in logit units) leave about 1e-7 relative error.
VCAL_TEST(hac_matches_an_independent_assembly) {
    const auto p = panel();
    const auto est = fit(p);
    VCAL_REQUIRE(!(est.flags & e::kHacSkipFlags));
    static const auto primary = vcal::quadrature::parity_rule();
    const auto g = grid();
    const double u0[2] = {vcal::grid::to_scaled(g.axis[0].scale, est.value[0]),
                          vcal::grid::to_scaled(g.axis[1].scale, est.value[1])};
    const auto lt = [&](const Obs& y, double a, double b) {
        const double v[2] = {vcal::grid::from_scaled(g.axis[0].scale, a), vcal::grid::from_scaled(g.axis[1].scale, b)};
        return Objective{}.log_contrib(y, Objective::theta(v), primary);
    };
    const double h = 1e-4;
    const auto T = static_cast<std::int64_t>(p.size());
    // score_t at (a, b), and the panel score as the sum.
    const auto score = [&](const Obs& y, double a, double b, double (&s)[2]) {
        s[0] = (lt(y, a + h, b) - lt(y, a - h, b)) / (2 * h);
        s[1] = (lt(y, a, b + h) - lt(y, a, b - h)) / (2 * h);
    };
    std::vector<double> s0(2 * p.size());
    for (std::size_t t = 0; t < p.size(); ++t) {
        double s[2];
        score(p[t], u0[0], u0[1], s);
        s0[2 * t] = s[0];
        s0[2 * t + 1] = s[1];
    }
    // H from the panel score's differences.
    double H[2][2];
    for (int b = 0; b < 2; ++b) {
        double up[2] = {0, 0}, dn[2] = {0, 0};
        for (const auto& y : p) {
            double s[2];
            score(y, u0[0] + (b == 0 ? h : 0.0), u0[1] + (b == 1 ? h : 0.0), s);
            up[0] += s[0];
            up[1] += s[1];
            score(y, u0[0] - (b == 0 ? h : 0.0), u0[1] - (b == 1 ? h : 0.0), s);
            dn[0] += s[0];
            dn[1] += s[1];
        }
        H[0][b] = (up[0] - dn[0]) / (2 * h);
        H[1][b] = (up[1] - dn[1]) / (2 * h);
    }
    const double h01 = 0.5 * (H[0][1] + H[1][0]);
    const double det = H[0][0] * H[1][1] - h01 * h01;
    const double A[2][2] = {{H[1][1] / det, -h01 / det}, {-h01 / det, H[0][0] / det}};
    for (const std::int32_t L : {0, 1, 2, 4}) {
        double mean[2] = {0, 0};
        for (std::int64_t t = 0; t < T; ++t) {
            mean[0] += s0[static_cast<std::size_t>(2 * t)] / static_cast<double>(T);
            mean[1] += s0[static_cast<std::size_t>(2 * t + 1)] / static_cast<double>(T);
        }
        double S[2][2] = {{0, 0}, {0, 0}};
        for (std::int64_t t = 0; t < T; ++t) {
            for (std::int64_t r = 0; r < T; ++r) {
                const std::int64_t m = t > r ? t - r : r - t;
                if (m > L) continue;
                const double k = 1.0 - static_cast<double>(m) / static_cast<double>(L + 1);
                for (int a = 0; a < 2; ++a) {
                    for (int b = 0; b < 2; ++b) {
                        S[a][b] += k * (s0[static_cast<std::size_t>(2 * t + a)] - mean[a]) *
                                   (s0[static_cast<std::size_t>(2 * r + b)] - mean[b]);
                    }
                }
            }
        }
        double V[2][2] = {{0, 0}, {0, 0}};
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b)
                for (int c = 0; c < 2; ++c)
                    for (int d = 0; d < 2; ++d) V[a][b] += A[a][c] * S[c][d] * A[d][b];
        const auto got = hac(p, est, L);
        VCAL_REQUIRE(got.flags == 0u);
        VCAL_CHECK_EQ(got.lag, L);
        const double rel0 = std::fabs(got.se_scaled[0] / std::sqrt(V[0][0]) - 1.0);
        const double rel1 = std::fabs(got.se_scaled[1] / std::sqrt(V[1][1]) - 1.0);
        const double dcorr = std::fabs(got.corr - V[0][1] / std::sqrt(V[0][0] * V[1][1]));
        vcal::test::note("L = " + std::to_string(L) + ": se_u " + describe(got.se_scaled[0]) + ", " +
                         describe(got.se_scaled[1]) + "; relative differences " + describe(rel0) + ", " +
                         describe(rel1) + "; corr difference " + describe(dcorr));
        VCAL_CHECK(rel0 <= 1e-6 && rel1 <= 1e-6 && dcorr <= 1e-6);
        for (int a = 0; a < 2; ++a) {
            const double j = vcal::grid::dvalue_dscaled(g.axis[a].scale, u0[a]);
            VCAL_CHECK_EQ(got.se[a], j * got.se_scaled[a]);
        }
    }
}

// Reversing the periods leaves S unchanged (each lag's sum is symmetric in the pair), so the SE agrees to
// rounding; any other reordering changes S once L > 0, and with L = 0 it only reorders the sums.
VCAL_TEST(hac_depends_on_the_order_only_through_the_lags) {
    const auto p = panel();
    auto rev = p;
    std::vector<Obs>(p.rbegin(), p.rend()).swap(rev);
    auto shuffled = p;
    for (std::size_t t = 0; t + 1 < shuffled.size(); t += 2) std::swap(shuffled[t], shuffled[t + 1]);
    const auto est = fit(p);  // one estimate for all three orders: only S is compared
    for (int a = 0; a < 2; ++a) {
        VCAL_CHECK_REL(hac(rev, est, 2).se[a], hac(p, est, 2).se[a], 1e-13);
        VCAL_CHECK_REL(hac(shuffled, est, 0).se[a], hac(p, est, 0).se[a], 1e-13);
        VCAL_CHECK(std::fabs(hac(shuffled, est, 2).se[a] / hac(p, est, 2).se[a] - 1.0) > 1e-3);
    }
}

// A fit without an iid Wald interval gets no HAC SE either: every value NaN, kIntervalNotComputed.
VCAL_TEST(hac_not_computed_on_flagged_fits) {
    std::vector<Obs> zeros(20, Obs{1000, 0});
    const auto est = fit(zeros);
    VCAL_REQUIRE((est.flags & e::kHacSkipFlags) != 0u);
    const auto got = hac(zeros, est, 2);
    VCAL_CHECK_EQ(got.flags, static_cast<std::uint32_t>(e::kIntervalNotComputed));
    VCAL_CHECK(std::isnan(got.se[0]) && std::isnan(got.se[1]) && std::isnan(got.corr));
    const auto w = e::hac_wald_interval(grid(), est, got, 0, 1.959963984540054);
    VCAL_CHECK(std::isnan(w.lo) && std::isnan(w.hi));
}

// Regression values on the six M1.7 validation panels (tests/golden/validation/panels.csv), the input of
// validation/scipy/hac_se.py: the estimate, the iid SEs, and the HAC SEs in logit units at L = 1, the rule
// (2 at T = 20) and twice the rule. Set VCAL_WRITE_REGRESSION_GOLDENS=1 to rewrite tests/golden/hac/
// validation.csv from this platform. Compared on every platform at TOL_XREF_CROSS_PLATFORM_REL (estimate)
// and TOL_XREF_CROSS_PLATFORM_SE_REL (SEs).
VCAL_TEST(hac_matches_regression_values) {
    const auto panels_csv = vcal::test::read_golden_csv("validation/panels.csv");
    std::vector<std::vector<Obs>> panels;
    for (const auto& r : panels_csv.rows) {
        const auto s = static_cast<std::size_t>(vcal::test::parse_int(r[panels_csv.column("scenario")]));
        if (panels.size() <= s) panels.resize(s + 1);
        panels[s].push_back({vcal::test::parse_int(r[panels_csv.column("n")]),
                             vcal::test::parse_int(r[panels_csv.column("d")])});
    }
    VCAL_REQUIRE(panels.size() == 6);
    const auto lags = [](std::int64_t T) {
        const std::int32_t rule = e::newey_west_lag(T);
        return std::vector<std::int32_t>{1, rule, 2 * rule};
    };
    const std::string path = vcal::test::golden_path("hac/validation.csv");
    if (vcal::test::env_or_empty("VCAL_WRITE_REGRESSION_GOLDENS") == "1") {
        std::ofstream out(path, std::ios::binary);
        out << "# Regression values: the HAC SE (engine/hac.hpp, D-209) on the six M1.7 validation panels, the input\n"
               "# of validation/scipy/hac_se.py. Recovery grid and parity rules; SEs in logit units. NOT a correctness\n"
               "# reference. Written by unit_hac hac_matches_regression_values with VCAL_WRITE_REGRESSION_GOLDENS=1.\n"
               "# Written by: " VCAL_TEST_COMPILER "\n"
               "scenario,lag,pd_hex,rho_hex,flags,iid_se_u_pd_hex,iid_se_u_rho_hex,iid_corr_hex,hac_se_u_pd_hex,"
               "hac_se_u_rho_hex,hac_corr_hex\n";
        for (std::size_t s = 0; s < panels.size(); ++s) {
            const auto est = fit(panels[s]);
            const auto g = grid();
            for (const std::int32_t L : lags(static_cast<std::int64_t>(panels[s].size()))) {
                const auto got = hac(panels[s], est, L);
                double iid[2];
                for (int a = 0; a < 2; ++a) {
                    iid[a] = est.se[a] / vcal::grid::dvalue_dscaled(
                                             g.axis[a].scale, vcal::grid::to_scaled(g.axis[a].scale, est.value[a]));
                }
                out << s << ',' << L << ',' << vcal::test::to_hex(est.value[0]) << ','
                    << vcal::test::to_hex(est.value[1]) << ',' << est.flags << ',' << vcal::test::to_hex(iid[0]) << ','
                    << vcal::test::to_hex(iid[1]) << ',' << vcal::test::to_hex(est.corr) << ','
                    << vcal::test::to_hex(got.se_scaled[0]) << ',' << vcal::test::to_hex(got.se_scaled[1]) << ','
                    << vcal::test::to_hex(got.corr) << '\n';
            }
        }
        vcal::test::note("wrote " + path);
        return;
    }
    const auto t = vcal::test::read_golden_csv("hac/validation.csv");
    std::size_t i = 0;
    double worst_est = 0.0, worst_se = 0.0;
    for (std::size_t s = 0; s < panels.size(); ++s) {
        const auto est = fit(panels[s]);
        for (const std::int32_t L : lags(static_cast<std::int64_t>(panels[s].size()))) {
            VCAL_REQUIRE(i < t.rows.size());
            const auto& r = t.rows[i++];
            VCAL_CHECK_EQ(vcal::test::parse_int(r[t.column("scenario")]), static_cast<std::int64_t>(s));
            VCAL_CHECK_EQ(vcal::test::parse_int(r[t.column("lag")]), static_cast<std::int64_t>(L));
            VCAL_CHECK_EQ(vcal::test::parse_int(r[t.column("flags")]), static_cast<std::int64_t>(est.flags));
            const auto got = hac(panels[s], est, L);
            for (int a = 0; a < 2; ++a) {
                const char* vc = a == 0 ? "pd_hex" : "rho_hex";
                const char* sc = a == 0 ? "hac_se_u_pd_hex" : "hac_se_u_rho_hex";
                worst_est = std::fmax(worst_est, std::fabs(est.value[a] / vcal::test::parse_double(r[t.column(vc)]) - 1.0));
                worst_se = std::fmax(worst_se, std::fabs(got.se_scaled[a] / vcal::test::parse_double(r[t.column(sc)]) - 1.0));
            }
            worst_se = std::fmax(worst_se, std::fabs(got.corr - vcal::test::parse_double(r[t.column("hac_corr_hex")])));
        }
    }
    VCAL_CHECK_EQ(i, t.rows.size());
    vcal::test::note("worst relative difference: estimate " + describe(worst_est) + ", HAC SE " + describe(worst_se));
    VCAL_CHECK(worst_est <= tol::TOL_XREF_CROSS_PLATFORM_REL);
    VCAL_CHECK(worst_se <= tol::TOL_XREF_CROSS_PLATFORM_SE_REL);
}
