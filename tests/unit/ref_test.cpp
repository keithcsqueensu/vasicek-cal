// SPDX-License-Identifier: Apache-2.0
// ref/ against the mpmath golden tables (M1.5; D-031). The goldens are data, not code, so
// measuring ref against them keeps ref independent of core while establishing its accuracy.
#include "ref/ref.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

using vcal::test::describe;
using vcal::test::parse_double;
using vcal::test::parse_int;
using vcal::test::UlpStats;
namespace tol = vcal::tol;

constexpr double kInf = std::numeric_limits<double>::infinity();

// Same metric as the quadrature tests: |a - e| / (eps * max(1, |e|)) on the log scale.
double eps_units(double actual, double expected) {
    return std::fabs(actual - expected) / (DBL_EPSILON * std::fmax(1.0, std::fabs(expected)));
}

}  // namespace

VCAL_TEST(ref_log_ncdf_matches_mpmath) {
    const auto t = vcal::test::read_golden_csv("special/log_phi.csv");
    const auto cx = t.column("x_hex");
    const auto ce = t.column("expected_hex");
    UlpStats cf, erfc_lower, upper;
    double subnormal_worst = 0.0;
    for (const auto& r : t.rows) {
        const double x = parse_double(r[cx]);
        const double e = parse_double(r[ce]);
        const double a = vcalref::log_ncdf(x);
        if (std::isinf(e)) {
            VCAL_CHECK_EQ(a, e);
        } else if (std::fabs(e) < DBL_MIN) {
            subnormal_worst = std::fmax(subnormal_worst, std::fabs(a - e));
        } else {
            (x <= -20 ? cf : x <= 0 ? erfc_lower : upper).add(x, a, e);
        }
    }
    VCAL_CHECK_ULP_STATS(cf, tol::TOL_REF_LOG_NCDF_ULP, "ref log_ncdf x <= -20 (continued fraction)");
    VCAL_CHECK_ULP_STATS(erfc_lower, tol::TOL_REF_LOG_NCDF_ULP, "ref log_ncdf -20 < x <= 0 (std::erfc)");
    VCAL_CHECK_ULP_STATS(upper, tol::TOL_REF_LOG_NCDF_ULP, "ref log_ncdf x > 0 (erfc to 2, then continued fraction)");
    VCAL_CHECK(subnormal_worst <= DBL_MIN);
}

VCAL_TEST(ref_ncdf_inv_matches_mpmath) {
    const auto t = vcal::test::read_golden_csv("special/probit.csv");
    const auto cp = t.column("p_hex");
    const auto ce = t.column("expected_hex");
    UlpStats s;
    for (const auto& r : t.rows) {
        const double p = parse_double(r[cp]);
        s.add(p, vcalref::ncdf_inv(p), parse_double(r[ce]));
    }
    VCAL_CHECK_ULP_STATS(s, tol::TOL_REF_NCDF_INV_ULP, "ref ncdf_inv");
    VCAL_CHECK_EQ(vcalref::ncdf_inv(0.0), -kInf);
    VCAL_CHECK_EQ(vcalref::ncdf_inv(1.0), kInf);
    VCAL_CHECK(std::isnan(vcalref::ncdf_inv(-0.5)));
}

// Rows with min(k, n-k) up to 1e8 (the direct sum's domain); huge-k rows are core-only.
VCAL_TEST(ref_lchoose_matches_mpmath) {
    const auto t = vcal::test::read_golden_csv("special/lbinom.csv");
    UlpStats s;
    for (const auto& r : t.rows) {
        const std::int64_t n = parse_int(r[t.column("n")]);
        const std::int64_t k = parse_int(r[t.column("k")]);
        if (std::min(k, n - k) > 100000000) {
            VCAL_CHECK(std::isnan(vcalref::lchoose(n, k)));
            continue;
        }
        s.add(static_cast<double>(n), vcalref::lchoose(n, k), parse_double(r[t.column("expected_hex")]));
    }
    VCAL_CHECK_ULP_STATS(s, tol::TOL_REF_LCHOOSE_ULP, "ref lchoose");
}

// Every golden mixture case, including the zero-default high-rho periods where core's
// adaptive GH is not at full precision (Q15): ref's adaptive panels must be.
VCAL_TEST(ref_log_mixture_matches_mpmath_everywhere) {
    const auto t = vcal::test::read_golden_csv("quadrature/mixture.csv");
    double worst = 0.0;
    std::string where;
    for (const auto& r : t.rows) {
        const double pd = parse_double(r[t.column("pd_hex")]);
        const double rho = parse_double(r[t.column("rho_hex")]);
        const std::int64_t n = parse_int(r[t.column("n")]);
        const std::int64_t d = parse_int(r[t.column("d")]);
        const double e = parse_double(r[t.column("log_integral_hex")]);
        const double err = eps_units(vcalref::log_mixture(pd, rho, n, d), e);
        if (!(err <= worst)) {
            worst = err;
            where = r[t.column("kind")] + " pd=" + describe(pd) + " rho=" + describe(rho) + " n=" + std::to_string(n) +
                    " d=" + std::to_string(d);
        }
    }
    vcal::test::note("ref log_mixture: worst " + describe(worst) + " eps at " + where);
    VCAL_CHECK(worst <= tol::TOL_REF_LOG_MIXTURE_EPS);
}

// The fit is a maximum: no nearby point in logit coordinates does better, the Hessian gives
// finite SEs, and an interior estimate is not flagged as on the boundary.
VCAL_TEST(ref_fit_finds_an_interior_maximum) {
    const std::vector<vcalref::Period> panel{{500, 1}, {500, 3}, {500, 0}, {500, 2}, {500, 6}, {500, 1}, {500, 0}, {500, 4}};
    const auto f = vcalref::fit(panel, 1e-4, 0.2, 1e-3, 0.5);
    vcal::test::note("ref fit: PD " + describe(f.pd) + " (se " + describe(f.se_pd) + "), rho " + describe(f.rho) +
                     " (se " + describe(f.se_rho) + "), corr " + describe(f.corr) + ", loglik " + describe(f.loglik));
    VCAL_CHECK(!f.on_boundary);
    VCAL_CHECK(std::isfinite(f.se_pd) && std::isfinite(f.se_rho) && std::isfinite(f.corr));
    const auto logit = [](double v) { return std::log(v / (1.0 - v)); };
    const auto inv = [](double u) { return 1.0 / (1.0 + std::exp(-u)); };
    for (const double du : {-1e-3, 1e-3}) {
        VCAL_CHECK(vcalref::loglik(inv(logit(f.pd) + du), f.rho, panel) <= f.loglik);
        VCAL_CHECK(vcalref::loglik(f.pd, inv(logit(f.rho) + du), panel) <= f.loglik);
    }
}

VCAL_TEST(ref_fit_reports_boundary_estimates) {
    const std::vector<vcalref::Period> no_defaults(4, vcalref::Period{500, 0});
    const auto f = vcalref::fit(no_defaults, 1e-4, 0.2, 1e-3, 0.5);
    VCAL_CHECK(f.on_boundary);
}
