// SPDX-License-Identifier: Apache-2.0
// core/special against mpmath golden tables (M1.2; D-033, D-067..D-071).
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

#include "core/special/erfcx.hpp"
#include "core/special/inverse_mills.hpp"
#include "core/special/lbinom.hpp"
#include "core/special/log_add_exp.hpp"
#include "core/special/log_phi.hpp"
#include "core/special/probit.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

using vcal::test::parse_double;
using vcal::test::parse_int;
using vcal::test::UlpStats;
namespace sp = vcal::special;
namespace tol = vcal::tol;

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

// --- erfcx --------------------------------------------------------------------------

VCAL_TEST(erfcx_matches_mpmath) {
    const auto t = vcal::test::read_golden_csv("special/erfcx.csv");
    const auto cx = t.column("x_hex");
    const auto ce = t.column("expected_hex");
    UlpStats negative, small, large;
    for (const auto& r : t.rows) {
        const double x = parse_double(r[cx]);
        const double e = parse_double(r[ce]);
        const double a = sp::erfcx(x);
        if (std::isinf(e)) {
            VCAL_CHECK_EQ(a, e);
            continue;
        }
        (x < 0 ? negative : x <= 8 ? small : large).add(x, a, e);
    }
    VCAL_CHECK_ULP_STATS(negative, tol::TOL_ERFCX_ULP, "erfcx x < 0");
    VCAL_CHECK_ULP_STATS(small, tol::TOL_ERFCX_ULP, "erfcx 0 <= x <= 8");
    VCAL_CHECK_ULP_STATS(large, tol::TOL_ERFCX_ULP, "erfcx x > 8");
}

VCAL_TEST(erfcx_special_values) {
    VCAL_CHECK(std::isnan(sp::erfcx(kNaN)));
    VCAL_CHECK_EQ(sp::erfcx(kInf), 0.0);
    VCAL_CHECK_EQ(sp::erfcx(-kInf), kInf);
    VCAL_CHECK_EQ(sp::erfcx(-0.0), sp::erfcx(0.0));
}

// --- log_phi ------------------------------------------------------------------------

VCAL_TEST(log_phi_matches_mpmath_relative_on_both_tails) {
    const auto t = vcal::test::read_golden_csv("special/log_phi.csv");
    const auto cx = t.column("x_hex");
    const auto ce = t.column("expected_hex");
    UlpStats lower, central, upper;
    double subnormal_worst = 0.0;
    for (const auto& r : t.rows) {
        const double x = parse_double(r[cx]);
        const double e = parse_double(r[ce]);
        const double a = sp::log_phi(x);
        if (std::isinf(e)) {
            VCAL_CHECK_EQ(a, e);
        } else if (std::fabs(e) < DBL_MIN) {  // log Phi(x) subnormal: absolute only (D-068)
            subnormal_worst = std::fmax(subnormal_worst, std::fabs(a - e));
        } else {
            (x < -5 ? lower : x > 5 ? upper : central).add(x, a, e);
        }
    }
    VCAL_CHECK_ULP_STATS(lower, tol::TOL_LOG_PHI_ULP, "log_phi x < -5 (relative)");
    VCAL_CHECK_ULP_STATS(central, tol::TOL_LOG_PHI_ULP, "log_phi -5 <= x <= 5 (relative)");
    VCAL_CHECK_ULP_STATS(upper, tol::TOL_LOG_PHI_ULP, "log_phi x > 5 (relative)");
    vcal::test::note("log_phi subnormal region: max abs error " + vcal::test::describe(subnormal_worst));
    VCAL_CHECK(subnormal_worst <= tol::TOL_LOG_PHI_SUBNORMAL_ABS);
}

VCAL_TEST(log_phi_special_values) {
    VCAL_CHECK(std::isnan(sp::log_phi(kNaN)));
    VCAL_CHECK_EQ(sp::log_phi(-kInf), -kInf);
    VCAL_CHECK_EQ(sp::log_phi(kInf), 0.0);
}

// --- inverse_mills ------------------------------------------------------------------

VCAL_TEST(inverse_mills_matches_mpmath) {
    const auto t = vcal::test::read_golden_csv("special/inverse_mills.csv");
    const auto cx = t.column("x_hex");
    const auto ce = t.column("expected_hex");
    UlpStats lower, upper;
    double underflow_worst = 0.0;
    for (const auto& r : t.rows) {
        const double x = parse_double(r[cx]);
        const double e = parse_double(r[ce]);
        const double a = sp::inverse_mills(x);
        if (std::fabs(e) < DBL_MIN) {  // phi(x) underflows for x above about 37.5
            underflow_worst = std::fmax(underflow_worst, std::fabs(a - e));
        } else {
            (x < 0 ? lower : upper).add(x, a, e);
        }
    }
    VCAL_CHECK_ULP_STATS(lower, tol::TOL_INVERSE_MILLS_ULP, "inverse_mills x < 0");
    VCAL_CHECK_ULP_STATS(upper, tol::TOL_INVERSE_MILLS_ULP, "inverse_mills x >= 0");
    VCAL_CHECK(underflow_worst <= DBL_MIN);
    VCAL_CHECK(std::isnan(sp::inverse_mills(kNaN)));
    VCAL_CHECK_EQ(sp::inverse_mills(-kInf), kInf);
    VCAL_CHECK_EQ(sp::inverse_mills(kInf), 0.0);
}

// --- probit -------------------------------------------------------------------------

VCAL_TEST(probit_matches_mpmath_near_0_and_1) {
    const auto t = vcal::test::read_golden_csv("special/probit.csv");
    const auto cp = t.column("p_hex");
    const auto ce = t.column("expected_hex");
    UlpStats near0, central, near1;
    for (const auto& r : t.rows) {
        const double p = parse_double(r[cp]);
        const double e = parse_double(r[ce]);
        (p < 1e-3 ? near0 : p > 1 - 1e-3 ? near1 : central).add(p, sp::probit(p), e);
    }
    VCAL_CHECK_ULP_STATS(near0, tol::TOL_PROBIT_ULP, "probit p < 1e-3");
    VCAL_CHECK_ULP_STATS(central, tol::TOL_PROBIT_ULP, "probit 1e-3 <= p <= 1-1e-3");
    VCAL_CHECK_ULP_STATS(near1, tol::TOL_PROBIT_ULP, "probit p > 1-1e-3");
}

// probit_upper(q) = Phi^-1(1 - q) must be accurate for q far below the spacing of doubles
// near 1, which no caller could reach by forming 1 - q (D-069).
VCAL_TEST(probit_upper_takes_the_complement) {
    const auto t = vcal::test::read_golden_csv("special/probit.csv");
    const auto cp = t.column("p_hex");
    const auto ce = t.column("expected_hex");
    UlpStats s;
    for (const auto& r : t.rows) {
        const double q = parse_double(r[cp]);
        if (q > 0.5) continue;
        const double a = sp::probit_upper(q);
        VCAL_CHECK_EQ(a, -sp::probit(q));
        s.add(q, a, -parse_double(r[ce]));  // Phi^-1(1 - q) = -Phi^-1(q)
    }
    VCAL_CHECK_ULP_STATS(s, tol::TOL_PROBIT_ULP, "probit_upper q <= 0.5");
}

VCAL_TEST(probit_special_values) {
    VCAL_CHECK_EQ(sp::probit(0.0), -kInf);
    VCAL_CHECK_EQ(sp::probit(1.0), kInf);
    VCAL_CHECK_EQ(sp::probit(0.5), 0.0);
    VCAL_CHECK(std::isnan(sp::probit(-0.1)));
    VCAL_CHECK(std::isnan(sp::probit(1.1)));
    VCAL_CHECK(std::isnan(sp::probit(kNaN)));
    VCAL_CHECK_EQ(sp::probit_upper(0.0), kInf);
}

// --- lbinom and the Stirling error ---------------------------------------------------

VCAL_TEST(lbinom_matches_mpmath) {
    const auto t = vcal::test::read_golden_csv("special/lbinom.csv");
    const auto cn = t.column("n");
    const auto ck = t.column("k");
    const auto ce = t.column("expected_hex");
    UlpStats exact_range, stirling_range;
    for (const auto& r : t.rows) {
        const std::int64_t n = parse_int(r[cn]);
        const std::int64_t k = parse_int(r[ck]);
        const double e = parse_double(r[ce]);
        const double a = sp::lbinom(n, k);
        if (k == 0 || k == n) {
            VCAL_CHECK_EQ(a, 0.0);  // exact, never via the general formula
            continue;
        }
        VCAL_CHECK_EQ(a, sp::lbinom(n, n - k));  // symmetric by construction
        (n <= sp::kLbinomExactMax ? exact_range : stirling_range).add(static_cast<double>(n), a, e);
    }
    VCAL_CHECK_ULP_STATS(exact_range, tol::TOL_LBINOM_ULP, "lbinom n <= 60 (exact integer)");
    VCAL_CHECK_ULP_STATS(stirling_range, tol::TOL_LBINOM_ULP, "lbinom n > 60 (Stirling differences)");
}

// Why D-070 does not use lgamma differences: at n = 1e6 they cancel to far worse
// accuracy than the Stirling form. Documents the choice; the bound is loose on purpose.
VCAL_TEST(lbinom_beats_lgamma_difference_at_large_n) {
    const std::int64_t n = 1000000;
    const std::int64_t k = 3;
    const auto t = vcal::test::read_golden_csv("special/lbinom.csv");
    double expected = kNaN;
    for (const auto& r : t.rows) {
        if (parse_int(r[t.column("n")]) == n && parse_int(r[t.column("k")]) == k) {
            expected = parse_double(r[t.column("expected_hex")]);
        }
    }
    VCAL_REQUIRE(!std::isnan(expected));
    const double via_lgamma = std::lgamma(1e6 + 1) - std::lgamma(4.0) - std::lgamma(1e6 - 2);
    const auto ours = vcal::test::ulp_distance(sp::lbinom(n, k), expected);
    const auto theirs = vcal::test::ulp_distance(via_lgamma, expected);
    vcal::test::note("lbinom(1e6, 3): Stirling form " + std::to_string(ours) + " ulp, lgamma difference " +
                     std::to_string(theirs) + " ulp");
    VCAL_CHECK(ours <= tol::TOL_LBINOM_ULP);
    VCAL_CHECK(theirs > 100 * tol::TOL_LBINOM_ULP);
}

VCAL_TEST(lbinom_invalid_arguments) {
    VCAL_CHECK(std::isnan(sp::lbinom(-1, 0)));
    VCAL_CHECK(std::isnan(sp::lbinom(5, -1)));
    VCAL_CHECK(std::isnan(sp::lbinom(5, 6)));
    VCAL_CHECK(std::isnan(sp::lbinom(sp::kLbinomMaxN + 1, 1)));
    VCAL_CHECK_EQ(sp::lbinom(0, 0), 0.0);
}

// Checks the generated small-n table against mpmath as well as the series (D-070).
VCAL_TEST(stirling_error_matches_mpmath) {
    const auto t = vcal::test::read_golden_csv("special/stirling_error.csv");
    const auto cn = t.column("n");
    const auto ce = t.column("expected_hex");
    UlpStats table, series;
    for (const auto& r : t.rows) {
        const std::int64_t n = parse_int(r[cn]);
        const double a = sp::stirling_error(n);
        (n <= vcal::special::generated::kStirlingTableMax ? table : series)
            .add(static_cast<double>(n), a, parse_double(r[ce]));
    }
    VCAL_CHECK_ULP_STATS(table, tol::TOL_STIRLING_ERROR_ULP, "stirling_error n <= 30 (table)");
    VCAL_CHECK_ULP_STATS(series, tol::TOL_STIRLING_ERROR_ULP, "stirling_error n > 30 (series)");
    VCAL_CHECK(std::isnan(sp::stirling_error(0)));
}

// --- log_add_exp (D-071) ---------------------------------------------------------------

VCAL_TEST(log_add_exp_infinities) {
    VCAL_CHECK_EQ(sp::log_add_exp(-kInf, -kInf), -kInf);
    VCAL_CHECK_EQ(sp::log_add_exp(-kInf, 3.0), 3.0);
    VCAL_CHECK_EQ(sp::log_add_exp(3.0, -kInf), 3.0);
    VCAL_CHECK_EQ(sp::log_add_exp(kInf, 3.0), kInf);
    VCAL_CHECK_EQ(sp::log_add_exp(kInf, -kInf), kInf);
    VCAL_CHECK(std::isnan(sp::log_add_exp(kNaN, 0.0)));
    VCAL_CHECK(std::isnan(sp::log_add_exp(0.0, kNaN)));
}

// The table's grid includes every equal pair (a (+) a = a + log 2) and log 2 (+) log 3.
VCAL_TEST(log_add_exp_matches_mpmath) {
    const auto t = vcal::test::read_golden_csv("special/log_add_exp.csv");
    const auto ca = t.column("a_hex");
    const auto cb = t.column("b_hex");
    const auto ce = t.column("expected_hex");
    UlpStats equal, distinct;
    for (const auto& r : t.rows) {
        const double a = parse_double(r[ca]);
        const double b = parse_double(r[cb]);
        const double got = sp::log_add_exp(a, b);
        VCAL_CHECK_EQ(got, sp::log_add_exp(b, a));  // symmetric bit for bit
        (a == b ? equal : distinct).add(a, got, parse_double(r[ce]));
    }
    VCAL_CHECK_ULP_STATS(equal, tol::TOL_LOG_ADD_EXP_ULP, "log_add_exp a == b");
    VCAL_CHECK_ULP_STATS(distinct, tol::TOL_LOG_ADD_EXP_ULP, "log_add_exp a != b");
    VCAL_CHECK_EQ(sp::log_add_exp(0.0, -800.0), 0.0);  // exp(-800) is below half an ulp of 1
}
