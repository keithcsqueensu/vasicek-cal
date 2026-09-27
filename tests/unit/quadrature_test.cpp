// SPDX-License-Identifier: Apache-2.0
// core/quadrature and the binomial-mixture integrand against mpmath (M1.3; D-035, D-037, D-080).
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "core/model/binomial_mixture.hpp"
#include "core/model/vasicek.hpp"
#include "core/quadrature/composite_legendre.hpp"
#include "core/quadrature/parity.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace q = vcal::quadrature;
namespace m = vcal::model;
namespace tol = vcal::tol;
using vcal::test::describe;
using vcal::test::parse_double;
using vcal::test::parse_int;
using Fixed = q::GaussHermiteFixed<vcal::PrecisionF64>;
using Adaptive = q::GaussHermiteAdaptive<vcal::PrecisionF64>;

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// |actual - expected| in units of DBL_EPSILON * max(1, |expected|): relative error of a
// log-scale quantity, and relative error of the integral itself when |expected| <= 1.
double eps_units(double actual, double expected) {
    return std::fabs(actual - expected) / (DBL_EPSILON * std::fmax(1.0, std::fabs(expected)));
}

struct MaxError {
    double value = 0.0;
    std::string where;
    void add(double v, const std::string& at) {
        if (!(v <= value)) {  // also catches NaN
            value = v;
            where = at;
        }
    }
};

struct MixtureCase {
    std::string kind;
    double pd;
    double rho;
    std::int64_t n;
    std::int64_t d;
    double expected;  // log E[p^d (1-p)^(n-d)]

    std::string label() const {
        return kind + " pd=" + describe(pd) + " rho=" + describe(rho) + " n=" + std::to_string(n) +
               " d=" + std::to_string(d);
    }
};

std::vector<MixtureCase> mixture_cases() {
    const auto t = vcal::test::read_golden_csv("quadrature/mixture.csv");
    std::vector<MixtureCase> out;
    for (const auto& r : t.rows) {
        out.push_back({r[t.column("kind")], parse_double(r[t.column("pd_hex")]),
                       parse_double(r[t.column("rho_hex")]), parse_int(r[t.column("n")]),
                       parse_int(r[t.column("d")]), parse_double(r[t.column("log_integral_hex")])});
    }
    return out;
}

template <class Integrator>
double log_mixture(const Integrator& integ, const MixtureCase& c) {
    const auto f = m::make_binomial_log_integrand(m::make_vasicek1f(c.pd, c.rho), c.n, c.d);
    return integ.log_integrate(f, m::binomial_mixture_hint(f));
}

}  // namespace

// --- the generated rules (D-080) ---------------------------------------------------------

VCAL_TEST(gauss_hermite_rule_lookup) {
    for (const int n : q::generated::kGaussHermiteSizes) {
        const auto r = q::gauss_hermite_rule(n);
        VCAL_CHECK_EQ(r.n, n);
        VCAL_CHECK(r.node != nullptr && r.log_weight != nullptr);
    }
    for (const int n : {0, 1, 7, 9, 100, 129, 255, 512}) {
        const auto r = q::gauss_hermite_rule(n);
        VCAL_CHECK_EQ(r.n, 0);
        VCAL_CHECK(r.node == nullptr && r.log_weight == nullptr);
    }
}

VCAL_TEST(gauss_hermite_tables_are_symmetric_and_sorted) {
    for (const int n : q::generated::kGaussHermiteSizes) {
        const auto r = q::gauss_hermite_rule(n);
        for (int i = 0; i < n; ++i) {
            VCAL_CHECK_EQ(r.node[i], -r.node[n - 1 - i]);
            VCAL_CHECK_EQ(r.log_weight[i], r.log_weight[n - 1 - i]);
            if (i > 0) VCAL_CHECK(r.node[i - 1] < r.node[i]);
        }
    }
}

// Every rule must integrate z^m exactly for m <= 2N-1 (E[Z^2k] = (2k-1)!!, odd moments 0).
// High moments are dominated by a few outer nodes and (2k-1)!! is huge, so the check is
// relative, with a tolerance proportional to the degree: each term z^m w carries about m
// rounding errors from the node. The sums run in log space (at N = 256 a degree-511 term is
// ~e^1260, far past double's range), scaled by the largest term and added smallest first.
VCAL_TEST(gauss_hermite_moments_are_exact_to_degree_2n_minus_1) {
    const auto t = vcal::test::read_golden_csv("quadrature/gh_moments.csv");
    std::vector<double> log_moment;
    for (const auto& r : t.rows) log_moment.push_back(parse_double(r[t.column("log_moment_hex")]));

    struct Term {
        double log_mag;
        bool negative;
    };
    MaxError even_per_degree, odd_rel;
    for (const int n : q::generated::kGaussHermiteSizes) {
        const auto r = q::gauss_hermite_rule(n);
        for (int k = 0; k < n; ++k) {
            for (const int deg : {2 * k, 2 * k + 1}) {
                std::vector<Term> terms;
                double log_max = -kInf;
                for (int i = 0; i < n; ++i) {
                    const double z = r.node[i];
                    const double log_mag = r.log_weight[i] + deg * std::log(std::fabs(z));
                    terms.push_back({log_mag, deg % 2 == 1 && z < 0});
                    log_max = std::fmax(log_max, log_mag);
                }
                std::sort(terms.begin(), terms.end(),
                          [](const Term& a, const Term& b) { return a.log_mag < b.log_mag; });
                double sum = 0.0;  // sum / e^log_max
                double abs_sum = 0.0;
                for (const Term& v : terms) {
                    const double scaled = std::exp(v.log_mag - log_max);
                    sum += v.negative ? -scaled : scaled;
                    abs_sum += scaled;
                }
                const std::string at = "N=" + std::to_string(n) + " degree " + std::to_string(deg);
                if (deg % 2 == 0) {
                    const double log_sum = log_max + std::log(sum);
                    const double rel = std::fabs(log_sum - log_moment[static_cast<std::size_t>(k)]);
                    even_per_degree.add(rel / (DBL_EPSILON * (deg + 1)), at);
                } else {
                    odd_rel.add(std::fabs(sum) / (DBL_EPSILON * abs_sum), at);
                }
            }
        }
    }
    vcal::test::note("even moments: max relative error " + describe(even_per_degree.value) +
                     " x eps x (degree+1), at " + even_per_degree.where);
    vcal::test::note("odd moments: max |sum| / sum|terms| " + describe(odd_rel.value) + " x eps, at " +
                     odd_rel.where);
    VCAL_CHECK(even_per_degree.value <= tol::TOL_GH_MOMENT_EPS_PER_DEGREE);
    VCAL_CHECK(odd_rel.value <= tol::TOL_GH_ODD_MOMENT_EPS);
}

// --- integrator mechanics ----------------------------------------------------------------

VCAL_TEST(online_log_sum_exp_edge_cases) {
    q::detail::OnlineLogSumExp all_neg_inf;
    all_neg_inf.add(-kInf);
    all_neg_inf.add(-kInf);
    VCAL_CHECK_EQ(all_neg_inf.result(), -kInf);

    q::detail::OnlineLogSumExp with_nan;
    with_nan.add(1.0);
    with_nan.add(kNaN);
    VCAL_CHECK(std::isnan(with_nan.result()));

    q::detail::OnlineLogSumExp with_inf;
    with_inf.add(1.0);
    with_inf.add(kInf);
    with_inf.add(kInf);
    VCAL_CHECK_EQ(with_inf.result(), kInf);

    q::detail::OnlineLogSumExp pair;  // log(e^0 + e^0) = log 2, independent of the -inf term
    pair.add(0.0);
    pair.add(-kInf);
    pair.add(0.0);
    VCAL_CHECK(vcal::test::ulp_distance(pair.result(), std::log(2.0)) <= 1);
}

// For g(z) = a + b z - c z^2 / 2, h = g - z^2/2 is quadratic, so the adaptive rule is exact at
// every N given the true mode b/(1+c) and scale 1/sqrt(1+c):
//     log I = a + b^2 / (2 (1 + c)) - log(1 + c) / 2.
VCAL_TEST(adaptive_rule_is_exact_for_a_gaussian_log_integrand) {
    MaxError err;
    for (const double b : {-40.0, -3.0, 0.0, 0.5, 12.0}) {
        for (const double c : {0.0, 0.5, 10.0, 1e4, 1e8}) {
            const double a = -7.25;
            const auto g = [=](double z) { return a + b * z - 0.5 * c * z * z; };
            const double expected = a + b * b / (2.0 * (1.0 + c)) - 0.5 * std::log1p(c);
            const double scale = 1.0 / std::sqrt(1.0 + c);
            const q::IntegrandHint hint{b / (1.0 + c), scale, b / (1.0 + c), scale, false};
            for (const int n : q::generated::kGaussHermiteSizes) {
                const double got = Adaptive{q::gauss_hermite_rule(n)}.log_integrate(g, hint);
                err.add(eps_units(got, expected),
                        "b=" + describe(b) + " c=" + describe(c) + " N=" + std::to_string(n));
            }
        }
    }
    vcal::test::note("gaussian log-integrand, adaptive: max " + describe(err.value) + " eps at " + err.where);
    VCAL_CHECK(err.value <= tol::TOL_AGH_GAUSSIAN_EPS);
}

VCAL_TEST(adaptive_rule_rejects_invalid_hints) {
    const Adaptive integ{q::gauss_hermite_rule(16)};
    const auto g = [](double) { return 0.0; };
    VCAL_CHECK(std::isnan(integ.log_integrate(g, {0.0, 0.0, 0.0, 0.0, false})));
    VCAL_CHECK(std::isnan(integ.log_integrate(g, {0.0, -1.0, 0.0, 0.0, false})));
    VCAL_CHECK(std::isnan(integ.log_integrate(g, {kNaN, 1.0, 0.0, 0.0, false})));
    VCAL_CHECK(std::isnan(integ.log_integrate(g, {0.0, kInf, 0.0, 0.0, false})));
    VCAL_CHECK_EQ(integ.log_integrate(g, {0.0, 1.0, 0.0, 0.0, false}), Fixed{q::gauss_hermite_rule(16)}.log_integrate(g, {}));
}

// --- binomial-mixture integrals against mpmath -------------------------------------------

// Convergence table: worst error per rule and N, over the closed forms (n <= 2: E[p] = PD,
// E[1-p], E[p^2] = Phi2, E[p(1-p)]) and the spikes (n up to 1e6).
VCAL_TEST(mixture_convergence_by_rule_and_n) {
    const auto cases = mixture_cases();
    for (const char* kind : {"closed", "spike"}) {
        for (const int n : q::generated::kGaussHermiteSizes) {
            MaxError fixed_err, adaptive_err;
            for (const auto& c : cases) {
                if (c.kind != kind) continue;
                const auto rule = q::gauss_hermite_rule(n);
                fixed_err.add(eps_units(log_mixture(Fixed{rule}, c), c.expected), c.label());
                adaptive_err.add(eps_units(log_mixture(Adaptive{rule}, c), c.expected), c.label());
            }
            vcal::test::note(std::string(kind) + " N=" + std::to_string(n) + ": fixed " +
                             describe(fixed_err.value) + " eps, adaptive " + describe(adaptive_err.value) +
                             " eps (worst: " + adaptive_err.where + ")");
        }
    }
    // Adaptive rule on each spike case.
    for (const auto& c : cases) {
        if (c.kind != "spike") continue;
        std::string line = c.label() + ":";
        for (const int n : q::generated::kGaussHermiteSizes) {
            const double e = eps_units(log_mixture(Adaptive{q::gauss_hermite_rule(n)}, c), c.expected);
            line += " N" + std::to_string(n) + "=" + describe(static_cast<float>(e));
        }
        vcal::test::note(line);
    }
    // Adaptive rule on the closed forms, broken down by rho and by n.
    for (const double rho : {1e-4, 0.01, 0.05, 0.12, 0.24, 0.5, 0.9, 0.95, 0.99}) {
        for (const std::int64_t nobs : {1, 2}) {
            std::string line = "closed rho=" + describe(rho) + " n=" + std::to_string(nobs) + ":";
            for (const int n : q::generated::kGaussHermiteSizes) {
                MaxError e;
                for (const auto& c : cases) {
                    if (c.kind == "closed" && c.rho == rho && c.n == nobs) {
                        e.add(eps_units(log_mixture(Adaptive{q::gauss_hermite_rule(n)}, c), c.expected), "");
                    }
                }
                line += " N" + std::to_string(n) + "=" + describe(static_cast<float>(e.value));
            }
            vcal::test::note(line);
        }
    }
}

// Worst error over the golden cases with rho <= max_rho and 0 < d < n, at every table size
// N >= min_n. Zero- and all-default periods (d in {0, n}) are reported separately, as absolute
// error in log L_t, and are not part of the precision claim (D-116).
static MaxError adaptive_regime_error(double max_rho, int min_n) {
    MaxError interior, monotone_abs;
    for (const auto& c : mixture_cases()) {
        if (c.rho > max_rho) continue;
        const bool monotone = c.d == 0 || c.d == c.n;
        for (const int n : q::generated::kGaussHermiteSizes) {
            if (n < min_n) continue;
            const double got = log_mixture(Adaptive{q::gauss_hermite_rule(n)}, c);
            const std::string at = c.label() + " N=" + std::to_string(n);
            if (monotone) {
                monotone_abs.add(std::fabs(got - c.expected), at);
            } else {
                interior.add(eps_units(got, c.expected), at);
            }
        }
    }
    vcal::test::note("adaptive, N >= " + std::to_string(min_n) + ", rho <= " + describe(max_rho) + ", 0 < d < n: max " +
                     describe(interior.value) + " eps at " + interior.where);
    vcal::test::note("  (not claimed) d in {0, n}: max abs error " + describe(monotone_abs.value) + " at " +
                     monotone_abs.where);
    return interior;
}

// Regime A: full precision for rho <= 0.24 at N >= 64, n from 2 to 1e6, periods with 0 < d < n.
VCAL_TEST(mixture_adaptive_full_precision_regime) {
    const auto err = adaptive_regime_error(tol::TOL_AGH_REGIME_A_MAX_RHO, tol::TOL_AGH_REGIME_A_MIN_N);
    VCAL_CHECK(err.value <= tol::TOL_AGH_REGIME_A_EPS);
}

// Regime B: rho <= 0.5 at N >= 128, periods with 0 < d < n. Zero- and all-default periods are
// not claimed at any rho (D-116): their integrand is the prior truncated by a steep survival
// factor, which polynomial rules resolve slowly, worst when few defaults are expected (tiny PD,
// large n). Parity integrates them with the composite rule instead (D-118).
VCAL_TEST(mixture_adaptive_rho_half_regime) {
    const auto err = adaptive_regime_error(tol::TOL_AGH_REGIME_B_MAX_RHO, tol::TOL_AGH_REGIME_B_MIN_N);
    VCAL_CHECK(err.value <= tol::TOL_AGH_REGIME_B_EPS);
}

// D-092, D-120: parity's per-run safeguard is the difference between the rule and its doubled
// variant. It is only a safeguard if it never misses: every golden period whose error in log L_t
// exceeds the flag threshold must also show a rule-vs-doubled difference above it.
//   - parity (D-118) vs doubled parity: reports the worst error; on the golden set nothing is
//     above the threshold, so this alone would not exercise detection.
//   - half parity (GH 64, K = 8 panels) vs parity: also accurate on the golden set.
//   - quarter parity (GH 32, K = 4) vs half parity as its doubled variant: the same family two
//     doublings coarser, inaccurate on many golden periods, so detection is actually tested.
VCAL_TEST(quadrature_check_detects_every_inaccurate_period) {
    const auto run = [](const auto& rule, const auto& doubled, const char* name, bool expect_flags) {
        MaxError worst_error, worst_check;
        double smallest_detection_ratio = kInf;
        std::string where;
        int flagged = 0;
        for (const auto& c : mixture_cases()) {
            const double l = log_mixture(rule, c);
            const double error = std::fabs(l - c.expected);
            const double check = std::fabs(l - log_mixture(doubled, c));
            worst_error.add(error, c.label());
            worst_check.add(check, c.label());
            if (error > tol::TOL_QUAD_CHECK_FLAG_ABS) {
                ++flagged;
                VCAL_CHECK(check > tol::TOL_QUAD_CHECK_FLAG_ABS);
                if (check / error < smallest_detection_ratio) {
                    smallest_detection_ratio = check / error;
                    where = c.label();
                }
            }
        }
        vcal::test::note(std::string(name) + ": worst |error in log L_t| " + describe(worst_error.value) + " at " +
                         worst_error.where + "; worst |rule - doubled| " + describe(worst_check.value) + "; " +
                         std::to_string(flagged) + " golden periods above the flag threshold" +
                         (flagged ? ", smallest check/error " + describe(smallest_detection_ratio) + " at " + where
                                  : std::string()));
        if (expect_flags) VCAL_CHECK(flagged > 0);
        VCAL_CHECK(flagged == 0 || smallest_detection_ratio >= tol::TOL_QUAD_CHECK_MIN_DETECTION_RATIO);
    };
    const auto parity = q::parity_rule();
    run(parity, q::parity_rule(true), "parity", false);
    using Split = q::SplitRule<vcal::PrecisionF64>;
    using Composite = q::CompositeLegendre<vcal::PrecisionF64>;
    const Split half{Adaptive{q::gauss_hermite_rule(64)}, Composite{q::gauss_legendre_rule(16), 8}};
    const Split quarter{Adaptive{q::gauss_hermite_rule(32)}, Composite{q::gauss_legendre_rule(16), 4}};
    run(half, parity, "half parity (GH 64, K = 8)", false);
    run(quarter, half, "quarter parity (GH 32, K = 4)", true);
}

// The parity rule against mpmath on every golden case: closed forms and spikes, interior and
// one-sided, rho up to 0.99 (D-118).
VCAL_TEST(parity_rule_matches_mpmath_on_every_golden_case) {
    const auto parity = q::parity_rule();
    MaxError interior, one_sided;
    for (const auto& c : mixture_cases()) {
        const double e = eps_units(log_mixture(parity, c), c.expected);
        (c.d == 0 || c.d == c.n ? one_sided : interior).add(e, c.label());
    }
    vcal::test::note("parity, 0 < d < n: worst " + describe(interior.value) + " eps at " + interior.where);
    vcal::test::note("parity, d in {0, n}: worst " + describe(one_sided.value) + " eps at " + one_sided.where);
    VCAL_CHECK(interior.value <= tol::TOL_PARITY_MIXTURE_EPS);
    VCAL_CHECK(one_sided.value <= tol::TOL_PARITY_MIXTURE_EPS);
}

// D-037: fixed-node GH is not a valid parity integrator at large n. Its error on an n = 1e6
// period is enormous even at N = 128 (documented, not toleranced; the bound is loose).
VCAL_TEST(fixed_rule_fails_on_large_n) {
    for (const auto& c : mixture_cases()) {
        if (c.n != 1000000 || c.d != 1000) continue;
        const double fixed = log_mixture(Fixed{q::gauss_hermite_rule(128)}, c);
        const double adaptive = log_mixture(Adaptive{q::gauss_hermite_rule(32)}, c);
        vcal::test::note(c.label() + ": expected " + describe(c.expected) + ", fixed N=128 " +
                         describe(fixed) + ", adaptive N=32 " + describe(adaptive));
        VCAL_CHECK(std::fabs(fixed - c.expected) > 1.0);
    }
}

// --- one-sided integrands: composite Gauss-Legendre (D-118) --------------------------------------

// Convergence by panel count K and points per panel M over the golden zero/all-default cases,
// for choosing the one fixed (K, M) that parity uses at every parameter point.
VCAL_TEST(composite_legendre_convergence_by_k_and_m) {
    using Composite = q::CompositeLegendre<vcal::PrecisionF64>;
    const auto cases = mixture_cases();
    for (const int m : {8, 12, 16, 20, 24, 32}) {
        std::string line = "M=" + std::to_string(m) + ":";
        for (const int k : {4, 8, 16, 24, 32}) {
            MaxError low, high;
            for (const auto& c : cases) {
                if (c.d != 0 && c.d != c.n) continue;
                const double got = log_mixture(Composite{q::gauss_legendre_rule(m), k}, c);
                (c.rho <= tol::TOL_AGH_REGIME_B_MAX_RHO ? low : high).add(eps_units(got, c.expected), c.label());
            }
            line += " K" + std::to_string(k) + "=" + describe(static_cast<float>(low.value)) + "/" +
                    describe(static_cast<float>(high.value));
        }
        vcal::test::note(line + "   (worst eps: rho <= 0.5 / rho > 0.5)");
    }
}
