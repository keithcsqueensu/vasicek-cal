// SPDX-License-Identifier: Apache-2.0
// Pluto and Tasche's most-prudent upper bounds (S-6; D-214, D-215; benchmark_models/pluto_tasche/bounds.hpp):
// the zero-default closed form; the one-factor bound tending to the independent one as rho_a -> 0; the solve's
// bracket; monotonicity in D, gamma and rho_a (P6's ordering); and the committed goldens, which
// validation/scipy/pluto_tasche.py replicates with plain scipy.
#include "benchmark_models/pluto_tasche/bounds.hpp"

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "core/quadrature/parity.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace pt = vcal::benchmark_models::pluto_tasche;
namespace tol = vcal::tol;
using vcal::test::describe;

const double kGammas[] = {0.5, 0.75, 0.9, 0.95, 0.99};

double logit(double p) { return std::log(p) - std::log1p(-p); }

// The solve's invariant: the bracket is at most the tolerance wide and holds the root (g(lo) >= 0 > g(hi)).
template <class LogCdf>
void check_bracket(const pt::Bound& b, double gamma, const LogCdf& log_cdf) {
    VCAL_CHECK(b.flags == 0);
    VCAL_CHECK(b.u_hi - b.u_lo <= tol::TOL_PT_SOLVE_LOGIT);
    VCAL_CHECK(b.u_lo <= b.u && b.u <= b.u_hi);
    VCAL_CHECK(log_cdf(b.u_lo) >= std::log1p(-gamma));
    VCAL_CHECK(log_cdf(b.u_hi) < std::log1p(-gamma));
}

}  // namespace

// (a) at D = 0 is 1 - (1 - gamma)^(1/N): within the solve's tolerance in logit, and the bracket holds the root.
VCAL_TEST(pt_independent_zero_defaults_is_the_closed_form) {
    double worst = 0.0;
    for (const std::int64_t n : {100, 1000, 10000}) {
        for (const std::int64_t T : {1, 20, 40, 100}) {
            for (const double g : kGammas) {
                const pt::Bound b = pt::independent_bound(n, T, 0, g);
                const double exact = -std::expm1(std::log1p(-g) / static_cast<double>(n * T));
                const double diff = std::fabs(b.u - logit(exact));
                worst = std::fmax(worst, diff);
                VCAL_CHECK(diff <= tol::TOL_PT_SOLVE_LOGIT);
                check_bracket(b, g, [&](double u) { return pt::log_cdf_independent(n * T, 0, u); });
            }
        }
    }
    vcal::test::note("worst |logit bound - logit closed form|: " + describe(worst));
}

// (b) -> (a) as rho_a -> 0: the gap in logit falls with rho_a, and at rho_a = 1e-8 it is within
// TOL_PT_RHO_LIMIT_LOGIT. Each (b) bound's bracket holds its root.
VCAL_TEST(pt_factor_tends_to_independent_as_rho_falls) {
    const auto rule = vcal::quadrature::parity_rule();
    pt::FactorWork w;
    double worst = 0.0;
    struct Case {
        std::int64_t n, T, D;
    };
    for (const Case c : {Case{100, 20, 0}, Case{100, 20, 3}, Case{1000, 40, 12}, Case{100, 100, 60}, Case{1000, 100, 150}}) {
        for (const double g : {0.5, 0.95, 0.99}) {
            const pt::Bound a = pt::independent_bound(c.n, c.T, c.D, g);
            double prev = 0.0;
            bool first = true;
            for (const double rho : {1e-2, 1e-4, 1e-6, 1e-8}) {
                const pt::Bound b = pt::factor_bound(c.n, c.T, c.D, rho, g, rule, w);
                const double gap = b.u - a.u;
                VCAL_CHECK(first || std::fabs(gap) < std::fabs(prev));
                prev = gap;
                first = false;
                if (rho == 1e-8) {
                    worst = std::fmax(worst, std::fabs(gap));
                    VCAL_CHECK(std::fabs(gap) <= tol::TOL_PT_RHO_LIMIT_LOGIT);
                    check_bracket(b, g, [&](double u) { return pt::log_cdf_factor(c.n, c.T, c.D, u, rho, rule, w); });
                }
            }
        }
    }
    vcal::test::note("worst |logit (b) - logit (a)| at rho_a = 1e-8: " + describe(worst));
}

// The bounds rise strictly with D and with gamma, and on every total and gamma follow P6's order:
// independent < (b, 0.02) < (b, 0.12) < (b, 0.24). Exact comparisons: the reference's closest approach is 4.3e-4
// in logit, far above the solve's tolerance.
VCAL_TEST(pt_bounds_are_monotone_in_d_gamma_and_rho) {
    const auto rule = vcal::quadrature::parity_rule();
    pt::FactorWork w;
    const double rhos[] = {0.02, 0.12, 0.24};
    struct Cell {
        std::int64_t n, T;
        std::vector<std::int64_t> totals;
    };
    const Cell cells[] = {{100, 20, {0, 1, 2, 5, 10}}, {1000, 40, {0, 10, 40, 80}}, {100, 100, {0, 50, 100, 197}}};
    std::int64_t compared = 0;
    double closest = 1e300;
    for (const Cell& c : cells) {
        std::vector<std::vector<double>> prev_d(4, std::vector<double>(5, -1e300));  // [law][gamma] at the previous D
        for (const std::int64_t D : c.totals) {
            for (std::size_t gi = 0; gi < 5; ++gi) {
                double u[4];
                u[0] = pt::independent_bound(c.n, c.T, D, kGammas[gi]).u;
                for (int k = 0; k < 3; ++k) u[k + 1] = pt::factor_bound(c.n, c.T, D, rhos[k], kGammas[gi], rule, w).u;
                for (int k = 0; k < 4; ++k) {
                    VCAL_CHECK(u[k] > prev_d[static_cast<std::size_t>(k)][gi]);  // rises with D
                    if (gi > 0) {
                        // rises with gamma: compare with the previous gamma's bound at this D
                        const double lower = k == 0 ? pt::independent_bound(c.n, c.T, D, kGammas[gi - 1]).u
                                                    : pt::factor_bound(c.n, c.T, D, rhos[k - 1], kGammas[gi - 1], rule, w).u;
                        VCAL_CHECK(u[k] > lower);
                    }
                    prev_d[static_cast<std::size_t>(k)][gi] = u[k];
                }
                for (int k = 0; k < 3; ++k) {
                    VCAL_CHECK(u[k] < u[k + 1]);  // P6
                    closest = std::fmin(closest, u[k + 1] - u[k]);
                    ++compared;
                }
            }
        }
    }
    vcal::test::note(std::to_string(compared) + " adjacent pairs ordered; closest " + describe(closest) + " in logit");
}

// Invalid input gives NaN, flagged; a total no p can reach in the bracket is reported at the bracket's end, flagged.
VCAL_TEST(pt_bounds_flag_invalid_input) {
    const auto rule = vcal::quadrature::parity_rule();
    pt::FactorWork w;
    for (const pt::Bound b : {pt::independent_bound(0, 20, 0, 0.95), pt::independent_bound(100, 20, -1, 0.95),
                              pt::independent_bound(100, 20, 2001, 0.95), pt::independent_bound(100, 20, 0, 1.0),
                              pt::factor_bound(100, 20, 0, 0.0, 0.95, rule, w),
                              pt::factor_bound(100, 20, 0, 0.12, 0.0, rule, w)}) {
        VCAL_CHECK(b.flags == pt::kBoundInvalidInput);
        VCAL_CHECK(std::isnan(b.pd));
    }
    // D = N: P(S <= N) = 1 for every p, so g >= 0 everywhere and the bound is the bracket's top.
    const pt::Bound top = pt::independent_bound(100, 2, 200, 0.95);
    VCAL_CHECK(top.flags == pt::kBoundAboveBracket);
    VCAL_CHECK(top.u == pt::kLogitHi);
}

// The committed goldens (tests/golden/pluto_tasche/bounds.csv, tests/studies/pluto_tasche_goldens.cpp): the engine
// reproduces them to TOL_PT_GOLDEN_LOGIT on every platform. validation/scipy/pluto_tasche.py replicates the same
// rows with plain scipy.
VCAL_TEST(pt_goldens_reproduce) {
    const auto t = vcal::test::read_golden_csv("pluto_tasche/bounds.csv");
    const auto rule = vcal::quadrature::parity_rule();
    pt::FactorWork w;
    double worst = 0.0;
    for (const auto& r : t.rows) {
        const std::int64_t n = std::stoll(r[t.column("n")]), T = std::stoll(r[t.column("periods")]),
                           D = std::stoll(r[t.column("d_total")]);
        const double rho = vcal::test::parse_double(r[t.column("rho_a_hex")]);
        const double g = vcal::test::parse_double(r[t.column("gamma_hex")]);
        const double want = vcal::test::parse_double(r[t.column("logit_bound_hex")]);
        const pt::Bound b = r[t.column("method")] == "independent" ? pt::independent_bound(n, T, D, g)
                                                                   : pt::factor_bound(n, T, D, rho, g, rule, w);
        const double diff = std::fabs(b.u - want);
        worst = std::fmax(worst, diff);
        VCAL_CHECK(b.flags == 0);
        VCAL_CHECK(diff <= tol::TOL_PT_GOLDEN_LOGIT);
    }
    vcal::test::note(std::to_string(t.rows.size()) + " rows; worst |logit engine - logit golden| " + describe(worst));
}
