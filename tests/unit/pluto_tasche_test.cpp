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

// The published examples: Pluto and Tasche (2005), arXiv:cond-mat/0411699v3, Tables 1-8 and the values in the
// text of Section 3 and footnote 4. Grades A, B and C have 100, 400 and 300 obligors, one period, and the most
// prudent bound for a grade pools it with every worse grade: N = 800, 700 and 300. Defaults: none (Tables 1-3, 7),
// or 0, 2 and 1 in A, B and C (Tables 4-6, 8: k = 3, 3, 1). Tables 7 and 8 are the one-factor bound at rho = 0.12.
// The tables print percentages to two decimals. Each value is reproduced within TOL_PT_PAPER_PCT percentage points:
// the rounding (0.005), plus the paper's own numerics (Table 8 matches its values rounded up; its 50% value for
// grade A, 0.72%, is printed 0.71% in Table 12). One printed value is not the bound it describes: Table 4's 0.65%
// at gamma = 75% (repeated in Table 11), where the exact Clopper-Pearson bound is 0.638% (scipy beta.ppf gives the
// same; at 0.65%, P(X <= 3) = 0.237 < 0.25, so the paper's value does not satisfy its own inequality (3.2)). That
// cell is checked against the exact value instead, and listed in docs/methodology/pluto_tasche.md.
VCAL_TEST(pt_reproduces_the_published_tables) {
    const auto rule = vcal::quadrature::parity_rule();
    pt::FactorWork w;
    const double g[6] = {0.5, 0.75, 0.9, 0.95, 0.99, 0.999};
    struct Row {
        const char* name;
        std::int64_t N, k;
        double rho;
        double pct[6];
    };
    const Row rows[] = {
        {"Table 1 A", 800, 0, 0.0, {0.09, 0.17, 0.29, 0.37, 0.57, 0.86}},
        {"Table 2 B", 700, 0, 0.0, {0.10, 0.20, 0.33, 0.43, 0.66, 0.98}},
        {"Table 3 C", 300, 0, 0.0, {0.23, 0.46, 0.76, 0.99, 1.52, 2.28}},
        {"Table 4 A", 800, 3, 0.0, {0.46, 0.65, 0.83, 0.97, 1.25, 1.62}},
        {"Table 5 B", 700, 3, 0.0, {0.52, 0.73, 0.95, 1.10, 1.43, 1.85}},
        {"Table 6 C", 300, 1, 0.0, {0.56, 0.90, 1.29, 1.57, 2.19, 3.04}},
        {"Table 7 A", 800, 0, 0.12, {0.15, 0.40, 0.86, 1.31, 2.65, 5.29}},
        {"Table 7 B", 700, 0, 0.12, {0.17, 0.45, 0.96, 1.45, 2.92, 5.77}},
        {"Table 7 C", 300, 0, 0.12, {0.37, 0.92, 1.89, 2.78, 5.30, 9.84}},
        {"Table 8 A", 800, 3, 0.12, {0.72, 1.42, 2.50, 3.42, 5.88, 10.08}},
        {"Table 8 B", 700, 3, 0.12, {0.81, 1.59, 2.77, 3.77, 6.43, 10.92}},
        {"Table 8 C", 300, 1, 0.12, {0.84, 1.76, 3.19, 4.41, 7.68, 13.14}},
    };
    double worst = 0.0;
    for (const Row& r : rows) {
        for (int i = 0; i < 6; ++i) {
            const pt::Bound b = r.rho == 0.0 ? pt::independent_bound(r.N, 1, r.k, g[i])
                                             : pt::factor_bound(r.N, 1, r.k, r.rho, g[i], rule, w);
            VCAL_CHECK(b.flags == 0);
            const double got = 100.0 * b.pd;
            if (std::string(r.name) == "Table 4 A" && i == 1) {
                // The erratum: the exact bound, 0.6378%, from the beta quantile (scipy beta.ppf(0.75, 4, 797)).
                VCAL_CHECK(std::fabs(got - 0.63783660637915540) <= 1e-9);
                continue;
            }
            worst = std::fmax(worst, std::fabs(got - r.pct[i]));
            VCAL_CHECK(std::fabs(got - r.pct[i]) <= tol::TOL_PT_PAPER_PCT);
        }
    }
    // Section 3's text: grade A alone (N = 100, no default) at 75% is 1.38%; footnote 4: grade B alone (N = 400,
    // two defaults) at 99.9% is 2.78%.
    VCAL_CHECK(std::fabs(100.0 * pt::independent_bound(100, 1, 0, 0.75).pd - 1.38) <= tol::TOL_PT_PAPER_PCT);
    VCAL_CHECK(std::fabs(100.0 * pt::independent_bound(400, 1, 2, 0.999).pd - 2.78) <= tol::TOL_PT_PAPER_PCT);
    vcal::test::note("71 table values (one erratum checked exactly) and 2 in the text; worst |engine - paper| " + describe(worst) + " percentage points");
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
