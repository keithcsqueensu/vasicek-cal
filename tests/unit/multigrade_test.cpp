// SPDX-License-Identifier: Apache-2.0
// The multi-grade binomial-mixture period log-likelihood (M6, D-151, D-411; core/objectives/multigrade.hpp):
//   - K = 1 is the single-grade likelihood bit for bit, value and derivatives, on every mixture golden case and with
//     empty grades inserted around it;
//   - two grades with equal PDs pool exactly (to rounding);
//   - rho -> 0 factorises into the grades' binomials;
//   - the derivatives against central differences of the value;
//   - the parity rule against mpmath on the K-grade goldens (tests/golden/multigrade/likelihood.csv), and against
//     ref/'s independent K-grade integral.
#include "core/objectives/multigrade.hpp"

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "core/special/lbinom.hpp"
#include "ref/ref.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace o = vcal::objectives;
namespace q = vcal::quadrature;
using vcal::test::parse_double;
using vcal::test::parse_int;
using Single = o::BinomialMixture<vcal::PrecisionF64>;
using o::GradeCount;

bool same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0 || (std::isnan(a) && std::isnan(b)); }

struct SingleCase {
    double pd, rho;
    std::int64_t n, d;
};

std::vector<SingleCase> single_cases() {
    const auto t = vcal::test::read_golden_csv("quadrature/mixture.csv");
    std::vector<SingleCase> out;
    for (const auto& r : t.rows) {
        out.push_back({parse_double(r[t.column("pd_hex")]), parse_double(r[t.column("rho_hex")]),
                       parse_int(r[t.column("n")]), parse_int(r[t.column("d")])});
    }
    // and a sweep across the interior and both one-sided shapes
    for (const double pd : {3e-4, 0.004, 0.03, 0.2}) {
        for (const double rho : {1e-3, 0.02, 0.12, 0.24, 0.5}) {
            for (const std::int64_t n : {1, 7, 100, 10000}) {
                for (const std::int64_t d : {std::int64_t{0}, std::int64_t{1}, n / 3, n}) out.push_back({pd, rho, n, d});
            }
        }
    }
    return out;
}

}  // namespace

// K = 1, alone and with empty grades around it, is BinomialMixture's value and derivatives bit for bit, under the
// parity rule and its doubled check.
VCAL_TEST(one_grade_is_the_single_grade_likelihood_bit_for_bit) {
    const auto parity = q::parity_rule();
    const auto check = q::parity_rule(true);
    std::size_t cases = 0, bad_value = 0, bad_derivs = 0;
    for (const auto& c : single_cases()) {
        const Single::Obs y{c.n, c.d};
        const Single::Theta th{c.pd, c.rho};
        for (const auto* integ : {&parity, &check}) {
            const double want = Single{}.log_contrib(y, th, *integ);
            double g1[2], h1[3];
            const double want_d = Single{}.log_contrib_derivs(y, th, *integ, g1, h1);
            // alone, then between empty grades
            const GradeCount alone[1] = {{c.n, c.d}};
            const double pd1[1] = {c.pd};
            const GradeCount padded[3] = {{0, 0}, {c.n, c.d}, {0, 0}};
            const double pd3[3] = {0.5, c.pd, 0.5};
            for (int variant = 0; variant < 2; ++variant) {
                const GradeCount* gy = variant == 0 ? alone : padded;
                const double* gp = variant == 0 ? pd1 : pd3;
                const int K = variant == 0 ? 1 : 3;
                const int at = variant == 0 ? 0 : 1;
                const double got = o::multigrade_log_contrib(gy, K, gp, c.rho, *integ);
                double grad[4], hess[16];
                const double got_d = o::multigrade_log_contrib_derivs(gy, K, gp, c.rho, *integ, grad, hess);
                const int n = K + 1;
                bad_value += same(got, want) && same(got_d, want_d) ? 0 : 1;
                bad_derivs += (same(grad[at], g1[0]) && same(grad[K], g1[1]) && same(hess[at * n + at], h1[0]) &&
                               same(hess[K * n + K], h1[1]) && same(hess[at * n + K], h1[2]) &&
                               same(hess[K * n + at], h1[2]))
                                  ? 0
                                  : 1;
                ++cases;
            }
        }
    }
    vcal::test::note(std::to_string(cases) + " cases; " + std::to_string(bad_value) + " values and " +
                     std::to_string(bad_derivs) + " derivative sets differ");
    VCAL_CHECK(bad_value == 0);
    VCAL_CHECK(bad_derivs == 0);
}

namespace {

struct GoldenCase {
    std::string kind;
    double rho;
    std::vector<GradeCount> grades;
    std::vector<double> pd;
    double expected;  // log E[prod_k p_k^d_k (1 - p_k)^(n_k - d_k)], without the binomial coefficients
    std::string label() const {
        std::string s = kind + " rho=" + std::to_string(rho) + " K=" + std::to_string(grades.size()) + " [";
        for (std::size_t k = 0; k < grades.size() && k < 4; ++k) {
            s += (k ? " " : "") + std::to_string(grades[k].n) + ":" + std::to_string(grades[k].d);
        }
        return s + (grades.size() > 4 ? " ...]" : "]");
    }
};

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::size_t a = 0;
    for (std::size_t b; (b = s.find(sep, a)) != std::string::npos; a = b + 1) out.push_back(s.substr(a, b - a));
    out.push_back(s.substr(a));
    return out;
}

std::vector<GoldenCase> golden_cases() {
    const auto t = vcal::test::read_golden_csv("multigrade/likelihood.csv");
    std::vector<GoldenCase> out;
    for (const auto& r : t.rows) {
        GoldenCase c;
        c.kind = r[t.column("kind")];
        c.rho = parse_double(r[t.column("rho_hex")]);
        for (const auto& g : split(r[t.column("grades")], ';')) {
            const auto nd = split(g, ':');
            c.grades.push_back({parse_int(nd[0]), parse_int(nd[1])});
        }
        for (const auto& p : split(r[t.column("pd_hex")], ';')) c.pd.push_back(parse_double(p));
        c.expected = parse_double(r[t.column("log_integral_hex")]);
        out.push_back(c);
    }
    return out;
}

double eps_units(double actual, double expected) {
    return std::fabs(actual - expected) / (DBL_EPSILON * std::fmax(1.0, std::fabs(expected)));
}

template <class I>
double core_log_integral(const GoldenCase& c, const I& integ) {
    const auto f = vcal::model::make_multigrade_log_integrand(c.grades.data(), c.pd.data(),
                                                              static_cast<int>(c.grades.size()), c.rho);
    return integ.log_integrate(f, vcal::model::multigrade_mixture_hint(f));
}

}  // namespace

// The parity rule against mpmath on every K-grade golden case (D-411 sec. 2), interior and one-sided reported
// separately, within TOL_PARITY_MIXTURE_EPS: the single-grade tolerance, unchanged.
VCAL_TEST(parity_rule_matches_mpmath_on_every_multigrade_golden_case) {
    const auto parity = q::parity_rule();
    double worst[2] = {0.0, 0.0};
    std::string where[2];
    int count[2] = {0, 0};
    for (const auto& c : golden_cases()) {
        const auto f = vcal::model::make_multigrade_log_integrand(c.grades.data(), c.pd.data(),
                                                                  static_cast<int>(c.grades.size()), c.rho);
        const int side = vcal::model::multigrade_is_one_sided(f) ? 1 : 0;
        const double e = eps_units(core_log_integral(c, parity), c.expected);
        ++count[side];
        if (!(e <= worst[side])) {
            worst[side] = e;
            where[side] = c.label();
        }
    }
    vcal::test::note("interior (" + std::to_string(count[0]) + " cases): worst " + std::to_string(worst[0]) +
                     " eps at " + where[0]);
    vcal::test::note("one-sided (" + std::to_string(count[1]) + " cases): worst " + std::to_string(worst[1]) +
                     " eps at " + where[1]);
    VCAL_CHECK(worst[0] <= vcal::tol::TOL_PARITY_MIXTURE_EPS);
    VCAL_CHECK(worst[1] <= vcal::tol::TOL_PARITY_MIXTURE_EPS);
}

namespace {

// Natural-scale derivatives mapped to logit coordinates u_i = logit(theta_i) (theta = PD_1..PD_K, rho):
// l_u_i = l_i v_i', l_u_i u_j = l_ij v_i' v_j' + [i == j] l_i v_i'', with v' = v(1 - v), v'' = v'(1 - 2v).
struct Logit {
    double l;
    std::vector<double> g, h;  // g[K + 1], h[(K + 1)^2]
};

Logit logit_derivs(const std::vector<GradeCount>& y, const std::vector<double>& u, const q::SplitRule<vcal::PrecisionF64>& rule) {
    const int K = static_cast<int>(y.size());
    const int n = K + 1;
    std::vector<double> v(n), v1(n), v2(n);
    for (int i = 0; i < n; ++i) {
        v[i] = 1.0 / (1.0 + std::exp(-u[i]));
        v1[i] = v[i] * (1.0 - v[i]);
        v2[i] = v1[i] * (1.0 - 2.0 * v[i]);
    }
    std::vector<double> g(n), h(n * n);
    Logit out{o::multigrade_log_contrib_derivs(y.data(), K, v.data(), v[K], rule, g.data(), h.data()), g, h};
    for (int i = 0; i < n; ++i) {
        out.g[i] = g[i] * v1[i];
        for (int j = 0; j < n; ++j) out.h[i * n + j] = h[i * n + j] * v1[i] * v1[j] + (i == j ? g[i] * v2[i] : 0.0);
    }
    return out;
}

double logit(double v) { return std::log(v) - std::log1p(-v); }

}  // namespace

// The analytic score and Hessian in all K + 1 coordinates against Richardson-extrapolated central differences in
// logit coordinates (as the single-grade test, TOL_SCORE_HESSIAN_FD_REL), on every K-grade golden case: each axis's
// step 0.01 / sqrt(max(1, |l_aa|)), relative to max(1, |analytic|); periods with a grade of n > 1e5 in the large-n
// tier. The step is half the single-grade test's 0.02: at 0.02 the 17-grade master-scale periods' rho-rho entry is
// limited by the differences' own O(h^4) truncation (2.1e-8; 3.4e-7 at 0.04, 1.3e-9 at 0.01, falling 16-fold per
// halving), not by the analytic formula. The value from the derivatives' fold equals multigrade_log_contrib's bit
// for bit.
VCAL_TEST(multigrade_score_and_hessian_match_central_differences) {
    const auto rule = q::parity_rule();
    double worst_g[2] = {0.0, 0.0}, worst_h[2] = {0.0, 0.0};
    std::string where_g[2], where_h[2];
    std::size_t value_bits = 0;
    for (const auto& c : golden_cases()) {
        const int K = static_cast<int>(c.grades.size());
        const int n = K + 1;
        std::vector<double> u(n);
        for (int k = 0; k < K; ++k) u[k] = logit(c.pd[k]);
        u[K] = logit(c.rho);
        const Logit a = logit_derivs(c.grades, u, rule);
        const double plain = o::multigrade_log_contrib(c.grades.data(), K, c.pd.data(), c.rho, rule);
        // the value via the logit round trip is at a slightly different point; compare at the exact inputs
        std::vector<double> gg(n), hh(n * n);
        const double via = o::multigrade_log_contrib_derivs(c.grades.data(), K, c.pd.data(), c.rho, rule, gg.data(), hh.data());
        value_bits += same(plain, via) ? 0 : 1;
        if (!std::isfinite(a.l)) continue;
        std::int64_t nmax = 0;
        for (const auto& gr : c.grades) nmax = std::max(nmax, gr.n);
        const int tier = nmax <= 100000 ? 0 : 1;
        std::vector<double> step(n);
        for (int i = 0; i < n; ++i) step[i] = 0.01 / std::sqrt(std::fmax(1.0, std::fabs(a.h[i * n + i])));
        auto diff = [&](double scale, std::vector<double>& dg, std::vector<double>& dh) {
            dg.assign(n, 0.0);
            dh.assign(n * n, 0.0);
            for (int j = 0; j < n; ++j) {
                std::vector<double> up = u, um = u;
                up[j] += scale * step[j];
                um[j] -= scale * step[j];
                const Logit p = logit_derivs(c.grades, up, rule), m = logit_derivs(c.grades, um, rule);
                dg[j] = (p.l - m.l) / (2 * scale * step[j]);
                for (int i = 0; i < n; ++i) dh[i * n + j] = (p.g[i] - m.g[i]) / (2 * scale * step[j]);
            }
        };
        std::vector<double> g1, h1, g2, h2;
        diff(1.0, g1, h1);
        diff(0.5, g2, h2);
        for (int i = 0; i < n; ++i) {
            const double fd = (4.0 * g2[i] - g1[i]) / 3.0;
            const double e = std::fabs(a.g[i] - fd) / std::fmax(1.0, std::fabs(a.g[i]));
            if (!(e <= worst_g[tier])) {
                worst_g[tier] = e;
                where_g[tier] = c.label() + " g" + std::to_string(i);
            }
            for (int j = 0; j < n; ++j) {
                const double fdh = (4.0 * h2[i * n + j] - h1[i * n + j]) / 3.0;
                const double eh = std::fabs(a.h[i * n + j] - fdh) / std::fmax(1.0, std::fabs(a.h[i * n + j]));
                if (!(eh <= worst_h[tier])) {
                    worst_h[tier] = eh;
                    where_h[tier] = c.label() + " h" + std::to_string(i) + "," + std::to_string(j);
                }
            }
        }
    }
    for (int tier = 0; tier < 2; ++tier) {
        vcal::test::note(std::string(tier == 0 ? "n <= 1e5" : "n > 1e5") + ": worst score " +
                         vcal::test::describe(worst_g[tier]) + " at " + where_g[tier] + "; Hessian " +
                         vcal::test::describe(worst_h[tier]) + " at " + where_h[tier]);
    }
    VCAL_CHECK(value_bits == 0);
    VCAL_CHECK(worst_g[0] <= vcal::tol::TOL_SCORE_HESSIAN_FD_REL && worst_h[0] <= vcal::tol::TOL_SCORE_HESSIAN_FD_REL);
    VCAL_CHECK(worst_g[1] <= vcal::tol::TOL_SCORE_HESSIAN_FD_REL_LARGE_N &&
               worst_h[1] <= vcal::tol::TOL_SCORE_HESSIAN_FD_REL_LARGE_N);
}

// Two grades with one PD pool exactly: their l_t is the pooled single grade's plus the log ratio of the binomial
// coefficients, log C(n1, d1) + log C(n2, d2) - log C(n1 + n2, d1 + d2). Different sums, so to rounding: in eps of
// the terms' size, within TOL_MULTIGRADE_POOLING_EPS, across PD, rho up to 0.9, and interior and one-sided pairs.
VCAL_TEST(two_grades_with_one_pd_pool_exactly) {
    const auto rule = q::parity_rule();
    double worst = 0.0;
    std::string where;
    int cases = 0;
    for (const double pd : {3e-4, 0.004, 0.03, 0.2}) {
        for (const double rho : {1e-3, 0.02, 0.12, 0.24, 0.5, 0.9}) {
            for (const auto& gg : std::vector<std::vector<GradeCount>>{{{400, 0}, {150, 0}}, {{400, 7}, {150, 4}},
                                                                       {{10000, 30}, {500, 0}}, {{1000000, 4000}, {200, 1}},
                                                                       {{50, 50}, {3, 3}}, {{1000, 0}, {20, 20}}}) {
                const double pds[2] = {pd, pd};
                const double two = o::multigrade_log_contrib(gg.data(), 2, pds, rho, rule);
                const Single::Obs pooled{gg[0].n + gg[1].n, gg[0].d + gg[1].d};
                const double c1 = vcal::special::lbinom(gg[0].n, gg[0].d), c2 = vcal::special::lbinom(gg[1].n, gg[1].d);
                const double cp = vcal::special::lbinom(pooled.n, pooled.d);
                const double one = Single{}.log_contrib(pooled, {pd, rho}, rule) + c1 + c2 - cp;
                const double scale = std::fabs(c1) + std::fabs(c2) + std::fabs(cp) + std::fabs(one);
                const double e = std::fabs(two - one) / (DBL_EPSILON * std::fmax(1.0, scale));
                ++cases;
                if (!(e <= worst)) {
                    worst = e;
                    where = "pd " + vcal::test::describe(pd) + " rho " + vcal::test::describe(rho) + " " +
                            std::to_string(gg[0].n) + ":" + std::to_string(gg[0].d);
                }
            }
        }
    }
    vcal::test::note(std::to_string(cases) + " pairs: worst " + vcal::test::describe(worst) + " eps at " + where);
    VCAL_CHECK(worst <= vcal::tol::TOL_MULTIGRADE_POOLING_EPS);
}

// rho -> 0 factorises: l_t tends to sum_k log Binom(d_k; n_k, PD_k). At rho = 1e-14 the gap is rho's own first-order
// effect (it grows 100-fold at rho = 1e-12), within TOL_MULTIGRADE_RHO_ZERO_ABS.
VCAL_TEST(small_rho_factorises_into_the_grades_binomials) {
    const auto rule = q::parity_rule();
    const std::vector<GradeCount> gg = {{1000, 3}, {300, 0}, {50, 7}, {20, 20}};
    const double pds[4] = {0.004, 1e-3, 0.1, 0.9};
    double gap[2];
    const double rhos[2] = {1e-14, 1e-12};
    for (int r = 0; r < 2; ++r) {
        const double l = o::multigrade_log_contrib(gg.data(), 4, pds, rhos[r], rule);
        double closed = 0.0;
        for (int k = 0; k < 4; ++k) {
            const double nn = static_cast<double>(gg[k].n), dd = static_cast<double>(gg[k].d);
            closed += vcal::special::lbinom(gg[k].n, gg[k].d) + (dd > 0 ? dd * std::log(pds[k]) : 0.0) +
                      (nn - dd > 0 ? (nn - dd) * std::log1p(-pds[k]) : 0.0);
        }
        gap[r] = std::fabs(l - closed);
    }
    vcal::test::note("|l - closed form|: " + vcal::test::describe(gap[0]) + " at rho 1e-14, " +
                     vcal::test::describe(gap[1]) + " at rho 1e-12");
    VCAL_CHECK(gap[0] <= vcal::tol::TOL_MULTIGRADE_RHO_ZERO_ABS);
    VCAL_CHECK(gap[1] > 10.0 * gap[0] && gap[1] < 1000.0 * gap[0]);  // first order in rho
}

// Core against ref/'s independent K-grade integral (its own bracket, golden-section and adaptive Gauss-Legendre;
// no shared quadrature code, D-031), on every K-grade golden case, and ref against mpmath: both in eps of the value,
// within TOL_XREF_PERIOD_LL_EPS and TOL_REF_LOG_MIXTURE_EPS, the single-grade tolerances.
VCAL_TEST(multigrade_core_and_ref_agree_on_every_golden_case) {
    const auto parity = q::parity_rule();
    double worst_core_ref = 0.0, worst_ref_mp = 0.0;
    std::string where_cr, where_rm;
    for (const auto& c : golden_cases()) {
        std::vector<vcalref::Period> gr;
        for (const auto& g : c.grades) gr.push_back({g.n, g.d});
        const double ref = vcalref::log_mixture_grades(c.pd, c.rho, gr);
        const double core = core_log_integral(c, parity);
        const double e1 = eps_units(core, ref), e2 = eps_units(ref, c.expected);
        if (!(e1 <= worst_core_ref)) {
            worst_core_ref = e1;
            where_cr = c.label();
        }
        if (!(e2 <= worst_ref_mp)) {
            worst_ref_mp = e2;
            where_rm = c.label();
        }
    }
    vcal::test::note("core vs ref: worst " + vcal::test::describe(worst_core_ref) + " eps at " + where_cr);
    vcal::test::note("ref vs mpmath: worst " + vcal::test::describe(worst_ref_mp) + " eps at " + where_rm);
    VCAL_CHECK(worst_core_ref <= vcal::tol::TOL_XREF_PERIOD_LL_EPS);
    VCAL_CHECK(worst_ref_mp <= vcal::tol::TOL_REF_LOG_MIXTURE_EPS);
}
