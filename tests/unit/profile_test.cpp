// SPDX-License-Identifier: Apache-2.0
// engine/profile.hpp (M2): profile-likelihood intervals (D-128..D-130); engine/conditional_pd.hpp (S-23).
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "engine/calibrate.hpp"
#include "engine/conditional_pd.hpp"
#include "engine/profile.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace e = vcal::engine;
namespace tol = vcal::tol;
using vcal::test::describe;
using Backend = vcal::backends::CpuBackend;
using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
using Adaptive = vcal::quadrature::GaussHermiteAdaptive<vcal::PrecisionF64>;
using Obs = Objective::Obs;

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

// The M1.4 test panel: 20 periods of 1000 obligors.
std::vector<Obs> panel() {
    const std::int64_t d[] = {2, 5, 1, 9, 3, 0, 4, 12, 6, 2, 1, 3, 7, 15, 2, 4, 0, 5, 8, 3};
    std::vector<Obs> p;
    for (const auto x : d) p.push_back({1000, x});
    return p;
}

vcal::Grid<2> grid() {
    return {{{2e-4, 0.1, 41, vcal::AxisScale::Logit}, {5e-3, vcal::kDefaultRhoUpper, 31, vcal::AxisScale::Logit}}};
}

struct Rules {
    Adaptive primary{vcal::quadrature::gauss_hermite_rule(64)};
    Adaptive check{vcal::quadrature::gauss_hermite_rule(128)};
};

struct Fitted {
    e::Estimate2 est{};
    std::vector<double> L;
    e::ProfileIntervals2 prof{};
};

Fitted fit(const std::vector<Obs>& p, const vcal::Grid<2>& g, int threads = 0) {
    const Rules r;
    Fitted f;
    VCAL_REQUIRE(e::calibrate(Backend{threads}, Objective{}, r.primary, r.check, p.data(),
                              static_cast<std::int64_t>(p.size()), g, f.L, f.est) == e::Status::Ok);
    f.prof = e::profile_intervals(Objective{}, r.primary, p.data(), static_cast<std::int64_t>(p.size()), g, f.L, f.est);
    return f;
}

// Independent of the engine's bracketing: golden-section over the whole box, very tight.
double golden_max(const std::function<double(double)>& f, double a, double b) {
    const double r = (std::sqrt(5.0) - 1.0) / 2.0;
    double x1 = b - r * (b - a), x2 = a + r * (b - a), f1 = f(x1), f2 = f(x2);
    while (b - a > 1e-11) {
        if (f1 >= f2) {
            b = x2; x2 = x1; f2 = f1; x1 = b - r * (b - a); f1 = f(x1);
        } else {
            a = x1; x1 = x2; f1 = f2; x2 = a + r * (b - a); f2 = f(x2);
        }
    }
    return f1 >= f2 ? f1 : f2;
}

}  // namespace

// Each solved endpoint satisfies P_a(u) = l_max - c, checked with an independent dense profile
// (golden-section over the whole nuisance range) and an independent maximum (nested
// golden-section over the whole box).
VCAL_TEST(profile_endpoints_solve_the_threshold) {
    const auto p = panel();
    const auto g = grid();
    const auto f = fit(p, g);
    const Rules r;
    const auto ll = [&](double u0, double u1) {
        const double v[2] = {vcal::grid::from_scaled(g.axis[0].scale, u0), vcal::grid::from_scaled(g.axis[1].scale, u1)};
        double s = 0.0;
        for (const auto& y : p) s += Objective{}.log_contrib(y, Objective::theta(v), r.primary);
        return s;
    };
    const double box[2][2] = {{g.axis[0].scaled_lo(), g.axis[0].scaled_hi()}, {g.axis[1].scaled_lo(), g.axis[1].scaled_hi()}};
    const auto profile = [&](int a, double u) {
        return golden_max([&](double w) { return a == 0 ? ll(u, w) : ll(w, u); }, box[1 - a][0], box[1 - a][1]);
    };
    const double l_max = golden_max([&](double u) { return profile(0, u); }, box[0][0], box[0][1]);
    vcal::test::note("l_max: engine " + describe(f.prof.loglik_max) + ", dense " + describe(l_max) + "; " +
                     std::to_string(f.prof.evaluations) + " panel evaluations; self-reported residual " +
                     describe(f.prof.residual_max));
    VCAL_CHECK(std::fabs(f.prof.loglik_max - l_max) <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
    double worst = 0.0;
    for (int a = 0; a < 2; ++a) {
        VCAL_CHECK_EQ(f.prof.flags[a], 0u);
        for (const double v : {f.prof.lo[a], f.prof.hi[a]}) {
            const double u = vcal::grid::to_scaled(g.axis[a].scale, v);
            worst = std::fmax(worst, std::fabs(profile(a, u) - (l_max - e::kProfileThreshold95)));
        }
        VCAL_CHECK(f.prof.lo[a] < f.est.value[a] && f.est.value[a] < f.prof.hi[a]);
    }
    vcal::test::note("PD [" + describe(f.prof.lo[0]) + ", " + describe(f.prof.hi[0]) + "], rho [" +
                     describe(f.prof.lo[1]) + ", " + describe(f.prof.hi[1]) +
                     "]; worst |P(endpoint) - (l_max - c)| against the dense profile " + describe(worst));
    VCAL_CHECK(worst <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
    VCAL_CHECK(f.prof.residual_max <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
}

// An estimate near the lower rho bound: the profile of rho stays above the threshold down to
// the bound, so that end is the bound, flagged, never extrapolated (D-130).
VCAL_TEST(profile_interval_is_truncated_at_a_bound) {
    auto g = grid();
    g.axis[1].lo = 0.03;  // rho-hat ~ 0.056 sits within 2 SEs of this bound (M1.7b near-bound test)
    const auto f = fit(panel(), g);
    vcal::test::note("rho [" + describe(f.prof.lo[1]) + ", " + describe(f.prof.hi[1]) + "], flags " +
                     std::to_string(f.prof.flags[1]));
    VCAL_CHECK(f.prof.flags[1] & e::kIntervalLowerTruncated);
    VCAL_CHECK(!(f.prof.flags[1] & e::kIntervalUpperTruncated));
    VCAL_CHECK(bits_equal(f.prof.lo[1], g.axis[1].lo));
    VCAL_CHECK(f.prof.hi[1] > f.est.value[1]);

    // No defaults at all: PD-hat is the lower bound (grid edge); its interval starts there.
    const auto z = fit(std::vector<Obs>(10, Obs{1000, 0}), grid());
    VCAL_CHECK(z.est.flags & e::kFlagGridEdge);
    VCAL_CHECK(z.prof.flags[0] & e::kIntervalLowerTruncated);
    VCAL_CHECK(bits_equal(z.prof.lo[0], grid().axis[0].lo));
    vcal::test::note("no defaults: PD [" + describe(z.prof.lo[0]) + ", " + describe(z.prof.hi[0]) + "], flags " +
                     std::to_string(z.prof.flags[0]));
}

// Flat or numeric fits get no interval.
VCAL_TEST(profile_interval_not_computed_for_flat_fits) {
    auto f = fit(panel(), grid());
    e::Estimate2 flat = f.est;
    flat.flags |= e::kFlagFlatSurface;
    const Rules r;
    const auto p = panel();
    const auto out = e::profile_intervals(Objective{}, r.primary, p.data(), static_cast<std::int64_t>(p.size()), grid(),
                                          f.L, flat);
    VCAL_CHECK_EQ(out.flags[0], e::kIntervalNotComputed);
    VCAL_CHECK_EQ(out.flags[1], e::kIntervalNotComputed);
    VCAL_CHECK(std::isnan(out.lo[0]) && std::isnan(out.hi[1]));
}

// The intervals do not depend on the number of threads used for the surface.
VCAL_TEST(profile_interval_is_bitwise_identical_for_any_thread_count) {
    const auto a = fit(panel(), grid(), 1);
    const auto b = fit(panel(), grid(), 3);
    for (int k = 0; k < 2; ++k) {
        VCAL_CHECK(bits_equal(a.prof.lo[k], b.prof.lo[k]));
        VCAL_CHECK(bits_equal(a.prof.hi[k], b.prof.hi[k]));
    }
    VCAL_CHECK(bits_equal(a.prof.loglik_max, b.prof.loglik_max));
}

// --- the conditional PD at the 99.9% adverse factor level (S-23) -----------------------------------

namespace {

// The dense profile of q, independent of the engine's bracketing: for c = logit^-1(s), golden-section
// over the whole feasible range of logit(rho), found by a dense scan, with Phi from std::erfc.
double dense_q_profile(const std::vector<Obs>& p, const vcal::Grid<2>& g, double s) {
    const Rules r;
    const double c = 1.0 / (1.0 + std::exp(-s));
    const double x = vcal::special::probit(c);
    const auto pd_at = [&](double w) {
        const double rho = 1.0 / (1.0 + std::exp(-w));
        return 0.5 * std::erfc(-(std::sqrt(1.0 - rho) * x - std::sqrt(rho) * e::kAdverseZ999) / std::sqrt(2.0));
    };
    const double w_lo = g.axis[1].scaled_lo(), w_hi = g.axis[1].scaled_hi();
    double a = 0.0, b = 0.0;
    bool any = false;
    for (int k = 0; k <= 4000; ++k) {
        const double w = w_lo + (w_hi - w_lo) * k / 4000.0;
        const double pd = pd_at(w);
        if (!(pd >= g.axis[0].lo && pd <= g.axis[0].hi)) continue;
        if (!any) a = w;
        b = w;
        any = true;
    }
    VCAL_REQUIRE(any);
    return golden_max(
        [&](double w) {
            const double v[2] = {pd_at(w), 1.0 / (1.0 + std::exp(-w))};
            double sum = 0.0;
            for (const auto& y : p) sum += Objective{}.log_contrib(y, Objective::theta(v), r.primary);
            return sum;
        },
        a, b);
}

e::ConditionalPdInterval q_interval(const std::vector<Obs>& p, const vcal::Grid<2>& g, const Fitted& f) {
    const Rules r;
    return e::conditional_pd_interval(Objective{}, r.primary, p.data(), static_cast<std::int64_t>(p.size()), g, f.L,
                                      f.est, f.prof);
}

}  // namespace

// q against mpmath (40 digits), including the corners of the recovery box; the logit form and
// the gradient against central differences.
VCAL_TEST(conditional_pd_matches_mpmath) {
    struct Case {
        double pd, rho, q, s;
    };
    const Case cases[] = {
        {0.001, 0.02, 0.0036795217954010484033, -5.6012861742597981117},
        {0.01, 0.12, 0.090325831326065272963, -2.3096629982791540639},
        {0.05, 0.24, 0.4402971526712318421, -0.23995616717710972256},
        {0.0001, 0.001, 0.0001455427998040301266, -8.8348948044874205319},
        {0.2, 0.5, 0.97128344958463191349, 3.5211447147349791015},
    };
    double worst = 0.0, worst_grad = 0.0;
    for (const auto& k : cases) {
        worst = std::fmax(worst, std::fabs(e::conditional_pd(k.pd, k.rho) / k.q - 1.0));
        const auto gr = e::conditional_pd_logit_gradient(k.pd, k.rho);
        worst = std::fmax(worst, std::fabs(gr.s / k.s - 1.0));
        const double u[2] = {std::log(k.pd / (1.0 - k.pd)), std::log(k.rho / (1.0 - k.rho))};
        const double step = 1e-5;
        for (int a = 0; a < 2; ++a) {
            const auto s_at = [&](double du) {
                double v[2] = {u[0], u[1]};
                v[a] += du;
                return e::conditional_pd_logit_gradient(1.0 / (1.0 + std::exp(-v[0])), 1.0 / (1.0 + std::exp(-v[1]))).s;
            };
            const double fd = (s_at(step) - s_at(-step)) / (2.0 * step);
            worst_grad = std::fmax(worst_grad, std::fabs(gr.ds_du[a] / fd - 1.0));
        }
    }
    vcal::test::note("q and logit(q) against mpmath: worst relative error " + describe(worst) +
                     "; gradient against central differences " + describe(worst_grad));
    VCAL_CHECK(worst <= tol::TOL_CONDITIONAL_PD_REL);
    VCAL_CHECK(worst_grad <= tol::TOL_CONDITIONAL_PD_GRADIENT_REL);
}

// Each end of q's interval satisfies P_q(c) = l_max - threshold against the dense profile; the
// interval contains q-hat and is not box-limited on this panel.
VCAL_TEST(conditional_pd_interval_endpoints_solve_the_threshold) {
    const auto p = panel();
    const auto g = grid();
    const auto f = fit(p, g);
    const auto q = q_interval(p, g, f);
    VCAL_CHECK_EQ(q.flags, 0u);
    VCAL_CHECK(q.lo < q.estimate && q.estimate < q.hi);
    double worst = 0.0;
    for (const double c : {q.lo, q.hi}) {
        const double s = std::log(c) - std::log1p(-c);
        worst = std::fmax(worst, std::fabs(dense_q_profile(p, g, s) - (f.prof.loglik_max - e::kProfileThreshold95)));
    }
    vcal::test::note("q-hat " + describe(q.estimate) + ", interval [" + describe(q.lo) + ", " + describe(q.hi) +
                     "]; " + std::to_string(q.evaluations) + " panel evaluations; self-reported residual " +
                     describe(q.residual_max) + "; against the dense profile " + describe(worst));
    VCAL_CHECK(worst <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
    VCAL_CHECK(q.residual_max <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
}

// No defaults at all: the lower end is the limit of q in the box, flagged truncated and
// box-limited, never extrapolated. With the rho axis starting near rho-hat, the lower end is held
// by that bound through the nuisance: box-limited, but solved, not truncated.
VCAL_TEST(conditional_pd_interval_at_the_box) {
    const auto z = std::vector<Obs>(10, Obs{1000, 0});
    const auto gz = grid();
    const auto fz = fit(z, gz);
    const auto qz = q_interval(z, gz, fz);
    const double q_min = e::conditional_pd(gz.axis[0].lo, gz.axis[1].lo);
    vcal::test::note("no defaults: q [" + describe(qz.lo) + ", " + describe(qz.hi) + "], flags " +
                     std::to_string(qz.flags) + "; q at the box's lower corner " + describe(q_min));
    VCAL_CHECK(qz.flags & e::kIntervalLowerTruncated);
    VCAL_CHECK(qz.flags & e::kIntervalLowerBoxLimited);
    VCAL_CHECK(!(qz.flags & e::kIntervalUpperTruncated));
    VCAL_CHECK_REL(qz.lo, q_min, tol::TOL_CONDITIONAL_PD_REL);

    auto g = grid();
    g.axis[1].lo = 0.03;
    const auto p = panel();
    const auto f = fit(p, g);
    const auto q = q_interval(p, g, f);
    vcal::test::note("rho from 0.03: q [" + describe(q.lo) + ", " + describe(q.hi) + "], flags " +
                     std::to_string(q.flags));
    VCAL_CHECK(q.flags & e::kIntervalLowerBoxLimited);
    VCAL_CHECK(!(q.flags & (e::kIntervalLowerTruncated | e::kIntervalUpperTruncated | e::kIntervalNotComputed)));
    VCAL_CHECK(q.lo < q.estimate && q.estimate < q.hi);
    VCAL_CHECK(q.residual_max <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
}

// Flat fits get no interval; the interval does not depend on the thread count.
VCAL_TEST(conditional_pd_interval_not_computed_for_flat_fits_and_thread_independent) {
    const auto p = panel();
    auto f = fit(p, grid());
    Fitted flat = f;
    flat.est.flags |= e::kFlagFlatSurface;
    const auto none = q_interval(p, grid(), flat);
    VCAL_CHECK_EQ(none.flags, e::kIntervalNotComputed);
    VCAL_CHECK(std::isnan(none.lo) && std::isnan(none.hi));
    const auto a = q_interval(p, grid(), fit(p, grid(), 1));
    const auto b = q_interval(p, grid(), fit(p, grid(), 3));
    VCAL_CHECK(bits_equal(a.lo, b.lo));
    VCAL_CHECK(bits_equal(a.hi, b.hi));
}
