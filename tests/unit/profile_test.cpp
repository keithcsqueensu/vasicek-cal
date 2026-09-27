// SPDX-License-Identifier: Apache-2.0
// engine/profile.hpp (M2): profile-likelihood intervals (D-128..D-130).
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
