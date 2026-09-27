// SPDX-License-Identifier: Apache-2.0
// engine/ + backends/cpu (M1.4): surface, weighted reduction, refinement, calibration and
// thread-count determinism (ARCHITECTURE.md §5.1, §6; D-038, D-089, D-092, D-095).
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "core/reducers/argmax.hpp"
#include "engine/calibrate.hpp"
#include "engine/refine.hpp"
#include "engine/surface.hpp"
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

static_assert(tol::TOL_QUAD_CHECK_FLAG_ABS == e::kQuadratureCheckFlagAbs,
              "the tolerance register and the engine must agree on the quadrature-check threshold");

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

bool same_estimate(const e::Estimate2& a, const e::Estimate2& b) {
    return bits_equal(a.value[0], b.value[0]) && bits_equal(a.value[1], b.value[1]) &&
           bits_equal(a.se[0], b.se[0]) && bits_equal(a.se[1], b.se[1]) && bits_equal(a.corr, b.corr) &&
           bits_equal(a.loglik, b.loglik) && bits_equal(a.quad_check_max, b.quad_check_max) &&
           bits_equal(a.quad_check_total, b.quad_check_total) && a.quad_check_flagged == b.quad_check_flagged &&
           a.grid_index == b.grid_index && a.nan_count == b.nan_count && a.flags == b.flags;
}

std::string flags_text(std::uint32_t f) {
    std::string s;
    if (f & e::kFlagGridEdge) s += "edge ";
    if (f & e::kFlagFlatSurface) s += "flat ";
    if (f & e::kFlagQuadratureUnconverged) s += "quad ";
    if (f & e::kFlagRefinementRejected) s += "rejected ";
    if (f & e::kFlagNumeric) s += "numeric ";
    if (f & e::kFlagNearBound) s += "near-bound ";
    return s.empty() ? "none" : s;
}

// Input panel (not an expected value): 20 periods of 1000 obligors, defaults varying enough
// to identify a correlation well inside the grid.
std::vector<Obs> panel() {
    const std::int64_t d[] = {2, 5, 1, 9, 3, 0, 4, 12, 6, 2, 1, 3, 7, 15, 2, 4, 0, 5, 8, 3};
    std::vector<Obs> p;
    for (const auto x : d) p.push_back({1000, x});
    return p;
}

// PD on a logit axis, rho on a logit axis capped at the default bound (D-089).
vcal::Grid<2> coarse_grid() {
    return {{{2e-4, 0.1, 41, vcal::AxisScale::Logit}, {5e-3, vcal::kDefaultRhoUpper, 31, vcal::AxisScale::Logit}}};
}

struct Rules {
    Adaptive primary{vcal::quadrature::gauss_hermite_rule(64)};
    Adaptive check{vcal::quadrature::gauss_hermite_rule(128)};
};

e::Estimate2 run(int threads, const std::vector<Obs>& p, const vcal::Grid<2>& g, std::vector<double>& L,
                 e::Status* status = nullptr) {
    const Rules r;
    e::Estimate2 est{};
    const auto st = e::calibrate(Backend{threads}, Objective{}, r.primary, r.check, p.data(),
                                 static_cast<std::int64_t>(p.size()), g, L, est);
    if (status) *status = st;
    return est;
}

// Synthetic surface f(u) = 3 - (u - mu)' A (u - mu) / 2 in scaled coordinates.
struct Quadratic {
    const vcal::Grid<2>* g;
    double mu[2];
    double a00, a01, a11;
    double operator()(std::int32_t i0, std::int32_t i1) const {
        const double d0 = g->axis[0].scaled_at(i0) - mu[0];
        const double d1 = g->axis[1].scaled_at(i1) - mu[1];
        return 3.0 - 0.5 * (a00 * d0 * d0 + 2.0 * a01 * d0 * d1 + a11 * d1 * d1);
    }
};

}  // namespace

// --- weighted sums and the tiled reduction ---------------------------------------------------

VCAL_TEST(weighted_sum_is_compensated_and_skips_zero_weights) {
    const double L1[] = {1e16, 1.0, -1e16};
    const double ones[] = {1.0, 1.0, 1.0};
    VCAL_CHECK_EQ(e::weighted_sum(L1, 3, 1, ones, 0), 1.0);  // naive summation gives 0

    const double L2[] = {1.5, -kInf};
    const double w2[] = {1.0, 0.0};  // zero weight: the -inf period is excluded, not 0 * -inf
    VCAL_CHECK_EQ(e::weighted_sum(L2, 2, 1, w2, 0), 1.5);
    const double w3[] = {1.0, 2.0};
    VCAL_CHECK_EQ(e::weighted_sum(L2, 2, 1, w3, 0), -kInf);
    const double L4[] = {1.0, kNaN};
    VCAL_CHECK(std::isnan(e::weighted_sum(L4, 2, 1, ones, 0)));
}

// Several weight rows and a K that straddles tile boundaries: every thread count must give
// the serial answer exactly.
VCAL_TEST(reduce_weighted_matches_serial_for_any_thread_count) {
    const std::int64_t T = 5;
    const std::int64_t K = 2 * e::kReduceTileSize + 37;
    const std::int64_t B = 3;
    std::vector<double> L(static_cast<std::size_t>(T * K));
    std::uint64_t x = 12345;
    for (auto& v : L) {
        x = x * 6364136223846793005ull + 1442695040888963407ull;
        v = static_cast<double>(x >> 40) * 1e-6 - 5.0;
    }
    const double W[] = {1, 1, 1, 1, 1, 0, 2, 0, 1, 3, 1, 0, 0, 0, 1};
    const vcal::reducers::ArgMax r;
    std::vector<vcal::reducers::ArgMax::State> serial(B);
    for (std::int64_t b = 0; b < B; ++b) {
        serial[static_cast<std::size_t>(b)] = r.init();
        for (std::int64_t k = 0; k < K; ++k) r.push(serial[static_cast<std::size_t>(b)], k, e::weighted_sum(L.data(), T, K, W + b * T, k));
    }
    for (const int threads : {1, 2, 3, 8}) {
        std::vector<vcal::reducers::ArgMax::State> out(B);
        e::reduce_weighted(Backend{threads}, r, L.data(), T, K, W, B, out.data());
        for (std::int64_t b = 0; b < B; ++b) {
            const auto& s = out[static_cast<std::size_t>(b)];
            const auto& t = serial[static_cast<std::size_t>(b)];
            VCAL_CHECK(bits_equal(s.best, t.best) && s.k == t.k && s.nan_count == t.nan_count);
        }
    }
}

// --- refinement (D-038, D-095) --------------------------------------------------------------

// On an exact quadratic the stencil recovers the vertex and the covariance A^-1 (in scaled
// coordinates) to rounding, including a strong cross term, and converts SEs to natural scale
// by the delta method on log and logit axes.
VCAL_TEST(refine_recovers_vertex_and_covariance_of_a_quadratic) {
    const vcal::Grid<2> g{{{1e-4, 0.2, 41, vcal::AxisScale::Log}, {1e-3, 0.5, 31, vcal::AxisScale::Logit}}};
    const std::int32_t i_star[2] = {20, 14};
    const double h0 = g.axis[0].step();
    const double h1 = g.axis[1].step();
    const Quadratic q{&g, {g.axis[0].scaled_at(20) + 0.3 * h0, g.axis[1].scaled_at(14) - 0.2 * h1}, 40.0, 25.0, 30.0};
    const auto r = e::refine_2d(g, g.flatten(i_star), q);
    VCAL_CHECK_EQ(r.flags, 0u);

    const double det = q.a00 * q.a11 - q.a01 * q.a01;
    const double cov_u[3] = {q.a11 / det, -q.a01 / det, q.a00 / det};  // A^-1
    const double j0 = vcal::grid::dvalue_dscaled(g.axis[0].scale, q.mu[0]);
    const double j1 = vcal::grid::dvalue_dscaled(g.axis[1].scale, q.mu[1]);
    const double se0 = j0 * std::sqrt(cov_u[0]);
    const double se1 = j1 * std::sqrt(cov_u[2]);
    const double corr = cov_u[1] / std::sqrt(cov_u[0] * cov_u[2]);
    vcal::test::note("quadratic: vertex error " + describe(std::fabs(r.scaled[0] - q.mu[0])) + ", " +
                     describe(std::fabs(r.scaled[1] - q.mu[1])) + " (scaled units); corr " + describe(r.corr) +
                     " vs " + describe(corr));
    const double rel[5] = {std::fabs(r.value[0] / vcal::grid::from_scaled(g.axis[0].scale, q.mu[0]) - 1),
                           std::fabs(r.value[1] / vcal::grid::from_scaled(g.axis[1].scale, q.mu[1]) - 1),
                           std::fabs(r.se[0] / se0 - 1), std::fabs(r.se[1] / se1 - 1), std::fabs(r.corr / corr - 1)};
    double worst_rel = 0.0;
    for (const double x : rel) worst_rel = std::fmax(worst_rel, x);
    vcal::test::note("quadratic: worst relative error of value/se/corr " + describe(worst_rel));
    VCAL_CHECK_NEAR(r.scaled[0], q.mu[0], tol::TOL_REFINE_QUADRATIC_ABS);
    VCAL_CHECK_NEAR(r.scaled[1], q.mu[1], tol::TOL_REFINE_QUADRATIC_ABS);
    VCAL_CHECK_REL(r.value[0], vcal::grid::from_scaled(g.axis[0].scale, q.mu[0]), tol::TOL_REFINE_QUADRATIC_REL);
    VCAL_CHECK_REL(r.value[1], vcal::grid::from_scaled(g.axis[1].scale, q.mu[1]), tol::TOL_REFINE_QUADRATIC_REL);
    VCAL_CHECK_REL(r.se[0], se0, tol::TOL_REFINE_QUADRATIC_REL);
    VCAL_CHECK_REL(r.se[1], se1, tol::TOL_REFINE_QUADRATIC_REL);
    VCAL_CHECK_REL(r.corr, corr, tol::TOL_REFINE_QUADRATIC_REL);
}

VCAL_TEST(refine_falls_back_on_flat_or_non_concave_surfaces) {
    const vcal::Grid<2> g{{{0.0, 1.0, 11, vcal::AxisScale::Linear}, {0.0, 1.0, 11, vcal::AxisScale::Linear}}};
    const std::int32_t c[2] = {5, 5};
    const auto flat = e::refine_2d(g, g.flatten(c), [](std::int32_t, std::int32_t) { return 1.0; });
    VCAL_CHECK_EQ(flat.flags, static_cast<std::uint32_t>(e::kFlagFlatSurface));
    VCAL_CHECK_EQ(flat.value[0], g.axis[0].value_at(5));  // grid point, not extrapolated
    VCAL_CHECK(std::isnan(flat.se[0]));

    const Quadratic saddle{&g, {0.5, 0.5}, 10.0, 0.0, -10.0};  // indefinite
    const auto s = e::refine_2d(g, g.flatten(c), saddle);
    VCAL_CHECK_EQ(s.flags, static_cast<std::uint32_t>(e::kFlagFlatSurface));

    const Quadratic ridge{&g, {0.5, 0.5}, 10.0, 9.9, 10.0};  // nearly singular: flat along a diagonal
    const std::int32_t off[2] = {3, 5};  // not the argmax: the vertex lies beyond the stencil
    const auto rj = e::refine_2d(g, g.flatten(off), ridge);
    VCAL_CHECK(rj.flags & e::kFlagRefinementRejected);
    VCAL_CHECK_EQ(rj.value[0], g.axis[0].value_at(3));
    VCAL_CHECK(std::isfinite(rj.se[0]));  // curvature is valid; only the vertex is rejected
}

VCAL_TEST(refine_flags_grid_edge) {
    const vcal::Grid<2> g{{{0.0, 1.0, 11, vcal::AxisScale::Linear}, {0.0, 1.0, 11, vcal::AxisScale::Linear}}};
    const std::int32_t edge[2] = {5, 10};
    const Quadratic q{&g, {0.5, 1.2}, 10.0, 0.0, 10.0};
    const auto r = e::refine_2d(g, g.flatten(edge), q);
    VCAL_CHECK_EQ(r.flags, static_cast<std::uint32_t>(e::kFlagGridEdge));
    VCAL_CHECK_EQ(r.value[1], 1.0);
    VCAL_CHECK(std::isnan(r.se[1]));
}

// --- calibration end to end ------------------------------------------------------------------

// M1.4 exit criterion: bitwise-identical surface and estimate for 1, 2, 3 and 8 threads.
VCAL_TEST(calibrate_is_bitwise_identical_for_any_thread_count) {
    const auto p = panel();
    const auto g = coarse_grid();
    std::vector<double> L1;
    const auto ref = run(1, p, g, L1);
    vcal::test::note("estimate PD " + describe(ref.value[0]) + " (se " + describe(ref.se[0]) + "), rho " +
                     describe(ref.value[1]) + " (se " + describe(ref.se[1]) + "), corr " + describe(ref.corr) +
                     ", loglik " + describe(ref.loglik) + ", flags " + flags_text(ref.flags) + ", OpenMP " +
                     (Backend::has_openmp() ? "on" : "off"));
    for (const int threads : {2, 3, 8}) {
        std::vector<double> L;
        const auto est = run(threads, p, g, L);
        VCAL_CHECK(L.size() == L1.size() && std::memcmp(L.data(), L1.data(), L.size() * sizeof(double)) == 0);
        VCAL_CHECK(same_estimate(est, ref));
    }
}

// A clean interior estimate: no flags, quadrature check far below its threshold, and loglik
// equal (bitwise) to the objective summed at the estimate.
VCAL_TEST(calibrate_interior_estimate) {
    const auto p = panel();
    std::vector<double> L;
    e::Status st{};
    const auto est = run(0, p, coarse_grid(), L, &st);
    VCAL_CHECK(st == e::Status::Ok);
    VCAL_CHECK_EQ(est.flags, 0u);
    VCAL_CHECK_EQ(est.quad_check_flagged, std::int64_t{0});
    VCAL_CHECK(est.quad_check_max <= e::kQuadratureCheckFlagAbs);
    vcal::test::note("quadrature check (N=64 vs 128): max " + describe(est.quad_check_max) + ", total " +
                     describe(est.quad_check_total));

    const Rules r;
    const double v[2] = {est.value[0], est.value[1]};
    std::vector<double> lt;
    for (const auto& o : p) lt.push_back(Objective{}.log_contrib(o, Objective::theta(v), r.primary));
    const std::vector<double> ones(lt.size(), 1.0);
    VCAL_CHECK(bits_equal(est.loglik, e::weighted_sum(lt.data(), static_cast<std::int64_t>(lt.size()), 1, ones.data(), 0)));
}

// Sub-grid accuracy: refining the coarse grid lands where a 10x finer grid around the estimate
// puts the maximum, to a small fraction of the standard error (measured in each axis's scaled
// coordinate), and the SEs agree. Both fits take SEs from the objective's Hessian at their own
// estimate (D-119), so the SE difference reflects only the difference between the estimates.
VCAL_TEST(calibrate_refinement_agrees_with_a_fine_grid) {
    const auto p = panel();
    const auto g = coarse_grid();
    std::vector<double> L;
    const auto coarse = run(0, p, g, L);
    vcal::Grid<2> fine = g;
    for (int a = 0; a < 2; ++a) {
        const auto s = g.axis[a].scale;
        const double u = vcal::grid::to_scaled(s, coarse.value[a]);
        const double h = g.axis[a].step();
        fine.axis[a] = {vcal::grid::from_scaled(s, u - 2.0 * h), vcal::grid::from_scaled(s, u + 2.0 * h), 41, s};
    }
    const auto fine_est = run(0, p, fine, L);
    // The fine box spans only +-2 coarse steps, so its estimate is near its bounds by design.
    VCAL_CHECK_EQ(fine_est.flags & ~static_cast<std::uint32_t>(e::kFlagNearBound), 0u);
    for (int a = 0; a < 2; ++a) {
        const auto s = g.axis[a].scale;
        const double du = std::fabs(vcal::grid::to_scaled(s, coarse.value[a]) - vcal::grid::to_scaled(s, fine_est.value[a]));
        const double se_u = fine_est.se[a] / vcal::grid::dvalue_dscaled(s, vcal::grid::to_scaled(s, fine_est.value[a]));
        vcal::test::note("axis " + std::to_string(a) + ": |coarse - fine| = " + describe(du / se_u) +
                         " SE; coarse SE / fine SE = " + describe(coarse.se[a] / fine_est.se[a]));
        VCAL_CHECK(du / se_u <= tol::TOL_REFINE_VS_FINE_GRID_SE);
        VCAL_CHECK_REL(coarse.se[a], fine_est.se[a], tol::TOL_SE_COARSE_VS_FINE_REL);
    }
}

// D-089: estimates on the grid's bounds are flagged, never silently clamped.
VCAL_TEST(calibrate_flags_estimates_on_the_grid_edge) {
    std::vector<double> L;
    std::vector<Obs> no_defaults(10, Obs{1000, 0});  // PD MLE below any positive grid point
    const auto a = run(0, no_defaults, coarse_grid(), L);
    VCAL_CHECK(a.flags & e::kFlagGridEdge);
    VCAL_CHECK_EQ(a.value[0], coarse_grid().axis[0].lo);

    std::vector<Obs> clustered;  // extreme clustering: rho MLE above the default bound 0.5
    for (int t = 0; t < 10; ++t) clustered.push_back({1000, t % 2 == 0 ? 0 : 400});
    const auto b = run(0, clustered, coarse_grid(), L);
    VCAL_CHECK(b.flags & e::kFlagGridEdge);
    VCAL_CHECK_EQ(b.value[1], vcal::kDefaultRhoUpper);
    vcal::test::note("edge cases: flags " + flags_text(a.flags) + "/ " + flags_text(b.flags));
}

VCAL_TEST(calibrate_rejects_invalid_input) {
    std::vector<double> L;
    e::Status st{};
    std::vector<Obs> bad{{1000, 3}, {10, 11}};
    run(0, bad, coarse_grid(), L, &st);
    VCAL_CHECK(st == e::Status::InvalidPanel);
    std::vector<Obs> zero_n{{0, 0}};
    run(0, zero_n, coarse_grid(), L, &st);
    VCAL_CHECK(st == e::Status::InvalidPanel);
    vcal::Grid<2> g = coarse_grid();
    g.axis[1].hi = 1.0;  // logit axis must stay below 1
    run(0, panel(), g, L, &st);
    VCAL_CHECK(st == e::Status::InvalidGrid);
}

// --- SEs from the objective's Hessian at the estimate (D-119) ------------------------------------

// The step is 0.15 x the stencil SE; 0.1 and 0.2 must give nearly the same SEs (neither too small,
// which amplifies rounding, nor too large, which reintroduces curvature bias).
VCAL_TEST(hessian_step_sensitivity) {
    const auto p = panel();
    const auto g = coarse_grid();
    const Rules r;
    std::vector<double> L;
    e::Estimate2 lo{}, mid{}, hi{};
    VCAL_REQUIRE(e::calibrate(Backend{}, Objective{}, r.primary, r.check, p.data(), 20, g, L, lo, 0.10) == e::Status::Ok);
    VCAL_REQUIRE(e::calibrate(Backend{}, Objective{}, r.primary, r.check, p.data(), 20, g, L, mid) == e::Status::Ok);
    VCAL_REQUIRE(e::calibrate(Backend{}, Objective{}, r.primary, r.check, p.data(), 20, g, L, hi, 0.20) == e::Status::Ok);
    double worst = 0.0;
    for (int a = 0; a < 2; ++a) {
        worst = std::fmax(worst, std::fabs(lo.se[a] / mid.se[a] - 1.0));
        worst = std::fmax(worst, std::fabs(hi.se[a] / mid.se[a] - 1.0));
    }
    vcal::test::note("SEs at step 0.10 / 0.15 / 0.20 x stencil SE: PD " + describe(lo.se[0]) + " / " + describe(mid.se[0]) +
                     " / " + describe(hi.se[0]) + "; rho " + describe(lo.se[1]) + " / " + describe(mid.se[1]) + " / " +
                     describe(hi.se[1]) + "; worst relative spread " + describe(worst));
    VCAL_CHECK(worst <= tol::TOL_SE_STEP_SENSITIVITY_REL);
}

// An interior estimate within kNearBoundSe SEs of an axis bound is flagged (Wald SEs are
// unreliable there); the same panel on the usual grid is not.
VCAL_TEST(calibrate_flags_estimates_near_a_bound) {
    const auto p = panel();
    std::vector<double> L;
    const auto far = run(0, p, coarse_grid(), L);
    VCAL_CHECK(!(far.flags & e::kFlagNearBound));
    vcal::Grid<2> tight = coarse_grid();
    tight.axis[1].lo = 0.03;  // rho-hat ~ 0.056 sits inside, but within 2 SEs of this bound
    const auto near = run(0, p, tight, L);
    vcal::test::note("near-bound case: rho " + describe(near.value[1]) + ", flags " + flags_text(near.flags));
    VCAL_CHECK(near.flags & e::kFlagNearBound);
    VCAL_CHECK(!(near.flags & e::kFlagGridEdge));
    VCAL_CHECK(std::isfinite(near.se[1]));  // reported, but flagged
}

// --- periods with equal observations share one evaluation (D-122) ----------------------------------

// The deduplicated surface is bitwise the surface with every period evaluated, for any thread
// count. The panel repeats observations (d = 0, 2, 3, 4 and 5 each occur more than once).
VCAL_TEST(surface_evaluates_repeated_observations_once_and_bitwise_identically) {
    const auto p = panel();
    const auto g = coarse_grid();
    const Rules r;
    const std::int64_t T = static_cast<std::int64_t>(p.size());
    const std::int64_t K = g.size();
    std::vector<double> naive(static_cast<std::size_t>(T * K));
    for (std::int64_t k = 0; k < K; ++k) {
        double v[2];
        g.values(k, v);
        for (std::int64_t t = 0; t < T; ++t) {
            naive[static_cast<std::size_t>(t * K + k)] = Objective{}.log_contrib(p[static_cast<std::size_t>(t)], Objective::theta(v), r.primary);
        }
    }
    int distinct = 0;
    for (std::int64_t t = 0; t < T; ++t) {
        bool seen = false;
        for (std::int64_t u = 0; u < t; ++u) seen = seen || p[static_cast<std::size_t>(u)] == p[static_cast<std::size_t>(t)];
        distinct += seen ? 0 : 1;
    }
    VCAL_REQUIRE(distinct < T);
    for (const int threads : {1, 3, 0}) {
        std::vector<double> L(naive.size());
        e::evaluate_surface(Backend{threads}, Objective{}, r.primary, p.data(), T, g, L.data());
        bool same = true;
        for (std::size_t i = 0; i < L.size(); ++i) same = same && bits_equal(L[i], naive[i]);
        VCAL_CHECK(same);
    }
    vcal::test::note(std::to_string(distinct) + " distinct observations over " + std::to_string(T) + " periods");
}
