// SPDX-License-Identifier: Apache-2.0
//
// Single calibration of a 2-parameter objective (ARCHITECTURE.md §5.1):
//   1. validate grid and panel (Objective::panel_error);
//   2. L = evaluate_surface (T x K);
//   3. W = one row of ones; reduce_weighted(ArgMax) -> grid argmax k*;
//   4. refine_2d on the 3x3 stencil of the same weighted sums -> estimate and flags;
//   4b. SEs from the observed information at the estimate (D-119): the Hessian of the objective
//      itself by central differences in the axes' scaled coordinates, step = hessian_step_fraction
//      (default 0.15) x the stencil SE, delta method to natural scale. Headline fits only;
//      resampling never calls this, so W x L reuse is unaffected. An estimate within
//      kNearBoundSe SEs of an axis bound is flagged kFlagNearBound: Wald SEs are unreliable
//      there, and profile-likelihood intervals (M2) are the right tool;
//   5. at the estimate, recompute each l_t with the primary integrator and with the check
//      integrator (parity: the rule of D-118 and its doubled variant). loglik is the primary sum. The check
//      reports max_t and sum_t |l_t(primary) - l_t(check)|, and sets
//      kFlagQuadratureUnconverged if any period differs by more than its threshold (D-092, D-120).
// The caller owns L (resized here) so resampling can reuse it (M2).
#pragma once

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <vector>

#include "core/grid.hpp"
#include "core/reducers/argmax.hpp"
#include "engine/refine.hpp"
#include "engine/surface.hpp"

namespace vcal::engine {

// A period counts as unconverged when |l_t(rule) - l_t(doubled rule)| exceeds
//     max(kQuadratureCheckFlagAbs, kQuadratureCheckRoundingUlps * eps * scale_t),
// scale_t = Objective::rounding_scale (the size of l_t's terms). Below that the difference is
// rounding, and no integration error could be detected there anyway (D-092, D-120).
inline constexpr double kQuadratureCheckFlagAbs = 1e-10;
inline constexpr double kQuadratureCheckRoundingUlps = 64.0;

// D-119: Hessian step as a fraction of the stencil SE, and the near-bound distance in SEs.
inline constexpr double kHessianStepFraction = 0.15;
inline constexpr double kNearBoundSe = 2.0;

enum class Status : std::int32_t { Ok = 0, InvalidGrid = 1, InvalidPanel = 2, SurfaceUndefined = 3 };

struct Estimate2 {
    double value[2];  // natural scale, grid axis order
    double se[2];
    double corr;
    double loglik;                // sum_t l_t at value, primary integrator
    double quad_check_max;        // max_t |l_t(primary) - l_t(check)| at value
    double quad_check_total;      // sum_t of the same
    std::int64_t quad_check_flagged;  // periods above the D-120 threshold (see above)
    std::int64_t grid_index;      // argmax k* before refinement
    std::int64_t nan_count;       // NaN surface values seen by the reducer
    std::uint32_t flags;
};

template <class Backend, class Objective, class Integrator>
Status calibrate(const Backend& backend, const Objective& objective, const Integrator& primary,
                 const Integrator& check, const typename Objective::Obs* obs, std::int64_t periods,
                 const Grid<2>& grid, std::vector<double>& L, Estimate2& out,
                 double hessian_step_fraction = kHessianStepFraction) {
    static_assert(Objective::n_params == 2, "calibrate handles 2-parameter objectives");
    if (grid_error(grid) != nullptr) return Status::InvalidGrid;
    if (Objective::panel_error(obs, periods) != nullptr) return Status::InvalidPanel;

    const std::int64_t K = grid.size();
    L.assign(static_cast<std::size_t>(periods * K), 0.0);
    evaluate_surface(backend, objective, primary, obs, periods, grid, L.data());

    const std::vector<double> ones(static_cast<std::size_t>(periods), 1.0);
    const reducers::ArgMax argmax;
    reducers::ArgMax::State best{};
    reduce_weighted(backend, argmax, L.data(), periods, K, ones.data(), 1, &best);
    if (best.k < 0) return Status::SurfaceUndefined;

    const auto value_at = [&](std::int32_t i0, std::int32_t i1) {
        const std::int32_t idx[2] = {i0, i1};
        return weighted_sum(L.data(), periods, K, ones.data(), grid.flatten(idx));
    };
    const Refined2 r = refine_2d(grid, best.k, value_at);

    out = Estimate2{};
    out.value[0] = r.value[0];
    out.value[1] = r.value[1];
    out.se[0] = r.se[0];
    out.se[1] = r.se[1];
    out.corr = r.corr;
    out.grid_index = best.k;
    out.nan_count = best.nan_count;
    out.flags = r.flags | (best.nan_count > 0 ? kFlagNumeric : 0u);

    const auto theta = Objective::theta(out.value);
    std::vector<double> l_primary(static_cast<std::size_t>(periods));
    for (std::int64_t t = 0; t < periods; ++t) {
        const double a = objective.log_contrib(obs[t], theta, primary);
        const double b = objective.log_contrib(obs[t], theta, check);
        const double diff = std::fabs(a - b);
        const double threshold = std::fmax(kQuadratureCheckFlagAbs, kQuadratureCheckRoundingUlps * DBL_EPSILON *
                                                                          Objective::rounding_scale(obs[t], a));
        l_primary[static_cast<std::size_t>(t)] = a;
        out.quad_check_max = std::fmax(out.quad_check_max, diff);
        out.quad_check_total += diff;
        if (!(diff <= threshold)) ++out.quad_check_flagged;  // NaN counts as flagged
    }
    out.loglik = weighted_sum(l_primary.data(), periods, 1, ones.data(), 0);
    if (out.quad_check_flagged > 0) out.flags |= kFlagQuadratureUnconverged;

    // D-119: observed information at the estimate from the objective itself.
    if (!(out.flags & (kFlagGridEdge | kFlagFlatSurface))) {
        const double s0 = hessian_step_fraction * r.se_scaled[0];
        const double s1 = hessian_step_fraction * r.se_scaled[1];
        std::vector<double> lt(static_cast<std::size_t>(periods));
        const auto ll = [&](double d0, double d1) {
            const double v[2] = {grid::from_scaled(grid.axis[0].scale, r.scaled[0] + d0),
                                 grid::from_scaled(grid.axis[1].scale, r.scaled[1] + d1)};
            const auto th = Objective::theta(v);
            for (std::int64_t t = 0; t < periods; ++t) lt[static_cast<std::size_t>(t)] = objective.log_contrib(obs[t], th, primary);
            return weighted_sum(lt.data(), periods, 1, ones.data(), 0);
        };
        const double f0 = out.loglik;
        const double h00 = (ll(s0, 0.0) - 2.0 * f0 + ll(-s0, 0.0)) / (s0 * s0);
        const double h11 = (ll(0.0, s1) - 2.0 * f0 + ll(0.0, -s1)) / (s1 * s1);
        const double h01 = (ll(s0, s1) - ll(s0, -s1) - ll(-s0, s1) + ll(-s0, -s1)) / (4.0 * s0 * s1);
        const double det = h00 * h11 - h01 * h01;
        if (h00 < 0.0 && det > 0.0) {
            const double su0 = std::sqrt(-h11 / det);
            const double su1 = std::sqrt(-h00 / det);
            const double j0 = grid::dvalue_dscaled(grid.axis[0].scale, r.scaled[0]);
            const double j1 = grid::dvalue_dscaled(grid.axis[1].scale, r.scaled[1]);
            out.se[0] = j0 * su0;
            out.se[1] = j1 * su1;
            out.corr = (h01 / det) / (su0 * su1);
            const double su[2] = {su0, su1};
            for (int a = 0; a < 2; ++a) {
                const double to_bound = std::fmin(r.scaled[a] - grid.axis[a].scaled_lo(), grid.axis[a].scaled_hi() - r.scaled[a]);
                if (to_bound < kNearBoundSe * su[a]) out.flags |= kFlagNearBound;
            }
        } else {
            out.se[0] = out.se[1] = out.corr = std::nan("");
            out.flags |= kFlagFlatSurface;
        }
    }
    return Status::Ok;
}

}  // namespace vcal::engine
