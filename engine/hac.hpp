// SPDX-License-Identifier: Apache-2.0
//
// The HAC (Newey-West) sandwich SE of a 2-parameter fit: a `native` option (S-48, D-209,
// docs/methodology/hac_se.md). It changes no parity result.
//
// calibrate's SE assumes independent periods: V = (-H)^-1, H the panel log-likelihood's Hessian at the
// estimate. When periods are dependent the scores are too, and V understates the estimate's variance
// (S-35: by f = (1 + phi) / (1 - phi) for PD under an AR(1) factor). The sandwich replaces the
// information by the long-run variance of the per-period scores:
//
//   u        the estimate in the axes' scaled coordinates (logit on the recovery grid), u_a =
//            to_scaled(value_a);
//   s_t      period t's score in u: d l_t / d u_a = (d l_t / d v_a) (dv_a / du_a), from the objective's
//            analytic derivatives (log_contrib_derivs, D-170), at the estimate;
//   c_t      s_t - s_bar, the scores centred on their mean over the T periods. At an exact maximum
//            s_bar = 0; calibrate's estimate is a refined grid point (refine_2d), whose score is small but
//            not zero, and centring keeps that from entering S;
//   S        sum_t c_t c_t' + sum_{l=1}^{L} w_l sum_{t=l}^{T-1} (c_t c_{t-l}' + c_{t-l} c_t'), with the
//            Bartlett weights w_l = 1 - l / (L + 1), which keep S positive semi-definite; plain sums in
//            period order, no small-sample factor;
//   H        the panel Hessian in u, the one calibrate's SE uses (profile_detail::panel_derivs);
//   V        H^-1 S H^-1, the sandwich in u. se_scaled_a = sqrt(V_aa); se_a = (dv_a / du_a) se_scaled_a,
//            the delta method to natural scale, as calibrate does.
//
// The lag L: newey_west_lag(T) = floor(4 (T / 100)^(2/9)), 2, 3 and 4 at T = 20, 40 and 100; a caller
// may pass another (S-48 reports L = 1 and twice the rule). L = 0 is the outer-product sandwich
// H^-1 (sum_t c_t c_t') H^-1; lags of T or more contribute nothing.
//
// Not computed (flags = kIntervalNotComputed, every value NaN) when the fit carries kFlagGridEdge,
// kFlagFlatSurface or kFlagNearBound, the fits that get no iid Wald interval either (so both Wald
// intervals are computed on the same fits); when H is not negative definite; or when V's diagonal is not
// positive and finite. The Wald interval is symmetric in u: u_a +/- z se_scaled_a (hac_wald_interval).
//
// Needs an objective with log_contrib_derivs. Serial and deterministic.
#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "core/grid.hpp"
#include "engine/calibrate.hpp"
#include "engine/derivs.hpp"
#include "engine/profile.hpp"

namespace vcal::engine {

// The fits that get no HAC SE, as they get no iid Wald interval.
inline constexpr std::uint32_t kHacSkipFlags = kFlagGridEdge | kFlagFlatSurface | kFlagNearBound;

// The fixed lag rule of S-48 (D-209): floor(4 (T / 100)^(2/9)).
inline std::int32_t newey_west_lag(std::int64_t periods) {
    return static_cast<std::int32_t>(std::floor(4.0 * std::pow(static_cast<double>(periods) / 100.0, 2.0 / 9.0)));
}

// S from per-period scores (row t = (s_t0, s_t1), T rows), centred and Bartlett-weighted as above.
// out = (S_00, S_11, S_01). lag < 0 is treated as 0.
inline void newey_west(const double* scores, std::int64_t periods, std::int32_t lag, double (&out)[3]) {
    out[0] = out[1] = out[2] = 0.0;
    if (periods <= 0) return;
    double mean[2] = {0.0, 0.0};
    for (std::int64_t t = 0; t < periods; ++t) {
        mean[0] += scores[2 * t];
        mean[1] += scores[2 * t + 1];
    }
    mean[0] /= static_cast<double>(periods);
    mean[1] /= static_cast<double>(periods);
    std::vector<double> c(static_cast<std::size_t>(2 * periods));
    for (std::int64_t t = 0; t < periods; ++t) {
        c[static_cast<std::size_t>(2 * t)] = scores[2 * t] - mean[0];
        c[static_cast<std::size_t>(2 * t + 1)] = scores[2 * t + 1] - mean[1];
    }
    const auto at = [&](std::int64_t t, int a) { return c[static_cast<std::size_t>(2 * t + a)]; };
    for (std::int64_t t = 0; t < periods; ++t) {
        out[0] += at(t, 0) * at(t, 0);
        out[1] += at(t, 1) * at(t, 1);
        out[2] += at(t, 0) * at(t, 1);
    }
    for (std::int64_t l = 1; l <= lag && l < periods; ++l) {
        const double w = 1.0 - static_cast<double>(l) / static_cast<double>(lag + 1);
        double g[3] = {0.0, 0.0, 0.0};
        for (std::int64_t t = l; t < periods; ++t) {
            g[0] += 2.0 * at(t, 0) * at(t - l, 0);
            g[1] += 2.0 * at(t, 1) * at(t - l, 1);
            g[2] += at(t, 0) * at(t - l, 1) + at(t - l, 0) * at(t, 1);
        }
        for (int i = 0; i < 3; ++i) out[i] += w * g[i];
    }
}

struct HacSe2 {
    double se[2];         // natural scale (delta method), grid axis order
    double se_scaled[2];  // in the axes' scaled coordinates
    double corr;
    std::int32_t lag;
    std::uint32_t flags;  // kIntervalNotComputed when none (see above)
};

template <class Objective, class Integrator>
HacSe2 hac_se(const Objective& objective, const Integrator& primary, const typename Objective::Obs* obs,
              std::int64_t periods, const Grid<2>& grid, const Estimate2& est, std::int32_t lag) {
    static_assert(Objective::n_params == 2, "hac_se handles 2-parameter objectives");
    static_assert(profile_detail::has_contrib_derivs<Objective, Integrator>::value,
                  "hac_se needs the objective's analytic derivatives (log_contrib_derivs)");
    const double nan = std::nan("");
    HacSe2 out{{nan, nan}, {nan, nan}, nan, lag < 0 ? 0 : lag, kIntervalNotComputed};
    if ((est.flags & kHacSkipFlags) || periods <= 0) return out;

    const double u[2] = {grid::to_scaled(grid.axis[0].scale, est.value[0]),
                         grid::to_scaled(grid.axis[1].scale, est.value[1])};
    const double j[2] = {grid::dvalue_dscaled(grid.axis[0].scale, u[0]),
                         grid::dvalue_dscaled(grid.axis[1].scale, u[1])};

    std::vector<std::int64_t> first, count;
    profile_detail::distinct_observations(obs, periods, first, count);
    const Axis* const ax[2] = {&grid.axis[0], &grid.axis[1]};
    const auto d = profile_detail::panel_derivs(objective, primary, obs, first, count, ax, u[0], u[1]);
    const double h00 = d.h[0], h11 = d.h[1], h01 = d.h[2];
    const double det = h00 * h11 - h01 * h01;
    if (!(h00 < 0.0 && det > 0.0)) return out;

    // Per-period scores in u: each distinct observation's score once, then in period order.
    const auto theta = Objective::theta(est.value);
    std::vector<double> distinct(2 * first.size());
    for (std::size_t k = 0; k < first.size(); ++k) {
        double g[2], h[3];
        (void)objective.log_contrib_derivs(obs[first[k]], theta, primary, g, h);
        distinct[2 * k] = g[0] * j[0];
        distinct[2 * k + 1] = g[1] * j[1];
    }
    std::vector<double> scores(static_cast<std::size_t>(2 * periods));
    for (std::int64_t t = 0; t < periods; ++t) {
        std::size_t k = 0;
        while (!(obs[first[k]] == obs[t])) ++k;
        scores[static_cast<std::size_t>(2 * t)] = distinct[2 * k];
        scores[static_cast<std::size_t>(2 * t + 1)] = distinct[2 * k + 1];
    }
    double s[3];
    newey_west(scores.data(), periods, out.lag, s);

    // V = H^-1 S H^-1, H^-1 = (1/det) [[h11, -h01], [-h01, h00]].
    const double a00 = h11 / det, a11 = h00 / det, a01 = -h01 / det;
    const double v00 = a00 * (a00 * s[0] + a01 * s[2]) + a01 * (a00 * s[2] + a01 * s[1]);
    const double v11 = a01 * (a01 * s[0] + a11 * s[2]) + a11 * (a01 * s[2] + a11 * s[1]);
    const double v01 = a00 * (a01 * s[0] + a11 * s[2]) + a01 * (a01 * s[2] + a11 * s[1]);
    if (!(v00 > 0.0 && v11 > 0.0 && std::isfinite(v00) && std::isfinite(v11))) return out;
    out.se_scaled[0] = std::sqrt(v00);
    out.se_scaled[1] = std::sqrt(v11);
    out.se[0] = j[0] * out.se_scaled[0];
    out.se[1] = j[1] * out.se_scaled[1];
    out.corr = v01 / (out.se_scaled[0] * out.se_scaled[1]);
    out.flags = 0u;
    return out;
}

// The HAC Wald interval for axis a, symmetric in the scaled coordinate: from_scaled(u_a -/+ z se_scaled_a).
// NaN ends when the SE was not computed.
struct WaldInterval {
    double lo, hi;
};

inline WaldInterval hac_wald_interval(const Grid<2>& grid, const Estimate2& est, const HacSe2& hac, int a,
                                      double z) {
    if (hac.flags & kIntervalNotComputed) return {std::nan(""), std::nan("")};
    const AxisScale sc = grid.axis[a].scale;
    const double u = grid::to_scaled(sc, est.value[a]);
    return {grid::from_scaled(sc, u - z * hac.se_scaled[a]), grid::from_scaled(sc, u + z * hac.se_scaled[a])};
}

}  // namespace vcal::engine
