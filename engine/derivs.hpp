// SPDX-License-Identifier: Apache-2.0
//
// The panel log-likelihood's analytic gradient and Hessian (D-170), shared by the profile solver
// (engine/profile.hpp) and calibrate's observed information (P-15, D-185). Moved here from
// profile.hpp unchanged, so that calibrate, which profile.hpp builds on, can use them.
#pragma once

#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>

#include "core/grid.hpp"

namespace vcal::engine::profile_detail {

// Distinct observations and their multiplicities, in order of first appearance.
template <class Obs>
void distinct_observations(const Obs* obs, std::int64_t periods, std::vector<std::int64_t>& first,
                           std::vector<std::int64_t>& count) {
    for (std::int64_t t = 0; t < periods; ++t) {
        std::size_t j = 0;
        while (j < first.size() && !(obs[first[j]] == obs[t])) ++j;
        if (j == first.size()) {
            first.push_back(t);
            count.push_back(1);
        } else {
            ++count[j];
        }
    }
}

// Does the objective give its contribution's score and Hessian (log_contrib_derivs, D-170) with
// this integrator?
template <class O, class I, class = void>
struct has_contrib_derivs : std::false_type {};
template <class O, class I>
struct has_contrib_derivs<
    O, I,
    std::void_t<decltype(std::declval<const O&>().log_contrib_derivs(
        std::declval<const typename O::Obs&>(), std::declval<const typename O::Theta&>(), std::declval<const I&>(),
        std::declval<double (&)[2]>(), std::declval<double (&)[3]>()))>> : std::true_type {};

// The panel log-likelihood with its gradient and Hessian in the axes' scaled coordinates.
// h = (h_00, h_11, h_01). l equals panel_loglik's bit for bit (same terms, same compensated sum).
struct PanelDerivs {
    double l;
    double g[2];
    double h[3];
};

// In the natural parameters v (grid axis order).
template <class Objective, class Integrator>
PanelDerivs panel_derivs_natural(const Objective& objective, const Integrator& primary,
                                 const typename Objective::Obs* obs, const std::vector<std::int64_t>& first,
                                 const std::vector<std::int64_t>& count, const double (&v)[2]) {
    const auto th = Objective::theta(v);
    double sum = 0.0, comp = 0.0;
    double gv[2] = {0.0, 0.0}, hv[3] = {0.0, 0.0, 0.0};
    for (std::size_t j = 0; j < first.size(); ++j) {
        double gj[2], hj[3];
        const double m = static_cast<double>(count[j]);
        const double x = m * objective.log_contrib_derivs(obs[first[j]], th, primary, gj, hj);
        const double s = sum + x;
        comp += std::fabs(sum) >= std::fabs(x) ? (sum - s) + x : (x - s) + sum;
        sum = s;
        for (int i = 0; i < 2; ++i) gv[i] += m * gj[i];
        for (int i = 0; i < 3; ++i) hv[i] += m * hj[i];
    }
    return {sum + comp, {gv[0], gv[1]}, {hv[0], hv[1], hv[2]}};
}

// In the axes' scaled coordinates (u0, u1): the natural derivatives by the chain rule, v = v(u)
// per axis.
template <class Objective, class Integrator>
PanelDerivs panel_derivs(const Objective& objective, const Integrator& primary, const typename Objective::Obs* obs,
                         const std::vector<std::int64_t>& first, const std::vector<std::int64_t>& count,
                         const Axis* const (&ax)[2], double u0, double u1) {
    const double u[2] = {u0, u1};
    const double v[2] = {grid::from_scaled(ax[0]->scale, u0), grid::from_scaled(ax[1]->scale, u1)};
    const PanelDerivs n = panel_derivs_natural(objective, primary, obs, first, count, v);
    double j1[2], j2[2];
    for (int a = 0; a < 2; ++a) {
        j1[a] = grid::dvalue_dscaled(ax[a]->scale, u[a]);
        j2[a] = grid::d2value_dscaled2(ax[a]->scale, u[a]);
    }
    PanelDerivs out{};
    out.l = n.l;
    out.g[0] = n.g[0] * j1[0];
    out.g[1] = n.g[1] * j1[1];
    out.h[0] = n.h[0] * j1[0] * j1[0] + n.g[0] * j2[0];
    out.h[1] = n.h[1] * j1[1] * j1[1] + n.g[1] * j2[1];
    out.h[2] = n.h[2] * j1[0] * j1[1];
    return out;
}

}  // namespace vcal::engine::profile_detail
