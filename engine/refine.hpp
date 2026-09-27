// SPDX-License-Identifier: Apache-2.0
//
// Sub-grid refinement and curvature standard errors for a 2-parameter surface (D-038, D-095).
//
// Everything happens in the axes' scaled coordinates u, where grid points are equally spaced.
// With f the surface on the 3x3 stencil around the grid argmax (steps h0, h1):
//     g_a  = (f(+e_a) - f(-e_a)) / (2 h_a)
//     H_aa = (f(+e_a) - 2 f(0) + f(-e_a)) / h_a^2
//     H_01 = (f(+,+) - f(+,-) - f(-,+) + f(-,-)) / (4 h0 h1)
// i.e. the exact quadratic through the stencil. The cross term is kept because PD and rho
// estimates are strongly correlated; refining each axis alone would bias the vertex.
//   - argmax on an axis end  -> grid point, kFlagGridEdge, SEs NaN (no interior stencil).
//   - H not negative definite -> grid point, kFlagFlatSurface, SEs NaN. Never extrapolate.
//   - vertex -H^-1 g outside the stencil (|delta_a| > h_a) -> grid point,
//     kFlagRefinementRejected. SEs are still reported, since H is valid.
//   - otherwise the vertex. Observed information is -H; Cov_u = (-H)^-1, converted to
//     natural scale by the delta method: Cov_v = J Cov_u J, J = diag(dv/du at the estimate).
// Curvature from a stencil of the grid's own step is O(h^2) accurate; the M1.4 tests
// measure it against a fine-step Hessian.
#pragma once

#include <cmath>
#include <cstdint>

#include "core/grid.hpp"

namespace vcal::engine {

enum : std::uint32_t {
    kFlagGridEdge = 1u << 0,
    kFlagFlatSurface = 1u << 1,
    kFlagQuadratureUnconverged = 1u << 2,
    kFlagRefinementRejected = 1u << 3,
    kFlagNumeric = 1u << 4,  // NaN somewhere in the surface
    kFlagNearBound = 1u << 5,  // within kNearBoundSe standard errors of an axis bound (D-119)
    kFlagRhoNotIdentified = 1u << 6,  // the data carry no information about rho (D-302); set by calibrate
};

struct Refined2 {
    double value[2];  // natural scale
    double scaled[2];
    double se[2];     // natural scale; NaN when curvature is unavailable
    double se_scaled[2];  // the same, in the axes' scaled coordinates (step size for D-119)
    double corr;
    std::uint32_t flags;
};

// f(i0, i1) returns the (weighted) surface value at grid indices (i0, i1).
template <class F>
Refined2 refine_2d(const Grid<2>& grid, std::int64_t k_star, const F& f) {
    const double nan = std::nan("");
    std::int32_t i[2];
    grid.unflatten(k_star, i);
    Refined2 r{{grid.axis[0].value_at(i[0]), grid.axis[1].value_at(i[1])},
               {grid.axis[0].scaled_at(i[0]), grid.axis[1].scaled_at(i[1])},
               {nan, nan},
               {nan, nan},
               nan,
               0u};
    for (int a = 0; a < 2; ++a) {
        if (i[a] == 0 || i[a] == grid.axis[a].n - 1) r.flags |= kFlagGridEdge;
    }
    if (r.flags & kFlagGridEdge) return r;

    const double h0 = grid.axis[0].step();
    const double h1 = grid.axis[1].step();
    auto at = [&](int d0, int d1) { return f(i[0] + d0, i[1] + d1); };
    const double f00 = at(0, 0);
    const double g0 = (at(1, 0) - at(-1, 0)) / (2.0 * h0);
    const double g1 = (at(0, 1) - at(0, -1)) / (2.0 * h1);
    const double H00 = (at(1, 0) - 2.0 * f00 + at(-1, 0)) / (h0 * h0);
    const double H11 = (at(0, 1) - 2.0 * f00 + at(0, -1)) / (h1 * h1);
    const double H01 = (at(1, 1) - at(1, -1) - at(-1, 1) + at(-1, -1)) / (4.0 * h0 * h1);
    const double det = H00 * H11 - H01 * H01;
    if (!(H00 < 0.0 && det > 0.0)) {  // not negative definite (also catches NaN)
        r.flags |= kFlagFlatSurface;
        return r;
    }

    const double delta0 = -(H11 * g0 - H01 * g1) / det;
    const double delta1 = -(H00 * g1 - H01 * g0) / det;
    if (std::fabs(delta0) <= h0 && std::fabs(delta1) <= h1) {
        r.scaled[0] += delta0;
        r.scaled[1] += delta1;
        r.value[0] = grid::from_scaled(grid.axis[0].scale, r.scaled[0]);
        r.value[1] = grid::from_scaled(grid.axis[1].scale, r.scaled[1]);
    } else {
        r.flags |= kFlagRefinementRejected;
    }

    // Cov_u = (-H)^-1 = [[-H11, H01], [H01, -H00]] / det
    const double j0 = grid::dvalue_dscaled(grid.axis[0].scale, r.scaled[0]);
    const double j1 = grid::dvalue_dscaled(grid.axis[1].scale, r.scaled[1]);
    const double var0 = j0 * j0 * (-H11 / det);
    const double var1 = j1 * j1 * (-H00 / det);
    const double cov01 = j0 * j1 * (H01 / det);
    r.se[0] = std::sqrt(var0);
    r.se[1] = std::sqrt(var1);
    r.se_scaled[0] = std::sqrt(-H11 / det);
    r.se_scaled[1] = std::sqrt(-H00 / det);
    r.corr = cov01 / (r.se[0] * r.se[1]);
    return r;
}

}  // namespace vcal::engine
