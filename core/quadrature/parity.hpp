// SPDX-License-Identifier: Apache-2.0
//
// The parity integration rule, in one place (D-035, D-088, D-118):
//   periods with 0 < d < n : adaptive Gauss-Hermite, N = 128;
//   periods with d in {0, n}: composite Gauss-Legendre, K = 16 panels of M = 16 points.
// K and M were chosen as the smallest pair at the rounding floor on every golden one-sided
// case up to rho = 0.99 (the composite convergence table in docs/methodology/tolerances.md),
// and they are fixed at every parameter point, so the surface is smooth in the parameters.
// The per-run convergence check (D-092) evaluates the same rule with twice the Gauss-Hermite
// nodes and twice the panels.
#pragma once

#include "core/precision.hpp"
#include "core/quadrature/composite_legendre.hpp"
#include "core/quadrature/gauss_hermite.hpp"

namespace vcal::quadrature {

inline constexpr int kParityGaussHermiteNodes = 128;
inline constexpr int kParityCompositePanels = 16;
inline constexpr int kParityCompositePoints = 16;

// Host-side: the rules point into generated host tables.
inline SplitRule<PrecisionF64> parity_rule(bool doubled = false) {
    const int f = doubled ? 2 : 1;
    return {GaussHermiteAdaptive<PrecisionF64>{gauss_hermite_rule(f * kParityGaussHermiteNodes)},
            CompositeLegendre<PrecisionF64>{gauss_legendre_rule(kParityCompositePoints), f * kParityCompositePanels}};
}

}  // namespace vcal::quadrature
