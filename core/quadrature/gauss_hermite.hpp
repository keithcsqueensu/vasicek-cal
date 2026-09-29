// SPDX-License-Identifier: Apache-2.0
//
// Gauss-Hermite integrators (ARCHITECTURE.md §3.4) over the mpmath-generated standard-normal
// rules in generated/gauss_hermite.hpp (D-080), for N in {8, 16, 24, 32, 48, 64, 96, 128}.
//
// GaussHermiteFixed:    log sum_i exp(log w_i + g(z_i)).
// GaussHermiteAdaptive: mode-centred, scaled rule (Liu & Pierce 1994). With z = mu + sigma u,
//     integral exp(g(z)) phi(z) dz = sigma * integral exp(g(mu + sigma u) + log phi(z) - log phi(u)) phi(u) du,
//   so log I ~ log sigma + log sum_i exp(log w_i + (u_i - z_i)(u_i + z_i)/2 + g(z_i)).
//   The exponent correction is written as a product so u^2 and z^2 never cancel. It is exact
//   (to rounding) whenever g(z) - z^2/2 is quadratic and the hint is its mode and scale;
//   fixed-node GH is not, and fails for the sharp integrands of large n (D-035, D-037).
//   A non-finite or non-positive hint returns NaN rather than guessing.
//
// Terms are folded in node order, so results do not depend on scheduling.
#pragma once

#include <cmath>

#include "core/precision.hpp"
#include "core/quadrature/generated/gauss_hermite.hpp"
#include "core/quadrature/integrator.hpp"

namespace vcal::quadrature {

template <class P>
struct GaussHermiteFixed {
    static_assert(check_precision<P>());
    GaussHermiteRule rule;

    template <class G>
    VCAL_HD typename P::accum_t log_integrate(const G& g, IntegrandHint hint) const {
        detail::OnlineLogSumExp acc;
        return log_integrate_with(g, hint, acc);
    }
    template <class G, class Acc>
    VCAL_HD typename P::accum_t log_integrate_with(const G& g, IntegrandHint /*unused*/, Acc& acc) const {
        for (int i = 0; i < rule.n; ++i) acc.add(rule.log_weight[i] + g(rule.node[i]), rule.node[i]);
        return acc.result();
    }
};

template <class P>
struct GaussHermiteAdaptive {
    static_assert(check_precision<P>());
    GaussHermiteRule rule;

    template <class G>
    VCAL_HD typename P::accum_t log_integrate(const G& g, IntegrandHint hint) const {
        detail::OnlineLogSumExp acc;
        return log_integrate_with(g, hint, acc);
    }
    template <class G, class Acc>
    VCAL_HD typename P::accum_t log_integrate_with(const G& g, IntegrandHint hint, Acc& acc) const {
        if (!(hint.scale > 0.0) || !std::isfinite(hint.scale) || !std::isfinite(hint.mode)) {
            return static_cast<double>(NAN);
        }
        for (int i = 0; i < rule.n; ++i) {
            const double u = rule.node[i];
            const double z = hint.mode + hint.scale * u;
            acc.add(rule.log_weight[i] + 0.5 * (u - z) * (u + z) + g(z), z);
        }
        return std::log(hint.scale) + acc.result();
    }
};

static_assert(check_integrator<GaussHermiteFixed<PrecisionF64>>());
static_assert(check_integrator<GaussHermiteAdaptive<PrecisionF64>>());

// Host-side lookup of a generated rule; n outside the supported set gives {nullptr, nullptr, 0}.
inline GaussHermiteRule gauss_hermite_rule(int n) {
    GaussHermiteRule r{nullptr, nullptr, 0};
    if (generated::gauss_hermite_lookup(n, &r.node, &r.log_weight)) r.n = n;
    return r;
}

}  // namespace vcal::quadrature
