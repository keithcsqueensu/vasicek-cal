// SPDX-License-Identifier: Apache-2.0
//
// Composite Gauss-Legendre rule for one-sided integrands (D-093, D-116, D-118), and the parity
// SplitRule that uses it for zero- and all-default periods.
//
// A one-sided integrand is the N(0,1) prior truncated by a steep survival factor: no interior
// peak, so polynomial rules centred on the mode resolve the cliff slowly (D-116). The rule:
//   1. Domain [z_L, z_R]: where h(z) = g(z) - z^2/2 has fallen by kDrop (50) from its value at
//      the hint's mode, on each side. h is concave, so each side has one crossing, found by
//      doubling a step then exactly 64 bisections: the result is the root to rounding, a smooth
//      function of the parameters.
//   2. Map z = z_c + w sinh(u), with z_c and w the survival factor's transition point and width
//      from the hint. Points concentrate at the cliff and spread geometrically towards the bulk.
//   3. K panels of M Gauss-Legendre points, equally spaced in u over [asinh((z_L-z_c)/w),
//      asinh((z_R-z_c)/w)]. K and M are fixed, the same at every parameter point, so the rule is
//      a smooth function of the parameters with no panel count to jump (owner's condition on Q16).
//   log I = h(mode) + log sum_k,i W_i (du/2) w cosh(u) exp(h(z) - h(mode)) - log sqrt(2 pi).
// z_c and w are clamped into a range around the domain. That only matters when the cliff is far
// outside the mass, where the rule is already exact to rounding.
#pragma once

#include <cmath>

#include "core/precision.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "core/quadrature/generated/gauss_legendre.hpp"
#include "core/quadrature/integrator.hpp"
#include "core/special/generated/constants.hpp"

namespace vcal::quadrature {

// Gauss-Legendre rule on [-1, 1], non-owning; n == 0 means "no rule".
struct GaussLegendreRule {
    const double* node;
    const double* weight;
    int n;
};

// Host-side lookup of a generated rule; n outside the supported set gives {nullptr, nullptr, 0}.
inline GaussLegendreRule gauss_legendre_rule(int n) {
    GaussLegendreRule r{nullptr, nullptr, 0};
    if (generated::gauss_legendre_lookup(n, &r.node, &r.weight)) r.n = n;
    return r;
}

inline constexpr double kCompositeDrop = 50.0;
inline constexpr int kCompositeBisections = 64;

template <class P>
struct CompositeLegendre {
    static_assert(check_precision<P>());
    GaussLegendreRule rule;
    int panels;

    template <class G>
    VCAL_HD typename P::accum_t log_integrate(const G& g, IntegrandHint hint) const {
        detail::OnlineLogSumExp acc;
        return log_integrate_with(g, hint, acc);
    }
    template <class G, class Acc>
    VCAL_HD typename P::accum_t log_integrate_with(const G& g, IntegrandHint hint, Acc& acc) const {
        if (!(hint.scale > 0.0) || !std::isfinite(hint.scale) || !std::isfinite(hint.mode) || panels < 1) {
            return static_cast<double>(NAN);
        }
        const double zm = hint.mode;
        const double hm = g(zm) - 0.5 * zm * zm;
        if (!std::isfinite(hm)) return hm == -HUGE_VAL ? -HUGE_VAL : static_cast<double>(NAN);
        const double target = hm - kCompositeDrop;
        const auto below = [&](double z) { return g(z) - 0.5 * z * z < target; };
        const auto extent = [&](double sign) {
            double inner = 0.0;
            double outer = hint.scale;
            for (int i = 0; i < 64 && !below(zm + sign * outer); ++i) {
                inner = outer;
                outer *= 2.0;
            }
            for (int i = 0; i < kCompositeBisections; ++i) {
                const double mid = 0.5 * (inner + outer);
                if (below(zm + sign * mid)) {
                    outer = mid;
                } else {
                    inner = mid;
                }
            }
            return outer;
        };
        const double z_lo = zm - extent(-1.0);
        const double z_hi = zm + extent(1.0);
        const double span = z_hi - z_lo;

        double z_c = hint.center;
        double w = hint.width;
        if (!std::isfinite(z_c)) z_c = zm;
        if (!(w > 0.0) || !std::isfinite(w)) w = hint.scale;
        z_c = std::fmin(std::fmax(z_c, z_lo - span), z_hi + span);
        w = std::fmin(std::fmax(w, 1e-12 * span), span);

        const double u_lo = std::asinh((z_lo - z_c) / w);
        const double u_hi = std::asinh((z_hi - z_c) / w);
        const double du = (u_hi - u_lo) / static_cast<double>(panels);
        const double half = 0.5 * du;
        for (int k = 0; k < panels; ++k) {
            const double mid = u_lo + (static_cast<double>(k) + 0.5) * du;
            for (int i = 0; i < rule.n; ++i) {
                const double u = mid + half * rule.node[i];
                const double z = z_c + w * std::sinh(u);
                const double h = g(z) - 0.5 * z * z;
                acc.add(std::log(rule.weight[i] * half * w * std::cosh(u)) + (h - hm), z);
            }
        }
        return hm + acc.result() - special::constants::kLnSqrt2Pi;
    }
};

// Parity integrator (D-118): adaptive Gauss-Hermite for integrands with an interior peak,
// composite Gauss-Legendre for one-sided ones; the hint says which.
template <class P>
struct SplitRule {
    static_assert(check_precision<P>());
    GaussHermiteAdaptive<P> interior;
    CompositeLegendre<P> one_sided;

    template <class G>
    VCAL_HD typename P::accum_t log_integrate(const G& g, IntegrandHint hint) const {
        return hint.one_sided ? one_sided.log_integrate(g, hint) : interior.log_integrate(g, hint);
    }
    template <class G, class Acc>
    VCAL_HD typename P::accum_t log_integrate_with(const G& g, IntegrandHint hint, Acc& acc) const {
        return hint.one_sided ? one_sided.log_integrate_with(g, hint, acc) : interior.log_integrate_with(g, hint, acc);
    }
};

static_assert(check_integrator<CompositeLegendre<PrecisionF64>>());
static_assert(check_integrator<SplitRule<PrecisionF64>>());

}  // namespace vcal::quadrature
