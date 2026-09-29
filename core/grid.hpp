// SPDX-License-Identifier: Apache-2.0
//
// Parameter grids (ARCHITECTURE.md §3.3). Each axis is uniform in a *scaled* coordinate
// u = S(v) and maps back to natural values v = S^-1(u):
//     Linear: u = v            Log:    u = log v
//     Logit:  u = log(v/(1-v)) Probit: u = Phi^-1(v)
// Refinement and curvature SEs are computed in u, where the points are equally spaced, and convert
// to natural scale with dv/du (the delta method; D-038, D-095).
//
// Axis point i has u_i = u_lo + i * step (i = 0 .. n-1), step = (u_hi - u_lo) / (n - 1); the
// endpoints are exactly lo and hi. Grids are row-major: axis 0 varies slowest.
#pragma once

#include <cmath>
#include <cstdint>

#include "core/precision.hpp"
#include "core/special/log_phi.hpp"
#include "core/special/probit.hpp"

namespace vcal {

enum class AxisScale : std::int32_t { Linear = 0, Log = 1, Logit = 2, Probit = 3 };

// Default upper bound of the rho axis (D-089): configurable per grid; estimates on it are
// flagged as grid-edge. Parity's precision claim covers rho <= 0.5 (D-092).
inline constexpr double kDefaultRhoUpper = 0.5;

namespace grid {

VCAL_HD double to_scaled(AxisScale s, double v) {
    switch (s) {
        case AxisScale::Log: return std::log(v);
        case AxisScale::Logit: return std::log(v) - std::log1p(-v);
        case AxisScale::Probit: return special::probit(v);
        case AxisScale::Linear: break;
    }
    return v;
}

VCAL_HD double from_scaled(AxisScale s, double u) {
    switch (s) {
        case AxisScale::Log: return std::exp(u);
        case AxisScale::Logit: {
            const double e = std::exp(-std::fabs(u));  // never overflows
            return u >= 0.0 ? 1.0 / (1.0 + e) : e / (1.0 + e);
        }
        case AxisScale::Probit: return std::exp(special::log_phi(u));
        case AxisScale::Linear: break;
    }
    return u;
}

// d2v/du2 at scaled coordinate u (the profile solves' chain rule, D-170).
VCAL_HD double d2value_dscaled2(AxisScale s, double u);

// dv/du at scaled coordinate u.
VCAL_HD double dvalue_dscaled(AxisScale s, double u) {
    switch (s) {
        case AxisScale::Log: return std::exp(u);
        case AxisScale::Logit: {
            const double v = from_scaled(s, u);
            return v * (1.0 - v);
        }
        case AxisScale::Probit: return std::exp(-0.5 * u * u) * special::constants::kInvSqrt2Pi;
        case AxisScale::Linear: break;
    }
    return 1.0;
}

VCAL_HD double d2value_dscaled2(AxisScale s, double u) {
    switch (s) {
        case AxisScale::Log: return std::exp(u);
        case AxisScale::Logit: {
            const double v = from_scaled(s, u);
            return v * (1.0 - v) * (1.0 - 2.0 * v);
        }
        case AxisScale::Probit: return -u * std::exp(-0.5 * u * u) * special::constants::kInvSqrt2Pi;
        case AxisScale::Linear: break;
    }
    return 0.0;
}

}  // namespace grid

struct Axis {
    double lo;
    double hi;
    std::int32_t n;
    AxisScale scale;

    VCAL_HD double scaled_lo() const { return grid::to_scaled(scale, lo); }
    VCAL_HD double scaled_hi() const { return grid::to_scaled(scale, hi); }
    VCAL_HD double step() const { return (scaled_hi() - scaled_lo()) / static_cast<double>(n - 1); }
    VCAL_HD double scaled_at(std::int32_t i) const {
        return i == n - 1 ? scaled_hi() : scaled_lo() + static_cast<double>(i) * step();
    }
    VCAL_HD double value_at(std::int32_t i) const {
        if (i == 0) return lo;
        if (i == n - 1) return hi;
        return grid::from_scaled(scale, scaled_at(i));
    }
};

// nullptr if the axis is usable; otherwise a static description of the problem.
inline const char* axis_error(const Axis& a) {
    if (a.n < 3) return "axis needs at least 3 points (refinement uses a 3-point stencil)";
    if (!(std::isfinite(a.lo) && std::isfinite(a.hi) && a.lo < a.hi)) return "axis bounds must be finite with lo < hi";
    switch (a.scale) {
        case AxisScale::Linear: return nullptr;
        case AxisScale::Log: return a.lo > 0.0 ? nullptr : "log axis needs lo > 0";
        case AxisScale::Logit:
        case AxisScale::Probit: return (a.lo > 0.0 && a.hi < 1.0) ? nullptr : "logit/probit axis needs 0 < lo < hi < 1";
    }
    return "unknown axis scale";
}

template <int D>
struct Grid {
    static_assert(D >= 1, "a grid needs at least one axis");
    Axis axis[D];

    VCAL_HD std::int64_t size() const {
        std::int64_t k = 1;
        for (int a = 0; a < D; ++a) k *= axis[a].n;
        return k;
    }
    VCAL_HD void unflatten(std::int64_t k, std::int32_t (&i)[D]) const {
        for (int a = D - 1; a >= 0; --a) {
            i[a] = static_cast<std::int32_t>(k % axis[a].n);
            k /= axis[a].n;
        }
    }
    VCAL_HD std::int64_t flatten(const std::int32_t (&i)[D]) const {
        std::int64_t k = 0;
        for (int a = 0; a < D; ++a) k = k * axis[a].n + i[a];
        return k;
    }
    VCAL_HD void values(std::int64_t k, double (&v)[D]) const {
        std::int32_t i[D];
        unflatten(k, i);
        for (int a = 0; a < D; ++a) v[a] = axis[a].value_at(i[a]);
    }
};

template <int D>
const char* grid_error(const Grid<D>& g) {
    for (int a = 0; a < D; ++a) {
        if (const char* e = axis_error(g.axis[a])) return e;
    }
    return nullptr;
}

}  // namespace vcal
