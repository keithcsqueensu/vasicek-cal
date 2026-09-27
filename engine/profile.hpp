// SPDX-License-Identifier: Apache-2.0
//
// Profile-likelihood intervals for a 2-parameter objective (M2, D-128..D-130).
//
// For parameter a with nuisance b, the profile is P_a(u) = max over u_b of l(u, u_b), in the
// axes' scaled (logit) coordinates, where l is the panel log-likelihood. The 95% interval is
// {u : P_a(u) >= l_max - c}, c = chi2_1(0.95) / 2 = 1.9207294103470630 (D-129: chi2_1 with
// truncation, also near rho = 0).
//
//   1. l_max: the maximum polished off-grid. The profile of PD is maximised by Brent's method
//      around the refined estimate, each value being an inner maximisation over rho.
//   2. Inner maximum: Brent's method over u_b, bracketed by the surface's own argmax over b in
//      the grid columns either side of u, +-2 steps. If the maximum lands on a bracket edge that
//      is not a bound of the box, the bracket widens by 2 steps and the search repeats.
//   3. Endpoint, each side: walk the grid outwards from the maximum. The grid profile
//      (max over grid points of the surface L) never exceeds the true profile, so a grid point
//      whose grid profile clears the threshold is inside the interval without further work. At
//      the first grid point where it does not, the true profile is evaluated. If that is below
//      the threshold the crossing is bracketed, and Brent's root finder solves
//      P_a(u) = l_max - c to kProfileRootTol in u. If no grid point falls below the threshold
//      before a bound of the box, the endpoint is that bound, flagged truncated: never
//      extrapolated (D-130).
//   Accuracy: each endpoint's residual |P_a(u) - (l_max - c)| is reported; the tests bound it
//   (TOL_PROFILE_ENDPOINT_RESIDUAL_LL) against an independent dense evaluation.
//
// Serial and deterministic: results do not depend on threads. Observations that repeat are
// evaluated once (as in the surface, D-122). Headline fits only; resampling does not call it.
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "core/grid.hpp"
#include "engine/calibrate.hpp"
#include "engine/surface.hpp"

namespace vcal::engine {

// chi2_1 0.95 quantile / 2 (mpmath, 20 digits: 1.9207294103470629792).
inline constexpr double kProfileThreshold95 = 1.9207294103470630;
// Brent tolerances, in scaled (logit) units: the inner maximisation and the root.
inline constexpr double kProfileMaxTol = 1e-9;
inline constexpr double kProfileRootTol = 1e-9;

enum : std::uint32_t {
    kIntervalLowerTruncated = 1u << 0,  // the profile never fell below the threshold: lower end = box bound
    kIntervalUpperTruncated = 1u << 1,
    kIntervalNotComputed = 1u << 2,  // the estimate carries flat / numeric flags, or the maximum is not finite
};

struct ProfileIntervals2 {
    double lo[2];  // natural scale, grid axis order
    double hi[2];
    std::uint32_t flags[2];
    double loglik_max;      // the polished maximum l_max
    double max_at[2];       // where it is attained (natural scale)
    double residual_max;    // max over solved (untruncated) endpoints of |P_a(u) - (l_max - c)|
    std::int64_t evaluations;  // panel log-likelihood evaluations used
};

namespace profile_detail {

struct Min1 {
    double x;
    double f;
};

// Brent's minimiser (Brent 1973, "localmin") of f on [a, b], absolute tolerance tol.
template <class F>
Min1 brent_min(const F& f, double a, double b, double tol) {
    const double c = 0.5 * (3.0 - std::sqrt(5.0));
    double x = a + c * (b - a), w = x, v = x;
    double fx = f(x), fw = fx, fv = fx;
    double d = 0.0, e = 0.0;
    for (int iter = 0; iter < 200; ++iter) {
        const double m = 0.5 * (a + b);
        const double t1 = 1e-12 * std::fabs(x) + tol / 3.0;
        const double t2 = 2.0 * t1;
        if (std::fabs(x - m) <= t2 - 0.5 * (b - a)) break;
        double p = 0.0, q = 0.0, r = 0.0;
        if (std::fabs(e) > t1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p; else q = -q;
            r = e;
            e = d;
        }
        if (std::fabs(p) < std::fabs(0.5 * q * r) && p > q * (a - x) && p < q * (b - x)) {
            d = p / q;  // parabolic step
            const double u = x + d;
            if (u - a < t2 || b - u < t2) d = x < m ? t1 : -t1;
        } else {
            e = (x < m ? b : a) - x;  // golden-section step
            d = c * e;
        }
        const double u = x + (std::fabs(d) >= t1 ? d : (d > 0.0 ? t1 : -t1));
        const double fu = f(u);
        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w; fv = fw; w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u; fv = fu;
            }
        }
    }
    return {x, fx};
}

// Brent's root finder (Brent 1973, "zero") on [a, b] with g(a) and g(b) of opposite sign.
template <class G>
double brent_root(const G& g, double a, double b, double ga, double gb, double tol) {
    double c = a, gc = ga, d = b - a, e = d;
    for (int iter = 0; iter < 200; ++iter) {
        if ((gb > 0.0) == (gc > 0.0)) { c = a; gc = ga; d = b - a; e = d; }
        if (std::fabs(gc) < std::fabs(gb)) { a = b; b = c; c = a; ga = gb; gb = gc; gc = ga; }
        const double t = 2e-16 * std::fabs(b) + 0.5 * tol;
        const double m = 0.5 * (c - b);
        if (std::fabs(m) <= t || gb == 0.0) break;
        if (std::fabs(e) >= t && std::fabs(ga) > std::fabs(gb)) {
            const double s = gb / ga;
            double p, q;
            if (a == c) {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                const double qq = ga / gc, r = gb / gc;
                p = s * (2.0 * m * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q; else p = -p;
            if (2.0 * p < 3.0 * m * q - std::fabs(t * q) && p < std::fabs(0.5 * e * q)) {
                e = d;
                d = p / q;
            } else {
                d = m;
                e = d;
            }
        } else {
            d = m;
            e = d;
        }
        a = b;
        ga = gb;
        b += std::fabs(d) > t ? d : (m > 0.0 ? t : -t);
        gb = g(b);
    }
    return b;
}

}  // namespace profile_detail

template <class Objective, class Integrator>
ProfileIntervals2 profile_intervals(const Objective& objective, const Integrator& primary,
                                    const typename Objective::Obs* obs, std::int64_t periods, const Grid<2>& grid,
                                    const std::vector<double>& L, const Estimate2& est,
                                    double threshold = kProfileThreshold95) {
    using profile_detail::brent_min;
    using profile_detail::brent_root;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    ProfileIntervals2 out{{nan, nan}, {nan, nan}, {0u, 0u}, nan, {nan, nan}, 0.0, 0};
    if ((est.flags & (kFlagFlatSurface | kFlagNumeric)) || !std::isfinite(est.value[0]) || !std::isfinite(est.value[1])) {
        out.flags[0] = out.flags[1] = kIntervalNotComputed;
        return out;
    }

    // Distinct observations and their multiplicities, in order of first appearance.
    std::vector<std::int64_t> first, count;
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
    const Axis* ax[2] = {&grid.axis[0], &grid.axis[1]};
    // Panel log-likelihood at scaled coordinates (u0, u1), compensated sum in a fixed order.
    const auto loglik = [&](double u0, double u1) {
        ++out.evaluations;
        const double v[2] = {grid::from_scaled(ax[0]->scale, u0), grid::from_scaled(ax[1]->scale, u1)};
        const auto th = Objective::theta(v);
        double sum = 0.0, comp = 0.0;
        for (std::size_t j = 0; j < first.size(); ++j) {
            const double x = static_cast<double>(count[j]) * objective.log_contrib(obs[first[j]], th, primary);
            const double s = sum + x;
            comp += std::fabs(sum) >= std::fabs(x) ? (sum - s) + x : (x - s) + sum;
            sum = s;
        }
        return sum + comp;
    };

    // The surface summed over periods, S[k], and each axis's grid profile.
    const std::int64_t K = grid.size();
    const std::vector<double> ones(static_cast<std::size_t>(periods), 1.0);
    std::vector<double> S(static_cast<std::size_t>(K));
    for (std::int64_t k = 0; k < K; ++k) S[static_cast<std::size_t>(k)] = weighted_sum(L.data(), periods, K, ones.data(), k);
    const auto surface = [&](int a, std::int32_t i, std::int32_t j) {  // i along axis a, j along the other
        const std::int32_t idx[2] = {a == 0 ? i : j, a == 0 ? j : i};
        return S[static_cast<std::size_t>(grid.flatten(idx))];
    };
    const auto grid_argmax_b = [&](int a, std::int32_t i) {
        const int b = 1 - a;
        std::int32_t best = 0;
        for (std::int32_t j = 1; j < ax[b]->n; ++j) {
            if (surface(a, i, j) > surface(a, i, best)) best = j;
        }
        return best;
    };
    const auto grid_profile = [&](int a, std::int32_t i) { return surface(a, i, grid_argmax_b(a, i)); };
    const auto clamp_index = [&](int a, double pos) {
        const double lim = static_cast<double>(ax[a]->n - 1);
        return static_cast<std::int32_t>(pos < 0.0 ? 0.0 : pos > lim ? lim : pos);
    };

    // P_a(u): maximise over the nuisance u_b; returns (argmax u_b, value).
    const auto inner = [&](int a, double u) {
        const int b = 1 - a;
        const double pos = (u - ax[a]->scaled_lo()) / ax[a]->step();
        const std::int32_t j0 = grid_argmax_b(a, clamp_index(a, std::floor(pos)));
        const std::int32_t j1 = grid_argmax_b(a, clamp_index(a, std::ceil(pos)));
        const double hb = ax[b]->step(), blo = ax[b]->scaled_lo(), bhi = ax[b]->scaled_hi();
        double lo = std::fmax(blo, ax[b]->scaled_at(j0 < j1 ? j0 : j1) - 2.0 * hb);
        double hi = std::fmin(bhi, ax[b]->scaled_at(j0 < j1 ? j1 : j0) + 2.0 * hb);
        const auto neg = [&](double w) { return a == 0 ? -loglik(u, w) : -loglik(w, u); };
        profile_detail::Min1 m{};
        for (int widen = 0; widen < 64; ++widen) {
            m = brent_min(neg, lo, hi, kProfileMaxTol);
            const double edge_tol = 1e3 * kProfileMaxTol;
            const bool at_lo = m.x - lo <= edge_tol && lo > blo;
            const bool at_hi = hi - m.x <= edge_tol && hi < bhi;
            if (!at_lo && !at_hi) break;
            if (at_lo) lo = std::fmax(blo, lo - 2.0 * hb);
            if (at_hi) hi = std::fmin(bhi, hi + 2.0 * hb);
        }
        return profile_detail::Min1{m.x, -m.f};
    };

    // 1. The polished maximum: maximise the PD profile around the estimate.
    const double u_hat[2] = {grid::to_scaled(ax[0]->scale, est.value[0]), grid::to_scaled(ax[1]->scale, est.value[1])};
    double lo0 = std::fmax(ax[0]->scaled_lo(), u_hat[0] - 2.0 * ax[0]->step());
    double hi0 = std::fmin(ax[0]->scaled_hi(), u_hat[0] + 2.0 * ax[0]->step());
    profile_detail::Min1 top{};
    for (int widen = 0; widen < 64; ++widen) {
        top = brent_min([&](double u) { return -inner(0, u).f; }, lo0, hi0, kProfileMaxTol);
        const bool at_lo = top.x - lo0 <= 1e3 * kProfileMaxTol && lo0 > ax[0]->scaled_lo();
        const bool at_hi = hi0 - top.x <= 1e3 * kProfileMaxTol && hi0 < ax[0]->scaled_hi();
        if (!at_lo && !at_hi) break;
        if (at_lo) lo0 = std::fmax(ax[0]->scaled_lo(), lo0 - 2.0 * ax[0]->step());
        if (at_hi) hi0 = std::fmin(ax[0]->scaled_hi(), hi0 + 2.0 * ax[0]->step());
    }
    const double u_max[2] = {top.x, inner(0, top.x).x};
    const double l_max = std::fmax(-top.f, est.loglik);
    if (!std::isfinite(l_max)) {
        out.flags[0] = out.flags[1] = kIntervalNotComputed;
        return out;
    }
    out.loglik_max = l_max;
    out.max_at[0] = grid::from_scaled(ax[0]->scale, u_max[0]);
    out.max_at[1] = grid::from_scaled(ax[1]->scale, u_max[1]);
    const double level = l_max - threshold;

    // 2. Each parameter, each side.
    for (int a = 0; a < 2; ++a) {
        const auto g = [&](double u) { return inner(a, u).f - level; };
        for (const int side : {-1, +1}) {
            const double bound = side < 0 ? ax[a]->scaled_lo() : ax[a]->scaled_hi();
            // Grid points strictly beyond the maximum, in order outwards.
            const double pos = (u_max[a] - ax[a]->scaled_lo()) / ax[a]->step();
            std::int32_t i = side < 0 ? static_cast<std::int32_t>(std::ceil(pos)) - 1
                                      : static_cast<std::int32_t>(std::floor(pos)) + 1;
            double inside = u_max[a];
            bool inside_exact = false;  // is g_inside g(inside) itself, or only a lower bound on it?
            double g_inside = 0.0;
            double endpoint = bound;
            bool truncated = true;
            bool failed = false;
            for (; i >= 0 && i < ax[a]->n; i += side) {
                const double u = ax[a]->scaled_at(i);
                if (grid_profile(a, i) >= level) {  // a lower bound on P_a(u): inside, no work needed
                    inside = u;
                    inside_exact = false;
                    continue;
                }
                const double g_u = g(u);
                if (!std::isfinite(g_u)) {
                    failed = true;
                    break;
                }
                if (g_u >= 0.0) {
                    inside = u;
                    g_inside = g_u;
                    inside_exact = true;
                    continue;
                }
                if (!inside_exact) g_inside = g(inside);  // Brent needs the value, not a bound
                if (!(g_inside >= 0.0)) {
                    failed = true;
                    break;
                }
                endpoint = brent_root(g, inside, u, g_inside, g_u, kProfileRootTol);
                truncated = false;
                break;
            }
            if (failed) {
                out.flags[a] |= kIntervalNotComputed;
                continue;
            }
            if (!truncated) out.residual_max = std::fmax(out.residual_max, std::fabs(g(endpoint)));
            const double v = truncated ? (side < 0 ? ax[a]->lo : ax[a]->hi) : grid::from_scaled(ax[a]->scale, endpoint);
            (side < 0 ? out.lo : out.hi)[a] = v;
            if (truncated) out.flags[a] |= side < 0 ? kIntervalLowerTruncated : kIntervalUpperTruncated;
        }
    }
    return out;
}

}  // namespace vcal::engine
