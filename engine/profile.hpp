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
// Solvers (D-170). The intervals are defined by what is computed, not by how: l_max is the
// maximum of the panel log-likelihood, and an endpoint is a point where P_a(u) = l_max - c to
// within the stated residual. Any solver meeting that is valid. Two are provided, with the same
// grid brackets and the same tolerances in u:
//   - profile_intervals_newton: safeguarded Newton, for objectives with analytic derivatives
//     (log_contrib_derivs). The score and Hessian come from posterior moments on the quadrature's
//     own nodes (the value is log_contrib's bit for bit); the outer derivatives from the envelope
//     theorem at the inner maximum. A Newton step that leaves the bracket, or is not an ascent
//     step, is replaced by bisection. About 20x fewer panel evaluations than nested Brent.
//   - profile_intervals_brent: nested Brent (steps 1-3 above), for objectives without them.
// profile_intervals picks Newton when it is available. The two agree to
// TOL_PROFILE_SOLVER_AGREEMENT_U (unit_profile: profile_newton_matches_nested_brent).
//
// Serial and deterministic: results do not depend on threads. Observations that repeat are
// evaluated once (as in the surface, D-122). Headline fits only; resampling does not call it.
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
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

// Panel log-likelihood at natural values v, compensated sum over the distinct observations in a
// fixed order.
template <class Objective, class Integrator>
double panel_loglik(const Objective& objective, const Integrator& primary, const typename Objective::Obs* obs,
                    const std::vector<std::int64_t>& first, const std::vector<std::int64_t>& count, const double (&v)[2]) {
    const auto th = Objective::theta(v);
    double sum = 0.0, comp = 0.0;
    for (std::size_t j = 0; j < first.size(); ++j) {
        const double x = static_cast<double>(count[j]) * objective.log_contrib(obs[first[j]], th, primary);
        const double s = sum + x;
        comp += std::fabs(sum) >= std::fabs(x) ? (sum - s) + x : (x - s) + sum;
        sum = s;
    }
    return sum + comp;
}

// The surface summed over periods, S[k].
inline std::vector<double> summed_surface(const Grid<2>& grid, const std::vector<double>& L, std::int64_t periods) {
    const std::int64_t K = grid.size();
    const std::vector<double> ones(static_cast<std::size_t>(periods), 1.0);
    std::vector<double> S(static_cast<std::size_t>(K));
    for (std::int64_t k = 0; k < K; ++k) S[static_cast<std::size_t>(k)] = weighted_sum(L.data(), periods, K, ones.data(), k);
    return S;
}

}  // namespace profile_detail

// Nested Brent (steps 1-3 above): for objectives without analytic derivatives.
template <class Objective, class Integrator>
ProfileIntervals2 profile_intervals_brent(const Objective& objective, const Integrator& primary,
                                          const typename Objective::Obs* obs, std::int64_t periods,
                                          const Grid<2>& grid, const std::vector<double>& L, const Estimate2& est,
                                          double threshold = kProfileThreshold95) {
    using profile_detail::brent_min;
    using profile_detail::brent_root;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    ProfileIntervals2 out{{nan, nan}, {nan, nan}, {0u, 0u}, nan, {nan, nan}, 0.0, 0};
    if ((est.flags & (kFlagFlatSurface | kFlagNumeric)) || !std::isfinite(est.value[0]) || !std::isfinite(est.value[1])) {
        out.flags[0] = out.flags[1] = kIntervalNotComputed;
        return out;
    }

    std::vector<std::int64_t> first, count;
    profile_detail::distinct_observations(obs, periods, first, count);
    const Axis* ax[2] = {&grid.axis[0], &grid.axis[1]};
    // Panel log-likelihood at scaled coordinates (u0, u1).
    const auto loglik = [&](double u0, double u1) {
        ++out.evaluations;
        const double v[2] = {grid::from_scaled(ax[0]->scale, u0), grid::from_scaled(ax[1]->scale, u1)};
        return profile_detail::panel_loglik(objective, primary, obs, first, count, v);
    };

    // The surface summed over periods, S[k], and each axis's grid profile.
    const std::vector<double> S = profile_detail::summed_surface(grid, L, periods);
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


namespace profile_detail {

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

// A point of a 1-D solve: position, value, first and second derivative, and for a profile the
// inner maximiser w there and its slope dw/du.
struct Point1 {
    double x;
    double f;
    double d1;
    double d2;
    double w;
    double dwdu;
};

// Safeguarded Newton maximisation of f on [lo, hi] from x0 (D-170). Every evaluation shrinks the
// bracket on its downhill side (d1 > 0: the maximum is to the right). The next point is the
// Newton step x - d1/d2 when f is concave there (d2 < 0) and the step lands strictly inside the
// bracket; otherwise the bracket's midpoint. Stops when the Newton step is at most tol, returning
// that evaluated point (its value is within d1^2 / (2|d2|) <= tol |d1| / 2 of the maximum), or
// when the bracket is narrower than tol, returning the best point evaluated, or the end of the
// original bracket it collapsed onto (evaluated exactly: a maximum on a bound of the box is the
// bound itself). A non-finite value is returned as is: the caller treats it as a failure.
template <class F>
Point1 newton_max(const F& f, double lo, double hi, double x0, double tol) {
    const double lo0 = lo, hi0 = hi;
    double x = std::fmin(std::fmax(x0, lo), hi);
    Point1 best{x, -HUGE_VAL, 0.0, 0.0, 0.0, 0.0};
    for (int iter = 0; iter < 200; ++iter) {
        const Point1 e = f(x);
        if (!std::isfinite(e.f) || !std::isfinite(e.d1)) return e;
        if (e.f > best.f) best = e;
        if (e.d1 > 0.0) {
            lo = x;
        } else if (e.d1 < 0.0) {
            hi = x;
        } else {
            return e;
        }
        double next = std::numeric_limits<double>::quiet_NaN();
        if (e.d2 < 0.0) {
            const double s = -e.d1 / e.d2;
            if (std::fabs(s) <= tol) return e;
            next = x + s;
        }
        if (hi - lo <= tol) {
            const double end = hi == hi0 ? hi0 : lo == lo0 ? lo0 : std::numeric_limits<double>::quiet_NaN();
            if (end == best.x || !std::isfinite(end)) return best;
            const Point1 at_end = f(end);
            return at_end.f > best.f ? at_end : best;
        }
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        x = next;
    }
    return best;
}

// Safeguarded Newton root of g on the bracket between `in` (g >= 0, known without evaluating it)
// and `out`, starting from the evaluated point at `out` (g < 0) (D-170). Convergence is by
// bracketing, never by the step alone, since a slope can be wrong (at an end held by a moving bound
// of the box, the envelope slope misses the bound's motion):
//   - a Newton step that leaves the bracket is replaced by its midpoint;
//   - a Newton step that does not halve |g| is followed by a bisection;
//   - a step shorter than tol / 2 is lengthened to tol / 2 towards the root, so its evaluation
//     confirms the root by a sign change;
//   - each evaluation moves the end with its sign, and the search stops when the bracket is no
//     wider than tol, returning the evaluated point with the smallest |g| (the residual reported).
// A non-finite value is returned as is.
template <class G>
Point1 newton_root(const G& g, double in, double out, Point1 at, double tol) {
    Point1 best = at;
    bool bisect = false;
    for (int iter = 0; iter < 200; ++iter) {
        const double a = std::fmin(in, out), b = std::fmax(in, out);
        if (b - a <= tol) break;
        double next = std::numeric_limits<double>::quiet_NaN();
        if (!bisect && at.d1 != 0.0 && std::isfinite(at.d1)) {
            double s = -at.f / at.d1;
            if (std::fabs(s) < 0.5 * tol) s = std::copysign(0.5 * tol, s);
            next = at.x + s;
        }
        if (!(next > a && next < b)) next = 0.5 * (in + out);
        const double before = std::fabs(at.f);
        at = g(next);
        if (!std::isfinite(at.f)) return at;
        if (std::fabs(at.f) < std::fabs(best.f)) best = at;
        if (at.f >= 0.0) {
            in = next;
        } else {
            out = next;
        }
        bisect = !(std::fabs(at.f) <= 0.5 * before);
    }
    return best;
}

}  // namespace profile_detail

// profile_intervals with analytic derivatives (D-170): the same intervals, defined by what is
// computed, solved by safeguarded Newton instead of nested Brent. The grid supplies every bracket
// (as in the Brent version); only the steps inside a bracket differ.
//   1. l_max: Newton on the PD profile P_0(u), whose derivatives come from the envelope theorem at
//      the inner maximum w*(u): P' = l_0, P'' = l_00 - l_01^2 / l_11. Bracket: the refined
//      estimate +-2 grid steps, widened by 2 steps while the maximum sits on a bracket edge that
//      is not a bound of the box.
//   2. Inner maximum over u_b: Newton on l_b, bracketed and widened exactly as the Brent
//      version's, started from the last inner maximum moved along its slope dw/du = -l_ab / l_bb.
//   3. Endpoints: the same grid walk. Between the last inside point and the first grid point
//      whose true profile is below the threshold, Newton on P_a(u) - (l_max - c) with P_a' = l_a
//      (envelope), from the outside point, converging by bracketing (newton_root). The endpoint is
//      the evaluated point with the smallest |P_a - (l_max - c)|, the residual reported: no extra
//      evaluation.
// Tolerances: kProfileMaxTol for the maxima and kProfileRootTol for the root, in scaled units, as
// the Brent version; the residual TOL_PROFILE_ENDPOINT_RESIDUAL_LL bounds the result.
template <class Objective, class Integrator>
ProfileIntervals2 profile_intervals_newton(const Objective& objective, const Integrator& primary,
                                           const typename Objective::Obs* obs, std::int64_t periods,
                                           const Grid<2>& grid, const std::vector<double>& L, const Estimate2& est,
                                           double threshold = kProfileThreshold95) {
    using profile_detail::Point1;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    ProfileIntervals2 out{{nan, nan}, {nan, nan}, {0u, 0u}, nan, {nan, nan}, 0.0, 0};
    if ((est.flags & (kFlagFlatSurface | kFlagNumeric)) || !std::isfinite(est.value[0]) || !std::isfinite(est.value[1])) {
        out.flags[0] = out.flags[1] = kIntervalNotComputed;
        return out;
    }

    std::vector<std::int64_t> first, count;
    profile_detail::distinct_observations(obs, periods, first, count);
    const Axis* const ax[2] = {&grid.axis[0], &grid.axis[1]};
    const auto derivs = [&](double u0, double u1) {
        ++out.evaluations;
        return profile_detail::panel_derivs(objective, primary, obs, first, count, ax, u0, u1);
    };

    const std::vector<double> S = profile_detail::summed_surface(grid, L, periods);
    const auto surface = [&](int a, std::int32_t i, std::int32_t j) {
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
    const double edge_tol = 1e3 * kProfileMaxTol;

    // Warm start of the inner maximisation along each axis: the last inner maximum and its slope.
    Point1 last[2] = {{nan, nan, nan, nan, nan, nan}, {nan, nan, nan, nan, nan, nan}};

    // P_a(u): maximise over the nuisance u_b. Returns P_a with its envelope derivatives, the
    // maximiser w and dw/du.
    const auto inner = [&](int a, double u) {
        const int b = 1 - a;
        const double pos = (u - ax[a]->scaled_lo()) / ax[a]->step();
        const std::int32_t j0 = grid_argmax_b(a, clamp_index(a, std::floor(pos)));
        const std::int32_t j1 = grid_argmax_b(a, clamp_index(a, std::ceil(pos)));
        const double hb = ax[b]->step(), blo = ax[b]->scaled_lo(), bhi = ax[b]->scaled_hi();
        double lo = std::fmax(blo, ax[b]->scaled_at(j0 < j1 ? j0 : j1) - 2.0 * hb);
        double hi = std::fmin(bhi, ax[b]->scaled_at(j0 < j1 ? j1 : j0) + 2.0 * hb);
        double w0 = 0.5 * (ax[b]->scaled_at(j0) + ax[b]->scaled_at(j1));
        if (std::isfinite(last[a].w) && std::isfinite(last[a].dwdu)) w0 = last[a].w + last[a].dwdu * (u - last[a].x);
        // newton_max returns its last or its best evaluation; keep the full derivatives of both.
        profile_detail::PanelDerivs at_last{}, at_best{};
        double w_last = nan, w_best = nan;
        const auto f = [&](double w) {
            at_last = a == 0 ? derivs(u, w) : derivs(w, u);
            w_last = w;
            if (!(at_last.l <= at_best.l) || std::isnan(w_best)) {
                at_best = at_last;
                w_best = w;
            }
            return Point1{w, at_last.l, at_last.g[b], at_last.h[b], 0.0, 0.0};
        };
        Point1 m{};
        profile_detail::PanelDerivs at_m{};
        for (int widen = 0; widen < 64; ++widen) {
            w_best = nan;
            m = profile_detail::newton_max(f, lo, hi, w0, kProfileMaxTol);
            at_m = m.x == w_last ? at_last : at_best;
            if (!std::isfinite(m.f)) break;
            const bool at_lo = m.x - lo <= edge_tol && lo > blo;
            const bool at_hi = hi - m.x <= edge_tol && hi < bhi;
            if (!at_lo && !at_hi) break;
            if (at_lo) lo = std::fmax(blo, lo - 2.0 * hb);
            if (at_hi) hi = std::fmin(bhi, hi + 2.0 * hb);
            w0 = m.x;
        }
        const double hab = at_m.h[2], hbb = at_m.h[b], haa = at_m.h[a];
        const bool curved = hbb < 0.0;
        Point1 r{u, m.f, at_m.g[a], curved ? haa - hab * hab / hbb : haa, m.x, curved ? -hab / hbb : 0.0};
        if (std::isfinite(r.f)) last[a] = r;
        return r;
    };

    // 1. The polished maximum: Newton on the PD profile around the estimate.
    const double u_hat[2] = {grid::to_scaled(ax[0]->scale, est.value[0]), grid::to_scaled(ax[1]->scale, est.value[1])};
    double lo0 = std::fmax(ax[0]->scaled_lo(), u_hat[0] - 2.0 * ax[0]->step());
    double hi0 = std::fmin(ax[0]->scaled_hi(), u_hat[0] + 2.0 * ax[0]->step());
    last[0] = Point1{u_hat[0], nan, nan, nan, u_hat[1], 0.0};
    Point1 top{};
    double x0 = u_hat[0];
    for (int widen = 0; widen < 64; ++widen) {
        top = profile_detail::newton_max([&](double u) { return inner(0, u); }, lo0, hi0, x0, kProfileMaxTol);
        if (!std::isfinite(top.f)) break;
        const bool at_lo = top.x - lo0 <= edge_tol && lo0 > ax[0]->scaled_lo();
        const bool at_hi = hi0 - top.x <= edge_tol && hi0 < ax[0]->scaled_hi();
        if (!at_lo && !at_hi) break;
        if (at_lo) lo0 = std::fmax(ax[0]->scaled_lo(), lo0 - 2.0 * ax[0]->step());
        if (at_hi) hi0 = std::fmin(ax[0]->scaled_hi(), hi0 + 2.0 * ax[0]->step());
        x0 = top.x;
    }
    const double l_max = std::fmax(top.f, est.loglik);
    if (!std::isfinite(l_max) || !std::isfinite(top.f)) {
        out.flags[0] = out.flags[1] = kIntervalNotComputed;
        return out;
    }
    const double u_max[2] = {top.x, top.w};
    out.loglik_max = l_max;
    out.max_at[0] = grid::from_scaled(ax[0]->scale, u_max[0]);
    out.max_at[1] = grid::from_scaled(ax[1]->scale, u_max[1]);
    const double level = l_max - threshold;

    // The warm start for axis 1 at the maximum: rho's profile is maximised over PD, whose
    // maximiser moves with u_1 at the slope -l_01 / l_00.
    {
        const auto at = derivs(u_max[0], u_max[1]);
        last[1] = Point1{u_max[1], nan, nan, nan, u_max[0], at.h[0] < 0.0 ? -at.h[2] / at.h[0] : 0.0};
    }
    const Point1 start0 = top;
    const Point1 start1 = last[1];

    // 2. Each parameter, each side.
    for (int a = 0; a < 2; ++a) {
        const auto g = [&](double u) {
            Point1 p = inner(a, u);
            p.f -= level;
            return p;
        };
        for (const int side : {-1, +1}) {
            last[a] = a == 0 ? start0 : start1;  // each side walks out from the maximum
            const double bound = side < 0 ? ax[a]->scaled_lo() : ax[a]->scaled_hi();
            const double pos = (u_max[a] - ax[a]->scaled_lo()) / ax[a]->step();
            std::int32_t i = side < 0 ? static_cast<std::int32_t>(std::ceil(pos)) - 1
                                      : static_cast<std::int32_t>(std::floor(pos)) + 1;
            double inside = u_max[a];
            double endpoint = bound;
            bool truncated = true;
            bool failed = false;
            double residual = 0.0;
            for (; i >= 0 && i < ax[a]->n; i += side) {
                const double u = ax[a]->scaled_at(i);
                if (grid_profile(a, i) >= level) {  // a lower bound on P_a(u): inside, no work needed
                    inside = u;
                    continue;
                }
                const Point1 g_u = g(u);
                if (!std::isfinite(g_u.f)) {
                    failed = true;
                    break;
                }
                if (g_u.f >= 0.0) {
                    inside = u;
                    continue;
                }
                const Point1 r = profile_detail::newton_root(g, inside, u, g_u, kProfileRootTol);
                if (!std::isfinite(r.f)) {
                    failed = true;
                    break;
                }
                endpoint = r.x;
                residual = std::fabs(r.f);
                truncated = false;
                break;
            }
            if (failed) {
                out.flags[a] |= kIntervalNotComputed;
                continue;
            }
            if (!truncated) out.residual_max = std::fmax(out.residual_max, residual);
            const double v = truncated ? (side < 0 ? ax[a]->lo : ax[a]->hi) : grid::from_scaled(ax[a]->scale, endpoint);
            (side < 0 ? out.lo : out.hi)[a] = v;
            if (truncated) out.flags[a] |= side < 0 ? kIntervalLowerTruncated : kIntervalUpperTruncated;
        }
    }
    return out;
}

// The profile-likelihood intervals: by safeguarded Newton when the objective gives analytic
// derivatives (D-170), otherwise by nested Brent. Both solve the same equations to the same
// tolerances; TOL_PROFILE_ENDPOINT_RESIDUAL_LL bounds either.
template <class Objective, class Integrator>
ProfileIntervals2 profile_intervals(const Objective& objective, const Integrator& primary,
                                    const typename Objective::Obs* obs, std::int64_t periods, const Grid<2>& grid,
                                    const std::vector<double>& L, const Estimate2& est,
                                    double threshold = kProfileThreshold95) {
    if constexpr (profile_detail::has_contrib_derivs<Objective, Integrator>::value) {
        return profile_intervals_newton(objective, primary, obs, periods, grid, L, est, threshold);
    } else {
        return profile_intervals_brent(objective, primary, obs, periods, grid, L, est, threshold);
    }
}

}  // namespace vcal::engine
