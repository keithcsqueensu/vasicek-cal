// SPDX-License-Identifier: Apache-2.0
//
// The conditional PD at an adverse factor quantile, and its profile-likelihood interval (S-23).
//
//     q(PD, rho) = Phi(x),   x = (Phi^-1(PD) + sqrt(rho) z_a) / sqrt(1 - rho),   z_a = Phi^-1(alpha),
//
// the PD given that the factor sits at its (1 - alpha) quantile: the adverse side under the
// engine's convention p(z) = Phi((Phi^-1(PD) - sqrt(rho) z) / sqrt(1 - rho)), where a higher Z means
// better conditions (R-2). alpha = 0.999 gives the quantity capital and stress calculations use.
// Everything is done in s = logit(q) = log Phi(x) - log Phi(-x), which keeps relative accuracy at
// both ends.
//
// Profile. On the curve q = c the PD is fixed by rho:
//     PD_c(rho) = Phi(sqrt(1 - rho) x_c - sqrt(rho) z_a),   x_c = Phi^-1(c),
// so P_q(c) = max over rho of l(PD_c(rho), rho), subject to (PD_c(rho), rho) lying in the box, is a
// one-dimensional maximisation. The 95% interval is {c : P_q(c) >= l_max - threshold}, with l_max
// the polished maximum of profile_intervals (D-129's threshold).
//   1. Inner maximum, over the rho axis's scaled coordinate w. A guide: on a fine scan of w (8
//      points per grid step) the surface, interpolated bilinearly in the axes' scaled coordinates at
//      (PD_c(rho), rho), at each point where PD_c(rho) lies in the box. The best point +-2 grid
//      steps, cut to the feasible stretch of w around it, brackets Brent's method, which widens by
//      2 steps (as in profile_intervals) until the maximum is interior or on a bound. A bound is an
//      end of the rho axis or a point where PD_c(rho) reaches an end of the PD axis.
//   2. End points. q is not a grid axis, so the crossing is bracketed by a walk in s from the
//      maximum outwards, in steps of the PD axis's step. Where the guide's maximum clears the
//      threshold the point is taken as inside without solving; the first point where it does not
//      is solved. The guide is not a bound on the profile, so the inside end of the bracket is
//      always solved too, and stepped back towards the maximum until it is inside. Brent's root
//      finder then solves P_q(c) = l_max - threshold to kProfileRootTol in s, and the residual is
//      reported (TOL_PROFILE_ENDPOINT_RESIDUAL_LL, as in D-128).
//   3. q attainable in the box runs from its minimum on the lower PD edge to its maximum on the
//      upper PD edge (q increases with PD). If the walk reaches that limit without leaving the
//      interval, the end point is the limit, flagged truncated: never extrapolated (D-130).
//   Box-limited: an end point is box-limited when the inner maximiser there lies on a bound of
//   the box, or the end point is the limit of q (truncated implies box-limited). Such an end can
//   be held by a bound through the nuisance without q itself reaching its range.
//
// Delta-method Wald (reported beside the profile): s-hat +- z_0.975 se_s, with se_s from the
// covariance of the axes' scaled coordinates at the estimate (D-119) and the gradient below.
//
// Serial and deterministic: results do not depend on threads. Headline fits only.
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "core/grid.hpp"
#include "core/special/generated/constants.hpp"
#include "core/special/log_phi.hpp"
#include "core/special/probit.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"

namespace vcal::engine {

// Phi^-1(0.999) (mpmath, 20 digits: 3.0902323061678135415).
inline constexpr double kAdverseZ999 = 3.0902323061678135;

// x(PD, rho), with q = Phi(x).
inline double conditional_pd_x(double pd, double rho, double z_a) {
    return (special::probit(pd) + std::sqrt(rho) * z_a) / std::sqrt(1.0 - rho);
}

inline double conditional_pd(double pd, double rho, double z_a = kAdverseZ999) {
    return std::exp(special::log_phi(conditional_pd_x(pd, rho, z_a)));
}

// s = logit(q), from x without cancellation.
inline double logit_phi(double x) { return special::log_phi(x) - special::log_phi(-x); }

// x = Phi^-1(c) for c = logit^-1(s), from whichever tail is small.
inline double probit_of_logit(double s) {
    return s >= 0.0 ? special::probit_upper(grid::from_scaled(AxisScale::Logit, -s))
                    : special::probit(grid::from_scaled(AxisScale::Logit, s));
}

struct ConditionalPdGradient {
    double s;         // logit(q)
    double ds_du[2];  // d s / d logit(PD), d s / d logit(rho)
};

// ds/dx = phi(x) / (q (1 - q)); dx/dlogit(PD) = PD (1 - PD) / (phi(Phi^-1(PD)) sqrt(1 - rho));
// dx/dlogit(rho) = rho (1 - rho) dx/drho = (z_a sqrt(rho) sqrt(1 - rho) + x rho) / 2.
inline ConditionalPdGradient conditional_pd_logit_gradient(double pd, double rho, double z_a = kAdverseZ999) {
    const double log_pdf_const = -special::constants::kLnSqrt2Pi;  // log phi(0)
    const double c = special::probit(pd);
    const double x = conditional_pd_x(pd, rho, z_a);
    const double ds_dx = std::exp(-0.5 * x * x + log_pdf_const - special::log_phi(x) - special::log_phi(-x));
    const double dx_dpd = std::exp(special::log_phi(c) + special::log_phi(-c) + 0.5 * c * c - log_pdf_const) /
                          std::sqrt(1.0 - rho);
    const double dx_drho = 0.5 * (z_a * std::sqrt(rho) * std::sqrt(1.0 - rho) + x * rho);
    return {logit_phi(x), {ds_dx * dx_dpd, ds_dx * dx_drho}};
}

enum : std::uint32_t {
    // With kIntervalLowerTruncated / kIntervalUpperTruncated (the end is the limit of q in the box)
    // and kIntervalNotComputed from profile.hpp.
    kIntervalLowerBoxLimited = 1u << 3,  // the inner maximiser at the lower end is on a bound of the box
    kIntervalUpperBoxLimited = 1u << 4,
};

struct ConditionalPdInterval {
    double estimate;  // q at the refined estimate
    double lo;        // natural scale
    double hi;
    std::uint32_t flags;
    double residual_max;       // max over solved end points of |P_q(c) - (l_max - threshold)|
    std::int64_t evaluations;  // panel log-likelihood evaluations used
};

template <class Objective, class Integrator>
ConditionalPdInterval conditional_pd_interval(const Objective& objective, const Integrator& primary,
                                              const typename Objective::Obs* obs, std::int64_t periods,
                                              const Grid<2>& grid, const std::vector<double>& L, const Estimate2& est,
                                              const ProfileIntervals2& prof, double z_a = kAdverseZ999,
                                              double threshold = kProfileThreshold95) {
    using profile_detail::brent_min;
    using profile_detail::brent_root;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    ConditionalPdInterval out{nan, nan, nan, 0u, 0.0, 0};
    if (std::isfinite(est.value[0]) && std::isfinite(est.value[1])) {
        out.estimate = conditional_pd(est.value[0], est.value[1], z_a);
    }
    if ((est.flags & (kFlagFlatSurface | kFlagNumeric)) || !std::isfinite(prof.loglik_max) ||
        !std::isfinite(out.estimate)) {
        out.flags = kIntervalNotComputed;
        return out;
    }

    std::vector<std::int64_t> first, count;
    profile_detail::distinct_observations(obs, periods, first, count);
    const Axis& apd = grid.axis[0];
    const Axis& arho = grid.axis[1];
    const auto loglik = [&](double pd, double rho) {
        ++out.evaluations;
        const double v[2] = {pd, rho};
        return profile_detail::panel_loglik(objective, primary, obs, first, count, v);
    };

    // The guide: the summed surface, bilinear in the scaled coordinates.
    const std::vector<double> S = profile_detail::summed_surface(grid, L, periods);
    const auto bilinear = [&](double u0, double u1) {
        const auto cell = [](const Axis& a, double u, std::int32_t& i) {
            const double pos = (u - a.scaled_lo()) / a.step();
            const double top = static_cast<double>(a.n - 2);
            const double f = std::floor(pos);
            i = static_cast<std::int32_t>(f < 0.0 ? 0.0 : f > top ? top : f);
            return pos - static_cast<double>(i);
        };
        std::int32_t i0 = 0, i1 = 0;
        const double t0 = cell(apd, u0, i0), t1 = cell(arho, u1, i1);
        const auto at = [&](std::int32_t a, std::int32_t b) {
            const std::int32_t idx[2] = {a, b};
            return S[static_cast<std::size_t>(grid.flatten(idx))];
        };
        return (1.0 - t0) * ((1.0 - t1) * at(i0, i1) + t1 * at(i0, i1 + 1)) +
               t0 * ((1.0 - t1) * at(i0 + 1, i1) + t1 * at(i0 + 1, i1 + 1));
    };

    // The rho axis's scaled coordinate w: its range, grid step and the fine scan of the guide.
    const double w_lo = arho.scaled_lo(), w_hi = arho.scaled_hi(), hw = arho.step();
    const std::int32_t fine = 8 * (arho.n - 1);
    const auto w_at = [&](std::int32_t k) { return k == fine ? w_hi : w_lo + static_cast<double>(k) * (hw / 8.0); };
    const double c_lo = special::probit(apd.lo), c_hi = special::probit(apd.hi);
    // PD on the curve with Phi^-1(q) = x at w; NaN outside the box's PD range.
    const auto pd_on_curve = [&](double x, double w) {
        const double rho = grid::from_scaled(arho.scale, w);
        const double y = std::sqrt(1.0 - rho) * x - std::sqrt(rho) * z_a;
        return (y >= c_lo && y <= c_hi) ? std::exp(special::log_phi(y)) : nan;
    };

    // The limits of q in the box: the minimum on the lower PD edge, the maximum on the upper.
    struct Limit {
        double s;
        double pd;
        double rho;
    };
    const auto edge_extreme = [&](double pd, double sign) {  // extreme of sign * s along w
        const auto f = [&](double w) {
            return -sign * logit_phi(conditional_pd_x(pd, grid::from_scaled(arho.scale, w), z_a));
        };
        std::int32_t best = 0;
        for (std::int32_t k = 1; k <= fine; ++k) {
            if (f(w_at(k)) < f(w_at(best))) best = k;
        }
        double w = w_at(best);
        if (best > 0 && best < fine) w = brent_min(f, w_at(best - 1), w_at(best + 1), kProfileMaxTol).x;
        const double rho = grid::from_scaled(arho.scale, w);
        return Limit{logit_phi(conditional_pd_x(pd, rho, z_a)), pd, rho};
    };
    const Limit limit[2] = {edge_extreme(apd.lo, -1.0), edge_extreme(apd.hi, +1.0)};

    // P_q at s: the inner maximum over w, and whether its maximiser is on a bound of the box.
    struct Inner {
        double value;
        bool on_bound;
    };
    const auto inner = [&](double s) -> Inner {
        const double x = probit_of_logit(s);
        const auto feasible = [&](double w) { return !std::isnan(pd_on_curve(x, w)); };
        // Seed: the best feasible point of the fine scan by the guide, else a limit's rho.
        std::int32_t best = -1;
        double best_g = -inf;
        for (std::int32_t k = 0; k <= fine; ++k) {
            const double w = w_at(k);
            const double pd = pd_on_curve(x, w);
            if (std::isnan(pd)) continue;
            const double gv = bilinear(grid::to_scaled(apd.scale, pd), w);
            if (best < 0 || gv > best_g) {
                best = k;
                best_g = gv;
            }
        }
        double w0 = nan;
        if (best >= 0) {
            w0 = w_at(best);
        } else {
            for (const Limit& lim : limit) {
                const double w = grid::to_scaled(arho.scale, lim.rho);
                if (feasible(w)) w0 = w;
            }
        }
        if (std::isnan(w0)) return {-inf, false};
        // The end of the feasible stretch from w0 towards `to` (a scan point beyond which, or at
        // which, feasibility is lost), by bisection; `bound` reports whether it is one.
        const auto cut = [&](double to, bool& bound) {
            bound = false;
            if (feasible(to)) return to;
            double in = w0, outp = to;
            // The first infeasible scan point between w0 and to, so the bisection brackets the
            // nearest change.
            const double step = (to > w0 ? 1.0 : -1.0) * hw / 8.0;
            for (double w = w0 + step; (to > w0) ? w < to : w > to; w += step) {
                if (!feasible(w)) {
                    outp = w;
                    break;
                }
                in = w;
            }
            for (int it = 0; it < 200 && std::fabs(outp - in) > 1e-13; ++it) {
                const double mid = 0.5 * (in + outp);
                (feasible(mid) ? in : outp) = mid;
            }
            bound = true;
            return in;
        };
        bool lo_bound = false, hi_bound = false;
        double lo = cut(std::fmax(w_lo, w0 - 2.0 * hw), lo_bound);
        double hi = cut(std::fmin(w_hi, w0 + 2.0 * hw), hi_bound);
        lo_bound = lo_bound || lo <= w_lo;
        hi_bound = hi_bound || hi >= w_hi;
        const auto neg = [&](double w) {
            const double rho = grid::from_scaled(arho.scale, w);
            const double y = std::sqrt(1.0 - rho) * x - std::sqrt(rho) * z_a;
            return -loglik(std::exp(special::log_phi(y)), rho);
        };
        profile_detail::Min1 m{};
        const double edge_tol = 1e3 * kProfileMaxTol;
        bool at_lo = false, at_hi = false;
        for (int widen = 0; widen < 64; ++widen) {
            m = brent_min(neg, lo, hi, kProfileMaxTol);
            at_lo = m.x - lo <= edge_tol;
            at_hi = hi - m.x <= edge_tol;
            const bool grow_lo = at_lo && !lo_bound, grow_hi = at_hi && !hi_bound;
            if (!grow_lo && !grow_hi) break;
            if (grow_lo) {
                lo = cut(std::fmax(w_lo, lo - 2.0 * hw), lo_bound);
                lo_bound = lo_bound || lo <= w_lo;
            }
            if (grow_hi) {
                hi = cut(std::fmin(w_hi, hi + 2.0 * hw), hi_bound);
                hi_bound = hi_bound || hi >= w_hi;
            }
        }
        return {-m.f, (at_lo && lo_bound) || (at_hi && hi_bound)};
    };
    // The guide's profile: the largest guide value along the curve (no likelihood evaluations).
    const auto guide = [&](double s) {
        const double x = probit_of_logit(s);
        double g = -inf;
        for (std::int32_t k = 0; k <= fine; ++k) {
            const double pd = pd_on_curve(x, w_at(k));
            if (!std::isnan(pd)) g = std::fmax(g, bilinear(grid::to_scaled(apd.scale, pd), w_at(k)));
        }
        return g;
    };

    const double level = prof.loglik_max - threshold;
    const double s_max = logit_phi(conditional_pd_x(prof.max_at[0], prof.max_at[1], z_a));
    const double h = apd.step();
    const auto g = [&](double s) { return inner(s).value - level; };
    for (const int side : {-1, +1}) {
        const Limit& lim = limit[side < 0 ? 0 : 1];
        const auto beyond = [&](double s) { return side < 0 ? s <= lim.s : s >= lim.s; };
        double inside = s_max;
        bool inside_exact = false;
        double g_inside = 0.0;
        double endpoint = lim.s;
        bool truncated = true, failed = false;
        for (std::int32_t k = 1;; ++k) {
            double s = s_max + static_cast<double>(side * k) * h;
            const bool at_limit = beyond(s);
            if (at_limit) s = lim.s;
            if (!at_limit && guide(s) >= level) {
                inside = s;
                inside_exact = false;
                continue;
            }
            const double g_s = at_limit ? loglik(lim.pd, lim.rho) - level : g(s);
            if (std::isnan(g_s)) {
                failed = true;
                break;
            }
            if (g_s >= 0.0) {
                inside = s;
                g_inside = g_s;
                inside_exact = true;
                if (at_limit) break;
                continue;
            }
            double outside = s, g_outside = g_s;
            if (!inside_exact) g_inside = g(inside);
            while (!(g_inside >= 0.0) && inside != s_max) {  // the guide was optimistic: step back
                outside = inside;
                g_outside = g_inside;
                inside -= static_cast<double>(side) * h;
                if (side * (inside - s_max) < 0.0) inside = s_max;
                g_inside = g(inside);
            }
            if (!(g_inside >= 0.0) || !(g_outside > -inf)) {
                failed = true;
                break;
            }
            endpoint = brent_root(g, inside, outside, g_inside, g_outside, kProfileRootTol);
            truncated = false;
            break;
        }
        if (failed) {
            out.flags |= kIntervalNotComputed;
            continue;
        }
        const std::uint32_t box = side < 0 ? kIntervalLowerBoxLimited : kIntervalUpperBoxLimited;
        if (truncated) {
            out.flags |= box | (side < 0 ? kIntervalLowerTruncated : kIntervalUpperTruncated);
        } else {
            const Inner at = inner(endpoint);
            out.residual_max = std::fmax(out.residual_max, std::fabs(at.value - level));
            if (at.on_bound) out.flags |= box;
        }
        (side < 0 ? out.lo : out.hi) = grid::from_scaled(AxisScale::Logit, endpoint);
    }
    if (out.flags & kIntervalNotComputed) out.lo = out.hi = nan;
    return out;
}

}  // namespace vcal::engine
