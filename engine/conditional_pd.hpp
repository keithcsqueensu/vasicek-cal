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
// Solvers (D-170): steps 1 and 2 solve by safeguarded Newton when the objective gives analytic
// derivatives (see conditional_pd_interval_solve), otherwise by Brent; everything else is shared.
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

// kNewton selects the solver at the two solve points (D-170): safeguarded Newton with analytic
// derivatives, or Brent. Everything else (guide, brackets, walk, limits) is shared.
//   Inner, along the curve at fixed x = Phi^-1(q): PD(rho) = Phi(y), y = sqrt(1 - rho) x - sqrt(rho) z_a,
//     y_rho = -x / (2 sqrt(1 - rho)) - z_a / (2 sqrt(rho)),
//     y_rho,rho = -x / (4 (1 - rho)^(3/2)) + z_a / (4 rho^(3/2)),
//     PD_rho = phi(y) y_rho,  PD_rho,rho = phi(y) (y_rho,rho - y y_rho^2),
//     dl/drho = l_PD PD_rho + l_rho,  d2l/drho2 = l_PD,PD PD_rho^2 + 2 l_PD,rho PD_rho + l_rho,rho + l_PD PD_rho,rho,
//   then to w by the axis's chain rule.
//   Outer: P_q'(s) = l_PD phi(y) sqrt(1 - rho) dx/ds at the inner maximum (envelope), with
//     dx/ds = q (1 - q) / phi(x).
// The root converges by bracketing (profile_detail::newton_root); the endpoint is its evaluated
// point with the smallest residual, reported without re-evaluation.
template <bool kNewton, class Objective, class Integrator>
ConditionalPdInterval conditional_pd_interval_solve(const Objective& objective, const Integrator& primary,
                                                    const typename Objective::Obs* obs, std::int64_t periods,
                                                    const Grid<2>& grid, const std::vector<double>& L,
                                                    const Estimate2& est, const ProfileIntervals2& prof,
                                                    double z_a = kAdverseZ999, double threshold = kProfileThreshold95) {
    using profile_detail::brent_min;
    using profile_detail::Point1;
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
        // At an end of the axis, the bound itself: a round trip through w could move it by an ulp,
        // and the limit must be the q of the box's corner exactly (an estimate there has that q).
        double rho = best == 0 ? arho.lo : arho.hi;
        if (best > 0 && best < fine) {
            rho = grid::from_scaled(arho.scale, brent_min(f, w_at(best - 1), w_at(best + 1), kProfileMaxTol).x);
        }
        return Limit{logit_phi(conditional_pd_x(pd, rho, z_a)), pd, rho};
    };
    const Limit limit[2] = {edge_extreme(apd.lo, -1.0), edge_extreme(apd.hi, +1.0)};

    // P_q at s: the inner maximum over w, whether its maximiser is on a bound of the box, and
    // (Newton) dP_q/ds there.
    struct Inner {
        double value;
        bool on_bound;
        double slope;
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
        if (std::isnan(w0)) return {-inf, false, nan};
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
        // Newton: l along the curve with its w-derivatives; dP_q/ds (envelope) rides in dwdu.
        const double dx_ds = std::exp(special::log_phi(x) + special::log_phi(-x) + 0.5 * x * x +
                                      special::constants::kLnSqrt2Pi);
        const auto along = [&](auto w) {  // generic: instantiated only by the Newton branch
            const double rho = grid::from_scaled(arho.scale, w);
            const double sr = std::sqrt(rho), s1 = std::sqrt(1.0 - rho);
            const double y = std::sqrt(1.0 - rho) * x - std::sqrt(rho) * z_a;
            const double v[2] = {std::exp(special::log_phi(y)), rho};
            ++out.evaluations;
            const auto d = profile_detail::panel_derivs_natural(objective, primary, obs, first, count, v);
            const double phi_y = special::constants::kInvSqrt2Pi * std::exp(-0.5 * y * y);
            const double y_r = -x / (2.0 * s1) - z_a / (2.0 * sr);
            const double y_rr = -x / (4.0 * s1 * s1 * s1) + z_a / (4.0 * sr * sr * sr);
            const double p_r = phi_y * y_r;
            const double p_rr = phi_y * (y_rr - y * y_r * y_r);
            const double l_r = d.g[0] * p_r + d.g[1];
            const double l_rr = d.h[0] * p_r * p_r + 2.0 * d.h[2] * p_r + d.h[1] + d.g[0] * p_rr;
            const double j1 = grid::dvalue_dscaled(arho.scale, w), j2 = grid::d2value_dscaled2(arho.scale, w);
            return Point1{w, d.l, l_r * j1, l_rr * j1 * j1 + l_r * j2, 0.0, d.g[0] * phi_y * s1 * dx_ds};
        };
        profile_detail::Min1 m{};
        double slope = nan;
        double start = w0;  // Newton's start: the guide's seed, then the last maximum when widening
        const double edge_tol = 1e3 * kProfileMaxTol;
        bool at_lo = false, at_hi = false;
        for (int widen = 0; widen < 64; ++widen) {
            if constexpr (kNewton) {
                const Point1 p = profile_detail::newton_max(along, lo, hi, start, kProfileMaxTol);
                m = {p.x, -p.f};
                slope = p.dwdu;
            } else {
                m = brent_min(neg, lo, hi, kProfileMaxTol);
            }
            start = m.x;
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
        return {-m.f, (at_lo && lo_bound) || (at_hi && hi_bound), slope};
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
    // g with its slope (Newton), and w = 1 when the inner maximiser is on a bound of the box.
    const auto gp = [&](double s) {
        const Inner in = inner(s);
        return Point1{s, in.value - level, in.slope, 0.0, in.on_bound ? 1.0 : 0.0, 0.0};
    };
    for (const int side : {-1, +1}) {
        const Limit& lim = limit[side < 0 ? 0 : 1];
        const auto beyond = [&](double s) { return side < 0 ? s <= lim.s : s >= lim.s; };
        double inside = s_max;
        bool inside_exact = false;
        double g_inside = 0.0;
        double endpoint = lim.s;
        bool truncated = true, failed = false;
        double residual = 0.0;
        bool end_on_bound = false;
        for (std::int32_t k = 1;; ++k) {
            double s = s_max + static_cast<double>(side * k) * h;
            const bool at_limit = beyond(s);
            if (at_limit) s = lim.s;
            if (!at_limit && guide(s) >= level) {
                inside = s;
                inside_exact = false;
                continue;
            }
            // At the limit there is no curve to maximise over: no slope.
            const Point1 p_s = at_limit ? Point1{s, loglik(lim.pd, lim.rho) - level, nan, 0.0, 1.0, 0.0} : gp(s);
            const double g_s = p_s.f;
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
            double outside = s;
            Point1 p_outside = p_s;
            Point1 p_inside{inside, g_inside, nan, 0.0, 0.0, 0.0};
            if (!inside_exact) p_inside = gp(inside);
            while (!(p_inside.f >= 0.0) && inside != s_max) {  // the guide was optimistic: step back
                outside = inside;
                p_outside = p_inside;
                inside -= static_cast<double>(side) * h;
                if (side * (inside - s_max) < 0.0) inside = s_max;
                p_inside = gp(inside);
            }
            if (!(p_inside.f >= 0.0) || !(p_outside.f > -inf)) {
                failed = true;
                break;
            }
            if constexpr (kNewton) {
                const Point1 r = profile_detail::newton_root(gp, inside, outside, p_outside, kProfileRootTol);
                if (!std::isfinite(r.f)) {
                    failed = true;
                    break;
                }
                endpoint = r.x;
                residual = std::fabs(r.f);
                end_on_bound = r.w != 0.0;
            } else {
                endpoint = brent_root(g, inside, outside, p_inside.f, p_outside.f, kProfileRootTol);
                const Inner at = inner(endpoint);
                residual = std::fabs(at.value - level);
                end_on_bound = at.on_bound;
            }
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
            out.residual_max = std::fmax(out.residual_max, residual);
            if (end_on_bound) out.flags |= box;
        }
        // A truncated end is q at the limit point, evaluated directly (not through logit and back).
        (side < 0 ? out.lo : out.hi) = truncated ? conditional_pd(lim.pd, lim.rho, z_a)
                                                 : grid::from_scaled(AxisScale::Logit, endpoint);
    }
    if (out.flags & kIntervalNotComputed) out.lo = out.hi = nan;
    return out;
}

// The q interval: Newton when the objective gives analytic derivatives (D-170), otherwise Brent.
template <class Objective, class Integrator>
ConditionalPdInterval conditional_pd_interval(const Objective& objective, const Integrator& primary,
                                              const typename Objective::Obs* obs, std::int64_t periods,
                                              const Grid<2>& grid, const std::vector<double>& L, const Estimate2& est,
                                              const ProfileIntervals2& prof, double z_a = kAdverseZ999,
                                              double threshold = kProfileThreshold95) {
    constexpr bool kNewton = profile_detail::has_contrib_derivs<Objective, Integrator>::value;
    return conditional_pd_interval_solve<kNewton>(objective, primary, obs, periods, grid, L, est, prof, z_a, threshold);
}

}  // namespace vcal::engine
