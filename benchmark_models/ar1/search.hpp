// SPDX-License-Identifier: Apache-2.0
//
// The AR(1)-aware binomial MLE (M6; D-195, D-196): maximises the filter's log-likelihood (filter.hpp) over
// (PD, rho, phi) in the box PD in [pd_lo, pd_hi], rho in [rho_lo, rho_hi], |phi| <= kPhiBound.
//
// Coordinates u = (logit PD, logit rho, phi). Every step is a safeguarded Newton step in u on the free
// coordinates: the gradient is the analytic score (Fisher's identity, filter.hpp), the Hessian central differences
// of it, damped towards gradient ascent (Levenberg) until it is negative definite, then a backtracking line search
// on the log-likelihood itself, projected onto the box. A coordinate at a bound whose gradient points out of the
// box is held there.
//
// The search, fixed and deterministic:
//   1. the scan: phi in {-0.6, -0.3, 0, 0.3, 0.6, 0.9}, a 2-D Newton in (PD, rho) at each, from the caller's start
//      at phi = 0 (the binomial fit) and warm-started outward from phi = 0 in each direction;
//   2. the polish: a 3-D Newton from each local maximum of the scan's profile values, the best kept.
// Reported for the owner's diagnostic: the number of local maxima along the scan and the log-likelihood gap to
// the second-best polished maximum.
//
// Checks at the estimate, as calibrate's (D-092, D-120): the doubled rule (parity_rule(true)) and the placement
// (Neighbours against Predictive), each flagged above max(kQuadratureCheckFlagAbs * T, kCheckUlps * eps * sum_t
// |l_t|).
#pragma once

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <vector>

#include "benchmark_models/ar1/filter.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "engine/refine.hpp"

namespace vcal::benchmark_models::ar1 {

inline constexpr std::uint32_t kFlagNotConverged = 1u << 7;    // a Newton search hit its iteration limit
inline constexpr std::uint32_t kFlagPlacementCheck = 1u << 8;  // the placement check exceeded its threshold
inline constexpr int kMaxNewtonIterations = 60;
inline constexpr double kStepTolerance = 1e-10;  // u units, the polish
inline constexpr double kScanTolerance = 1e-6;   // u units, the scan (each local maximum is polished)
inline constexpr double kHessianStep = 1e-4;     // u units, central differences of the score
inline constexpr double kCheckUlps = 256.0;
inline constexpr double kScanPhi[] = {-0.6, -0.3, 0.0, 0.3, 0.6, 0.9};

struct Box {
    double pd_lo, pd_hi, rho_lo, rho_hi;
};

struct Estimate3 {
    double value[3];  // PD, rho, phi
    double se[3];     // from the observed information at the estimate (NaN if not negative definite)
    double loglik;
    double doubled_check;    // |l(parity) - l(parity doubled)| at the estimate
    double placement_check;  // |Neighbours - Predictive| at the estimate
    int scan_local_maxima;
    int scan_unconverged;       // scan points that hit the iteration limit (a diagnostic; the polish decides)
    double second_maximum_gap;  // loglik minus the second-best polished maximum (inf if only one)
    double phi0_value[2];       // the scan's maximum over (PD, rho) at phi = 0: the binomial-mixture MLE
    double phi0_loglik;         // the log-likelihood there (S-35 M6's checks before scoring, PREDICTION.md)
    std::int64_t evaluations;   // log-likelihood evaluations (the three passes each)
    std::int64_t scores;        // score evaluations (a backward pass with pairs, after an evaluation)
    std::uint32_t flags;        // engine kFlag* (kFlagGridEdge: at a bound) and the ones above
};

namespace detail {

inline double logit(double p) { return std::log(p / (1.0 - p)); }
inline double inv_logit(double u) { return 1.0 / (1.0 + std::exp(-u)); }

class Search {
public:
    Search(const Obs* obs, std::int64_t periods, const Box& box)
        : obs_(obs), T_(periods), lo_{logit(box.pd_lo), logit(box.rho_lo), -kPhiBound},
          hi_{logit(box.pd_hi), logit(box.rho_hi), kPhiBound} {}

    static Theta theta(const double (&u)[3]) { return {inv_logit(u[0]), inv_logit(u[1]), u[2]}; }

    // The last evaluation and the search's current point, each with its per-period terms, so fit() can take the
    // estimate's log-likelihood and placement check from the evaluation the search already made there (P-18): the
    // same theta through the same code, so the same numbers bit for bit.
    struct Point {
        double u[3];
        LogLikelihood ll;
        std::vector<double> contrib;
        bool valid = false;
    };

    void record(const double (&u)[3], const LogLikelihood& ll) {
        for (int i = 0; i < 3; ++i) last_.u[i] = u[i];
        last_.ll = ll;
        last_.contrib = work_.contrib;
        last_.valid = true;
    }

    double value(const double (&u)[3]) {
        ++evaluations;
        const auto ll = log_likelihood(obs_, T_, theta(u), rule_, work_);
        record(u, ll);
        return ll.value;
    }

    // The log-likelihood and its gradient in u.
    double value_and_gradient(const double (&u)[3], double (&g)[3]) {
        const Theta th = theta(u);
        ++evaluations;
        const auto ll = log_likelihood(obs_, T_, th, rule_, work_);
        record(u, ll);
        const double l = ll.value;
        if (!std::isfinite(l)) {
            g[0] = g[1] = g[2] = static_cast<double>(NAN);
            return l;
        }
        ++scores;
        double s[3];
        score(obs_, T_, th, work_, s);
        g[0] = s[0] * th.pd * (1.0 - th.pd);
        g[1] = s[1] * th.rho * (1.0 - th.rho);
        g[2] = s[2];
        return l;
    }

    // The gradient in u at u, where the last log-likelihood evaluation (value or value_and_gradient) was made.
    void gradient_of_last(const double (&u)[3], double (&g)[3]) {
        const Theta th = theta(u);
        ++scores;
        double s[3];
        score(obs_, T_, th, work_, s);
        g[0] = s[0] * th.pd * (1.0 - th.pd);
        g[1] = s[1] * th.rho * (1.0 - th.rho);
        g[2] = s[2];
    }

    // Quasi-Newton on the coordinates with free[i]; u is updated in place. Returns the log-likelihood at u. The
    // Hessian is central differences of the score at the start, whenever the active set changes and whenever the
    // last step needed damping, and BFGS updates from the analytic score in between. BFGS keeps a negative definite
    // matrix negative definite (an update with y's <= 0 is skipped) but cannot make an indefinite one definite: far
    // from the maximum the differenced Hessian can be indefinite, and carried forward it would leave every later
    // step damped towards gradient ascent, converging linearly. So a damped step is followed by a fresh Hessian.
    double newton(double (&u)[3], const bool (&free)[3], double tolerance, bool& converged, bool warm = false) {
        double g[3];
        double l = value_and_gradient(u, g);
        current_ = last_;
        converged = false;
        double H[3][3] = {};
        bool have_h = false;
        bool prev_act[3] = {false, false, false};
        if (warm && warm_valid_) {  // the previous solve's Hessian, for the same active set only
            for (int i = 0; i < 3; ++i) {
                prev_act[i] = warm_act_[i];
                for (int j = 0; j < 3; ++j) H[i][j] = warm_h_[i][j];
            }
            have_h = true;
        }
        int backtracked = 0;  // consecutive steps the line search had to shorten
        bool was_damped = false;  // the last step needed lambda > 0: H was not negative definite
        for (int it = 0; it < kMaxNewtonIterations && std::isfinite(l); ++it) {
            // Active set: a free coordinate at a bound with its gradient pointing out is held.
            bool act[3];
            for (int i = 0; i < 3; ++i) {
                act[i] = free[i] && !((u[i] <= lo_[i] && g[i] < 0.0) || (u[i] >= hi_[i] && g[i] > 0.0));
            }
            if (!have_h || was_damped || backtracked >= 2 || act[0] != prev_act[0] || act[1] != prev_act[1] ||
                act[2] != prev_act[2]) {
                backtracked = 0;
                hessian(u, act, H);
                have_h = true;
                for (int i = 0; i < 3; ++i) prev_act[i] = act[i];
            }
            double s[3] = {0.0, 0.0, 0.0};
            if (!damped_step(H, g, act, s, was_damped)) break;
            // A coordinate at a bound whose step points out of the box is held too, and the step recomputed
            // (projected Newton); otherwise the projection returns it to the bound and the search cycles.
            for (int pass = 0; pass < 3; ++pass) {
                bool changed = false;
                for (int i = 0; i < 3; ++i) {
                    if (act[i] && ((u[i] <= lo_[i] && s[i] < 0.0) || (u[i] >= hi_[i] && s[i] > 0.0))) {
                        act[i] = false;
                        changed = true;
                    }
                }
                if (!changed) break;
                bool damped_again = false;
                if (!damped_step(H, g, act, s, damped_again)) {
                    s[0] = s[1] = s[2] = 0.0;
                    break;
                }
                was_damped = damped_again;  // the step taken is the one on the reduced set
            }
            // Backtracking line search, projected onto the box.
            double slope = 0.0;
            for (int i = 0; i < 3; ++i) slope += g[i] * s[i];
            double alpha = 1.0;
            double un[3];
            double ln = -HUGE_VAL;
            bool accepted = false;
            for (int k = 0; k < 40; ++k) {
                for (int i = 0; i < 3; ++i) un[i] = std::fmin(hi_[i], std::fmax(lo_[i], u[i] + alpha * s[i]));
                ln = value(un);
                if (std::isfinite(ln) && ln >= l + 1e-4 * alpha * slope) {
                    accepted = true;
                    accepted_ = last_;  // un's evaluation, before the Hessian's or gradient's probes
                    break;
                }
                alpha *= 0.5;
            }
            double moved = 0.0;
            for (int i = 0; i < 3; ++i) moved = std::fmax(moved, std::fabs(un[i] - u[i]));
            if (!accepted || moved < tolerance) {
                converged = true;  // no ascent left at this resolution
                break;
            }
            double gn[3];
            const double lnew = ln;
            gradient_of_last(un, gn);  // the filter's state is still un's, from the accepted line-search point
            // BFGS on A = -H with step d = un - u and y = -(gn - g), on the active coordinates.
            double d[3], y[3];
            for (int i = 0; i < 3; ++i) {
                d[i] = act[i] ? un[i] - u[i] : 0.0;
                y[i] = act[i] ? -(gn[i] - g[i]) : 0.0;
            }
            double yd = 0.0, Ad[3] = {0.0, 0.0, 0.0}, dAd = 0.0;
            for (int i = 0; i < 3; ++i) {
                yd += y[i] * d[i];
                for (int j = 0; j < 3; ++j) Ad[i] += -H[i][j] * d[j];
            }
            for (int i = 0; i < 3; ++i) dAd += d[i] * Ad[i];
            if (yd > 1e-12 * std::sqrt(dAd * dAd + 1e-300) && dAd > 0.0) {
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) H[i][j] -= -Ad[i] * Ad[j] / dAd + y[i] * y[j] / yd;
                }
            }
            for (int i = 0; i < 3; ++i) {
                u[i] = un[i];
                g[i] = gn[i];
            }
            l = lnew;
            current_ = accepted_;
            backtracked = alpha < 1.0 ? backtracked + 1 : 0;
            if (moved < 10.0 * tolerance) {
                converged = true;
                break;
            }
        }
        for (int i = 0; i < 3; ++i) last_g_[i] = g[i];
        if (have_h) {
            warm_valid_ = true;
            for (int i = 0; i < 3; ++i) {
                warm_act_[i] = prev_act[i];
                for (int j = 0; j < 3; ++j) warm_h_[i][j] = H[i][j];
            }
        }
        return l;
    }

    // Central differences of the gradient in u on the active coordinates (inactive rows and columns zero).
    void hessian(const double (&u)[3], const bool (&act)[3], double (&H)[3][3]) {
        for (auto& row : H) row[0] = row[1] = row[2] = 0.0;
        for (int j = 0; j < 3; ++j) {
            if (!act[j]) continue;
            double up[3] = {u[0], u[1], u[2]}, dn[3] = {u[0], u[1], u[2]};
            up[j] += kHessianStep;
            dn[j] -= kHessianStep;
            double gp[3], gm[3];
            value_and_gradient(up, gp);
            value_and_gradient(dn, gm);
            for (int i = 0; i < 3; ++i) H[i][j] = (gp[i] - gm[i]) / (2.0 * kHessianStep);
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < i; ++j) {
                const double m = 0.5 * (H[i][j] + H[j][i]);
                H[i][j] = H[j][i] = m;
            }
        }
    }

    // Solves (-H + lambda I) s = g on the active coordinates, raising lambda from 0 until -H + lambda I is
    // positive definite (Cholesky); damped is whether lambda > 0 was needed. False if no lambda works (non-finite
    // input).
    static bool damped_step(const double (&H)[3][3], const double (&g)[3], const bool (&act)[3], double (&s)[3],
                            bool& damped) {
        damped = false;
        int idx[3];
        int m = 0;
        for (int i = 0; i < 3; ++i) {
            if (act[i]) idx[m++] = i;
        }
        if (m == 0) return false;
        double scale = 0.0;
        for (int a = 0; a < m; ++a) scale = std::fmax(scale, std::fabs(H[idx[a]][idx[a]]));
        if (!(scale > 0.0) || !std::isfinite(scale)) scale = 1.0;
        double lambda = 0.0;
        for (int attempt = 0; attempt < 60; ++attempt) {
            double A[3][3], L[3][3] = {};
            for (int a = 0; a < m; ++a) {
                for (int b = 0; b < m; ++b) A[a][b] = -H[idx[a]][idx[b]] + (a == b ? lambda : 0.0);
            }
            bool ok = true;
            for (int a = 0; a < m && ok; ++a) {
                for (int b = 0; b <= a; ++b) {
                    double sum = A[a][b];
                    for (int k = 0; k < b; ++k) sum -= L[a][k] * L[b][k];
                    if (a == b) {
                        if (!(sum > 0.0)) {
                            ok = false;
                            break;
                        }
                        L[a][a] = std::sqrt(sum);
                    } else {
                        L[a][b] = sum / L[b][b];
                    }
                }
            }
            if (ok) {
                double y[3] = {0.0, 0.0, 0.0}, x[3] = {0.0, 0.0, 0.0};  // every entry read is written first; GCC cannot see it
                for (int a = 0; a < m; ++a) {
                    double sum = g[idx[a]];
                    for (int k = 0; k < a; ++k) sum -= L[a][k] * y[k];
                    y[a] = sum / L[a][a];
                }
                for (int a = m - 1; a >= 0; --a) {
                    double sum = y[a];
                    for (int k = a + 1; k < m; ++k) sum -= L[k][a] * x[k];
                    x[a] = sum / L[a][a];
                }
                s[0] = s[1] = s[2] = 0.0;
                for (int a = 0; a < m; ++a) s[idx[a]] = x[a];
                damped = lambda > 0.0;
                return std::isfinite(x[0]);
            }
            lambda = lambda == 0.0 ? 1e-6 * scale : 4.0 * lambda;
        }
        return false;
    }

    const Obs* obs_;
    std::int64_t T_;
    double lo_[3], hi_[3];
    quadrature::SplitRule<PrecisionF64> rule_ = quadrature::parity_rule();
    FilterWork work_;
    std::int64_t evaluations = 0;
    std::int64_t scores = 0;
    double last_g_[3] = {0.0, 0.0, 0.0};  // the gradient in u at newton's result
    Point last_, accepted_, current_;      // the last evaluation; the line search's accepted one; newton's result
    bool warm_valid_ = false;
    bool warm_act_[3] = {false, false, false};
    double warm_h_[3][3] = {};
};

}  // namespace detail

// The MLE. start_pd and start_rho: where the scan starts at phi = 0 (the binomial fit's estimate; any point in the
// box works). Returns false if the log-likelihood is not finite at the start.
inline bool fit(const Obs* obs, std::int64_t periods, const Box& box, double start_pd, double start_rho, Estimate3& out) {
    detail::Search s(obs, periods, box);
    out = Estimate3{};
    constexpr int K = static_cast<int>(sizeof(kScanPhi) / sizeof(kScanPhi[0]));
    double us[K][3];
    double ls[K];
    bool all_converged = true;
    int k0 = 0;
    while (kScanPhi[k0] != 0.0) ++k0;
    const bool pd_rho[3] = {true, true, false};
    {
        double u[3] = {std::fmin(s.hi_[0], std::fmax(s.lo_[0], detail::logit(start_pd))),
                       std::fmin(s.hi_[1], std::fmax(s.lo_[1], detail::logit(start_rho))), 0.0};
        bool conv;
        ls[k0] = s.newton(u, pd_rho, kScanTolerance, conv);
        if (!conv) ++out.scan_unconverged;
        for (int i = 0; i < 3; ++i) us[k0][i] = u[i];
        out.phi0_value[0] = detail::inv_logit(u[0]);
        out.phi0_value[1] = detail::inv_logit(u[1]);
        out.phi0_loglik = ls[k0];
        if (!std::isfinite(ls[k0])) return false;
    }
    for (int dir : {1, -1}) {
        for (int k = k0 + dir; k >= 0 && k < K; k += dir) {
            double u[3] = {us[k - dir][0], us[k - dir][1], kScanPhi[k]};
            bool conv;
            ls[k] = s.newton(u, pd_rho, kScanTolerance, conv, true);
            if (!conv) ++out.scan_unconverged;
            for (int i = 0; i < 3; ++i) us[k][i] = u[i];
        }
    }
    // Local maxima of the scan's profile (in phi order), each polished in 3-D.
    const bool all[3] = {true, true, true};
    double best_l = -HUGE_VAL, second_l = -HUGE_VAL;
    double best_u[3] = {us[k0][0], us[k0][1], us[k0][2]};
    detail::Search::Point best_point;
    for (int k = 0; k < K; ++k) {
        const bool left = k == 0 || ls[k] >= ls[k - 1];
        const bool right = k == K - 1 || ls[k] >= ls[k + 1];
        if (!(left && right && std::isfinite(ls[k]))) continue;
        ++out.scan_local_maxima;
        double u[3] = {us[k][0], us[k][1], us[k][2]};
        bool conv;
        const double l = s.newton(u, all, kStepTolerance, conv);
        all_converged = all_converged && conv;
        if (l > best_l) {
            second_l = best_l;
            best_l = l;
            for (int i = 0; i < 3; ++i) best_u[i] = u[i];
            best_point = s.current_;
        } else if (l > second_l) {
            second_l = l;
        }
    }
    out.second_maximum_gap = second_l == -HUGE_VAL ? HUGE_VAL : best_l - second_l;

    // The estimate, its checks and its observed information.
    const Theta th = detail::Search::theta(best_u);
    out.value[0] = th.pd;
    out.value[1] = th.rho;
    out.value[2] = th.phi;
    // The search evaluated the estimate already (P-18): take its numbers, else evaluate it once.
    const bool cached = best_point.valid && best_point.u[0] == best_u[0] && best_point.u[1] == best_u[1] &&
                        best_point.u[2] == best_u[2];
    if (!cached) {
        best_point.ll = log_likelihood(obs, periods, th, quadrature::parity_rule(), s.work_);
        best_point.contrib = s.work_.contrib;
    }
    const auto est = best_point.ll;
    out.loglik = est.value;
    out.placement_check = est.placement;
    double scale = 0.0;
    for (const double l : best_point.contrib) scale += std::fabs(l);
    detail::Search d(obs, periods, box);
    d.rule_ = quadrature::parity_rule(true);
    const double doubled = log_likelihood(obs, periods, th, d.rule_, d.work_).value;
    out.doubled_check = std::fabs(doubled - est.value);
    const double threshold =
        std::fmax(engine::kQuadratureCheckFlagAbs * static_cast<double>(periods), kCheckUlps * DBL_EPSILON * scale);
    if (!(out.doubled_check <= threshold)) out.flags |= engine::kFlagQuadratureUnconverged;
    if (!(out.placement_check <= threshold)) out.flags |= kFlagPlacementCheck;
    if (!all_converged) out.flags |= kFlagNotConverged;
    bool at_bound = false;
    for (int i = 0; i < 3; ++i) {
        if (best_u[i] <= s.lo_[i] || best_u[i] >= s.hi_[i]) at_bound = true;
    }
    if (at_bound) out.flags |= engine::kFlagGridEdge;

    // SEs: the inverse of the negative Hessian in u, by the delta method to the natural scale.
    double H[3][3];
    s.hessian(best_u, all, H);
    double A[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) A[i][j] = -H[i][j];
    }
    const double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                       A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                       A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    const double inv_diag[3] = {(A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det,
                                (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det,
                                (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det};
    const double jac[3] = {th.pd * (1.0 - th.pd), th.rho * (1.0 - th.rho), 1.0};
    const bool pos_def = A[0][0] > 0.0 && (A[0][0] * A[1][1] - A[0][1] * A[1][0]) > 0.0 && det > 0.0;
    for (int i = 0; i < 3; ++i) out.se[i] = pos_def ? jac[i] * std::sqrt(inv_diag[i]) : static_cast<double>(NAN);
    out.evaluations = s.evaluations;
    out.scores = s.scores;
    return std::isfinite(out.loglik);
}

// Profile-likelihood intervals, 95% (D-129's threshold, engine::kProfileThreshold95), for the parameters with
// which[i] (0 PD, 1 rho, 2 phi). As the engine's (engine/profile.hpp) but in three dimensions: the profile of
// parameter a is P_a(v) = max over the other two of l, in u coordinates, each value a warm-started 2-D Newton to
// kStepTolerance. dP_a/dv is the score's a-component at the inner maximum (the envelope theorem). Each side: from
// the Wald end (1.96 SE in u), Newton extrapolation outward until P_a falls below l_max - c, then safeguarded
// Newton on P_a(v) = l_max - c within the bracket, to engine::kProfileRootTol in u. If P_a
// stays above the threshold up to a bound, the end is that bound, flagged truncated, never extrapolated (D-130).
// residual_max is the largest |P_a - (l_max - c)| at a solved end.
struct ProfileIntervals3 {
    double lo[3], hi[3];  // natural scale
    std::uint32_t flags[3];
    double residual_max;
    std::int64_t evaluations;
    std::int64_t scores;
};

inline void profile_intervals(const Obs* obs, std::int64_t periods, const Box& box, const Estimate3& est,
                              const bool (&which)[3], ProfileIntervals3& out,
                              double threshold = engine::kProfileThreshold95) {
    out = ProfileIntervals3{};
    for (int a = 0; a < 3; ++a) {
        out.lo[a] = out.hi[a] = static_cast<double>(NAN);
        out.flags[a] = which[a] ? 0u : engine::kIntervalNotComputed;
    }
    detail::Search s(obs, periods, box);
    const double u_hat[3] = {detail::logit(est.value[0]), detail::logit(est.value[1]), est.value[2]};
    const double target = est.loglik - threshold;
    if (!std::isfinite(est.loglik)) {
        for (int a = 0; a < 3; ++a) out.flags[a] |= engine::kIntervalNotComputed;
        return;
    }
    const double jac[3] = {est.value[0] * (1.0 - est.value[0]), est.value[1] * (1.0 - est.value[1]), 1.0};
    for (int a = 0; a < 3; ++a) {
        if (!which[a]) continue;
        bool others[3] = {true, true, true};
        others[a] = false;
        const double se_u = est.se[a] / jac[a];
        const double step0 = std::isfinite(se_u) && se_u > 0.0 ? se_u : 0.25;
        for (int side : {-1, 1}) {
            double warm[3] = {u_hat[0], u_hat[1], u_hat[2]};
            // P_a(v) and its derivative dP_a/dv, the score's a-component at the inner maximum (envelope theorem).
            const auto profile = [&](double v, double& dp) {
                double u[3] = {warm[0], warm[1], warm[2]};
                u[a] = v;
                bool conv;
                const double l = s.newton(u, others, kStepTolerance, conv, true);
                for (int i = 0; i < 3; ++i) warm[i] = u[i];
                dp = s.last_g_[a];
                return l;
            };
            const double bound = side < 0 ? s.lo_[a] : s.hi_[a];
            const auto clamp = [&](double v) { return side < 0 ? std::fmax(v, bound) : std::fmin(v, bound); };
            double v_in = u_hat[a];  // inside: f >= 0
            double v_out = v_in, f_out = 0.0;
            bool bracketed = false;
            // First trial: the Wald end; then Newton extrapolation on f (or doubling), until f < 0 or the bound.
            double v = clamp(u_hat[a] + side * 1.959963984540054 * step0);
            double dp = 0.0;
            for (int k = 0; k < 40; ++k) {
                const double f = profile(v, dp) - target;
                if (f < 0.0) {
                    v_out = v;
                    f_out = f;
                    bracketed = true;
                    break;
                }
                const double prev = v_in;
                v_in = v;
                if (v == bound) break;
                double next = dp * side < 0.0 ? v - f / dp : v + 2.0 * (v - prev);
                if (!(side * (next - v) > 0.0) || !std::isfinite(next)) next = v + 2.0 * (v - prev);
                // Overshoot a little, so the crossing is bracketed rather than approached from inside.
                v = clamp(v + 1.25 * (next - v));
            }
            double end;
            if (!bracketed) {
                end = bound;
                out.flags[a] |= side < 0 ? engine::kIntervalLowerTruncated : engine::kIntervalUpperTruncated;
            } else {
                // Safeguarded Newton on f within [v_in, v_out]: bisection whenever the step leaves the bracket.
                double vc = v_out, fc = f_out, dc = dp;
                for (int k = 0; k < 60 && std::fabs(v_out - v_in) > engine::kProfileRootTol; ++k) {
                    double next = dc != 0.0 ? vc - fc / dc : 0.5 * (v_in + v_out);
                    const double lo_b = std::fmin(v_in, v_out), hi_b = std::fmax(v_in, v_out);
                    if (!(next > lo_b && next < hi_b)) next = 0.5 * (v_in + v_out);
                    vc = next;
                    fc = profile(vc, dc) - target;
                    if (fc >= 0.0) {
                        v_in = vc;
                    } else {
                        v_out = vc;
                    }
                    if (std::fabs(fc) < 1e-12) break;
                }
                end = vc;
                out.residual_max = std::fmax(out.residual_max, std::fabs(fc));
            }
            const double nat = a < 2 ? detail::inv_logit(end) : end;
            (side < 0 ? out.lo[a] : out.hi[a]) = nat;
        }
    }
    out.evaluations = s.evaluations;
    out.scores = s.scores;
}

}  // namespace vcal::benchmark_models::ar1
