// SPDX-License-Identifier: Apache-2.0
//
// The multi-grade binomial MLE with a shared rho (M6; D-151, D-411, D-413; docs/methodology/multigrade_mle.md §3,
// §4): maximises sum_t l_t(PD_1, ..., PD_K, rho), l_t the K-grade period log-likelihood of
// core/objectives/multigrade.hpp under the parity rule, over the box PD_k in [pd_lo, pd_hi], rho in [rho_lo, rho_hi].
//
// Coordinates u_k = logit PD_k and v = logit rho. The search, fixed and deterministic:
//   1. Outer: v on kRhoGridPoints points uniform in v, endpoints exactly the box's. At each the inner solve gives the
//      profile l_p(v) = max over the PDs of l, and its derivative dl_p/dv, the v-component of the score at the
//      inner maximum (the envelope theorem).
//   2. Inner: the PDs at fixed v by safeguarded Newton in u with the analytic, dense K x K Hessian (the grades share
//      the factor, so the cross terms are not zero): the step solves (-H + lambda I) s = g with lambda raised from 0
//      until Cholesky succeeds, then a backtracking line search on l, projected onto the box. A grade at a bound
//      whose gradient, or whose step, points out of the box is held there (an active set). The first grid point
//      starts from each grade's pooled rate (sum d + 1/2) / (sum n + 1), clamped into the box; each later one from
//      the previous one's solution.
//   3. Refinement: from the grid argmax, safeguarded Newton on l_p(v) inside the bracket of its two neighbours, with
//      dl_p/dv from the envelope theorem and d2l_p/dv2 the Schur complement H_vv - H_vu H_uu^-1 H_uv on the free
//      grades; a bisection step on the sign of dl_p/dv whenever Newton's leaves the bracket or the curvature is not
//      negative. At a box end whose dl_p/dv points out, rho is that end and the fit is flagged (kFlagGridEdge).
//      Then one joint Newton step in every free coordinate, kept only if it raises l.
//
// At the estimate:
//   - SEs: the inverse of -H on the free coordinates (the grades off their bounds, and v unless rho is on a
//     bound), by Cholesky, then the delta method to the natural scale (D-185's convention). A coordinate on a
//     bound has SE NaN; if -H is not positive definite every SE is NaN and the fit is flagged kFlagFlatSurface.
//   - The rho interval: the profile-likelihood interval at engine::kProfileThreshold95, each end bracketed by the
//     grid's own profile values and solved by safeguarded Newton on l_p - (l^ - c) to engine::kProfileRootTol in
//     v; if the profile stays above the threshold up to the box, that end is the bound, flagged truncated
//     (D-128 to D-130).
//   - Each PD_k: the logit-scale Wald interval, u_k +- 1.96 SE_u. NaN for a grade on a bound.
//   - Checks, as calibrate's (D-120, D-178): the doubled rule (parity_rule(true)) at the estimate, flagged above
//     max(engine::kQuadratureCheckFlagAbs * T, kCheckUlps * eps * sum_t |l_t|).
//
// A grade with n = 0 in every period is empty: not estimated (PD, SE and interval NaN, kGradeEmpty). A grade with
// no defaults in any period has its MLE at PD = 0; it lands on pd_lo with kGradeAtLower and the fit's
// kFlagGradeEdge, its SE NaN, and rho is still estimated from every grade (owner's choice, 2026-11-10).
#pragma once

#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "core/grid.hpp"
#include "core/objectives/multigrade.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "engine/refine.hpp"

namespace vcal::benchmark_models::multigrade {

using model::GradeCount;
using model::kMaxGrades;

// Fit flags: engine kFlag* (kFlagGridEdge: rho on a box end; kFlagFlatSurface: -H not positive definite at the
// estimate; kFlagQuadratureUnconverged: the doubled-rule check; kFlagNumeric: a non-finite profile value) and these.
inline constexpr std::uint32_t kFlagGradeEdge = 1u << 7;     // some grade's PD is on a box bound
inline constexpr std::uint32_t kFlagNotConverged = 1u << 8;  // an inner solve or the refinement hit its limit
// Per-grade flags.
inline constexpr std::uint32_t kGradeAtLower = 1u << 0;
inline constexpr std::uint32_t kGradeAtUpper = 1u << 1;
inline constexpr std::uint32_t kGradeEmpty = 1u << 2;  // n = 0 in every period: not estimated

inline constexpr int kRhoGridPoints = 41;           // D-151: the recovery grid's rho axis
inline constexpr int kMaxInnerIterations = 100;
inline constexpr int kMaxRefineIterations = 60;
inline constexpr double kInnerTolerance = 1e-10;    // u units: the largest coordinate move that ends a solve
inline constexpr double kRefineTolerance = 1e-10;   // v units
inline constexpr double kCheckUlps = 256.0;
inline constexpr double kWaldZ95 = 1.959963984540054;  // Phi^-1(0.975)

struct Box {
    double pd_lo = 1e-5;  // D-411: below the 17-grade master scale's 0.03% top grade
    double pd_hi = 0.5;
    double rho_lo = 1e-3;
    double rho_hi = kDefaultRhoUpper;
};

enum class Status : std::int32_t { Ok = 0, InvalidPanel = 1, InvalidBox = 2, NotFinite = 3 };

struct Fit {
    int K = 0;
    std::vector<double> pd, pd_se, pd_lo, pd_hi;  // natural scale; the Wald interval
    std::vector<std::uint32_t> grade_flags;
    double rho = NAN, rho_se = NAN;
    double rho_lo = NAN, rho_hi = NAN;  // the profile-likelihood interval
    std::uint32_t rho_interval_flags = 0;  // engine kInterval*
    double loglik = NAN;
    double profile_curvature = NAN;  // d2 l_p / dv2 at the estimate (the Schur complement), NaN on a bound
    double doubled_check = NAN;      // |l(parity) - l(parity doubled)| at the estimate
    double interval_residual = 0.0;  // the largest |l_p - (l^ - c)| at a solved interval end
    int grid_index = -1;             // the profile's grid argmax
    std::vector<double> grid_rho, grid_profile;  // the outer grid and l_p on it
    int inner_unconverged = 0;       // inner solves that hit kMaxInnerIterations
    bool joint_step_taken = false;
    std::int64_t evaluations = 0;    // panel log-likelihood evaluations
    std::int64_t derivative_evaluations = 0;  // panel score-and-Hessian evaluations (each also gives the value)
    std::uint32_t flags = 0;
};

namespace detail {

inline double logit(double p) { return grid::to_scaled(AxisScale::Logit, p); }
inline double inv_logit(double u) { return grid::from_scaled(AxisScale::Logit, u); }

constexpr int kMaxDim = kMaxGrades + 1;

// Solves (-H + lambda I) s = g on the coordinates with act[i] (H is n x n row-major), raising lambda from 0 until
// Cholesky succeeds; s is zero off the active set. False if no active coordinate or no lambda works.
inline bool damped_step(const double* H, const double* g, const bool* act, int n, double* s, bool& damped) {
    damped = false;
    int idx[kMaxDim];
    int m = 0;
    for (int i = 0; i < n; ++i) {
        if (act[i]) idx[m++] = i;
    }
    for (int i = 0; i < n; ++i) s[i] = 0.0;
    if (m == 0) return false;
    double scale = 0.0;
    for (int a = 0; a < m; ++a) scale = std::fmax(scale, std::fabs(H[idx[a] * n + idx[a]]));
    if (!(scale > 0.0) || !std::isfinite(scale)) scale = 1.0;
    double lambda = 0.0;
    double L[kMaxDim * kMaxDim];
    for (int attempt = 0; attempt < 80; ++attempt) {
        bool ok = true;
        for (int a = 0; a < m && ok; ++a) {
            for (int b = 0; b <= a; ++b) {
                double sum = -H[idx[a] * n + idx[b]] + (a == b ? lambda : 0.0);
                for (int k = 0; k < b; ++k) sum -= L[a * kMaxDim + k] * L[b * kMaxDim + k];
                if (a == b) {
                    if (!(sum > 0.0)) {
                        ok = false;
                        break;
                    }
                    L[a * kMaxDim + a] = std::sqrt(sum);
                } else {
                    L[a * kMaxDim + b] = sum / L[b * kMaxDim + b];
                }
            }
        }
        if (ok) {
            double y[kMaxDim], x[kMaxDim];
            for (int a = 0; a < m; ++a) {
                double sum = g[idx[a]];
                for (int k = 0; k < a; ++k) sum -= L[a * kMaxDim + k] * y[k];
                y[a] = sum / L[a * kMaxDim + a];
            }
            bool finite = true;
            for (int a = m - 1; a >= 0; --a) {
                double sum = y[a];
                for (int k = a + 1; k < m; ++k) sum -= L[k * kMaxDim + a] * x[k];
                x[a] = sum / L[a * kMaxDim + a];
                finite = finite && std::isfinite(x[a]);
            }
            for (int a = 0; a < m; ++a) s[idx[a]] = x[a];
            damped = lambda > 0.0;
            return finite;
        }
        lambda = lambda == 0.0 ? 1e-8 * scale : 4.0 * lambda;
    }
    return false;
}

// The Cholesky factor of A (m x m, row-major, dense) into L; false if A is not positive definite.
inline bool cholesky(const double* A, int m, double* L) {
    for (int a = 0; a < m; ++a) {
        for (int b = 0; b <= a; ++b) {
            double sum = A[a * m + b];
            for (int k = 0; k < b; ++k) sum -= L[a * m + k] * L[b * m + k];
            if (a == b) {
                if (!(sum > 0.0) || !std::isfinite(sum)) return false;
                L[a * m + a] = std::sqrt(sum);
            } else {
                L[a * m + b] = sum / L[b * m + b];
            }
        }
        for (int b = a + 1; b < m; ++b) L[a * m + b] = 0.0;
    }
    return true;
}

// Solves L L' x = b in place.
inline void cholesky_solve(const double* L, int m, double* b) {
    for (int a = 0; a < m; ++a) {
        double sum = b[a];
        for (int k = 0; k < a; ++k) sum -= L[a * m + k] * b[k];
        b[a] = sum / L[a * m + a];
    }
    for (int a = m - 1; a >= 0; --a) {
        double sum = b[a];
        for (int k = a + 1; k < m; ++k) sum -= L[k * m + a] * b[k];
        b[a] = sum / L[a * m + a];
    }
}

class Search {
public:
    Search(const GradeCount* y, std::int64_t periods, int K, const Box& box)
        : y_(y), T_(periods), K_(K), n_(K + 1), v_lo_(logit(box.rho_lo)), v_hi_(logit(box.rho_hi)),
          u_lo_(logit(box.pd_lo)), u_hi_(logit(box.pd_hi)), pd_(static_cast<std::size_t>(K)),
          grad_(static_cast<std::size_t>(n_)), hess_(static_cast<std::size_t>(n_ * n_)),
          free_(static_cast<std::size_t>(K), true) {
        for (int k = 0; k < K_; ++k) {
            bool any = false;
            for (std::int64_t t = 0; t < T_; ++t) any = any || at(t, k).n > 0;
            free_[static_cast<std::size_t>(k)] = any;
        }
    }

    const GradeCount& at(std::int64_t t, int k) const {
        return y_[static_cast<std::size_t>(t) * static_cast<std::size_t>(K_) + static_cast<std::size_t>(k)];
    }
    const GradeCount* row(std::int64_t t) const {
        return y_ + static_cast<std::size_t>(t) * static_cast<std::size_t>(K_);
    }

    // The PDs at u; an empty grade's entry is any valid PD (the objective ignores grades with n = 0).
    void set_pd(const double* u) {
        for (int k = 0; k < K_; ++k) {
            pd_[static_cast<std::size_t>(k)] = free_[static_cast<std::size_t>(k)] ? inv_logit(u[k]) : 0.5;
        }
    }

    // The panel log-likelihood at (u, v) under rule; scale, if given, gets sum_t |l_t|.
    double value(const double* u, double v, const quadrature::SplitRule<PrecisionF64>& rule, double* scale = nullptr) {
        ++evaluations;
        set_pd(u);
        const double rho = inv_logit(v);
        double l = 0.0, s = 0.0;
        for (std::int64_t t = 0; t < T_; ++t) {
            const double lt = objectives::multigrade_log_contrib(row(t), K_, pd_.data(), rho, rule);
            l += lt;
            s += std::fabs(lt);
        }
        if (scale != nullptr) *scale = s;
        return l;
    }
    double value(const double* u, double v) { return value(u, v, rule_); }

    // The panel log-likelihood with its gradient g (n_) and Hessian H (n_ x n_, row-major) in (u, v). Empty grades
    // have zero rows and columns.
    double derivs(const double* u, double v, double* g, double* H) {
        ++derivative_evaluations;
        set_pd(u);
        const double rho = inv_logit(v);
        const int n = n_;
        for (int i = 0; i < n; ++i) g[i] = 0.0;
        for (int i = 0; i < n * n; ++i) H[i] = 0.0;
        double l = 0.0;
        for (std::int64_t t = 0; t < T_; ++t) {
            l += objectives::multigrade_log_contrib_derivs(row(t), K_, pd_.data(), rho, rule_, grad_.data(),
                                                           hess_.data());
            for (int i = 0; i < n; ++i) g[i] += grad_[static_cast<std::size_t>(i)];
            for (int i = 0; i < n * n; ++i) H[i] += hess_[static_cast<std::size_t>(i)];
        }
        // The chain rule to (u, v): theta_i = sigma(w_i), J_i = theta_i (1 - theta_i), J2_i = J_i (1 - 2 theta_i).
        double J[kMaxDim], J2[kMaxDim];
        for (int i = 0; i < n; ++i) {
            const double th = i < K_ ? pd_[static_cast<std::size_t>(i)] : rho;
            J[i] = th * (1.0 - th);
            J2[i] = J[i] * (1.0 - 2.0 * th);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) H[i * n + j] *= J[i] * J[j];
            H[i * n + i] += g[i] * J2[i];
        }
        for (int i = 0; i < n; ++i) g[i] *= J[i];
        return l;
    }

    // Safeguarded Newton on the coordinates with use[i] of w = (u, v); v is coordinate K_. Updates u (and v if
    // use[K_]) in place; g and H are the derivatives at the result.
    double newton(double* u, double v_in, const bool* use, int max_iterations, double tolerance, double* g, double* H,
                  bool& converged) {
        converged = false;
        double w[kMaxDim], lo[kMaxDim], hi[kMaxDim];
        for (int i = 0; i < K_; ++i) {
            w[i] = u[i];
            lo[i] = u_lo_;
            hi[i] = u_hi_;
        }
        w[K_] = v_in;
        lo[K_] = v_lo_;
        hi[K_] = v_hi_;
        double l = derivs(w, w[K_], g, H);
        for (int it = 0; it < max_iterations && std::isfinite(l); ++it) {
            bool act[kMaxDim];
            for (int i = 0; i < n_; ++i) {
                act[i] = use[i] && !((w[i] <= lo[i] && g[i] < 0.0) || (w[i] >= hi[i] && g[i] > 0.0));
            }
            double s[kMaxDim];
            bool damped = false;
            if (!damped_step(H, g, act, n_, s, damped)) {
                converged = true;  // nothing left to move
                break;
            }
            // A coordinate on a bound whose step points out is held as well, and the step recomputed.
            for (int pass = 0; pass < n_; ++pass) {
                bool changed = false;
                for (int i = 0; i < n_; ++i) {
                    if (act[i] && ((w[i] <= lo[i] && s[i] < 0.0) || (w[i] >= hi[i] && s[i] > 0.0))) {
                        act[i] = false;
                        changed = true;
                    }
                }
                if (!changed) break;
                if (!damped_step(H, g, act, n_, s, damped)) {
                    for (int i = 0; i < n_; ++i) s[i] = 0.0;
                    break;
                }
            }
            double slope = 0.0;
            for (int i = 0; i < n_; ++i) slope += g[i] * s[i];
            double alpha = 1.0, wn[kMaxDim];
            bool accepted = false;
            for (int k = 0; k < 60; ++k) {
                for (int i = 0; i < n_; ++i) wn[i] = std::fmin(hi[i], std::fmax(lo[i], w[i] + alpha * s[i]));
                const double ln = value(wn, wn[K_]);
                if (std::isfinite(ln) && ln >= l + 1e-4 * alpha * slope) {
                    accepted = true;
                    break;
                }
                alpha *= 0.5;
            }
            double moved = 0.0;
            for (int i = 0; i < n_; ++i) moved = std::fmax(moved, std::fabs(wn[i] - w[i]));
            if (!accepted || moved < tolerance) {
                converged = true;  // no ascent left at this resolution
                break;
            }
            for (int i = 0; i < n_; ++i) w[i] = wn[i];
            l = derivs(w, w[K_], g, H);
            if (moved < 10.0 * tolerance && !damped) {
                converged = true;
                break;
            }
        }
        for (int i = 0; i < K_; ++i) u[i] = w[i];
        last_v = w[K_];
        return l;
    }

    // The inner solve at fixed v: u updated in place. Returns l; g and H (in (u, v), all n_ coordinates) are the
    // derivatives at the returned u. converged: the last move was below kInnerTolerance.
    double newton_at_fixed_v(double* u, double v, double* g, double* H, bool& converged) {
        bool act[kMaxDim];
        for (int i = 0; i < n_; ++i) act[i] = i < K_ && free_[static_cast<std::size_t>(i)];
        return newton(u, v, act, kMaxInnerIterations, kInnerTolerance, g, H, converged);
    }

    // d2 l_p / dv2: H_vv - H_vu H_uu^-1 H_uv on the free grades off their bounds. NaN if -H_uu is not positive
    // definite.
    double schur(const double* u, const double* H) const {
        int idx[kMaxDim];
        int m = 0;
        for (int k = 0; k < K_; ++k) {
            if (free_[static_cast<std::size_t>(k)] && u[k] > u_lo_ && u[k] < u_hi_) idx[m++] = k;
        }
        const int n = n_;
        double hvv = H[K_ * n + K_];
        if (m == 0) return hvv;
        double A[kMaxDim * kMaxDim] = {}, L[kMaxDim * kMaxDim] = {}, x[kMaxDim] = {};  // GCC cannot see every read written
        for (int a = 0; a < m; ++a) {
            for (int b = 0; b < m; ++b) A[a * m + b] = -H[idx[a] * n + idx[b]];
            x[a] = H[idx[a] * n + K_];
        }
        if (!cholesky(A, m, L)) return static_cast<double>(NAN);
        cholesky_solve(L, m, x);  // x = (-H_uu)^-1 H_uv
        for (int a = 0; a < m; ++a) hvv += H[K_ * n + idx[a]] * x[a];
        return hvv;
    }

    const GradeCount* y_;
    std::int64_t T_;
    int K_, n_;
    double v_lo_, v_hi_, u_lo_, u_hi_;
    quadrature::SplitRule<PrecisionF64> rule_ = quadrature::parity_rule();
    std::vector<double> pd_, grad_, hess_;
    std::vector<bool> free_;  // non-empty grades
    double last_v = NAN;  // v at newton's result
    std::int64_t evaluations = 0;
    std::int64_t derivative_evaluations = 0;
};

}  // namespace detail

// nullptr if 0 < pd_lo < pd_hi < 1 and 0 < rho_lo < rho_hi < 1.
inline const char* box_error(const Box& box) {
    if (!(box.pd_lo > 0.0 && box.pd_lo < box.pd_hi && box.pd_hi < 1.0)) return "invalid PD box";
    if (!(box.rho_lo > 0.0 && box.rho_lo < box.rho_hi && box.rho_hi < 1.0)) return "invalid rho box";
    return nullptr;
}

// nullptr if the panel is valid: 1 <= K <= kMaxGrades, T >= 1, every period valid for multigrade_period_error.
inline const char* panel_error(const GradeCount* y, std::int64_t periods, int K) {
    if (!(K >= 1 && K <= kMaxGrades)) return "K must be in [1, kMaxGrades]";
    if (y == nullptr || periods < 1) return "the panel needs at least one period";
    for (std::int64_t t = 0; t < periods; ++t) {
        const char* e =
            objectives::multigrade_period_error(y + static_cast<std::size_t>(t) * static_cast<std::size_t>(K), K);
        if (e != nullptr) return e;
    }
    return nullptr;
}

// The MLE. y: periods x K counts, row-major (period t's grades at y[t K .. t K + K - 1]).
inline Status fit(const GradeCount* y, std::int64_t periods, int K, const Box& box, Fit& out) {
    out = Fit{};
    out.K = K;
    if (panel_error(y, periods, K) != nullptr) return Status::InvalidPanel;
    if (box_error(box) != nullptr) return Status::InvalidBox;
    using detail::inv_logit;
    detail::Search s(y, periods, K, box);
    const int n = K + 1;
    const auto Ks = static_cast<std::size_t>(K);
    out.pd.assign(Ks, NAN);
    out.pd_se.assign(Ks, NAN);
    out.pd_lo.assign(Ks, NAN);
    out.pd_hi.assign(Ks, NAN);
    out.grade_flags.assign(Ks, 0u);
    for (int k = 0; k < K; ++k) {
        if (!s.free_[static_cast<std::size_t>(k)]) out.grade_flags[static_cast<std::size_t>(k)] |= kGradeEmpty;
    }

    // The start: each grade's pooled rate, clamped into the box.
    std::vector<double> u(Ks, 0.0);
    for (int k = 0; k < K; ++k) {
        double nn = 0.0, dd = 0.0;
        for (std::int64_t t = 0; t < periods; ++t) {
            nn += static_cast<double>(s.at(t, k).n);
            dd += static_cast<double>(s.at(t, k).d);
        }
        u[static_cast<std::size_t>(k)] = std::fmin(s.u_hi_, std::fmax(s.u_lo_, detail::logit((dd + 0.5) / (nn + 1.0))));
    }

    // 1. The outer grid.
    const int G = kRhoGridPoints;
    std::vector<double> gv(static_cast<std::size_t>(G)), grid_v(static_cast<std::size_t>(G));
    std::vector<std::vector<double>> grid_u(static_cast<std::size_t>(G));
    out.grid_rho.resize(static_cast<std::size_t>(G));
    out.grid_profile.resize(static_cast<std::size_t>(G));
    std::vector<double> g(static_cast<std::size_t>(n)), H(static_cast<std::size_t>(n * n));
    bool all_converged = true;
    const double step = (s.v_hi_ - s.v_lo_) / static_cast<double>(G - 1);
    for (int j = 0; j < G; ++j) {
        const auto js = static_cast<std::size_t>(j);
        const double v = j == G - 1 ? s.v_hi_ : s.v_lo_ + static_cast<double>(j) * step;
        bool conv = false;
        const double l = s.newton_at_fixed_v(u.data(), v, g.data(), H.data(), conv);
        if (!conv) {
            ++out.inner_unconverged;
            all_converged = false;
        }
        grid_v[js] = v;
        grid_u[js] = u;
        gv[js] = g[static_cast<std::size_t>(K)];
        out.grid_rho[js] = inv_logit(v);
        out.grid_profile[js] = l;
        if (!std::isfinite(l)) out.flags |= engine::kFlagNumeric;
    }
    if (out.flags & engine::kFlagNumeric) {
        out.evaluations = s.evaluations;
        out.derivative_evaluations = s.derivative_evaluations;
        return Status::NotFinite;
    }
    int jstar = 0;
    for (int j = 1; j < G; ++j) {
        if (out.grid_profile[static_cast<std::size_t>(j)] > out.grid_profile[static_cast<std::size_t>(jstar)]) jstar = j;
    }
    out.grid_index = jstar;

    // 2. The refinement of v.
    const auto js = static_cast<std::size_t>(jstar);
    double v_hat = grid_v[js];
    u = grid_u[js];
    double l_hat = 0.0;
    bool rho_on_bound = false;
    {
        bool conv = false;
        l_hat = s.newton_at_fixed_v(u.data(), v_hat, g.data(), H.data(), conv);  // g and H at the grid argmax
        all_converged = all_converged && conv;
    }
    const auto& gk = g[static_cast<std::size_t>(K)];
    if ((jstar == 0 && gk <= 0.0) || (jstar == G - 1 && gk >= 0.0)) {
        rho_on_bound = true;
    } else {
        double a = grid_v[static_cast<std::size_t>(jstar > 0 ? jstar - 1 : 0)];
        double b = grid_v[static_cast<std::size_t>(jstar < G - 1 ? jstar + 1 : G - 1)];
        if (gk > 0.0) {
            a = v_hat;
        } else {
            b = v_hat;
        }
        double v = v_hat;
        bool refined = false;
        for (int it = 0; it < kMaxRefineIterations; ++it) {
            const double d1 = g[static_cast<std::size_t>(K)];
            const double d2 = s.schur(u.data(), H.data());
            double next = std::isfinite(d2) && d2 < 0.0 ? v - d1 / d2 : NAN;
            if (!(next > a && next < b)) next = 0.5 * (a + b);
            if (std::fabs(next - v) < kRefineTolerance || b - a < kRefineTolerance) {
                refined = true;
                break;
            }
            v = next;
            bool conv = false;
            const double l = s.newton_at_fixed_v(u.data(), v, g.data(), H.data(), conv);
            all_converged = all_converged && conv;
            if (g[static_cast<std::size_t>(K)] > 0.0) {
                a = v;
            } else {
                b = v;
            }
            v_hat = v;
            l_hat = l;
        }
        if (!refined) all_converged = false;
    }

    // 3. One joint Newton step in every free coordinate, kept only if it raises l.
    {
        bool use[detail::kMaxDim];
        for (int i = 0; i < n; ++i) use[i] = i < K ? s.free_[static_cast<std::size_t>(i)] : !rho_on_bound;
        std::vector<double> uj = u, gj(g.size()), Hj(H.size());
        bool conv = false;
        const double lj = s.newton(uj.data(), v_hat, use, 1, kInnerTolerance, gj.data(), Hj.data(), conv);
        if (lj > l_hat && std::isfinite(lj)) {
            u = uj;
            v_hat = s.last_v;
            l_hat = lj;
            g = gj;
            H = Hj;
            out.joint_step_taken = true;
        }
    }
    rho_on_bound = v_hat <= s.v_lo_ || v_hat >= s.v_hi_;
    if (rho_on_bound) out.flags |= engine::kFlagGridEdge;
    if (!all_converged) out.flags |= kFlagNotConverged;

    // The estimate and its flags.
    // A coordinate on a bound reports the box's own value, not inv_logit(logit(bound)).
    out.rho = v_hat <= s.v_lo_ ? box.rho_lo : v_hat >= s.v_hi_ ? box.rho_hi : inv_logit(v_hat);
    out.loglik = l_hat;
    for (int k = 0; k < K; ++k) {
        const auto ks = static_cast<std::size_t>(k);
        if (!s.free_[ks]) continue;
        out.pd[ks] = u[ks] <= s.u_lo_ ? box.pd_lo : u[ks] >= s.u_hi_ ? box.pd_hi : inv_logit(u[ks]);
        if (u[ks] <= s.u_lo_) out.grade_flags[ks] |= kGradeAtLower;
        if (u[ks] >= s.u_hi_) out.grade_flags[ks] |= kGradeAtUpper;
        if (out.grade_flags[ks] & (kGradeAtLower | kGradeAtUpper)) out.flags |= kFlagGradeEdge;
    }

    // The doubled-rule check.
    {
        double scale = 0.0;
        const double primary = s.value(u.data(), v_hat, s.rule_, &scale);
        const auto doubled_rule = quadrature::parity_rule(true);
        const double doubled = s.value(u.data(), v_hat, doubled_rule);
        out.doubled_check = std::fabs(doubled - primary);
        const double threshold = std::fmax(engine::kQuadratureCheckFlagAbs * static_cast<double>(periods),
                                           kCheckUlps * DBL_EPSILON * scale);
        if (!(out.doubled_check <= threshold)) out.flags |= engine::kFlagQuadratureUnconverged;
    }

    // SEs from -H on the free coordinates.
    {
        int idx[detail::kMaxDim];
        int m = 0;
        for (int k = 0; k < K; ++k) {
            const auto ks = static_cast<std::size_t>(k);
            if (s.free_[ks] && !(out.grade_flags[ks] & (kGradeAtLower | kGradeAtUpper))) idx[m++] = k;
        }
        if (!rho_on_bound) idx[m++] = K;
        std::vector<double> A(static_cast<std::size_t>(m * m)), L(A.size());
        for (int a = 0; a < m; ++a) {
            for (int b = 0; b < m; ++b) A[static_cast<std::size_t>(a * m + b)] = -H[static_cast<std::size_t>(idx[a] * n + idx[b])];
        }
        if (m > 0 && !detail::cholesky(A.data(), m, L.data())) {
            out.flags |= engine::kFlagFlatSurface;
        } else {
            std::vector<double> e(static_cast<std::size_t>(m));
            for (int a = 0; a < m; ++a) {
                for (int b = 0; b < m; ++b) e[static_cast<std::size_t>(b)] = a == b ? 1.0 : 0.0;
                detail::cholesky_solve(L.data(), m, e.data());
                const double var_w = e[static_cast<std::size_t>(a)];
                const double se_w = std::sqrt(var_w);
                const int i = idx[a];
                if (i < K) {
                    const auto is = static_cast<std::size_t>(i);
                    const double p = out.pd[is];
                    out.pd_se[is] = p * (1.0 - p) * se_w;
                    out.pd_lo[is] = inv_logit(u[is] - kWaldZ95 * se_w);
                    out.pd_hi[is] = inv_logit(u[is] + kWaldZ95 * se_w);
                } else {
                    out.rho_se = out.rho * (1.0 - out.rho) * se_w;
                }
            }
        }
        if (!rho_on_bound) out.profile_curvature = s.schur(u.data(), H.data());
    }

    // The rho profile interval.
    const double target = l_hat - engine::kProfileThreshold95;
    for (const int side : {-1, 1}) {
        const double bound = side < 0 ? s.v_lo_ : s.v_hi_;
        // The first grid point beyond v_hat on this side whose profile is below the target, and the last above.
        double v_in = v_hat;
        int j_out = -1;
        for (int j = 0; j < G; ++j) {
            const auto jj = static_cast<std::size_t>(side < 0 ? G - 1 - j : j);
            if (side * (grid_v[jj] - v_hat) <= 0.0) continue;
            if (out.grid_profile[jj] < target) {
                j_out = static_cast<int>(jj);
                break;
            }
            v_in = grid_v[jj];
        }
        double end;
        if (j_out < 0) {
            end = bound;
            out.rho_interval_flags |= side < 0 ? engine::kIntervalLowerTruncated : engine::kIntervalUpperTruncated;
        } else {
            double v_out = grid_v[static_cast<std::size_t>(j_out)];
            double vc = v_out, fc = out.grid_profile[static_cast<std::size_t>(j_out)] - target;
            double dc = gv[static_cast<std::size_t>(j_out)];
            std::vector<double> uc = grid_u[static_cast<std::size_t>(j_out)];
            for (int k = 0; k < 60 && std::fabs(v_out - v_in) > engine::kProfileRootTol; ++k) {
                double next = dc != 0.0 ? vc - fc / dc : 0.5 * (v_in + v_out);
                const double lo_b = std::fmin(v_in, v_out), hi_b = std::fmax(v_in, v_out);
                if (!(next > lo_b && next < hi_b)) next = 0.5 * (v_in + v_out);
                vc = next;
                bool conv = false;
                fc = s.newton_at_fixed_v(uc.data(), vc, g.data(), H.data(), conv) - target;
                dc = g[static_cast<std::size_t>(K)];
                if (!conv) out.flags |= kFlagNotConverged;
                if (fc >= 0.0) {
                    v_in = vc;
                } else {
                    v_out = vc;
                }
                if (std::fabs(fc) < 1e-12) break;
            }
            end = vc;
            out.interval_residual = std::fmax(out.interval_residual, std::fabs(fc));
        }
        (side < 0 ? out.rho_lo : out.rho_hi) =
            end <= s.v_lo_ ? box.rho_lo : end >= s.v_hi_ ? box.rho_hi : inv_logit(end);
    }
    out.evaluations = s.evaluations;
    out.derivative_evaluations = s.derivative_evaluations;
    return std::isfinite(out.loglik) ? Status::Ok : Status::NotFinite;
}

}  // namespace vcal::benchmark_models::multigrade
