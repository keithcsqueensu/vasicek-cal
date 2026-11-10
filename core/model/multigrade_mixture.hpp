// SPDX-License-Identifier: Apache-2.0
//
// The multi-grade binomial-mixture log-integrand and its quadrature hint (M6, D-151, D-411;
// docs/methodology/multigrade_mle.md). K grades share one factor Z and one rho; grade k has its own PD_k:
//
//   g(z) = sum_k [ d_k log Phi(x_k(z)) + (n_k - d_k) log Phi(-x_k(z)) ],
//   x_k(z) = (Phi^-1(PD_k) - sqrt(rho) z) / sqrt(1 - rho),
//
// so log of the integral of exp(g(z)) phi(z) dz is the period's log-likelihood without its binomial coefficients.
// Each grade's term is concave in z (log Phi of an affine function, twice), so h(z) = g(z) - z^2/2 is strictly concave
// and has exactly one mode: the single-grade hint's Newton-with-bisection applies to the summed slopes.
//
// K = 1 is the single-grade integrand (core/model/binomial_mixture.hpp) bit for bit: each grade's terms are that
// integrand's expressions in its order, the grades summed in index order from 0, and with one non-empty grade the
// hint is binomial_mixture_hint itself. Grades with n = 0 are dropped when the integrand is made, so they change
// nothing.
//
// The rule follows from the counts alone (D-411): interior (adaptive Gauss-Hermite) when some grade has 0 < d < n, or
// one grade has d = 0 and another d = n; one-sided (the composite rule) when every grade has d = 0, or every grade
// d = n. For a one-sided period the transition point z_c is where the product of the grades' survival factors is
// 1/2, and the width is 1 / |d log S / dz| there.
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "core/model/binomial_mixture.hpp"
#include "core/model/vasicek.hpp"
#include "core/precision.hpp"
#include "core/quadrature/integrator.hpp"
#include "core/special/erfcx_batch.hpp"
#include "core/special/generated/constants.hpp"
#include "core/special/inverse_mills.hpp"
#include "core/special/log_phi.hpp"
#include "core/special/probit.hpp"

namespace vcal::model {

inline constexpr int kMaxGrades = 32;

struct MultiGradeLogIntegrand {
    int grades = 0;              // non-empty grades only
    double c[kMaxGrades];        // Phi^-1(PD_k)
    double defaults[kMaxGrades];   // d_k
    double survivors[kMaxGrades];  // n_k - d_k
    double sqrt_rho;
    double sqrt_one_minus_rho;

    VCAL_HD double threshold(int k, double z) const { return (c[k] - sqrt_rho * z) / sqrt_one_minus_rho; }
    VCAL_HD double beta() const { return sqrt_rho / sqrt_one_minus_rho; }

    // One grade's term, BinomialLogIntegrand::operator()'s expressions, added to v.
    VCAL_HD void add_grade(int k, double z, double& v) const {
        const double x = threshold(k, z);
        if (defaults[k] > 0.0 && survivors[k] > 0.0) {
            const special::LogPhiPair lp = special::log_phi_pair(x);
            v += defaults[k] * lp.at_x;
            v += survivors[k] * lp.at_neg_x;
            return;
        }
        if (defaults[k] > 0.0) v += defaults[k] * special::log_phi(x);
        if (survivors[k] > 0.0) v += survivors[k] * special::log_phi(-x);
    }

    VCAL_HD double operator()(double z) const {
        double v = 0.0;
        for (int k = 0; k < grades; ++k) add_grade(k, z, v);
        return v;
    }

    // out[i] = (*this)(z[i]) bit for bit, for n <= kBatchMax: each grade's erfcx(|x| / sqrt 2) at every node from the
    // vector kernels (P-7, D-212), then each node finished grade by grade with operator()'s expressions. A node where
    // any grade's |x| is not finite, and every node if no kernel is selected, takes operator() itself. Host only.
    static constexpr int kBatchMax = 256;
    void batch(const double* z, double* out, int n) const {
        bool redo[kBatchMax];
        for (int i = 0; i < n; ++i) {
            out[i] = 0.0;
            redo[i] = false;
        }
        double x[kBatchMax];
        double a[kBatchMax] = {};
        double e[kBatchMax];
        bool finite[kBatchMax];
        for (int k = 0; k < grades; ++k) {
            for (int i = 0; i < n; ++i) {
                x[i] = threshold(k, z[i]);
                const double ai = std::fabs(x[i]) * special::constants::kInvSqrt2;
                finite[i] = std::isfinite(ai);
                a[i] = finite[i] ? ai : 0.0;
            }
            const bool ok = special::erfcx_batch::nonnegative(a, e, static_cast<std::size_t>(n));
            for (int i = 0; i < n; ++i) {
                if (!ok || !finite[i]) {
                    redo[i] = true;
                    continue;
                }
                if (redo[i]) continue;
                double& v = out[i];
                if (defaults[k] > 0.0 && survivors[k] > 0.0) {
                    const special::LogPhiPair lp = special::log_phi_pair_given_erfcx(x[i], e[i]);
                    v += defaults[k] * lp.at_x;
                    v += survivors[k] * lp.at_neg_x;
                    continue;
                }
                if (defaults[k] > 0.0) v += defaults[k] * special::log_phi_given_erfcx(x[i], e[i]);
                if (survivors[k] > 0.0) v += survivors[k] * special::log_phi_given_erfcx(-x[i], e[i]);
            }
        }
        for (int i = 0; i < n; ++i) {
            if (redo[i]) out[i] = (*this)(z[i]);
        }
    }

    // The grade as a single-grade integrand (for the hint at one grade, and for tests).
    VCAL_HD BinomialLogIntegrand grade(int k) const {
        return {Vasicek1F{c[k], sqrt_rho, sqrt_one_minus_rho}, defaults[k], survivors[k]};
    }
};

// One grade's counts.
struct GradeCount {
    std::int64_t n;
    std::int64_t d;
};

// The integrand for grades 0..K-1 with PDs pd[0..K-1] and a shared rho; grades with n = 0 are dropped. Requires
// 0 < pd_k < 1 for every non-empty grade, 0 < rho < 1 and K <= kMaxGrades; otherwise grades = -1 (invalid).
VCAL_HD MultiGradeLogIntegrand make_multigrade_log_integrand(const GradeCount* y, const double* pd, int K, double rho) {
    MultiGradeLogIntegrand f;
    f.grades = 0;
    if (!(K >= 0 && K <= kMaxGrades && rho > 0.0 && rho < 1.0)) {
        f.grades = -1;
        return f;
    }
    f.sqrt_rho = std::sqrt(rho);
    f.sqrt_one_minus_rho = std::sqrt(1.0 - rho);
    for (int k = 0; k < K; ++k) {
        if (y[k].n == 0) continue;
        if (!(pd[k] > 0.0 && pd[k] < 1.0)) {
            f.grades = -1;
            return f;
        }
        const int j = f.grades++;
        f.c[j] = special::probit(pd[k]);
        f.defaults[j] = static_cast<double>(y[k].d);
        f.survivors[j] = static_cast<double>(y[k].n - y[k].d);
    }
    return f;
}

namespace detail {

// h'(z) and h''(z) for the K-grade integrand: mixture_slopes's per-grade sums, accumulated over the grades in index
// order, then the same final expressions.
VCAL_HD Slopes multigrade_slopes(const MultiGradeLogIntegrand& f, double z) {
    const double b = f.beta();
    double g1 = 0.0;
    double g2 = 0.0;
    for (int k = 0; k < f.grades; ++k) {
        const double x = f.threshold(k, z);
        const double d = f.defaults[k], s = f.survivors[k];
        double lam_x = 0.0, lam_neg_x = 0.0;
        if (d > 0.0 && s > 0.0) {
            const special::InverseMillsPair p = special::inverse_mills_pair(x);
            lam_x = p.at_x;
            lam_neg_x = p.at_neg_x;
        } else if (d > 0.0) {
            lam_x = special::inverse_mills(x);
        } else if (s > 0.0) {
            lam_neg_x = special::inverse_mills(-x);
        }
        if (d > 0.0) {
            g1 += d * lam_x;
            g2 += d * inverse_mills_slope_given(x, lam_x);
        }
        if (s > 0.0) {
            g1 -= s * lam_neg_x;
            g2 += s * inverse_mills_slope_given(-x, lam_neg_x);
        }
    }
    return {-b * g1 - z, b * b * g2 - 1.0};
}

// log S(z) for a one-sided period: sum_k n_k log Phi(-x_k(z)) if every grade has d = 0 (increasing in z), or
// sum_k n_k log Phi(x_k(z)) if every grade has d = n (decreasing in z); with its z-derivative, sum_k n_k beta
// lambda(-x_k) or -sum_k n_k beta lambda(x_k).
struct LogSurvival {
    double value;
    double slope;
};

VCAL_HD LogSurvival multigrade_log_survival(const MultiGradeLogIntegrand& f, double z, bool zero_defaults) {
    const double b = f.beta();
    double v = 0.0, s = 0.0;
    for (int k = 0; k < f.grades; ++k) {
        const double n = f.defaults[k] + f.survivors[k];
        const double x = f.threshold(k, z);
        if (zero_defaults) {
            v += n * special::log_phi(-x);
            s += n * b * special::inverse_mills(-x);
        } else {
            v += n * special::log_phi(x);
            s -= n * b * special::inverse_mills(x);
        }
    }
    return {v, s};
}

}  // namespace detail

// How a period is integrated (D-411): interior if some grade has 0 < d < n, or grades with d = 0 and with d = n are
// both present; one-sided otherwise.
VCAL_HD bool multigrade_is_one_sided(const MultiGradeLogIntegrand& f) {
    bool any_zero = false, any_full = false;
    for (int k = 0; k < f.grades; ++k) {
        if (f.defaults[k] > 0.0 && f.survivors[k] > 0.0) return false;
        if (f.defaults[k] == 0.0) any_zero = true;
        if (f.survivors[k] == 0.0) any_full = true;
    }
    return !(any_zero && any_full);
}

inline constexpr int kMultigradeTransitionIterations = 200;

VCAL_HD quadrature::IntegrandHint multigrade_mixture_hint(const MultiGradeLogIntegrand& f) {
    if (f.grades <= 0) return {0.0, 1.0, 0.0, 1.0, false};  // no observations: the integrand is the prior itself
    if (f.grades == 1) return binomial_mixture_hint(f.grade(0));  // K = 1: the single-grade hint, bit for bit

    // Start: the precision-weighted mean of the grades' own large-n starts, shrunk toward the prior mean as the
    // single-grade hint does (tau z_L / (tau + 1), tau the summed precisions).
    double tau = 0.0, tz = 0.0;
    for (int k = 0; k < f.grades; ++k) {
        const BinomialLogIntegrand gk = f.grade(k);
        const double n = gk.defaults + gk.survivors;
        const double x_l = special::probit((gk.defaults + 0.5) / (n + 1.0));
        const double z_l = (gk.model.c - gk.model.sqrt_one_minus_rho * x_l) / gk.model.sqrt_rho;
        const double tau_k = -(detail::mixture_slopes(gk, z_l).d2 + 1.0);
        tau += tau_k;
        tz += tau_k * z_l;
    }
    double z = tz / (tau + 1.0);
    const double step0 = 1.0 / std::sqrt(tau + 1.0);
    if (!std::isfinite(z)) z = 0.0;

    // Bracket the root of h', then Newton with a bisection fallback: binomial_mixture_hint's loops on the summed
    // slopes.
    detail::Slopes s = detail::multigrade_slopes(f, z);
    double lo = z;
    double hi = z;
    double step = std::isfinite(step0) && step0 > 0.0 ? step0 : 1.0;
    for (int i = 0; i < kHintMaxIterations && s.d1 != 0.0; ++i) {
        if (s.d1 > 0.0) {
            lo = hi;
            hi = z + step;
            if (detail::multigrade_slopes(f, hi).d1 < 0.0) break;
        } else {
            hi = lo;
            lo = z - step;
            if (detail::multigrade_slopes(f, lo).d1 > 0.0) break;
        }
        step *= 2.0;
    }
    for (int i = 0; i < kHintMaxIterations && s.d1 != 0.0; ++i) {
        if (s.d1 > 0.0) {
            lo = z;
        } else {
            hi = z;
        }
        double next = z - s.d1 / s.d2;
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        const bool done = std::fabs(next - z) <= kHintRelativeTolerance * (1.0 + std::fabs(z));
        z = next;
        s = detail::multigrade_slopes(f, z);
        if (done) break;
    }
    const double scale = 1.0 / std::sqrt(-s.d2);
    if (!multigrade_is_one_sided(f)) return {z, scale, z, scale, false};

    // One-sided: the transition point, where log S = -ln 2, by bracketing from the mode and bisection (log S is
    // monotone in z), and the width 1 / |d log S / dz| there.
    bool zero_defaults = true;
    for (int k = 0; k < f.grades; ++k) zero_defaults = zero_defaults && f.defaults[k] == 0.0;
    const double target = -special::constants::kLn2;
    // r(z) = log S(z) - target is increasing in z for d = 0, decreasing for d = n; orient it to increase.
    const auto r = [&](double zz) {
        const double v = detail::multigrade_log_survival(f, zz, zero_defaults).value - target;
        return zero_defaults ? v : -v;
    };
    double a = z, b = z;
    double width_step = scale;
    if (r(z) < 0.0) {
        for (int i = 0; i < 200 && r(b) < 0.0; ++i) {
            a = b;
            b += width_step;
            width_step *= 2.0;
        }
    } else {
        for (int i = 0; i < 200 && r(a) >= 0.0; ++i) {
            b = a;
            a -= width_step;
            width_step *= 2.0;
        }
    }
    for (int i = 0; i < kMultigradeTransitionIterations; ++i) {
        const double mid = 0.5 * (a + b);
        if (mid == a || mid == b) break;
        (r(mid) < 0.0 ? a : b) = mid;
    }
    const double z_c = 0.5 * (a + b);
    const double slope = std::fabs(detail::multigrade_log_survival(f, z_c, zero_defaults).slope);
    return {z, scale, z_c, 1.0 / slope, true};
}

// Parameter derivatives of the K-grade log-integrand (D-170, extended): with theta = (PD_1, ..., PD_K, rho) and
// E the expectation under exp(g) phi / I,
//     d log I / d theta_i = E[g_i],   d2 log I / d theta_i d theta_j = E[g_ij + g_i g_j] - E[g_i] E[g_j].
// Grade k's term depends on (PD_k, rho) only, through x_k, so g_{PD_k PD_l} = 0 for k != l, and every per-grade
// quantity is BinomialLogIntegrandDerivs's expression for that grade; the rho terms are summed over the grades.
// operator() writes M = (K + 1) + (K + 1)(K + 2)/2 node functions: the K + 1 first derivatives, then, for i <= j in
// row order over (PD_1, ..., PD_K, rho), g_ij + g_i g_j. At K = 1 they are BinomialLogIntegrandDerivs's five, bit for
// bit: (g_PD, g_rho, then (PD, PD), (PD, rho), (rho, rho)); the caller maps them.
struct MultiGradeLogIntegrandDerivs {
    int grades = 0;
    double c[kMaxGrades], c1[kMaxGrades], c2[kMaxGrades];
    double defaults[kMaxGrades], survivors[kMaxGrades];
    double a, b, rho;

    static VCAL_HD int moments(int K) { return (K + 1) + (K + 1) * (K + 2) / 2; }
    // Index of the second moment (i, j), i <= j < K + 1.
    static VCAL_HD int second(int K, int i, int j) {
        const int n = K + 1;
        return n + i * n - i * (i - 1) / 2 + (j - i);
    }

    VCAL_HD void operator()(double z, double* f) const {
        const int K = grades;
        double g1[kMaxGrades], gx_[kMaxGrades], gxx_[kMaxGrades], x2_[kMaxGrades], x22_[kMaxGrades];
        double gr = 0.0;   // g_rho
        double grr = 0.0;  // g_rho,rho
        const double b2 = b * b;
        const double ab = a * b;
        for (int k = 0; k < K; ++k) {
            const double x = (c[k] - a * z) / b;
            double lam_x = 0.0, lam_neg_x = 0.0;
            const double d = defaults[k], s = survivors[k];
            if (d > 0.0 && s > 0.0) {
                const special::InverseMillsPair p = special::inverse_mills_pair(x);
                lam_x = p.at_x;
                lam_neg_x = p.at_neg_x;
            } else if (d > 0.0) {
                lam_x = special::inverse_mills(x);
            } else if (s > 0.0) {
                lam_neg_x = special::inverse_mills(-x);
            }
            double gx = 0.0, gxx = 0.0;
            if (d > 0.0) {
                const double lam = lam_x;
                gx += d * lam;
                gxx += d * (-lam * (x + lam));
            }
            if (s > 0.0) {
                const double lam = lam_neg_x;
                gx -= s * lam;
                gxx += s * (-lam * (lam - x));
            }
            const double x1 = c1[k] / b;
            const double x2 = -z / (2.0 * a * b) + x / (2.0 * b2);
            const double x22 = z * (1.0 - 2.0 * rho) / (4.0 * ab * ab * ab) + x2 / (2.0 * b2) + x / (2.0 * b2 * b2);
            g1[k] = gx * x1;
            gx_[k] = gx;
            gxx_[k] = gxx;
            x2_[k] = x2;
            x22_[k] = x22;
            gr += gx * x2;
            grr += gxx * x2 * x2 + gx * x22;
        }
        for (int k = 0; k < K; ++k) f[k] = g1[k];
        f[K] = gr;
        for (int i = 0; i < K; ++i) {
            const double x1 = c1[i] / b;
            const double x11 = c2[i] / b;
            const double x12 = c1[i] / (2.0 * b2 * b);
            f[second(K, i, i)] = gxx_[i] * x1 * x1 + gx_[i] * x11 + g1[i] * g1[i];
            for (int j = i + 1; j < K; ++j) f[second(K, i, j)] = g1[i] * g1[j];
            f[second(K, i, K)] = gxx_[i] * x1 * x2_[i] + gx_[i] * x12 + g1[i] * gr;
        }
        f[second(K, K, K)] = grr + gr * gr;
        (void)x22_;
    }
};

// Requires the same as make_multigrade_log_integrand; otherwise grades = -1.
VCAL_HD MultiGradeLogIntegrandDerivs make_multigrade_log_integrand_derivs(const GradeCount* y, const double* pd, int K,
                                                                           double rho) {
    MultiGradeLogIntegrandDerivs df;
    if (!(K >= 0 && K <= kMaxGrades && rho > 0.0 && rho < 1.0)) {
        df.grades = -1;
        return df;
    }
    df.a = std::sqrt(rho);
    df.b = std::sqrt(1.0 - rho);
    df.rho = rho;
    for (int k = 0; k < K; ++k) {
        if (y[k].n == 0) continue;
        if (!(pd[k] > 0.0 && pd[k] < 1.0)) {
            df.grades = -1;
            return df;
        }
        const int j = df.grades++;
        const double ck = special::probit(pd[k]);
        df.c[j] = ck;
        df.c1[j] = 1.0 / (special::constants::kInvSqrt2Pi * std::exp(-0.5 * ck * ck));
        df.c2[j] = ck * df.c1[j] * df.c1[j];
        df.defaults[j] = static_cast<double>(y[k].d);
        df.survivors[j] = static_cast<double>(y[k].n - y[k].d);
    }
    return df;
}

}  // namespace vcal::model
