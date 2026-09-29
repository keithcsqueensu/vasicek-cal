// SPDX-License-Identifier: Apache-2.0
//
// The binomial-mixture log-integrand and its adaptive-quadrature hint (ARCHITECTURE.md §1, §3.5).
//
//   g(z) = d log Phi(x(z)) + (n - d) log Phi(-x(z)),     so that
//   log E[p(Z)^d (1 - p(Z))^(n-d)] = log of the integral of exp(g(z)) phi(z) dz,
// the per-period log-likelihood without log C(n, d). A zero count's term is skipped, not
// multiplied by 0, so it can never produce 0 * inf.
//
// Hint: the mode and scale of h(z) = g(z) - z^2/2. log Phi is concave, so h'' <= -1: h is
// strictly concave and h' is strictly decreasing with exactly one root. With lambda the inverse
// Mills ratio, lambda'(t) = -lambda(t)(t + lambda(t)) in [-1, 0], and s = n - d:
//   h'(z)  = -beta (d lambda(x) - s lambda(-x)) - z
//   h''(z) =  beta^2 (d lambda'(x) + s lambda'(-x)) - 1
// Start: z_L, where p(z_L) = (d + 1/2)/(n + 1) (the large-n mode), pulled toward 0 by the
// N(0,1) prior in proportion to the precisions (tau_L z_L / (tau_L + 1)). Then bracket the
// root and run Newton with a bisection fallback, for at most kHintMaxIterations steps.
// Every step is plain arithmetic in a fixed order, so the hint is deterministic.
// scale = 1/sqrt(-h''(mode)) <= 1.
#pragma once

#include <cmath>
#include <cstdint>

#include "core/model/vasicek.hpp"
#include "core/precision.hpp"
#include "core/quadrature/integrator.hpp"
#include "core/special/inverse_mills.hpp"
#include "core/special/log_phi.hpp"
#include "core/special/generated/constants.hpp"
#include "core/special/probit.hpp"

namespace vcal::model {

struct BinomialLogIntegrand {
    Vasicek1F model;
    double defaults;   // d, exact as a double for counts up to 2^53
    double survivors;  // n - d

    VCAL_HD double operator()(double z) const {
        const double x = model.threshold(z);
        double v = 0.0;
        if (defaults > 0.0) v += defaults * special::log_phi(x);
        if (survivors > 0.0) v += survivors * special::log_phi(-x);
        return v;
    }
};

VCAL_HD BinomialLogIntegrand make_binomial_log_integrand(const Vasicek1F& m, std::int64_t n, std::int64_t d) {
    return {m, static_cast<double>(d), static_cast<double>(n - d)};
}

// Parameter derivatives of the log-integrand, for the score and Hessian of log I (D-170). With
// theta = (PD, rho), log I = log of the integral of exp(g) phi, and E the expectation under
// exp(g) phi / I (the nodes are fixed in z, so differentiation passes under the integral):
//     d log I / d theta_i          = E[g_i]
//     d2 log I / d theta_i d theta_j = E[g_ij + g_i g_j] - E[g_i] E[g_j].
// g depends on theta only through x = (c - sqrt(rho) z) / sqrt(1 - rho), c = Phi^-1(PD), so with
// a = sqrt(rho), b = sqrt(1 - rho) and lambda the inverse Mills ratio:
//     g_x  = d lambda(x) - s lambda(-x),   g_xx = d lambda'(x) + s lambda'(-x)
//     x_PD = c'/b, x_PD,PD = c''/b, where c' = 1/phi(c), c'' = c c'^2
//     x_rho = -z/(2ab) + x/(2b^2),  x_PD,rho = c'/(2b^3)
//     x_rho,rho = z(1 - 2rho)/(4a^3 b^3) + x_rho/(2b^2) + x/(2b^4)
//     g_i = g_x x_i,  g_ij = g_xx x_i x_j + g_x x_ij.
// operator() writes the five node functions whose expectations give the score and Hessian:
// (g_PD, g_rho, g_PD,PD + g_PD^2, g_rho,rho + g_rho^2, g_PD,rho + g_PD g_rho).
struct BinomialLogIntegrandDerivs {
    double c, a, b, rho;
    double c1, c2;  // dc/dPD, d2c/dPD2
    double defaults, survivors;

    VCAL_HD void operator()(double z, double (&f)[5]) const {
        const double x = (c - a * z) / b;
        double gx = 0.0, gxx = 0.0;
        if (defaults > 0.0) {
            const double lam = special::inverse_mills(x);
            gx += defaults * lam;
            gxx += defaults * (-lam * (x + lam));
        }
        if (survivors > 0.0) {
            const double lam = special::inverse_mills(-x);
            gx -= survivors * lam;
            gxx += survivors * (-lam * (lam - x));
        }
        const double b2 = b * b;
        const double x1 = c1 / b;
        const double x2 = -z / (2.0 * a * b) + x / (2.0 * b2);
        const double x11 = c2 / b;
        const double x12 = c1 / (2.0 * b2 * b);
        const double ab = a * b;
        const double x22 = z * (1.0 - 2.0 * rho) / (4.0 * ab * ab * ab) + x2 / (2.0 * b2) + x / (2.0 * b2 * b2);
        const double g1 = gx * x1, g2 = gx * x2;
        f[0] = g1;
        f[1] = g2;
        f[2] = gxx * x1 * x1 + gx * x11 + g1 * g1;
        f[3] = gxx * x2 * x2 + gx * x22 + g2 * g2;
        f[4] = gxx * x1 * x2 + gx * x12 + g1 * g2;
    }
};

// Requires 0 < pd < 1 and 0 < rho < 1 (as make_vasicek1f); otherwise the fields are NaN.
VCAL_HD BinomialLogIntegrandDerivs make_binomial_log_integrand_derivs(double pd, double rho, std::int64_t n,
                                                                       std::int64_t d) {
    const Vasicek1F m = make_vasicek1f(pd, rho);
    const double c1 = 1.0 / (special::constants::kInvSqrt2Pi * std::exp(-0.5 * m.c * m.c));
    return {m.c, m.sqrt_rho, m.sqrt_one_minus_rho, rho, c1, m.c * c1 * c1, static_cast<double>(d),
            static_cast<double>(n - d)};
}

namespace detail {

// lambda'(t) = -lambda(t)(t + lambda(t)), clamped to its true range [-1, 0]: for t << 0 the sum
// t + lambda(t) cancels, and only the sign and rough size matter for the hint.
VCAL_HD double inverse_mills_slope(double t) {
    const double lam = special::inverse_mills(t);
    const double v = -lam * (t + lam);
    return v < -1.0 ? -1.0 : (v > 0.0 ? 0.0 : v);
}

struct Slopes {
    double d1;  // h'(z)
    double d2;  // h''(z)
};

VCAL_HD Slopes mixture_slopes(const BinomialLogIntegrand& f, double z) {
    const double x = f.model.threshold(z);
    const double b = f.model.beta();
    double g1 = 0.0;
    double g2 = 0.0;
    if (f.defaults > 0.0) {
        g1 += f.defaults * special::inverse_mills(x);
        g2 += f.defaults * inverse_mills_slope(x);
    }
    if (f.survivors > 0.0) {
        g1 -= f.survivors * special::inverse_mills(-x);
        g2 += f.survivors * inverse_mills_slope(-x);
    }
    return {-b * g1 - z, b * b * g2 - 1.0};
}

}  // namespace detail

inline constexpr int kHintMaxIterations = 100;
inline constexpr double kHintRelativeTolerance = 1e-12;

VCAL_HD quadrature::IntegrandHint binomial_mixture_hint(const BinomialLogIntegrand& f) {
    const double n = f.defaults + f.survivors;
    if (!(n > 0.0)) return {0.0, 1.0, 0.0, 1.0, false};  // no observations: the integrand is the prior itself

    // Large-n starting point, shrunk toward the prior mean.
    const double x_l = special::probit((f.defaults + 0.5) / (n + 1.0));
    const double z_l = (f.model.c - f.model.sqrt_one_minus_rho * x_l) / f.model.sqrt_rho;
    const double tau_l = -(detail::mixture_slopes(f, z_l).d2 + 1.0);  // -g''(z_l) >= 0
    double z = tau_l * z_l / (tau_l + 1.0);
    const double step0 = 1.0 / std::sqrt(tau_l + 1.0);

    // Bracket the root of h': lo has h' > 0, hi has h' < 0.
    detail::Slopes s = detail::mixture_slopes(f, z);
    double lo = z;
    double hi = z;
    double step = step0;
    for (int i = 0; i < kHintMaxIterations && s.d1 != 0.0; ++i) {
        if (s.d1 > 0.0) {
            lo = hi;
            hi = z + step;
            if (detail::mixture_slopes(f, hi).d1 < 0.0) break;
        } else {
            hi = lo;
            lo = z - step;
            if (detail::mixture_slopes(f, lo).d1 > 0.0) break;
        }
        step *= 2.0;
    }

    // Newton, falling back to bisection whenever the step would leave the bracket.
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
        s = detail::mixture_slopes(f, z);
        if (done) break;
    }
    const double scale = 1.0 / std::sqrt(-s.d2);
    if (f.defaults > 0.0 && f.survivors > 0.0) return {z, scale, z, scale, false};

    // One-sided (d = 0 or d = n, D-118): the transition of the survival factor S, where S = 1/2.
    // d = 0: S = (1 - p)^n = 1/2 at p = q = 1 - 2^(-1/n); d = n: S = p^n = 1/2 at p = 1 - q.
    // Width: 1 / |d log S / dz| there, with d log S / dz = n beta lambda(-/+x).
    const double q = -std::expm1(-special::constants::kLn2 / n);
    const double b = f.model.beta();
    double x_c;
    double width;
    if (f.defaults == 0.0) {
        x_c = special::probit(q);
        width = 1.0 / (n * b * special::inverse_mills(-x_c));
    } else {
        x_c = special::probit_upper(q);
        width = 1.0 / (n * b * special::inverse_mills(x_c));
    }
    const double z_c = (f.model.c - f.model.sqrt_one_minus_rho * x_c) / f.model.sqrt_rho;
    return {z, scale, z_c, width, true};
}

}  // namespace vcal::model
