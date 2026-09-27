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
