// SPDX-License-Identifier: Apache-2.0
//
// Independent FP64 reference (D-012, D-031). Host-only, C++ standard library only, and no
// code or algorithm shared with core/ (the layering check enforces the code part, D-050):
//
//   piece          core/                                   ref/
//   log Phi        erfcx Chebyshev series, split x^2       std::erfc; Lentz continued fraction
//                                                          for the Mills ratio below x = -20
//   Phi^-1         Wichura AS241                           safeguarded Newton on ref's log Phi
//   log C(n,k)     Stirling differences                    compensated sum of log((n-k+i)/i)
//   mixture        adaptive Gauss-Hermite, stored tables   adaptive Gauss-Legendre 10/20 on
//   integral       (mpmath), Laplace mode/scale            panels around the mode, sized by the
//                                                          one-sided half-drop widths; nodes
//                                                          computed at start-up
//   estimate       grid argmax + 3x3 quadratic stencil     nested golden-section search
//                                                          (profile likelihood) in logit
//                                                          coordinates; finite-difference Hessian
//
// Speed is not a goal. Used to cross-check core (M1.7); its own accuracy is measured against
// the mpmath golden tables (tests/unit/ref_test.cpp).
#pragma once

#include <cstdint>
#include <vector>

namespace vcalref {

double log_ncdf(double x);                        // log Phi(x)
double ncdf_inv(double p);                        // Phi^-1(p); p = 0 -> -inf, 1 -> +inf, else outside [0,1] -> NaN
double lchoose(std::int64_t n, std::int64_t k);   // log C(n, k); NaN outside 0 <= k <= n or if min(k, n-k) > 1e8

struct Period {
    std::int64_t n;
    std::int64_t d;
};

// log E[p(Z)^d (1 - p(Z))^(n-d)] for the one-factor Vasicek p(z), 0 < pd < 1, 0 < rho < 1.
double log_mixture(double pd, double rho, std::int64_t n, std::int64_t d);

double period_loglik(double pd, double rho, const Period& y);  // lchoose(n, d) + log_mixture
double loglik(double pd, double rho, const std::vector<Period>& panel);

struct Fit {
    double pd;
    double rho;
    double loglik;
    double se_pd;  // from the observed information; NaN if the Hessian is not negative definite
    double se_rho;
    double corr;
    bool on_boundary;  // an estimate within the search tolerance of a bound
};

// Maximum-likelihood estimate over [pd_lo, pd_hi] x [rho_lo, rho_hi] (bounds in (0, 1)).
Fit fit(const std::vector<Period>& panel, double pd_lo, double pd_hi, double rho_lo, double rho_hi);

struct ProfileInterval {
    double lo;  // natural scale
    double hi;
    bool lo_at_bound;  // the profile stays above the threshold all the way to the bound
    bool hi_at_bound;
};

// Profile-likelihood interval for one parameter (0 = PD, 1 = rho) over the same box as fit():
// {v : max over the other parameter of loglik >= fit.loglik - threshold}. Independent of core's
// method: the profile maximises by golden-section over the whole range of the other parameter,
// and each end is found by bisection between the estimate and the bound, in logit coordinates.
ProfileInterval profile_interval(const std::vector<Period>& panel, const Fit& fit, int param, double pd_lo,
                                 double pd_hi, double rho_lo, double rho_hi, double threshold);

// The log-likelihood under an AR(1) factor, Z_1 ~ N(0, 1), Z_t | Z_{t-1} ~ N(phi Z_{t-1}, 1 - phi^2) (M6, D-195).
// Independent of core's AR(1) filter (benchmark_models/ar1/): a dense filter on one fixed uniform grid in z over
// [-10, 10], trapezoid weights, the transition density built directly, and ref's own log Phi and log C(n, k). The
// grid depends only on the panel (step = min(0.02, 0.4 / sqrt(max n)) unless given), never on the parameters, so
// the result is a smooth function of (PD, rho, phi): the oracle for core's goldens and its smoothness test. Valid
// where the conditional likelihood of every period is resolved by that step (rho <= 0.5 for n up to 10^4) and
// carries no mass beyond |z| = 10; |phi| < 1. Slow by design: O(T G^2) with G = 20 / step.
double ar1_loglik(double pd, double rho, double phi, const std::vector<Period>& panel, double step = 0.0);

}  // namespace vcalref
