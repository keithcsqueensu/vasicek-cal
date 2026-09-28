// SPDX-License-Identifier: Apache-2.0
//
// Method of moments, joint-default-probability form (M3; D-042; studies/mle-vs-mom/PREDICTION.md).
//
// Two moments of the one-factor model: the PD, and the joint default probability of two obligors
// in the same period,
//     PD_2(PD, rho) = P(both default) = E[p(Z)^2] = Phi_2(c, c; rho),   c = Phi^-1(PD),
// the bivariate normal CDF at (c, c) with correlation rho. From a panel of counts both are estimated
// exactly in finite n by pooled ratios of per-period sufficient statistics,
//     PD-hat = sum_t d_t / sum_t n_t,   PD_2-hat = sum_t d_t (d_t - 1) / sum_t n_t (n_t - 1)
// (d(d - 1) / (n(n - 1)) is unbiased for p(z)^2 given z), and rho-hat solves PD_2(PD-hat, rho) =
// PD_2-hat. The statistics are sums over periods, so a resampling weight row W reweights them
// (D-042): the weighted sums replace the plain ones.
//
// E[p(Z)^2] is the binomial-mixture integral of a period with n = 2 and d = 2 (log C(2, 2) = 0), so
// it is computed by the same log integrand and the same integrator as the binomial MLE, in log
// space. It increases with rho, so rho-hat is the root of log PD_2(PD-hat, rho) - log PD_2-hat,
// found by Brent's method in logit(rho) over [rho_lo, rho_hi] (the box). Edge cases are flagged,
// never extrapolated:
//   - no defaults (PD-hat = 0): refused, no estimate;
//   - PD_2-hat at or below PD_2(PD-hat, rho_lo) (including PD_2-hat <= PD-hat^2, no dependence
//     seen): rho-hat = rho_lo, flagged kMomRhoAtFloor;
//   - PD_2-hat at or above PD_2(PD-hat, rho_hi): rho-hat = rho_hi, flagged kMomRhoAtCap;
//   - PD-hat outside [pd_lo, pd_hi]: flagged kMomPdOutsideBox, rho still solved at PD-hat.
//
// For a series of rates without counts (D-046), mom_from_rates uses PD-hat = mean r_t and
// PD_2-hat = mean r_t^2 (the moment form without the binomial correction).
#pragma once

#include <cmath>
#include <cstdint>

#include "core/model/binomial_mixture.hpp"
#include "core/model/vasicek.hpp"
#include "core/special/probit.hpp"
#include "engine/profile.hpp"

namespace vcal::engine {

enum : std::uint32_t {
    kMomRefused = 1u << 0,       // no defaults: no estimate
    kMomRhoAtFloor = 1u << 1,    // PD_2-hat at or below the floor's joint default probability
    kMomRhoAtCap = 1u << 2,      // PD_2-hat at or above the cap's
    kMomPdOutsideBox = 1u << 3,  // PD-hat outside the PD box
};

// Brent tolerance for rho-hat, in logit(rho).
inline constexpr double kMomRootTolerance = 1e-12;

struct MomentSums {
    double defaults = 0.0;       // sum d
    double pairs = 0.0;          // sum d (d - 1)
    double obligors = 0.0;       // sum n
    double obligor_pairs = 0.0;  // sum n (n - 1)
};

struct MomEstimate {
    double pd;
    double rho;
    double pd2;  // the joint default probability the estimate matches
    std::uint32_t flags;
};

// Sums over the periods of a count panel, each weighted by w[t] (w = nullptr: all 1). Accumulated
// in period order; the counts are exact in double up to 2^53.
template <class Obs>
MomentSums moment_sums(const Obs* obs, std::int64_t periods, const double* w = nullptr) {
    MomentSums s;
    for (std::int64_t t = 0; t < periods; ++t) {
        const double wt = w ? w[t] : 1.0;
        const double n = static_cast<double>(obs[t].n), d = static_cast<double>(obs[t].d);
        s.defaults += wt * d;
        s.pairs += wt * d * (d - 1.0);
        s.obligors += wt * n;
        s.obligor_pairs += wt * n * (n - 1.0);
    }
    return s;
}

// log PD_2(pd, rho) = log E[p(Z)^2], by the integrator.
template <class Integrator>
double log_joint_default(const Integrator& integrator, double pd, double rho) {
    const auto f = model::make_binomial_log_integrand(model::make_vasicek1f(pd, rho), 2, 2);
    return integrator.log_integrate(f, model::binomial_mixture_hint(f));
}

// rho-hat from PD-hat and PD_2-hat (both > 0), within [rho_lo, rho_hi].
template <class Integrator>
MomEstimate mom_from_moments(const Integrator& integrator, double pd, double pd2, double pd_lo, double pd_hi,
                             double rho_lo, double rho_hi) {
    MomEstimate out{pd, std::nan(""), pd2, 0u};
    if (!(pd > 0.0)) {
        out.flags |= kMomRefused;
        return out;
    }
    if (pd < pd_lo || pd > pd_hi) out.flags |= kMomPdOutsideBox;
    const double target = pd2 > 0.0 ? std::log(pd2) : -HUGE_VAL;
    const auto logit = [](double v) { return std::log(v) - std::log1p(-v); };
    const auto g = [&](double u) { return log_joint_default(integrator, pd, 1.0 / (1.0 + std::exp(-u))) - target; };
    const double a = logit(rho_lo), b = logit(rho_hi);
    const double ga = g(a), gb = g(b);
    if (!(ga < 0.0)) {
        out.rho = rho_lo;
        out.flags |= kMomRhoAtFloor;
    } else if (!(gb > 0.0)) {
        out.rho = rho_hi;
        out.flags |= kMomRhoAtCap;
    } else {
        const double u = profile_detail::brent_root(g, a, b, ga, gb, kMomRootTolerance);
        out.rho = 1.0 / (1.0 + std::exp(-u));
    }
    return out;
}

// The joint-default-probability MoM on a count panel, optionally weighted (a W row).
template <class Integrator, class Obs>
MomEstimate mom_from_counts(const Integrator& integrator, const Obs* obs, std::int64_t periods, double pd_lo,
                            double pd_hi, double rho_lo, double rho_hi, const double* w = nullptr) {
    const MomentSums s = moment_sums(obs, periods, w);
    if (!(s.defaults > 0.0) || !(s.obligor_pairs > 0.0)) {
        MomEstimate out{s.obligors > 0.0 ? s.defaults / s.obligors : std::nan(""), std::nan(""), std::nan(""), kMomRefused};
        return out;
    }
    return mom_from_moments(integrator, s.defaults / s.obligors, s.pairs / s.obligor_pairs, pd_lo, pd_hi, rho_lo, rho_hi);
}

// The moment form for a series of rates (no counts): PD-hat = mean r, PD_2-hat = mean r^2.
template <class Integrator>
MomEstimate mom_from_rates(const Integrator& integrator, const double* rates, std::int64_t periods, double pd_lo,
                           double pd_hi, double rho_lo, double rho_hi) {
    double s1 = 0.0, s2 = 0.0;
    for (std::int64_t t = 0; t < periods; ++t) {
        s1 += rates[t];
        s2 += rates[t] * rates[t];
    }
    const double T = static_cast<double>(periods);
    return mom_from_moments(integrator, s1 / T, s2 / T, pd_lo, pd_hi, rho_lo, rho_hi);
}

}  // namespace vcal::engine
