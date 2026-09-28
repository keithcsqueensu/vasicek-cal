// SPDX-License-Identifier: Apache-2.0
//
// Vasicek-rate objective (ARCHITECTURE.md §3.5, M3; studies/mle-vs-mom/PREDICTION.md): the
// log-likelihood of an observed default rate r_t under the Vasicek rate distribution, closed form,
// no integral. With x = Phi^-1(r) and c = Phi^-1(PD),
//     l_t(PD, rho) = 1/2 log(1 - rho) - 1/2 log(rho) + x^2 / 2 - (sqrt(1 - rho) x - c)^2 / (2 rho),
// the log density of r = Phi((c + sqrt(rho) Z) / sqrt(1 - rho)), Z ~ N(0, 1) (the same model as the
// binomial mixture's p(z), in the engine's convention; the rate is its limit as n -> infinity).
//
// Rates of exactly 0 or 1 have no density. What happens to them is the zero-rate treatment, a
// template parameter so that panel_error can see it (D-044):
//   Refuse     (parity): the panel is refused; zero_rate_periods names the periods.
//   Censor     (native): a censored likelihood. A zero rate is a rate below the detection limit
//              `detect` and contributes log P(R <= detect) = log Phi((sqrt(1 - rho) Phi^-1(detect) - c)
//              / sqrt(rho)); a rate of 1 contributes log P(R >= 1 - detect).
//   Substitute (native): a continuity correction, a data edit stated as such. A zero rate enters
//              the density at `detect` (a rate of 1 at 1 - detect).
// Dropping such periods (the fourth treatment of S-28) is also a data edit, done by the caller
// before the fit (drop_zero_rate_periods); it is not an objective.
//
// For count data (n_t, d_t), rate_obs gives r = d/n and detect = 1/(2n), half a default.
#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "core/precision.hpp"
#include "core/special/log_phi.hpp"
#include "core/special/probit.hpp"

namespace vcal::objectives {

enum class ZeroRates { Refuse, Censor, Substitute };

struct RateObs {
    double rate;    // the observed default rate, 0 <= rate <= 1
    double detect;  // detection limit for rates of 0 and 1 (count data: 1/(2n)); 0 < detect < 1/2
    friend VCAL_HD constexpr bool operator==(const RateObs& a, const RateObs& b) {
        return a.rate == b.rate && a.detect == b.detect;
    }
};

// Count data: r = d/n, detect = 1/(2n). Requires n > 0 and 0 <= d <= n.
inline RateObs rate_obs(std::int64_t n, std::int64_t d) {
    return {static_cast<double>(d) / static_cast<double>(n), 0.5 / static_cast<double>(n)};
}

// The parameters at one grid point, shared by every treatment.
struct RateTheta {
    double pd;
    double rho;
};

inline bool is_boundary_rate(const RateObs& y) { return y.rate == 0.0 || y.rate == 1.0; }

// The periods (0-based) with a rate of 0 or 1, the ones parity refuses (D-044 names them).
inline std::vector<std::int64_t> zero_rate_periods(const RateObs* obs, std::int64_t periods) {
    std::vector<std::int64_t> out;
    for (std::int64_t t = 0; t < periods; ++t) {
        if (is_boundary_rate(obs[t])) out.push_back(t);
    }
    return out;
}

// The drop treatment (S-28): the panel without its boundary-rate periods. An explicit data edit.
inline std::vector<RateObs> drop_zero_rate_periods(const RateObs* obs, std::int64_t periods) {
    std::vector<RateObs> out;
    for (std::int64_t t = 0; t < periods; ++t) {
        if (!is_boundary_rate(obs[t])) out.push_back(obs[t]);
    }
    return out;
}

namespace detail {

// Phi^-1(r) for 0 < r < 1, taking the complement above 1/2 (1 - r is exact there, Sterbenz), so a
// rate near 1 keeps its digits (D-069).
VCAL_HD double rate_probit(double r) { return r < 0.5 ? special::probit(r) : special::probit_upper(1.0 - r); }

// The log density at x = Phi^-1(r), given c = Phi^-1(PD).
VCAL_HD double rate_log_density(double x, double c, double rho) {
    const double u = std::sqrt(1.0 - rho) * x - c;
    return 0.5 * std::log1p(-rho) - 0.5 * std::log(rho) + 0.5 * x * x - u * u / (2.0 * rho);
}

}  // namespace detail

template <class P, ZeroRates Z = ZeroRates::Refuse>
struct VasicekRate {
    static_assert(check_precision<P>());

    using Obs = RateObs;
    using Theta = RateTheta;
    static constexpr int n_params = 2;  // grid axis 0 = PD, axis 1 = rho
    static constexpr ZeroRates zero_rates = Z;

    static VCAL_HD Theta theta(const double (&v)[2]) { return {v[0], v[1]}; }

    static const char* panel_error(const Obs* obs, std::int64_t periods) {
        if (periods < 1) return "panel has no periods";
        for (std::int64_t t = 0; t < periods; ++t) {
            if (!(obs[t].rate >= 0.0 && obs[t].rate <= 1.0)) return "each period needs 0 <= rate <= 1";
            if (!(obs[t].detect > 0.0 && obs[t].detect < 0.5)) return "each period needs 0 < detect < 1/2";
            if (Z == ZeroRates::Refuse && is_boundary_rate(obs[t])) {
                return "a period has a default rate of 0 or 1: the parity rate likelihood refuses it "
                       "(zero_rate_periods names them; D-044)";
            }
        }
        return nullptr;
    }

    // The engine's quadrature check compares two integrators; this objective has no integral, so
    // the two agree exactly and the scale only has to be positive.
    static VCAL_HD double rounding_scale(const Obs&, double l) { return std::fabs(l) + 1.0; }

    template <class I>
    VCAL_HD typename P::accum_t log_contrib(const Obs& y, const Theta& th, const I&) const {
        const double nan = static_cast<double>(NAN);
        if (!(th.pd > 0.0 && th.pd < 1.0 && th.rho > 0.0 && th.rho < 1.0)) return nan;
        const double c = special::probit(th.pd);
        if (!is_boundary_rate(y)) return detail::rate_log_density(detail::rate_probit(y.rate), c, th.rho);
        const double xd = special::probit(y.detect);  // < 0
        const bool zero = y.rate == 0.0;
        if (Z == ZeroRates::Substitute) return detail::rate_log_density(zero ? xd : -xd, c, th.rho);
        if (Z == ZeroRates::Censor) {
            const double s1 = std::sqrt(1.0 - th.rho), s = std::sqrt(th.rho);
            // P(R <= detect) = Phi((s1 xd - c)/s); P(R >= 1 - detect) = Phi((s1 xd + c)/s).
            return special::log_phi(zero ? (s1 * xd - c) / s : (s1 * xd + c) / s);
        }
        return nan;  // Refuse: panel_error has already refused the panel
    }
};

// The closed-form maximum for a panel without boundary rates, over the whole parameter space
// (not the box): the MLE is the mean mu and 1/T variance v of x_t = Phi^-1(r_t), with
// rho-hat = v / (1 + v) and PD-hat = Phi(mu / sqrt(1 + v)). An independent check of the grid path;
// returns false if any rate is 0 or 1 or T < 2.
struct RateClosedForm {
    double pd;
    double rho;
};

inline bool vasicek_rate_closed_form(const RateObs* obs, std::int64_t periods, RateClosedForm& out) {
    if (periods < 2) return false;
    double sum = 0.0;
    for (std::int64_t t = 0; t < periods; ++t) {
        if (is_boundary_rate(obs[t]) || !(obs[t].rate > 0.0 && obs[t].rate < 1.0)) return false;
        sum += detail::rate_probit(obs[t].rate);
    }
    const double mu = sum / static_cast<double>(periods);
    double ss = 0.0;
    for (std::int64_t t = 0; t < periods; ++t) {
        const double e = detail::rate_probit(obs[t].rate) - mu;
        ss += e * e;
    }
    const double v = ss / static_cast<double>(periods);
    out.rho = v / (1.0 + v);
    out.pd = 0.5 * std::erfc(-(mu / std::sqrt(1.0 + v)) / std::sqrt(2.0));
    return true;
}

}  // namespace vcal::objectives
