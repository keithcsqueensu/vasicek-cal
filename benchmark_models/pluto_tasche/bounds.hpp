// SPDX-License-Identifier: Apache-2.0
//
// Pluto and Tasche's most-prudent upper bounds for a low-default portfolio (S-6; D-214, D-215;
// docs/methodology/pluto_tasche.md). A `native` benchmark: not parity, not in the C ABI.
//
// A panel of T periods with n obligors each and D = sum_t d_t defaults in total. At confidence gamma the
// bound is the PD p at which the observed total is just at the 1 - gamma quantile of its assumed law:
//     P(S <= D | p) = 1 - gamma.
//   - independent:  S ~ Bin(N, p), N = n T: the one-sided Clopper-Pearson bound. At D = 0 it is
//                   1 - (1 - gamma)^(1/N).
//   - one factor:   S = sum_t S_t, the S_t independent Vasicek mixed binomials at (p, rho_a): each period's
//                   PMF on 0..D is exp(log_contrib) of the binomial-mixture objective (lbinom plus the
//                   parity quadrature), and the T PMFs are convolved on 0..D only, by repeated squaring.
// Everything is in log space: each convolution term is a log-sum-exp over its pairs, shifted by its own
// maximum, so no term underflows however far the trial p is from the root.
//
// The solve: P(S <= D | p) falls as p rises (the conditional PD rises with p at every factor value), so the
// root is unique and the bound rises with D. It is found by bisection on u = logit p, from the bracket
// [kLogitLo, kLogitHi], until the bracket is at most kSolveTolLogit wide; the bound is the bracket's midpoint.
// The bracket's ends keep g(u) = log P(S <= D | p(u)) - log(1 - gamma) >= 0 at the lower end and < 0 at the
// upper, so the root lies inside it. If g does not change sign over the starting bracket, the bound is that
// end of the bracket, flagged.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "core/objectives/binomial_mixture.hpp"
#include "core/precision.hpp"
#include "core/special/lbinom.hpp"

namespace vcal::benchmark_models::pluto_tasche {

using Objective = objectives::BinomialMixture<PrecisionF64>;

// The bisection's stopping width in logit p (TOL_PT_SOLVE_LOGIT, docs/methodology/tolerances.md). The
// registration's reference finds no true PD within 8.4e-5 in logit of a bound at the total where coverage
// switches, so solver error cannot flip a coverage flag.
inline constexpr double kSolveTolLogit = 1e-10;
// The starting bracket in logit p: p from 9.4e-14 to 1 - 3.1e-7.
inline constexpr double kLogitLo = -30.0;
inline constexpr double kLogitHi = 15.0;

enum : std::uint32_t {
    kBoundBelowBracket = 1u << 0,  // g < 0 already at kLogitLo: the bound is reported as the bracket's low end
    kBoundAboveBracket = 1u << 1,  // g >= 0 still at kLogitHi: the bound is reported as the bracket's high end
    kBoundInvalidInput = 1u << 2,  // n < 1, T < 1, D < 0, D > n T, or gamma outside (0, 1): the bound is NaN
};

struct Bound {
    double pd;        // the bound, natural scale
    double u;         // its logit
    double u_lo;      // the final bracket in logit p: g(u_lo) >= 0 > g(u_hi)
    double u_hi;
    std::int32_t iterations;
    std::uint32_t flags;
};

namespace detail {

// log p and log(1 - p) from u = logit p, without forming p.
inline double log_p(double u) { return -std::log1p(std::exp(-u)); }
inline double log_q(double u) { return -std::log1p(std::exp(u)); }
inline double inv_logit(double u) { return 1.0 / (1.0 + std::exp(-u)); }

// log sum_i exp(x_i), the maximum shifted out; -inf if every term is -inf.
inline double log_sum_exp(const double* x, std::size_t m) {
    double mx = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < m; ++i) mx = std::fmax(mx, x[i]);
    if (!(mx > -std::numeric_limits<double>::infinity())) return mx;
    double s = 0.0;
    for (std::size_t i = 0; i < m; ++i) s += std::exp(x[i] - mx);
    return mx + std::log(s);
}

// c[k] = log sum_{j=0..k} exp(a[j] + b[k - j]), k = 0..D: the log PMF of the sum, truncated at D.
inline void log_convolve(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c,
                         std::vector<double>& scratch) {
    const std::size_t m = a.size();
    c.assign(m, 0.0);
    scratch.resize(m);
    for (std::size_t k = 0; k < m; ++k) {
        for (std::size_t j = 0; j <= k; ++j) scratch[j] = a[j] + b[k - j];
        c[k] = log_sum_exp(scratch.data(), k + 1);
    }
}

inline bool valid(std::int64_t n, std::int64_t T, std::int64_t D, double gamma) {
    return n >= 1 && T >= 1 && D >= 0 && D <= n * T && gamma > 0.0 && gamma < 1.0;
}

// Bisection on u for g(u) = log_cdf(u) - log(1 - gamma), g decreasing.
template <class LogCdf>
Bound solve(const LogCdf& log_cdf, double gamma) {
    const double target = std::log1p(-gamma);
    Bound out{};
    double lo = kLogitLo, hi = kLogitHi;
    if (!(log_cdf(lo) >= target)) {
        out.flags |= kBoundBelowBracket;
        hi = lo;
    } else if (log_cdf(hi) >= target) {
        out.flags |= kBoundAboveBracket;
        lo = hi;
    }
    while (hi - lo > kSolveTolLogit) {
        const double mid = 0.5 * (lo + hi);
        if (log_cdf(mid) >= target) {
            lo = mid;
        } else {
            hi = mid;
        }
        ++out.iterations;
    }
    out.u_lo = lo;
    out.u_hi = hi;
    out.u = 0.5 * (lo + hi);
    out.pd = inv_logit(out.u);
    return out;
}

inline Bound invalid() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return {nan, nan, nan, nan, 0, kBoundInvalidInput};
}

}  // namespace detail

// log P(Bin(N, p) <= D), p = logit^-1(u).
inline double log_cdf_independent(std::int64_t N, std::int64_t D, double u) {
    const double lp = detail::log_p(u), lq = detail::log_q(u);
    std::vector<double> terms(static_cast<std::size_t>(D + 1));
    for (std::int64_t k = 0; k <= D; ++k) {
        terms[static_cast<std::size_t>(k)] =
            special::lbinom(N, k) + static_cast<double>(k) * lp + static_cast<double>(N - k) * lq;
    }
    return detail::log_sum_exp(terms.data(), terms.size());
}

// Work space for log_cdf_factor, reused across calls.
struct FactorWork {
    std::vector<double> period, power, result, tmp, scratch;
};

// log P(S <= D), S the sum of T independent Vasicek mixed binomials with n obligors at (p, rho),
// p = logit^-1(u). Each period's log PMF on 0..min(D, n) is the objective's log_contrib; counts above n
// have probability 0.
template <class Integrator>
double log_cdf_factor(std::int64_t n, std::int64_t T, std::int64_t D, double u, double rho,
                      const Integrator& integrator, FactorWork& w) {
    const Objective obj{};
    const Objective::Theta th{detail::inv_logit(u), rho};
    const auto m = static_cast<std::size_t>(D + 1);
    w.period.assign(m, -std::numeric_limits<double>::infinity());
    for (std::int64_t d = 0; d <= std::min(D, n); ++d) {
        w.period[static_cast<std::size_t>(d)] = obj.log_contrib(Objective::Obs{n, d}, th, integrator);
    }
    // result = period^T by repeated squaring, every product truncated at D.
    w.power = w.period;
    bool have = false;
    for (std::int64_t e = T;;) {
        if (e & 1) {
            if (have) {
                detail::log_convolve(w.result, w.power, w.tmp, w.scratch);
                w.result.swap(w.tmp);
            } else {
                w.result = w.power;
                have = true;
            }
        }
        e >>= 1;
        if (e == 0) break;
        detail::log_convolve(w.power, w.power, w.tmp, w.scratch);
        w.power.swap(w.tmp);
    }
    return detail::log_sum_exp(w.result.data(), w.result.size());
}

// (a) The independent bound: P(Bin(n T, p) <= D) = 1 - gamma.
inline Bound independent_bound(std::int64_t n, std::int64_t T, std::int64_t D, double gamma) {
    if (!detail::valid(n, T, D, gamma)) return detail::invalid();
    const std::int64_t N = n * T;
    return detail::solve([&](double u) { return log_cdf_independent(N, D, u); }, gamma);
}

// (b) The one-factor bound at an assumed rho_a: P(S <= D | p, rho_a) = 1 - gamma.
template <class Integrator>
Bound factor_bound(std::int64_t n, std::int64_t T, std::int64_t D, double rho_a, double gamma,
                   const Integrator& integrator, FactorWork& w) {
    if (!detail::valid(n, T, D, gamma) || !(rho_a > 0.0 && rho_a < 1.0)) return detail::invalid();
    return detail::solve([&](double u) { return log_cdf_factor(n, T, D, u, rho_a, integrator, w); }, gamma);
}

}  // namespace vcal::benchmark_models::pluto_tasche
