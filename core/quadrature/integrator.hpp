// SPDX-License-Identifier: Apache-2.0
//
// Integrator contract (ARCHITECTURE.md §3.4). An integrator computes
//     log of the integral of exp(g(z)) phi(z) dz
// for a log-integrand g, entirely in log space. It exposes
//     accum_t log_integrate(const G& g, IntegrandHint hint) const
// where g is callable as double(double). C++17 trait + per-clause static_asserts (D-063).
#pragma once

#include <cmath>
#include <type_traits>
#include <utility>

#include "core/precision.hpp"

namespace vcal::quadrature {

// Where exp(g(z)) phi(z) peaks: the mode of h(z) = g(z) - z^2/2 and 1/sqrt(-h''(mode)).
// Fixed rules ignore it; adaptive rules centre and scale on it. For one-sided integrands (the
// prior truncated by a steep survival factor, D-118) `center` and `width` locate the factor's
// transition, and `one_sided` selects the composite rule in a SplitRule. Members left out of
// an initialiser are zero: center/width unused, one_sided false.
struct IntegrandHint {
    double mode;
    double scale;
    double center;
    double width;
    bool one_sided;
};

// Standard-normal rule: E[f(Z)] ~ sum_i exp(log_weight[i]) f(node[i]). Non-owning;
// n == 0 means "no rule".
struct GaussHermiteRule {
    const double* node;
    const double* log_weight;
    int n;
};

namespace detail {

// One-pass log-sum-exp, terms folded in call order (so results are deterministic).
// -inf terms contribute nothing; any NaN term gives NaN; any +inf term gives +inf.
class OnlineLogSumExp {
public:
    VCAL_HD void add(double v) {
        if (std::isnan(v)) {
            nan_ = true;
        } else if (v == HUGE_VAL) {
            inf_ = true;
        } else if (v == -HUGE_VAL) {
            // exp(-inf) = 0
        } else if (v > max_) {
            sum_ = sum_ * std::exp(max_ - v) + 1.0;
            max_ = v;
        } else {
            sum_ += std::exp(v - max_);
        }
    }

    VCAL_HD double result() const {
        if (nan_) return static_cast<double>(NAN);
        if (inf_) return HUGE_VAL;
        if (max_ == -HUGE_VAL) return -HUGE_VAL;
        return max_ + std::log(sum_);
    }

private:
    double max_ = -HUGE_VAL;
    double sum_ = 0.0;
    bool nan_ = false;
    bool inf_ = false;
};

struct LogIntegrandArchetype {
    VCAL_HD double operator()(double) const { return 0.0; }
};

template <class I, class = void>
struct has_log_integrate : std::false_type {};
template <class I>
struct has_log_integrate<I, std::void_t<decltype(std::declval<const I&>().log_integrate(
                                std::declval<const LogIntegrandArchetype&>(), std::declval<IntegrandHint>()))>>
    : std::true_type {};

template <class I>
using log_integrate_result_t = decltype(std::declval<const I&>().log_integrate(
    std::declval<const LogIntegrandArchetype&>(), std::declval<IntegrandHint>()));

}  // namespace detail

template <class I, class = void>
struct is_integrator : std::false_type {};
template <class I>
struct is_integrator<I, std::enable_if_t<detail::has_log_integrate<I>::value>>
    : std::is_same<detail::log_integrate_result_t<I>, double> {};
template <class I>
inline constexpr bool is_integrator_v = is_integrator<I>::value;

template <class I>
constexpr bool check_integrator() {
    static_assert(detail::has_log_integrate<I>::value,
                  "vcal: an Integrator must provide log_integrate(g, IntegrandHint) const, "
                  "callable with a double(double) log-integrand g");
    if constexpr (detail::has_log_integrate<I>::value) {
        static_assert(std::is_same_v<detail::log_integrate_result_t<I>, double>,
                      "vcal: Integrator::log_integrate must return accum_t (double, D-039)");
    }
    return true;
}

}  // namespace vcal::quadrature
