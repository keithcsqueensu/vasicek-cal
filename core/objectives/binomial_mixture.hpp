// SPDX-License-Identifier: Apache-2.0
//
// Binomial-mixture objective (ARCHITECTURE.md §1, §3.5): the per-period log-likelihood
//     l_t(PD, rho) = log C(n_t, d_t) + log E[p(Z)^d_t (1 - p(Z))^(n_t - d_t)]
// with the expectation from the integrator (parity: the split rule of D-118).
// log C(n, d) is constant in the parameters but kept, so absolute log-likelihoods, AIC/BIC
// and goldens are the true values (D-070).
//
// Observations must satisfy 0 < n <= 2^53 and 0 <= d <= n (panel_error checks a panel).
#pragma once

#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "core/model/binomial_mixture.hpp"
#include "core/model/vasicek.hpp"
#include "core/precision.hpp"
#include "core/quadrature/integrator.hpp"
#include "core/special/lbinom.hpp"

namespace vcal::objectives {

// nullptr if every period has n > 0 and 0 <= d <= n (and n <= 2^53, the lbinom domain).
template <class Obs>
const char* binomial_panel_error(const Obs* obs, std::int64_t periods) {
    if (periods < 1) return "panel has no periods";
    for (std::int64_t t = 0; t < periods; ++t) {
        if (!(obs[t].n > 0 && obs[t].n <= special::kLbinomMaxN)) return "each period needs 0 < n <= 2^53";
        if (!(obs[t].d >= 0 && obs[t].d <= obs[t].n)) return "each period needs 0 <= d <= n";
    }
    return nullptr;
}

template <class P>
struct BinomialMixture {
    static_assert(check_precision<P>());

    struct Obs {
        std::int64_t n;
        std::int64_t d;
        friend VCAL_HD constexpr bool operator==(const Obs& a, const Obs& b) { return a.n == b.n && a.d == b.d; }
    };
    struct Theta {
        double pd;
        double rho;
    };
    static constexpr int n_params = 2;  // grid axis 0 = PD, axis 1 = rho

    static VCAL_HD Theta theta(const double (&v)[2]) { return {v[0], v[1]}; }
    static const char* panel_error(const Obs* obs, std::int64_t periods) { return binomial_panel_error(obs, periods); }

    // Size of the terms of l_t = log C(n, d) + log I: the scale of its rounding error. The two
    // terms can nearly cancel (n = 1e6, d = 3e5: each ~6e5, sum ~ -15), so |l_t| understates it.
    // Used by the engine's quadrature check (D-120).
    static VCAL_HD double rounding_scale(const Obs& y, double l) {
        const double lc = special::lbinom(y.n, y.d);
        return std::fabs(lc) + std::fabs(l - lc);
    }

    template <class I>
    VCAL_HD typename P::accum_t log_contrib(const Obs& y, const Theta& th, const I& integrator) const {
        const auto f = model::make_binomial_log_integrand(model::make_vasicek1f(th.pd, th.rho), y.n, y.d);
        return special::lbinom(y.n, y.d) + integrator.log_integrate(f, model::binomial_mixture_hint(f));
    }
};

// --- contract (C++17 trait + per-clause static_asserts, D-063) -------------------------------

namespace detail {

struct IntegratorArchetype {
    template <class G>
    VCAL_HD double log_integrate(const G&, quadrature::IntegrandHint) const {
        return 0.0;
    }
};

template <class O, class = void>
struct has_objective_members : std::false_type {};
template <class O>
struct has_objective_members<
    O, std::void_t<typename O::Obs, typename O::Theta, decltype(O::n_params),
                   decltype(std::declval<const O&>().log_contrib(std::declval<const typename O::Obs&>(),
                                                                 std::declval<const typename O::Theta&>(),
                                                                 std::declval<const IntegratorArchetype&>()))>>
    : std::true_type {};

// Periods with equal observations share one surface evaluation (D-122), so Obs needs ==.
template <class O, class = void>
struct has_obs_equality : std::false_type {};
template <class O>
struct has_obs_equality<O, std::void_t<decltype(std::declval<const typename O::Obs&>() ==
                                                std::declval<const typename O::Obs&>())>> : std::true_type {};

}  // namespace detail

template <class O>
constexpr bool check_objective() {
    static_assert(detail::has_objective_members<O>::value,
                  "vcal: an Objective must declare Obs, Theta, n_params and "
                  "log_contrib(const Obs&, const Theta&, const Integrator&) const");
    if constexpr (detail::has_objective_members<O>::value) {
        static_assert(std::is_same_v<decltype(std::declval<const O&>().log_contrib(
                                         std::declval<const typename O::Obs&>(),
                                         std::declval<const typename O::Theta&>(),
                                         std::declval<const detail::IntegratorArchetype&>())),
                                     double>,
                      "vcal: Objective::log_contrib must return accum_t (double, D-039)");
        static_assert(detail::has_obs_equality<O>::value,
                      "vcal: Objective::Obs must be equality-comparable (bool operator==), D-122");
    }
    return true;
}

static_assert(check_objective<BinomialMixture<PrecisionF64>>());

}  // namespace vcal::objectives
