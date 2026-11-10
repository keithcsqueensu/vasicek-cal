// SPDX-License-Identifier: Apache-2.0
//
// The multi-grade binomial-mixture period log-likelihood (M6, D-151, D-411; docs/methodology/multigrade_mle.md):
//
//   l_t(PD_1, ..., PD_K, rho) = sum_k log C(n_k, d_k) + log E[ prod_k p_k(Z)^d_k (1 - p_k(Z))^(n_k - d_k) ],
//
// one integral per period over the product of the grades' binomials, with the integrator's rule chosen by the
// period's counts (core/model/multigrade_mixture.hpp). The binomial coefficients are kept (D-070). A grade with n = 0
// contributes nothing. At K = 1 every function here is BinomialMixture's bit for bit.
//
// The estimator over a panel (the rho grid, the inner Newton over the PDs) is PR 2's (D-411); this header is the
// likelihood and its derivatives only.
#pragma once

#include <cmath>
#include <cstdint>

#include "core/model/multigrade_mixture.hpp"
#include "core/precision.hpp"
#include "core/quadrature/integrator.hpp"
#include "core/special/lbinom.hpp"

namespace vcal::objectives {

using model::GradeCount;
using model::kMaxGrades;

// nullptr if 1 <= K <= kMaxGrades and every grade has 0 <= n <= 2^53 and 0 <= d <= n, with at least one n > 0.
inline const char* multigrade_period_error(const GradeCount* y, int K) {
    if (!(K >= 1 && K <= kMaxGrades)) return "a period needs 1 <= K <= kMaxGrades grades";
    bool any = false;
    for (int k = 0; k < K; ++k) {
        if (!(y[k].n >= 0 && y[k].n <= special::kLbinomMaxN)) return "each grade needs 0 <= n <= 2^53";
        if (!(y[k].d >= 0 && y[k].d <= y[k].n)) return "each grade needs 0 <= d <= n";
        any = any || y[k].n > 0;
    }
    return any ? nullptr : "a period needs at least one grade with n > 0";
}

// sum_k log C(n_k, d_k), over the non-empty grades in index order.
VCAL_HD double multigrade_lbinom(const GradeCount* y, int K) {
    double lc = 0.0;
    for (int k = 0; k < K; ++k) {
        if (y[k].n == 0) continue;
        lc += special::lbinom(y[k].n, y[k].d);
    }
    return lc;
}

// The period's log-likelihood. NaN if a non-empty grade's PD or rho is outside (0, 1), or K is out of range.
template <class I>
VCAL_HD double multigrade_log_contrib(const GradeCount* y, int K, const double* pd, double rho, const I& integrator) {
    const auto f = model::make_multigrade_log_integrand(y, pd, K, rho);
    if (f.grades < 0) return static_cast<double>(NAN);
    return multigrade_lbinom(y, K) + integrator.log_integrate(f, model::multigrade_mixture_hint(f));
}

// The contribution with its score and Hessian in theta = (PD_1, ..., PD_K, rho): grad[0..K], and hess, the full
// (K + 1) x (K + 1) symmetric matrix row-major. From the same nodes and fold as multigrade_log_contrib, so the value
// equals it bit for bit, and at K = 1 grad and hess equal BinomialMixture::log_contrib_derivs's bit for bit. A grade
// with n = 0 has zero score and zero Hessian row and column. Needs an integrator with log_integrate_with.
template <class I>
VCAL_HD double multigrade_log_contrib_derivs(const GradeCount* y, int K, const double* pd, double rho,
                                             const I& integrator, double* grad, double* hess) {
    const int n = K + 1;
    for (int i = 0; i < n; ++i) grad[i] = 0.0;
    for (int i = 0; i < n * n; ++i) hess[i] = 0.0;
    const auto f = model::make_multigrade_log_integrand(y, pd, K, rho);
    const auto df = model::make_multigrade_log_integrand_derivs(y, pd, K, rho);
    if (f.grades < 0 || df.grades < 0) {
        for (int i = 0; i < n; ++i) grad[i] = static_cast<double>(NAN);
        for (int i = 0; i < n * n; ++i) hess[i] = static_cast<double>(NAN);
        return static_cast<double>(NAN);
    }
    constexpr int kMaxMoments = (kMaxGrades + 1) + (kMaxGrades + 1) * (kMaxGrades + 2) / 2;
    const int Kd = df.grades;  // the non-empty grades
    quadrature::detail::OnlineLogSumExpMomentsN<model::MultiGradeLogIntegrandDerivs, kMaxMoments> acc(
        df, model::MultiGradeLogIntegrandDerivs::moments(Kd));
    const double l = multigrade_lbinom(y, K) + integrator.log_integrate_with(f, model::multigrade_mixture_hint(f), acc);

    // Map the non-empty grades back to their indices; rho is index K.
    int idx[kMaxGrades + 1];
    int j = 0;
    for (int k = 0; k < K; ++k) {
        if (y[k].n != 0) idx[j++] = k;
    }
    idx[Kd] = K;
    double m[kMaxGrades + 1];
    for (int a = 0; a <= Kd; ++a) {
        m[a] = acc.mean(a);
        grad[idx[a]] = m[a];
    }
    for (int a = 0; a <= Kd; ++a) {
        for (int b = a; b <= Kd; ++b) {
            const double h = acc.mean(model::MultiGradeLogIntegrandDerivs::second(Kd, a, b)) - m[a] * m[b];
            hess[idx[a] * n + idx[b]] = h;
            hess[idx[b] * n + idx[a]] = h;
        }
    }
    return l;
}

}  // namespace vcal::objectives
