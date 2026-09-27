// SPDX-License-Identifier: Apache-2.0
//
// Jackknife-based estimates and intervals for the shared jackknife run (S-3, S-5, S-21; D-155).
// `native` options, reported beside the parity results: nothing here changes a parity estimate or
// interval. The delete-one estimates come from W x L with jackknife_weights (weights.hpp) and
// replicate_estimates (bootstrap.hpp), i.e. the grid refinement on each reweighted surface.
//
// delete_two_weights(T)        the T(T - 1)/2 leave-two-out rows, pairs (s, t) with s < t in
//                              lexicographic order: W[row][u] = 1 - [u == s] - [u == t].
// jackknife_bias_corrected     Quenouille's correction T theta_hat - (T - 1) mean_t theta_(-t).
// bca_interval                 Efron's bias-corrected and accelerated interval (Efron 1987; Efron
//                              and Tibshirani 1993, ch. 14), in whatever coordinate the inputs use:
//     z0 = Phi^-1((#{b : x_b < theta_hat} + #{b : x_b = theta_hat} / 2) / B'), B' finite replicates;
//     a  = sum_t d_t^3 / (6 (sum_t d_t^2)^(3/2)),  d_t = mean_s theta_(-s) - theta_(-t);
//     alpha_1,2 = Phi(z0 + (z0 + z) / (1 - a (z0 + z))),  z = -+ Phi^-1((1 + level) / 2);
//     the ends are type-7 quantiles (as percentile_interval) of the finite replicates at alpha_1,2.
//   Not computed (computed = false, NaN ends) when z0 is infinite (every replicate on one side of
//   theta_hat), the jackknife has no spread, or 1 - a (z0 + z) <= 0 for either end.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "core/special/log_phi.hpp"
#include "core/special/probit.hpp"
#include "resample/bootstrap.hpp"

namespace vcal::resample {

inline std::vector<double> delete_two_weights(std::int64_t periods) {
    if (periods < 3) throw std::invalid_argument("the leave-two-out jackknife needs at least 3 periods");
    std::vector<double> w;
    w.reserve(static_cast<std::size_t>(periods * (periods - 1) / 2 * periods));
    for (std::int64_t s = 0; s < periods; ++s) {
        for (std::int64_t t = s + 1; t < periods; ++t) {
            for (std::int64_t u = 0; u < periods; ++u) w.push_back(u == s || u == t ? 0.0 : 1.0);
        }
    }
    return w;
}

inline double jackknife_bias_corrected(double theta_hat, const double* theta_minus, std::int64_t periods) {
    double sum = 0.0;
    for (std::int64_t t = 0; t < periods; ++t) sum += theta_minus[t];
    const double T = static_cast<double>(periods);
    return T * theta_hat - (T - 1.0) * (sum / T);
}

struct BcaInterval {
    double lo;
    double hi;
    double z0;
    double acceleration;
    bool computed;
};

inline BcaInterval bca_interval(const double* replicates, std::int64_t count, double theta_hat,
                                const double* theta_minus, std::int64_t periods, double level = 0.95) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    BcaInterval out{nan, nan, nan, nan, false};
    std::vector<double> x;
    x.reserve(static_cast<std::size_t>(count));
    double below = 0.0;
    for (std::int64_t b = 0; b < count; ++b) {
        const double v = replicates[b];
        if (!std::isfinite(v)) continue;
        x.push_back(v);
        below += v < theta_hat ? 1.0 : v == theta_hat ? 0.5 : 0.0;
    }
    if (x.empty() || !std::isfinite(theta_hat)) return out;
    out.z0 = special::probit(below / static_cast<double>(x.size()));
    double mean = 0.0;
    for (std::int64_t t = 0; t < periods; ++t) mean += theta_minus[t];
    mean /= static_cast<double>(periods);
    double s2 = 0.0, s3 = 0.0;
    for (std::int64_t t = 0; t < periods; ++t) {
        const double d = mean - theta_minus[t];
        s2 += d * d;
        s3 += d * d * d;
    }
    out.acceleration = s3 / (6.0 * std::pow(s2, 1.5));
    if (!std::isfinite(out.z0) || !(s2 > 0.0) || !std::isfinite(out.acceleration)) return out;
    const double z = special::probit(0.5 * (1.0 + level));
    double alpha[2];
    for (int k = 0; k < 2; ++k) {
        const double zk = out.z0 + (k == 0 ? -z : z);
        const double denom = 1.0 - out.acceleration * zk;
        if (!(denom > 0.0)) return out;
        alpha[k] = std::exp(special::log_phi(out.z0 + zk / denom));
    }
    std::sort(x.begin(), x.end());
    out.lo = quantile_type7(x, alpha[0]);
    out.hi = quantile_type7(x, alpha[1]);
    out.computed = true;
    return out;
}

}  // namespace vcal::resample
