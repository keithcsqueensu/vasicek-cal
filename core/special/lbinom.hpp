// SPDX-License-Identifier: Apache-2.0
//
// log C(n, k) and the Stirling error it is built on (D-070).
//
// stirling_error(n) = log n! - [(n + 1/2) log n - n + log sqrt(2 pi)], n >= 1 (Loader 2000):
//   mpmath-generated table for n <= 30, series beyond (generated/stirling_error.hpp).
//
// lbinom(n, k) for 0 <= k <= n <= 2^53 (NaN otherwise):
//   k = 0 or k = n : exactly 0 (the general form would evaluate 0 * inf).
//   n <= 60        : log of the exact integer C(n, k), which fits in uint64.
//   otherwise, with k <= n/2 by symmetry and m = n - k:
//     log C(n, k) = [delta(n) - delta(k) - delta(m)] + 1/2 log(n / (k m)) - log sqrt(2 pi)
//                   + k log(n / k) - m log1p(-k / n)
//   No step subtracts quantities of size log n!, so the error scales with the result.
//   lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) does not: at n ~ 1e6 its absolute error is
//   ~1e-9 whatever the result. (std::lgamma is also not thread-safe in glibc, and differs
//   between host and CUDA.)
// Accuracy: TOL_LBINOM_ULP, TOL_STIRLING_ERROR_ULP.
#pragma once

#include <cmath>
#include <cstdint>

#include "core/precision.hpp"
#include "core/special/generated/constants.hpp"
#include "core/special/generated/stirling_error.hpp"

namespace vcal::special {

inline constexpr std::int64_t kLbinomExactMax = 60;
inline constexpr std::int64_t kLbinomMaxN = std::int64_t{1} << 53;

VCAL_HD double stirling_error(std::int64_t n) {
    if (n < 1) return static_cast<double>(NAN);
    if (n <= generated::kStirlingTableMax) return generated::stirling_error_table(static_cast<int>(n));
    return generated::stirling_error_series(1.0 / static_cast<double>(n));
}

VCAL_HD double lbinom(std::int64_t n, std::int64_t k) {
    if (n < 0 || k < 0 || k > n || n > kLbinomMaxN) return static_cast<double>(NAN);
    if (n - k < k) k = n - k;
    if (k == 0) return 0.0;

    if (n <= kLbinomExactMax) {
        // c_i = C(n - k + i, i); the product c_{i-1} * (n - k + i) divides exactly by i. Its
        // maximum over n <= 60, k <= n/2 is 3.55e18 (n = 60, k = i = 30; ~19% of 2^64), by
        // exhaustive enumeration (DECISIONS.md D-087).
        std::uint64_t c = 1;
        for (std::int64_t i = 1; i <= k; ++i) {
            c = c * static_cast<std::uint64_t>(n - k + i) / static_cast<std::uint64_t>(i);
        }
        return std::log(static_cast<double>(c));
    }

    const std::int64_t m = n - k;
    const double nd = static_cast<double>(n);
    const double kd = static_cast<double>(k);
    const double md = static_cast<double>(m);
    const double corrections = (stirling_error(n) - stirling_error(k) - stirling_error(m)) +
                               (0.5 * std::log(nd / (kd * md)) - constants::kLnSqrt2Pi);
    return corrections + (kd * std::log(nd / kd) - md * std::log1p(-kd / nd));
}

}  // namespace vcal::special
