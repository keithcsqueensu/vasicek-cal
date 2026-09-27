// SPDX-License-Identifier: Apache-2.0
//
// log Phi(x), the log of the standard normal CDF, on the whole real line. This is the
// core primitive: log p(z) = log_phi(x) and log(1 - p(z)) = log_phi(-x) (ARCHITECTURE.md §1).
//
//   x < 0  : log Phi(x) = log erfcx(-x/sqrt2) - log 2 - x^2/2. All three terms have the same
//            sign, so nothing cancels; x^2 is split so its large part is exact.
//   x >= 0 : log Phi(x) = log1p(-Phi(-x)), Phi(-x) = erfcx(x/sqrt2) exp(-x^2/2) / 2 (D-068).
//            This keeps relative accuracy as log Phi(x) -> 0, where log(Phi(x)) would not.
//            Once the result is subnormal (x above about 37.5) only absolute accuracy is
//            meaningful, and it reaches -0.0 near x = 38.5.
// Accuracy: TOL_LOG_PHI_ULP, TOL_LOG_PHI_SUBNORMAL_ABS.
#pragma once

#include <cmath>

#include "core/precision.hpp"
#include "core/special/erfcx.hpp"
#include "core/special/generated/constants.hpp"
#include "core/special/square_split.hpp"

namespace vcal::special {

VCAL_HD double log_phi(double x) {
    if (std::isnan(x)) return x;
    if (x < 0.0) {
        const detail::Square s = detail::split_square(x);
        return ((std::log(erfcx(-x * constants::kInvSqrt2)) - constants::kLn2) - 0.5 * s.lo) -
               0.5 * s.hi;
    }
    return std::log1p(-0.5 * erfcx(x * constants::kInvSqrt2) * detail::exp_neg_half_square(x));
}

}  // namespace vcal::special
