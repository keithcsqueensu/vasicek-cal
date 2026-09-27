// SPDX-License-Identifier: Apache-2.0
//
// Inverse Mills ratio lambda(x) = phi(x) / Phi(x) = d/dx log Phi(x), used for the slope and
// curvature of binomial-mixture log-integrands (core/model/binomial_mixture.hpp).
//   x < 0  : lambda(x) = sqrt(2/pi) / erfcx(-x/sqrt2). Exact identity, so no cancellation
//            between phi and Phi however deep the tail; lambda(x) ~ -x as x -> -inf.
//   x >= 0 : phi(x) / (1 - Phi(-x)), where Phi(x) >= 1/2 so the division is benign.
// Accuracy: TOL_INVERSE_MILLS_ULP.
#pragma once

#include <cmath>

#include "core/precision.hpp"
#include "core/special/erfcx.hpp"
#include "core/special/generated/constants.hpp"
#include "core/special/square_split.hpp"

namespace vcal::special {

VCAL_HD double inverse_mills(double x) {
    if (std::isnan(x)) return x;
    if (x < 0.0) return constants::kSqrt2OverPi / erfcx(-x * constants::kInvSqrt2);
    const double e = detail::exp_neg_half_square(x);
    return constants::kInvSqrt2Pi * e / (1.0 - 0.5 * erfcx(x * constants::kInvSqrt2) * e);
}

}  // namespace vcal::special
