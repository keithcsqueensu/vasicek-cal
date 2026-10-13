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

// lambda(x) given e = erfcx(|x| / sqrt 2) and E = exp(-x^2 / 2) as exp_neg_half_square computes it (E is used only
// for x >= 0), for x not NaN (P-21, D-400): inverse_mills's own expression for x's sign, so a caller that already has
// e and E (the value's node batch, or the pair below) gets inverse_mills(x) bit for bit. For x < 0 inverse_mills
// takes erfcx(-x / sqrt 2), and -x = |x| exactly; for x >= 0 (-0 included) erfcx(x / sqrt 2), equal at -0.
// exp_neg_half_square squares |x|, so the same E serves x and -x.
VCAL_HD double inverse_mills_given_erfcx(double x, double e, double E) {
    if (x < 0.0) return constants::kSqrt2OverPi / e;
    return constants::kInvSqrt2Pi * E / (1.0 - 0.5 * e * E);
}

// lambda(x) and lambda(-x) together (P-21, D-400): bit for bit the two inverse_mills calls, from one erfcx and one
// exp pair. NaN gives {x, -x}, as the two calls do.
struct InverseMillsPair {
    double at_x;      // lambda(x)
    double at_neg_x;  // lambda(-x)
};

VCAL_HD InverseMillsPair inverse_mills_pair(double x) {
    if (std::isnan(x)) return {x, -x};
    const double e = erfcx(std::fabs(x) * constants::kInvSqrt2);
    const double E = detail::exp_neg_half_square(x);
    return {inverse_mills_given_erfcx(x, e, E), inverse_mills_given_erfcx(-x, e, E)};
}

}  // namespace vcal::special
