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

// log Phi(x) and log Phi(-x) together (P-7, D-199): bit for bit the two log_phi calls, with erfcx(|x|/sqrt2) and the
// square split computed once. Both calls need the same two (for x < 0, -x/sqrt2 = |x|/sqrt2 exactly, and
// split_square takes |x|), and each value is then finished with log_phi's own expression for its sign. At x = +-0
// both arguments are >= 0 in log_phi's sense, so both take the x >= 0 branch, as here. NaN gives {x, -x}, as the
// two calls do.
struct LogPhiPair {
    double at_x;      // log Phi(x)
    double at_neg_x;  // log Phi(-x)
};

VCAL_HD LogPhiPair log_phi_pair(double x) {
    if (std::isnan(x)) return {x, -x};
    const double ax = std::fabs(x);
    const double e = erfcx(ax * constants::kInvSqrt2);
    const detail::Square s = detail::split_square(x);
    const double upper = std::log1p(-0.5 * e * (std::exp(-0.5 * s.hi) * std::exp(-0.5 * s.lo)));  // log Phi(|x|)
    if (x == 0.0) return {upper, upper};
    const double lower = ((std::log(e) - constants::kLn2) - 0.5 * s.lo) - 0.5 * s.hi;  // log Phi(-|x|)
    return x < 0.0 ? LogPhiPair{lower, upper} : LogPhiPair{upper, lower};
}

}  // namespace vcal::special
