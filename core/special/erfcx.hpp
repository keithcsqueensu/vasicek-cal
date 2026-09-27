// SPDX-License-Identifier: Apache-2.0
//
// Scaled complementary error function erfcx(x) = exp(x^2) erfc(x). In-house because the
// host standard library has none (ARCHITECTURE.md §3.2).
//   x >= 0 : Chebyshev series fitted in mpmath (generated/erfcx_chebyshev.hpp); four
//            pieces, the last in 8/x so the series stays well conditioned as x -> inf.
//   x <  0 : erfcx(x) = 2 exp(x^2) - erfcx(-x); overflows to +inf below about -26.6.
// Accuracy: TOL_ERFCX_ULP (tests/tolerances.hpp, docs/methodology/tolerances.md).
#pragma once

#include <cmath>

#include "core/precision.hpp"
#include "core/special/generated/erfcx_chebyshev.hpp"
#include "core/special/square_split.hpp"

namespace vcal::special {

VCAL_HD double erfcx(double x) {
    if (std::isnan(x)) return x;
    if (x < 0.0) return 2.0 * detail::exp_square(x) - generated::erfcx_nonnegative(-x);
    if (x == HUGE_VAL) return 0.0;
    return generated::erfcx_nonnegative(x);
}

}  // namespace vcal::special
