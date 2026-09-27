// SPDX-License-Identifier: Apache-2.0
//
// Accurate x^2 for the Gaussian exponent. For |x| < 2^16, x^2 = hi + lo where
// xh = |x| truncated to 10 fraction bits has at most 26 significant bits, so
// hi = xh*xh is exact and lo = (|x| - xh)(|x| + xh) is small (|x| - xh is exact).
// Rounding therefore touches only the small part, which keeps exp(+-x^2) and
// log Phi accurate to a few ulp where a plain x*x would cost ~x^2 ulp.
// Beyond 2^16 the exponentials over/underflow anyway, so hi = x*x, lo = 0.
#pragma once

#include <cmath>

#include "core/precision.hpp"

namespace vcal::special::detail {

struct Square {
    double hi;
    double lo;
};

VCAL_HD Square split_square(double x) {
    const double ax = std::fabs(x);
    if (!(ax < 65536.0)) return {ax * ax, 0.0};  // also NaN and inf
    const double xh = std::trunc(ax * 1024.0) / 1024.0;
    return {xh * xh, (ax - xh) * (ax + xh)};
}

// exp(x^2)
VCAL_HD double exp_square(double x) {
    const Square s = split_square(x);
    return std::exp(s.hi) * std::exp(s.lo);
}

// exp(-x^2 / 2)
VCAL_HD double exp_neg_half_square(double x) {
    const Square s = split_square(x);
    return std::exp(-0.5 * s.hi) * std::exp(-0.5 * s.lo);
}

}  // namespace vcal::special::detail
