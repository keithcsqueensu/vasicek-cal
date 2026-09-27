// SPDX-License-Identifier: Apache-2.0
//
// log(exp(a) + exp(b)) without overflow or underflow. Zero-probability terms arrive as
// -inf routinely, and -inf - (-inf) is NaN, so the infinities are handled first (D-071):
// either input -inf returns the other (both -inf -> -inf); either +inf -> +inf.
// Equal arguments give a + log 2.
// Accuracy: TOL_LOG_ADD_EXP_ULP.
#pragma once

#include <cmath>

#include "core/precision.hpp"

namespace vcal::special {

VCAL_HD double log_add_exp(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return a + b;
    const double hi = a > b ? a : b;
    const double lo = a > b ? b : a;
    if (lo == -HUGE_VAL || hi == HUGE_VAL) return hi;
    return hi + std::log1p(std::exp(lo - hi));
}

}  // namespace vcal::special
