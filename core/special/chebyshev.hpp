// SPDX-License-Identifier: Apache-2.0
//
// Clenshaw evaluation of f(t) = sum_k c[k] T_k(t) on [-1, 1], with c[0] already halved.
// tools/gen_special_tables.py measures fit error with this exact operation order.
#pragma once

#include <cstddef>

#include "core/precision.hpp"

namespace vcal::special {

template <std::size_t N>
VCAL_HD double chebyshev_eval(const double (&c)[N], double t) {
    static_assert(N >= 2, "chebyshev_eval needs at least two coefficients");
    double b1 = 0.0;
    double b2 = 0.0;
    for (std::size_t k = N - 1; k > 0; --k) {
        const double b0 = c[k] + 2.0 * t * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return c[0] + t * b1 - b2;
}

}  // namespace vcal::special
