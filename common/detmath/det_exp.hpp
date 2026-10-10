// SPDX-License-Identifier: Apache-2.0
//
// det_exp: e^x from operations IEEE 754 makes exact or correctly rounded (+ - * /, floor, ldexp), so it is
// identical on every platform and to the Python mirror in tools/gen_dgp_tables.py (D-052, D-053). Within 2 ulp of
// mpmath over its range (TOL_DGP_EXP_ULP; tests/golden/dgp/exp.csv).
//
// Shared, not copied (D-199): the DGP draws with it (dgp/det_math.cpp) and the AR(1) filter's transition
// accumulates with it (benchmark_models/ar1/), so one definition serves both. common/ includes nothing else in the
// repository, and dgp/, core/ and benchmark_models/ may include it (cmake/LayeringCheck.cmake). The vector
// versions (benchmark_models/ar1/transition_*.cpp) use the same data and the same operations, lane by lane, and are
// tested bit for bit against this one.
//
//   e^x = 2^k e^r,  k = floor(x / ln 2 + 1/2),  r = (x - k ln2_hi) - k ln2_lo  (Cody-Waite),
//   e^r by its Taylor series to r^16 (|r| <= 0.347), by Horner's rule, scaled by ldexp (correctly rounded).
//
// Compiled with FP contraction off (cmake/VcalTargetDefaults.cmake): a fused multiply-add would round differently.
#pragma once

#include <cmath>
#include <limits>

#include "common/detmath/generated/det_exp_constants.hpp"

namespace vcal::detmath {

inline double det_exp(double x) {
    namespace gc = generated;
    if (std::isnan(x)) return x;
    if (x > 800.0) return std::numeric_limits<double>::infinity();
    if (x < -800.0) return 0.0;
    const double kd = std::floor(x * gc::kInvLn2 + 0.5);
    const double r = (x - kd * gc::kLn2Hi) - kd * gc::kLn2Lo;
    double p = gc::kExpCoef[gc::kExpTerms];
    for (int n = gc::kExpTerms - 1; n >= 0; --n) p = p * r + gc::kExpCoef[n];
    return std::ldexp(p, static_cast<int>(kd));
}

}  // namespace vcal::detmath
