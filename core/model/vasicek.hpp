// SPDX-License-Identifier: Apache-2.0
//
// One-factor Vasicek/ASRF model (ARCHITECTURE.md §1):
//     p(z) = Phi(x(z)),   x(z) = (Phi^-1(PD) - sqrt(rho) z) / sqrt(1 - rho),   dx/dz = -beta.
// Precomputes Phi^-1(PD), sqrt(rho) and sqrt(1 - rho) once per parameter point.
#pragma once

#include <cmath>

#include "core/precision.hpp"
#include "core/special/probit.hpp"

namespace vcal::model {

struct Vasicek1F {
    double c;                   // Phi^-1(PD)
    double sqrt_rho;            // sqrt(rho)
    double sqrt_one_minus_rho;  // sqrt(1 - rho); 1 - rho is exact for rho >= 1/2

    VCAL_HD double threshold(double z) const { return (c - sqrt_rho * z) / sqrt_one_minus_rho; }
    VCAL_HD double beta() const { return sqrt_rho / sqrt_one_minus_rho; }
};

// Requires 0 < pd < 1 and 0 < rho < 1; otherwise every field is NaN.
VCAL_HD Vasicek1F make_vasicek1f(double pd, double rho) {
    if (!(pd > 0.0 && pd < 1.0 && rho > 0.0 && rho < 1.0)) {
        const double nan = static_cast<double>(NAN);
        return {nan, nan, nan};
    }
    return {special::probit(pd), std::sqrt(rho), std::sqrt(1.0 - rho)};
}

}  // namespace vcal::model
