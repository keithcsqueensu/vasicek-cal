// SPDX-License-Identifier: Apache-2.0
//
// The compensated sum over periods behind every grid reduction (engine/surface.hpp, D-008, D-009), and the same
// sum on a run of adjacent grid points at once (P-6, D-404).
//
// weighted_sum: sum_t w[t] * L[t * K + k] in ascending t, Neumaier-compensated; periods with w[t] == 0 are skipped,
// so a zero weight never meets a non-finite L. The single definition every reducer and refinement uses, so they all
// see the same values bit for bit. (Moved here from engine/surface.hpp by P-6, unchanged; engine:: keeps the name.)
//
// weighted_sums: out[i] = weighted_sum(L, T, K, w, k0 + i) for i < n, bit for bit. Outside device code it runs the
// vector kernels of sum_batch.hpp (one lane per grid point, each lane the scalar operations in their order), and
// otherwise, or if no kernel is selected, the scalar loop.
#pragma once

#include <cmath>
#include <cstdint>

#include "core/precision.hpp"
#include "core/reducers/sum_batch.hpp"

namespace vcal::reducers {

VCAL_HD double weighted_sum(const double* L, std::int64_t T, std::int64_t K, const double* w, std::int64_t k) {
    double sum = 0.0;
    double comp = 0.0;
    double nonfinite = 0.0;  // +-inf / NaN terms bypass compensation, which they would poison
    bool any_nonfinite = false;
    for (std::int64_t t = 0; t < T; ++t) {
        if (w[t] == 0.0) continue;
        const double x = w[t] * L[t * K + k];
        if (!std::isfinite(x)) {
            nonfinite += x;
            any_nonfinite = true;
            continue;
        }
        const double s = sum + x;
        comp += std::fabs(sum) >= std::fabs(x) ? (sum - s) + x : (x - s) + sum;
        sum = s;
    }
    return any_nonfinite ? nonfinite + (sum + comp) : sum + comp;
}

VCAL_HD void weighted_sums(const double* L, std::int64_t T, std::int64_t K, const double* w, std::int64_t k0,
                           std::int64_t n, double* out) {
#if !defined(__CUDA_ARCH__)
    if (sum_batch::weighted_sums(L, T, K, w, k0, n, out)) return;
#endif
    for (std::int64_t i = 0; i < n; ++i) out[i] = weighted_sum(L, T, K, w, k0 + i);
}

}  // namespace vcal::reducers
