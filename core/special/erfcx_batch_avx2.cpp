// SPDX-License-Identifier: Apache-2.0
//
// The AVX2 erfcx kernel (P-7, D-212): erfcx_lanes.hpp on four lanes. This file alone is compiled for AVX2, and
// erfcx_batch.cpp calls it only after checking at run time that the CPU and the OS support AVX2. It therefore
// includes nothing with external-linkage inline functions (no standard library beyond <cstddef>); see
// erfcx_lanes.hpp. Whole blocks of four only: erfcx_batch.cpp finishes any remainder with the scalar kernel.
#include <cstddef>

#if defined(VCAL_ERFCX_AVX2)

#include <immintrin.h>

#include "core/special/erfcx_lanes.hpp"

namespace vcal::special::erfcx_batch {
namespace {

struct V256 {
    using T = __m256d;
    using M = __m256d;
    static T set1(double v) { return _mm256_set1_pd(v); }
    static T add(T a, T b) { return _mm256_add_pd(a, b); }
    static T sub(T a, T b) { return _mm256_sub_pd(a, b); }
    static T mul(T a, T b) { return _mm256_mul_pd(a, b); }
    static T div(T a, T b) { return _mm256_div_pd(a, b); }
    static M le(T a, T b) { return _mm256_cmp_pd(a, b, _CMP_LE_OQ); }
    static T blend(M m, T a, T b) { return _mm256_blendv_pd(b, a, m); }  // m ? a : b
    static bool all(M m) { return _mm256_movemask_pd(m) == 15; }
    static bool none(M m) { return _mm256_movemask_pd(m) == 0; }
};

}  // namespace

std::size_t nonnegative_avx2(const double* x, double* out, std::size_t n) {
    constexpr int kB = 4;  // four vectors of four interleaved (erfcx_lanes.hpp)
    std::size_t i = 0;
    for (; i + 4 * kB <= n; i += 4 * kB) {
        __m256d in[kB], res[kB];
        for (int j = 0; j < kB; ++j) in[j] = _mm256_loadu_pd(x + i + 4 * j);
        lanes_erfcx_nonnegative<V256, kB>(in, res);
        for (int j = 0; j < kB; ++j) _mm256_storeu_pd(out + i + 4 * j, res[j]);
    }
    for (; i + 4 <= n; i += 4) {
        __m256d in[1] = {_mm256_loadu_pd(x + i)}, res[1];
        lanes_erfcx_nonnegative<V256, 1>(in, res);
        _mm256_storeu_pd(out + i, res[0]);
    }
    return i;
}

}  // namespace vcal::special::erfcx_batch

#endif
