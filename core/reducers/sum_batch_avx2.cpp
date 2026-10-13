// SPDX-License-Identifier: Apache-2.0
//
// The AVX2 kernels of the grid reduction's sums (P-6, D-404): sum_lanes.hpp on four lanes. This file alone is
// compiled for AVX2, and sum_batch.cpp calls it only after checking at run time that the CPU and the OS support
// AVX2. It therefore includes nothing with external-linkage inline functions (no standard library); see
// sum_lanes.hpp. Whole blocks of four only: sum_batch.cpp finishes any remainder with the scalar code.
#if defined(VCAL_REDUCE_AVX2)

#include <immintrin.h>

#include "core/reducers/sum_lanes.hpp"

namespace vcal::reducers::sum_batch {
namespace {

struct V256 {
    using T = __m256d;
    using M = __m256d;
    static constexpr int width = 4;
    static T set1(double v) { return _mm256_set1_pd(v); }
    static T load(const double* p) { return _mm256_loadu_pd(p); }
    static void store(double* p, T v) { _mm256_storeu_pd(p, v); }
    static T add(T a, T b) { return _mm256_add_pd(a, b); }
    static T sub(T a, T b) { return _mm256_sub_pd(a, b); }
    static T mul(T a, T b) { return _mm256_mul_pd(a, b); }
    static T abs(T a) { return _mm256_andnot_pd(_mm256_set1_pd(-0.0), a); }  // clears the sign bit, as fabs
    static M ge(T a, T b) { return _mm256_cmp_pd(a, b, _CMP_GE_OQ); }
    static M isfinite(T a) {  // |a| < +inf: false for +-inf and NaN (+inf by its bit pattern: no <cmath> here)
        return _mm256_cmp_pd(abs(a), _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FF0000000000000LL)), _CMP_LT_OQ);
    }
    static M not_(M m) { return _mm256_xor_pd(m, _mm256_castsi256_pd(_mm256_set1_epi64x(-1))); }
    static M or_(M a, M b) { return _mm256_or_pd(a, b); }
    static M none_mask() { return _mm256_setzero_pd(); }
    static bool none(M m) { return _mm256_movemask_pd(m) == 0; }
    static T blend(M m, T a, T b) { return _mm256_blendv_pd(b, a, m); }  // m ? a : b
};

}  // namespace

long long weighted_sums_avx2(const double* L, long long T, long long K, const double* w, long long n, double* out,
                             unsigned char* dirty) {
    constexpr int kB = 4;  // four vectors of four interleaved
    long long i = 0;
    for (; i + 4 * kB <= n; i += 4 * kB) {
        if (!lanes_weighted_sums<V256, kB>(L + i, T, K, w, out + i)) {
            for (int j = 0; j < kB; ++j) dirty[i / 4 + j] = 1;
        }
    }
    for (; i + 4 <= n; i += 4) {
        if (!lanes_weighted_sums<V256, 1>(L + i, T, K, w, out + i)) dirty[i / 4] = 1;
    }
    return i;
}

long long tile_bounds_avx2(const double* U, long long periods, long long tiles, const double* w, double margin,
                           double* bound) {
    constexpr int kB = 4;
    long long i = 0;
    for (; i + 4 * kB <= tiles; i += 4 * kB) lanes_tile_bounds<V256, kB>(U + i, periods, tiles, w, margin, bound + i);
    for (; i + 4 <= tiles; i += 4) lanes_tile_bounds<V256, 1>(U + i, periods, tiles, w, margin, bound + i);
    return i;
}

}  // namespace vcal::reducers::sum_batch

#endif
