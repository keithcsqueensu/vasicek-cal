// SPDX-License-Identifier: Apache-2.0
//
// The AVX2 derivative node-function kernel (P-22, D-402): derivs_lanes.hpp on four lanes. This file alone is compiled
// for AVX2, and derivs_batch.cpp calls it only after checking at run time that the CPU and the OS support AVX2. It
// therefore includes nothing with external-linkage inline functions (no standard library beyond <cstddef>); see
// derivs_lanes.hpp. Whole blocks of four only: derivs_batch.cpp finishes any remainder with the scalar kernel.
#include <cstddef>

#if defined(VCAL_DERIVS_AVX2)

#include <immintrin.h>

#include "core/model/derivs_lanes.hpp"

namespace vcal::model::derivs_batch {
namespace {

struct V256 {
    using T = __m256d;
    using M = __m256d;
    static T set1(double v) { return _mm256_set1_pd(v); }
    static T add(T a, T b) { return _mm256_add_pd(a, b); }
    static T sub(T a, T b) { return _mm256_sub_pd(a, b); }
    static T mul(T a, T b) { return _mm256_mul_pd(a, b); }
    static T div(T a, T b) { return _mm256_div_pd(a, b); }
    static T neg(T a) { return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); }  // the sign bit, as unary minus
    static M lt(T a, T b) { return _mm256_cmp_pd(a, b, _CMP_LT_OQ); }
    static M gt(T a, T b) { return _mm256_cmp_pd(a, b, _CMP_GT_OQ); }
    static T blend(M m, T a, T b) { return _mm256_blendv_pd(b, a, m); }  // m ? a : b
};

template <int kB>
void block(const Params& p, const PeriodConstants& k, const double* z, const double* e, const double* E,
           double* f, std::size_t i) {
    __m256d zv[kB], ev[kB], Ev[kB], out[5][kB];
    for (int j = 0; j < kB; ++j) {
        zv[j] = _mm256_loadu_pd(z + i + 4 * j);
        ev[j] = _mm256_loadu_pd(e + i + 4 * j);
        Ev[j] = _mm256_loadu_pd(E + i + 4 * j);
    }
    lanes_node_functions<V256, kB>(p, k, zv, ev, Ev, out);
    for (int j = 0; j < kB; ++j) {  // node-major: each node's five values together
        double lane[5][4];
        for (int m = 0; m < 5; ++m) _mm256_storeu_pd(lane[m], out[m][j]);
        double* q = f + 5 * (i + 4 * j);
        for (int l = 0; l < 4; ++l) {
            for (int m = 0; m < 5; ++m) q[5 * l + m] = lane[m][l];
        }
    }
}

}  // namespace

std::size_t node_functions_avx2(const Params& p, const PeriodConstants& k, const double* z, const double* e,
                                const double* E, double* f, std::size_t n) {
    constexpr int kB = 4;  // four vectors of four interleaved
    std::size_t i = 0;
    for (; i + 4 * kB <= n; i += 4 * kB) block<kB>(p, k, z, e, E, f, i);
    for (; i + 4 <= n; i += 4) block<1>(p, k, z, e, E, f, i);
    return i;
}

}  // namespace vcal::model::derivs_batch

#endif
