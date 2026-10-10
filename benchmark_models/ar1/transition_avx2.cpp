// SPDX-License-Identifier: Apache-2.0
//
// The AVX2 transition kernels (P-18, D-199): transition_lanes.hpp on four lanes. This file alone is compiled for
// AVX2, and transition.cpp calls it only after checking at run time that the CPU and the OS support AVX2. It
// therefore includes nothing with external-linkage inline functions (no standard library beyond <cstddef>), and
// builds its infinities from bit patterns; see transition_lanes.hpp. Whole blocks of four only: transition.cpp
// finishes any remainder with the scalar lane type.
#include <cstddef>

#if defined(VCAL_AR1_AVX2)

#include <immintrin.h>

#include "benchmark_models/ar1/transition_lanes.hpp"

namespace vcal::benchmark_models::ar1::kernel {
namespace {

struct V256 {
    using T = __m256d;
    using M = __m256d;
    static constexpr int kLanes = 4;
    static constexpr bool kCounts = false;
    static T set1(double v) { return _mm256_set1_pd(v); }
    static T load(const double* p) { return _mm256_loadu_pd(p); }
    static void store(double* p, T v) { _mm256_storeu_pd(p, v); }
    static T add(T a, T b) { return _mm256_add_pd(a, b); }
    static T sub(T a, T b) { return _mm256_sub_pd(a, b); }
    static T mul(T a, T b) { return _mm256_mul_pd(a, b); }
    static T vmin(T x, T c) { return _mm256_min_pd(x, c); }
    static T vmax(T x, T c) { return _mm256_max_pd(x, c); }
    static T floor(T x) { return _mm256_round_pd(x, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC); }
    static T pow2(T k) {
        const __m128i e = _mm_add_epi32(_mm256_cvttpd_epi32(k), _mm_set1_epi32(1023));
        return _mm256_castsi256_pd(_mm256_slli_epi64(_mm256_cvtepu32_epi64(e), 52));
    }
    static T neg_abs(T x) { return _mm256_or_pd(x, _mm256_set1_pd(-0.0)); }
    static T pinf() { return _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FF0000000000000LL)); }
    static T ninf() { return _mm256_castsi256_pd(_mm256_set1_epi64x(static_cast<long long>(0xFFF0000000000000ULL))); }
    static M gt(T a, T b) { return _mm256_cmp_pd(a, b, _CMP_GT_OQ); }
    static M lt(T a, T b) { return _mm256_cmp_pd(a, b, _CMP_LT_OQ); }
    static M eq(T a, T b) { return _mm256_cmp_pd(a, b, _CMP_EQ_OQ); }
    static M unord(T a, T b) { return _mm256_cmp_pd(a, b, _CMP_UNORD_Q); }
    static M none() { return _mm256_setzero_pd(); }
    static M bor(M a, M b) { return _mm256_or_pd(a, b); }
    static M band(M a, M b) { return _mm256_and_pd(a, b); }
    static M bnot(M a) { return _mm256_xor_pd(a, _mm256_castsi256_pd(_mm256_set1_epi64x(-1))); }
    static bool all(M a) { return _mm256_movemask_pd(a) == 15; }
    static unsigned bits(M a) { return static_cast<unsigned>(_mm256_movemask_pd(a)); }
    static T blend(M m, T a, T b) { return _mm256_blendv_pd(b, a, m); }
};

}  // namespace

std::size_t forward_avx2(const double* pu, const double* lhat, std::size_t np, const double* z, const double* a,
                         std::size_t nd, double c, double hlv, double* mx, double* sum, unsigned char* flags) {
    std::size_t k = 0, terms = 0, below = 0;
    for (; k + 4 <= nd; k += 4) {
        forward_block<V256>(pu, lhat, np, z + k, a + k, c, hlv, mx + k, sum + k, flags + k, &terms, &below);
    }
    return k;
}

std::size_t backward_avx2(const double* pu, std::size_t nc, const double* z, const double* a, const double* v,
                          const double* b, std::size_t nn, double c, double hlv, double* mx, double* sum,
                          unsigned char* flags, double* e) {
    std::size_t j = 0, terms = 0, below = 0;
    for (; j + 4 <= nc; j += 4) {
        backward_block<V256>(pu + j, z, a, v, b, nn, c, hlv, mx + j, sum + j, flags + j,
                             e != nullptr ? e + j * nn : nullptr, &terms, &below);
    }
    return j;
}

std::size_t exp_shifted_avx2(const double* in, std::size_t n, double shift, double* out) {
    std::size_t i = 0, terms = 0, below = 0;
    for (; i + 4 <= n; i += 4) exp_shifted_block<V256>(in + i, shift, out + i, &terms, &below);
    return i;
}

std::size_t exp_avx2(const double* in, std::size_t n, double* out) {
    std::size_t i = 0;
    for (; i + 4 <= n; i += 4) exp_block<V256>(in + i, out + i);
    return i;
}

}  // namespace vcal::benchmark_models::ar1::kernel

#endif  // VCAL_AR1_AVX2
