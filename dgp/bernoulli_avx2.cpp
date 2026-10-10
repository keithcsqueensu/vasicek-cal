// SPDX-License-Identifier: Apache-2.0
//
// The AVX2 Bernoulli kernel (P-17, D-194): dgp/philox_lanes.hpp on 8 lanes, 2 registers interleaved.
// This file alone is compiled for AVX2, and dgp/bernoulli.cpp calls it only after checking at run time
// that the CPU and the OS support AVX2. It therefore includes nothing with external-linkage inline
// functions (no standard library beyond <cstdint>); see dgp/philox_lanes.hpp.
#include <cstdint>

#if defined(VCAL_DGP_AVX2)

#include <immintrin.h>

#include "dgp/philox_lanes.hpp"

namespace vcal::dgp {
namespace {

struct V256 {
    using T = __m256i;
    static constexpr int kLanes = 8;
    static T set1(std::uint32_t v) { return _mm256_set1_epi32(static_cast<int>(v)); }
    static T iota(std::uint32_t v) {
        return _mm256_add_epi32(_mm256_set1_epi32(static_cast<int>(v)), _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
    }
    static T even_mask() { return _mm256_set_epi32(0, -1, 0, -1, 0, -1, 0, -1); }
    static T sub(T a, T b) { return _mm256_sub_epi32(a, b); }
    static T x(T a, T b) { return _mm256_xor_si256(a, b); }
    static T band(T a, T b) { return _mm256_and_si256(a, b); }
    static T bor(T a, T b) { return _mm256_or_si256(a, b); }
    static T mul_even(T a, T b) { return _mm256_mul_epu32(a, b); }
    static T srli64_32(T a) { return _mm256_srli_epi64(a, 32); }
    static T slli64_32(T a) { return _mm256_slli_epi64(a, 32); }
    static T srli32_6(T a) { return _mm256_srli_epi32(a, 6); }
    static T lt(T a, T b) { return _mm256_cmpgt_epi32(b, a); }
    static T eq(T a, T b) { return _mm256_cmpeq_epi32(a, b); }
    static std::uint32_t hsum(T a) {
        const __m128i s = _mm_add_epi32(_mm256_castsi256_si128(a), _mm256_extracti128_si256(a, 1));
        const __m128i t = _mm_add_epi32(s, _mm_shuffle_epi32(s, 0x4E));
        const __m128i u = _mm_add_epi32(t, _mm_shuffle_epi32(t, 0xB1));
        return static_cast<std::uint32_t>(_mm_cvtsi128_si32(u));
    }
};

}  // namespace

// Declared in dgp/bernoulli.cpp. nblocks is a multiple of kAvx2Step = 16.
std::uint64_t count_blocks_avx2(std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                                std::uint32_t block0, std::uint32_t key0, std::uint32_t key1, std::uint64_t nblocks,
                                std::uint64_t threshold) {
    const std::uint64_t d =
        count_blocks_lanes<V256, 2>(scenario, replicate, period, block0, key0, key1, nblocks, threshold);
    _mm256_zeroupper();
    return d;
}

}  // namespace vcal::dgp

#endif  // VCAL_DGP_AVX2
