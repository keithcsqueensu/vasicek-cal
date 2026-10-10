// SPDX-License-Identifier: Apache-2.0
//
// The vectorised Bernoulli count (P-17, D-194): Philox4x32-10 on consecutive Bernoulli blocks, one
// block per 32-bit lane, and each obligor's default decided by comparing integer keys, with no
// floating point. Included only by dgp/bernoulli.cpp (SSE2) and dgp/bernoulli_avx2.cpp (AVX2),
// each with its own vector type V.
//
// Everything here has internal linkage on purpose. bernoulli_avx2.cpp is compiled for AVX2 and runs
// only when the CPU has it; an inline function with external linkage instantiated there could be
// merged by the linker with the SSE2 file's copy and run on a CPU without AVX2.
//
// Bit for bit the scalar loop in dgp/dgp.cpp's original form, by construction:
//   - Philox: each round's 64-bit products M0 * c0 and M1 * c2 come from mul_epu32, which
//     multiplies the even 32-bit lanes; the odd lanes are shifted down first. The high and low
//     halves are then put back into their lanes, so every lane computes exactly philox4x32.
//   - The comparison: obligor i defaults iff U_i < p, with U_i = (k + 1/2) 2^-52 and
//     k = hi * 2^26 + lo, hi = a >> 6, lo = b >> 6 (dgp/dgp.hpp). With K the number of integers
//     k in [0, 2^52) satisfying (k + 1/2) 2^-52 < p (dgp/bernoulli.cpp, threshold), U_i < p iff
//     k < K, iff hi < K_hi or (hi == K_hi and lo < K_lo), K_hi = K >> 26, K_lo = K mod 2^26.
//     Every operand is below 2^31, so the signed 32-bit comparisons SSE2 has are exact.
#pragma once

#include <cstdint>

namespace vcal::dgp {
namespace {

// Counts the defaults among the 2 * nblocks obligors of Bernoulli blocks block0 .. block0 + nblocks - 1
// of one period, nblocks a multiple of V::kLanes * R. R registers are interleaved for instruction-level
// parallelism. V wraps one vector width: kLanes 32-bit lanes and the handful of operations used here.
template <class V, int R>
std::uint64_t count_blocks_lanes(std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                                 std::uint32_t block0, std::uint32_t key0, std::uint32_t key1,
                                 std::uint64_t nblocks, std::uint64_t threshold) {
    using T = typename V::T;
    constexpr std::uint32_t kStep = static_cast<std::uint32_t>(V::kLanes * R);
    // Each lane of the accumulator gains at most 2 R per step; it is emptied every kFlush steps,
    // long before it could wrap.
    constexpr std::uint64_t kFlush = std::uint64_t{1} << 24;
    const T m0 = V::set1(0xD2511F53u);
    const T m1 = V::set1(0xCD9E8D57u);
    const T even = V::even_mask();
    const T odd = V::x(even, V::set1(0xFFFFFFFFu));
    const T k_hi = V::set1(static_cast<std::uint32_t>(threshold >> 26));
    const T k_lo = V::set1(static_cast<std::uint32_t>(threshold & ((std::uint64_t{1} << 26) - 1)));
    const T vs = V::set1(scenario);
    const T vr = V::set1(replicate);
    const T vt = V::set1(period);
    std::uint64_t total = 0;
    T acc = V::set1(0);
    std::uint64_t since_flush = 0;
    for (std::uint64_t j = 0; j < nblocks; j += kStep) {
        T c0[R], c1[R], c2[R], c3[R];
        for (int q = 0; q < R; ++q) {
            c0[q] = vs;
            c1[q] = vr;
            c2[q] = vt;
            c3[q] = V::iota(block0 + static_cast<std::uint32_t>(j) + static_cast<std::uint32_t>(q * V::kLanes));
        }
        std::uint32_t k0 = key0;
        std::uint32_t k1 = key1;
        for (int round = 0; round < 10; ++round) {
            if (round > 0) {
                k0 += 0x9E3779B9u;
                k1 += 0xBB67AE85u;
            }
            const T vk0 = V::set1(k0);
            const T vk1 = V::set1(k1);
            for (int q = 0; q < R; ++q) {
                const T e0 = V::mul_even(c0[q], m0);
                const T o0 = V::mul_even(V::srli64_32(c0[q]), m0);
                const T e1 = V::mul_even(c2[q], m1);
                const T o1 = V::mul_even(V::srli64_32(c2[q]), m1);
                const T hi0 = V::bor(V::srli64_32(e0), V::band(o0, odd));
                const T lo0 = V::bor(V::band(e0, even), V::slli64_32(o0));
                const T hi1 = V::bor(V::srli64_32(e1), V::band(o1, odd));
                const T lo1 = V::bor(V::band(e1, even), V::slli64_32(o1));
                const T n0 = V::x(V::x(hi1, c1[q]), vk0);
                const T n2 = V::x(V::x(hi0, c3[q]), vk1);
                c1[q] = lo1;
                c3[q] = lo0;
                c0[q] = n0;
                c2[q] = n2;
            }
        }
        for (int q = 0; q < R; ++q) {
            // Words (0, 1) are obligor 2j, words (2, 3) obligor 2j + 1. A true comparison is all ones,
            // that is -1, so subtracting it counts one.
            const T a0 = V::srli32_6(c0[q]);
            const T b0 = V::srli32_6(c1[q]);
            const T a1 = V::srli32_6(c2[q]);
            const T b1 = V::srli32_6(c3[q]);
            const T below0 = V::bor(V::lt(a0, k_hi), V::band(V::eq(a0, k_hi), V::lt(b0, k_lo)));
            const T below1 = V::bor(V::lt(a1, k_hi), V::band(V::eq(a1, k_hi), V::lt(b1, k_lo)));
            acc = V::sub(V::sub(acc, below0), below1);
        }
        if (++since_flush == kFlush) {
            total += V::hsum(acc);
            acc = V::set1(0);
            since_flush = 0;
        }
    }
    return total + V::hsum(acc);
}

}  // namespace
}  // namespace vcal::dgp
