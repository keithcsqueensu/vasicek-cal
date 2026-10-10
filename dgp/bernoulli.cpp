// SPDX-License-Identifier: Apache-2.0
//
// The Bernoulli count d_t of one period (P-17, D-194; dgp/dgp.hpp, count_defaults). Three kernels give
// the same count bit for bit: the scalar loop (the specification, operation for operation; any CPU),
// SSE2 (every x86-64 CPU) and AVX2 (chosen at run time when the CPU and the OS support it). The vector
// kernels compare integer keys against a threshold instead of uniforms against p; dgp/philox_lanes.hpp
// gives the argument that they equal the scalar loop, and tests/unit/dgp_test.cpp checks it.
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>

#include "dgp/dgp.hpp"
#include "dgp/philox.hpp"

#if defined(__x86_64__) || defined(_M_X64)
#define VCAL_DGP_X86_64 1
#include <emmintrin.h>
#if defined(_MSC_VER)
#include <intrin.h>
#endif
#include "dgp/philox_lanes.hpp"
#endif

namespace vcal::dgp {

#if defined(VCAL_DGP_AVX2)
// dgp/bernoulli_avx2.cpp. nblocks is a multiple of kAvx2Step.
std::uint64_t count_blocks_avx2(std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                                std::uint32_t block0, std::uint32_t key0, std::uint32_t key1, std::uint64_t nblocks,
                                std::uint64_t threshold);
#endif

namespace {

#if defined(VCAL_DGP_X86_64)
constexpr std::uint64_t kSse2Step = 4 * 4;  // 4 lanes, 4 registers
#endif
#if defined(VCAL_DGP_AVX2)
constexpr std::uint64_t kAvx2Step = 8 * 2;  // 8 lanes, 2 registers
#endif

// The specification's loop (dgp/dgp.hpp), operation for operation: the scalar kernel, and the reference
// the vector kernels are tested against.
std::int64_t count_defaults_spec(std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                                 const PhiloxKey& key, std::int64_t n, double p) {
    std::int64_t d = 0;
    PhiloxCounter w{};
    for (std::int64_t i = 0; i < n; ++i) {
        if (i % 2 == 0) {
            const auto block = static_cast<std::uint32_t>(kFactorBlocks + static_cast<std::uint64_t>(i / 2));
            w = philox4x32({scenario, replicate, period, block}, key);
        }
        const double u = i % 2 == 0 ? uniform_from_words(w[0], w[1]) : uniform_from_words(w[2], w[3]);
        if (u < p) ++d;
    }
    return d;
}

// The obligor's integer key k, U = (k + 1/2) 2^-52 (uniform_from_words).
std::uint64_t key52(std::uint32_t a, std::uint32_t b) {
    return (static_cast<std::uint64_t>(a >> 6) << 26) | static_cast<std::uint64_t>(b >> 6);
}

// K, the number of k in [0, 2^52) with (k + 1/2) 2^-52 < p, so that U < p iff k < K. Both floating-point
// steps are exact for 0 <= p <= 1: scaling by 2^52 is exact, and x = p 2^52 - 1/2 is a multiple of
// ulp(p 2^52) <= 1/2 no larger in magnitude than p 2^52. k + 1/2 < p 2^52 iff k < x iff k < ceil(x).
// A NaN p counts no defaults, as U < NaN is false.
std::uint64_t threshold(double p) {
    const double x = p * 0x1p52 - 0.5;
    if (!(x > 0.0)) return 0;
    if (x >= 0x1p52) return std::uint64_t{1} << 52;
    return static_cast<std::uint64_t>(std::ceil(x));
}

std::uint64_t count_blocks_scalar(std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                                  const PhiloxKey& key, std::uint64_t first, std::uint64_t last, std::uint64_t k) {
    std::uint64_t d = 0;
    for (std::uint64_t j = first; j < last; ++j) {
        const PhiloxCounter w = philox4x32({scenario, replicate, period, static_cast<std::uint32_t>(j)}, key);
        d += key52(w[0], w[1]) < k ? 1u : 0u;
        d += key52(w[2], w[3]) < k ? 1u : 0u;
    }
    return d;
}

#if defined(VCAL_DGP_X86_64)
struct V128 {
    using T = __m128i;
    static constexpr int kLanes = 4;
    static T set1(std::uint32_t v) { return _mm_set1_epi32(static_cast<int>(v)); }
    static T iota(std::uint32_t v) {
        return _mm_add_epi32(_mm_set1_epi32(static_cast<int>(v)), _mm_setr_epi32(0, 1, 2, 3));
    }
    static T even_mask() { return _mm_set_epi32(0, -1, 0, -1); }
    static T sub(T a, T b) { return _mm_sub_epi32(a, b); }
    static T x(T a, T b) { return _mm_xor_si128(a, b); }
    static T band(T a, T b) { return _mm_and_si128(a, b); }
    static T bor(T a, T b) { return _mm_or_si128(a, b); }
    static T mul_even(T a, T b) { return _mm_mul_epu32(a, b); }
    static T srli64_32(T a) { return _mm_srli_epi64(a, 32); }
    static T slli64_32(T a) { return _mm_slli_epi64(a, 32); }
    static T srli32_6(T a) { return _mm_srli_epi32(a, 6); }
    static T lt(T a, T b) { return _mm_cmplt_epi32(a, b); }
    static T eq(T a, T b) { return _mm_cmpeq_epi32(a, b); }
    static std::uint32_t hsum(T a) {
        const T s = _mm_add_epi32(a, _mm_shuffle_epi32(a, 0x4E));
        const T t = _mm_add_epi32(s, _mm_shuffle_epi32(s, 0xB1));
        return static_cast<std::uint32_t>(_mm_cvtsi128_si32(t));
    }
};
#endif

#if defined(VCAL_DGP_AVX2)
bool cpu_has_avx2() {
#if defined(_MSC_VER)
    // CPUID leaf 7 EBX bit 5 is AVX2; the OS must also save the YMM state (OSXSAVE, then XCR0 bits 1-2).
    int r[4];
    __cpuid(r, 0);
    if (r[0] < 7) return false;
    __cpuid(r, 1);
    const bool osxsave = (r[2] & (1 << 27)) != 0;
    const bool avx = (r[2] & (1 << 28)) != 0;
    if (!osxsave || !avx || (_xgetbv(0) & 6u) != 6u) return false;
    __cpuidex(r, 7, 0);
    return (r[1] & (1 << 5)) != 0;
#else
    // libgcc and compiler-rt check the OS's XCR0 for AVX features as well.
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") != 0;
#endif
}
#endif

}  // namespace

bool bernoulli_kernel_available(BernoulliKernel kernel) {
    switch (kernel) {
        case BernoulliKernel::Scalar:
            return true;
        case BernoulliKernel::Sse2:
#if defined(VCAL_DGP_X86_64)
            return true;
#else
            return false;
#endif
        case BernoulliKernel::Avx2: {
#if defined(VCAL_DGP_AVX2)
            static const bool has = cpu_has_avx2();
            return has;
#else
            return false;
#endif
        }
        case BernoulliKernel::Unavailable:
            break;
    }
    return false;
}

const char* bernoulli_kernel_name(BernoulliKernel kernel) {
    switch (kernel) {
        case BernoulliKernel::Scalar: return "scalar";
        case BernoulliKernel::Sse2: return "sse2";
        case BernoulliKernel::Avx2: return "avx2";
        case BernoulliKernel::Unavailable: break;
    }
    return "unavailable";
}

namespace {

BernoulliKernel fastest_available() {
    if (bernoulli_kernel_available(BernoulliKernel::Avx2)) return BernoulliKernel::Avx2;
    if (bernoulli_kernel_available(BernoulliKernel::Sse2)) return BernoulliKernel::Sse2;
    return BernoulliKernel::Scalar;
}

// VCAL_DGP_KERNEL if set (Unavailable for a value that is not a kernel this CPU and build have), else the
// fastest available.
BernoulliKernel initial_selection() {
    std::string value;
#if defined(_MSC_VER)
    char* raw = nullptr;
    std::size_t length = 0;
    if (_dupenv_s(&raw, &length, "VCAL_DGP_KERNEL") == 0 && raw != nullptr) {
        value = raw;
        std::free(raw);
    }
#else
    if (const char* raw = std::getenv("VCAL_DGP_KERNEL")) value = raw;
#endif
    if (value.empty()) return fastest_available();
    for (const auto k : {BernoulliKernel::Scalar, BernoulliKernel::Sse2, BernoulliKernel::Avx2}) {
        if (value == bernoulli_kernel_name(k)) return bernoulli_kernel_available(k) ? k : BernoulliKernel::Unavailable;
    }
    return BernoulliKernel::Unavailable;
}

std::atomic<std::int32_t>& selection() {
    static std::atomic<std::int32_t> s{static_cast<std::int32_t>(initial_selection())};
    return s;
}

}  // namespace

BernoulliKernel bernoulli_kernel_selected() {
    return static_cast<BernoulliKernel>(selection().load(std::memory_order_relaxed));
}

bool select_bernoulli_kernel(BernoulliKernel kernel) {
    if (kernel != BernoulliKernel::Unavailable && !bernoulli_kernel_available(kernel)) return false;
    selection().store(static_cast<std::int32_t>(kernel), std::memory_order_relaxed);
    return true;
}

std::uint64_t bernoulli_threshold(double p) { return threshold(p); }

std::int64_t count_defaults(std::uint64_t seed, std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                            std::int64_t n, double p, BernoulliKernel kernel) {
    if (!(n >= 0 && n <= kMaxObligors) || !bernoulli_kernel_available(kernel)) return -1;
    const PhiloxKey key{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)};
    if (kernel == BernoulliKernel::Scalar) return count_defaults_spec(scenario, replicate, period, key, n, p);
    const std::uint64_t k = threshold(p);
    const auto un = static_cast<std::uint64_t>(n);
    const std::uint64_t pairs = un / 2;  // blocks kFactorBlocks .. kFactorBlocks + pairs - 1 hold two obligors
    std::uint64_t bulk = 0;              // blocks handled by a vector kernel, from kFactorBlocks
    std::uint64_t d = 0;
#if defined(VCAL_DGP_AVX2)
    if (kernel == BernoulliKernel::Avx2) {
        bulk = pairs / kAvx2Step * kAvx2Step;
        d += count_blocks_avx2(scenario, replicate, period, kFactorBlocks, key[0], key[1], bulk, k);
    }
#endif
#if defined(VCAL_DGP_X86_64)
    if (kernel == BernoulliKernel::Sse2) {
        bulk = pairs / kSse2Step * kSse2Step;
        d += count_blocks_lanes<V128, 4>(scenario, replicate, period, kFactorBlocks, key[0], key[1], bulk, k);
    }
#endif
    d += count_blocks_scalar(scenario, replicate, period, key, kFactorBlocks + bulk, kFactorBlocks + pairs, k);
    if (un % 2 != 0) {  // the last obligor, words (0, 1) of block kFactorBlocks + pairs
        const auto block = static_cast<std::uint32_t>(kFactorBlocks + pairs);
        const PhiloxCounter w = philox4x32({scenario, replicate, period, block}, key);
        d += key52(w[0], w[1]) < k ? 1u : 0u;
    }
    return static_cast<std::int64_t>(d);
}

}  // namespace vcal::dgp
