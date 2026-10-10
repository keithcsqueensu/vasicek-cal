// SPDX-License-Identifier: Apache-2.0
//
// The AR(1) transition kernels (P-18, D-199; transition.hpp): the scalar and SSE2 lane types, the run-time
// selection, and the dispatch. Each call runs the selected kernel on as many whole blocks of its width as fit and
// finishes the remainder with the scalar lane type, which does the same operations: the results do not depend on
// the kernel or on where a block boundary falls.
#include "benchmark_models/ar1/transition.hpp"

#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>

#include "common/cpu/avx2.hpp"

#if defined(__x86_64__) || defined(_M_X64)
#define VCAL_AR1_X86_64 1
#include <emmintrin.h>
#endif

#include "benchmark_models/ar1/transition_lanes.hpp"

namespace vcal::benchmark_models::ar1::kernel {

#if defined(VCAL_AR1_AVX2)
// transition_avx2.cpp: whole blocks of four only; each returns how many it did.
std::size_t forward_avx2(const double* pu, const double* lhat, std::size_t np, const double* z, const double* a,
                         std::size_t nd, double c, double hlv, double* mx, double* sum, unsigned char* flags);
std::size_t backward_avx2(const double* pu, std::size_t nc, const double* z, const double* a, const double* v,
                          const double* b, std::size_t nn, double c, double hlv, double* mx, double* sum,
                          unsigned char* flags, double* e);
std::size_t exp_shifted_avx2(const double* in, std::size_t n, double shift, double* out);
std::size_t exp_avx2(const double* in, std::size_t n, double* out);
#endif

namespace {

// One double per lane: the reference every vector kernel must match.
struct V1 {
    using T = double;
    using M = bool;
    static constexpr int kLanes = 1;
    static constexpr bool kCounts = true;
    static T set1(double v) { return v; }
    static T load(const double* p) { return *p; }
    static void store(double* p, T v) { *p = v; }
    static T add(T a, T b) { return a + b; }
    static T sub(T a, T b) { return a - b; }
    static T mul(T a, T b) { return a * b; }
    static T vmin(T x, T c) { return x < c ? x : c; }  // c when x is NaN, as minpd(x, c)
    static T vmax(T x, T c) { return x > c ? x : c; }  // c when x is NaN, as maxpd(x, c)
    static T floor(T x) { return std::floor(x); }
    static T pow2(T k) {  // 2^k for an integral k in [-1022, 1023]
        const auto bits = static_cast<std::uint64_t>(static_cast<std::int64_t>(k) + 1023) << 52;
        double out;
        std::memcpy(&out, &bits, sizeof out);
        return out;
    }
    static T neg_abs(T x) { return -std::fabs(x); }
    static T pinf() { return HUGE_VAL; }
    static T ninf() { return -HUGE_VAL; }
    static M gt(T a, T b) { return a > b; }
    static M lt(T a, T b) { return a < b; }
    static M eq(T a, T b) { return a == b; }
    static M unord(T a, T b) { return a != a || b != b; }
    static M none() { return false; }
    static M bor(M a, M b) { return a || b; }
    static M band(M a, M b) { return a && b; }
    static M bnot(M a) { return !a; }
    static bool all(M a) { return a; }
    static unsigned bits(M a) { return a ? 1u : 0u; }
    static T blend(M m, T a, T b) { return m ? a : b; }
};

#if defined(VCAL_AR1_X86_64)
// Two doubles per lane (SSE2, every x86-64 CPU). floor by truncation and a correction, exact for |x| < 2^31,
// which covers every argument here (the exponent's k and k/2, at most 1,156 in magnitude).
struct V128 {
    using T = __m128d;
    using M = __m128d;
    static constexpr int kLanes = 2;
    static constexpr bool kCounts = false;
    static T set1(double v) { return _mm_set1_pd(v); }
    static T load(const double* p) { return _mm_loadu_pd(p); }
    static void store(double* p, T v) { _mm_storeu_pd(p, v); }
    static T add(T a, T b) { return _mm_add_pd(a, b); }
    static T sub(T a, T b) { return _mm_sub_pd(a, b); }
    static T mul(T a, T b) { return _mm_mul_pd(a, b); }
    static T vmin(T x, T c) { return _mm_min_pd(x, c); }
    static T vmax(T x, T c) { return _mm_max_pd(x, c); }
    static T floor(T x) {
        const T t = _mm_cvtepi32_pd(_mm_cvttpd_epi32(x));
        return _mm_sub_pd(t, _mm_and_pd(_mm_cmpgt_pd(t, x), _mm_set1_pd(1.0)));
    }
    static T pow2(T k) {
        const __m128i e = _mm_add_epi32(_mm_cvttpd_epi32(k), _mm_set1_epi32(1023));
        return _mm_castsi128_pd(_mm_slli_epi64(_mm_unpacklo_epi32(e, _mm_setzero_si128()), 52));
    }
    static T neg_abs(T x) { return _mm_or_pd(x, _mm_set1_pd(-0.0)); }
    static T pinf() { return _mm_set1_pd(HUGE_VAL); }
    static T ninf() { return _mm_set1_pd(-HUGE_VAL); }
    static M gt(T a, T b) { return _mm_cmpgt_pd(a, b); }
    static M lt(T a, T b) { return _mm_cmplt_pd(a, b); }
    static M eq(T a, T b) { return _mm_cmpeq_pd(a, b); }
    static M unord(T a, T b) { return _mm_cmpunord_pd(a, b); }
    static M none() { return _mm_setzero_pd(); }
    static M bor(M a, M b) { return _mm_or_pd(a, b); }
    static M band(M a, M b) { return _mm_and_pd(a, b); }
    static M bnot(M a) { return _mm_xor_pd(a, _mm_castsi128_pd(_mm_set1_epi32(-1))); }
    static bool all(M a) { return _mm_movemask_pd(a) == 3; }
    static unsigned bits(M a) { return static_cast<unsigned>(_mm_movemask_pd(a)); }
    static T blend(M m, T a, T b) { return _mm_or_pd(_mm_and_pd(m, a), _mm_andnot_pd(m, b)); }
};
#endif

Kernel fastest_available() {
    if (available(Kernel::Avx2)) return Kernel::Avx2;
    if (available(Kernel::Sse2)) return Kernel::Sse2;
    return Kernel::Scalar;
}

// VCAL_AR1_KERNEL if set (Unavailable for a value that is not a kernel this CPU and build have), else the fastest.
Kernel initial_selection() {
    std::string value;
#if defined(_MSC_VER)
    char* raw = nullptr;
    std::size_t length = 0;
    if (_dupenv_s(&raw, &length, "VCAL_AR1_KERNEL") == 0 && raw != nullptr) {
        value = raw;
        std::free(raw);
    }
#else
    if (const char* raw = std::getenv("VCAL_AR1_KERNEL")) value = raw;
#endif
    if (value.empty()) return fastest_available();
    for (const auto k : {Kernel::Scalar, Kernel::Sse2, Kernel::Avx2}) {
        if (value == name(k)) return available(k) ? k : Kernel::Unavailable;
    }
    return Kernel::Unavailable;
}

std::atomic<std::int32_t>& selection() {
    static std::atomic<std::int32_t> s{static_cast<std::int32_t>(initial_selection())};
    return s;
}

}  // namespace

bool available(Kernel k) {
    switch (k) {
        case Kernel::Scalar:
            return true;
        case Kernel::Sse2:
#if defined(VCAL_AR1_X86_64)
            return true;
#else
            return false;
#endif
        case Kernel::Avx2: {
#if defined(VCAL_AR1_AVX2)
            static const bool has = cpu::has_avx2();
            return has;
#else
            return false;
#endif
        }
        case Kernel::Unavailable:
            break;
    }
    return false;
}

const char* name(Kernel k) {
    switch (k) {
        case Kernel::Scalar: return "scalar";
        case Kernel::Sse2: return "sse2";
        case Kernel::Avx2: return "avx2";
        case Kernel::Unavailable: break;
    }
    return "unavailable";
}

Kernel selected() { return static_cast<Kernel>(selection().load(std::memory_order_relaxed)); }

bool select(Kernel k) {
    if (k != Kernel::Unavailable && !available(k)) return false;
    selection().store(static_cast<std::int32_t>(k), std::memory_order_relaxed);
    return true;
}

bool forward(const double* pu, const double* lhat, std::size_t np, const double* z, const double* a, std::size_t nd,
             double c, double hlv, double* mx, double* sum, unsigned char* flags, Counts* counts) {
    const Kernel kern = selected();
    if (kern == Kernel::Unavailable) return false;
    std::size_t k = 0, terms = 0, below = 0;
#if defined(VCAL_AR1_AVX2)
    if (kern == Kernel::Avx2) k = forward_avx2(pu, lhat, np, z, a, nd, c, hlv, mx, sum, flags);
#endif
#if defined(VCAL_AR1_X86_64)
    if (kern == Kernel::Sse2) {
        for (; k + 2 <= nd; k += 2) {
            forward_block<V128>(pu, lhat, np, z + k, a + k, c, hlv, mx + k, sum + k, flags + k, &terms, &below);
        }
    }
#endif
    for (; k < nd; ++k) forward_block<V1>(pu, lhat, np, z + k, a + k, c, hlv, mx + k, sum + k, flags + k, &terms, &below);
    if (counts != nullptr) {
        counts->terms += terms;
        counts->below += below;
    }
    return true;
}

bool backward(const double* pu, std::size_t nc, const double* z, const double* a, const double* v, const double* b,
              std::size_t nn, double c, double hlv, double* mx, double* sum, unsigned char* flags, double* e,
              Counts* counts) {
    const Kernel kern = selected();
    if (kern == Kernel::Unavailable) return false;
    std::size_t j = 0, terms = 0, below = 0;
#if defined(VCAL_AR1_AVX2)
    if (kern == Kernel::Avx2) j = backward_avx2(pu, nc, z, a, v, b, nn, c, hlv, mx, sum, flags, e);
#endif
#if defined(VCAL_AR1_X86_64)
    if (kern == Kernel::Sse2) {
        for (; j + 2 <= nc; j += 2) {
            backward_block<V128>(pu + j, z, a, v, b, nn, c, hlv, mx + j, sum + j, flags + j,
                                 e != nullptr ? e + j * nn : nullptr, &terms, &below);
        }
    }
#endif
    for (; j < nc; ++j) {
        backward_block<V1>(pu + j, z, a, v, b, nn, c, hlv, mx + j, sum + j, flags + j, e != nullptr ? e + j * nn : nullptr,
                           &terms, &below);
    }
    if (counts != nullptr) {
        counts->terms += terms;
        counts->below += below;
    }
    return true;
}

bool exp_shifted(const double* in, std::size_t n, double shift, double* out, Counts* counts) {
    const Kernel kern = selected();
    if (kern == Kernel::Unavailable) return false;
    std::size_t i = 0, terms = 0, below = 0;
#if defined(VCAL_AR1_AVX2)
    if (kern == Kernel::Avx2) i = exp_shifted_avx2(in, n, shift, out);
#endif
#if defined(VCAL_AR1_X86_64)
    if (kern == Kernel::Sse2) {
        for (; i + 2 <= n; i += 2) exp_shifted_block<V128>(in + i, shift, out + i, &terms, &below);
    }
#endif
    for (; i < n; ++i) exp_shifted_block<V1>(in + i, shift, out + i, &terms, &below);
    if (counts != nullptr) {
        counts->terms += terms;
        counts->below += below;
    }
    return true;
}

bool exp_lanes(const double* in, std::size_t n, double* out) {
    const Kernel kern = selected();
    if (kern == Kernel::Unavailable) return false;
    std::size_t i = 0;
#if defined(VCAL_AR1_AVX2)
    if (kern == Kernel::Avx2) i = exp_avx2(in, n, out);
#endif
#if defined(VCAL_AR1_X86_64)
    if (kern == Kernel::Sse2) {
        for (; i + 2 <= n; i += 2) exp_block<V128>(in + i, out + i);
    }
#endif
    for (; i < n; ++i) exp_block<V1>(in + i, out + i);
    return true;
}

}  // namespace vcal::benchmark_models::ar1::kernel
