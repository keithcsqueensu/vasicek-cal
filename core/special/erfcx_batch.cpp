// SPDX-License-Identifier: Apache-2.0
//
// The erfcx lane kernels (P-7, D-210; erfcx_batch.hpp): the scalar kernel, the SSE2 lane type, the run-time
// selection and the dispatch. A call runs the selected kernel on as many whole blocks of its width as fit and
// finishes the remainder with the scalar kernel, which gives the same values: the results do not depend on the
// kernel or on where a block boundary falls.
#include "core/special/erfcx_batch.hpp"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <string>

#include "common/cpu/avx2.hpp"
#include "core/special/generated/erfcx_chebyshev.hpp"

#if defined(__x86_64__) || defined(_M_X64)
#define VCAL_ERFCX_X86_64 1
#include <emmintrin.h>

#include "core/special/erfcx_lanes.hpp"
#endif

namespace vcal::special::erfcx_batch {

#if defined(VCAL_ERFCX_AVX2)
// erfcx_batch_avx2.cpp: whole blocks of four only; returns how many it did.
std::size_t nonnegative_avx2(const double* x, double* out, std::size_t n);
#endif

namespace {

#if defined(VCAL_ERFCX_X86_64)
struct V128 {
    using T = __m128d;
    using M = __m128d;
    static T set1(double v) { return _mm_set1_pd(v); }
    static T add(T a, T b) { return _mm_add_pd(a, b); }
    static T sub(T a, T b) { return _mm_sub_pd(a, b); }
    static T mul(T a, T b) { return _mm_mul_pd(a, b); }
    static T div(T a, T b) { return _mm_div_pd(a, b); }
    static M le(T a, T b) { return _mm_cmple_pd(a, b); }
    static T blend(M m, T a, T b) { return _mm_or_pd(_mm_and_pd(m, a), _mm_andnot_pd(m, b)); }  // m ? a : b
    static bool all(M m) { return _mm_movemask_pd(m) == 3; }
    static bool none(M m) { return _mm_movemask_pd(m) == 0; }
};

std::size_t nonnegative_sse2(const double* x, double* out, std::size_t n) {
    constexpr int kB = 4;  // four vectors of two interleaved (erfcx_lanes.hpp)
    std::size_t i = 0;
    for (; i + 2 * kB <= n; i += 2 * kB) {
        __m128d in[kB], res[kB];
        for (int j = 0; j < kB; ++j) in[j] = _mm_loadu_pd(x + i + 2 * j);
        lanes_erfcx_nonnegative<V128, kB>(in, res);
        for (int j = 0; j < kB; ++j) _mm_storeu_pd(out + i + 2 * j, res[j]);
    }
    for (; i + 2 <= n; i += 2) {
        __m128d in[1] = {_mm_loadu_pd(x + i)}, res[1];
        lanes_erfcx_nonnegative<V128, 1>(in, res);
        _mm_storeu_pd(out + i, res[0]);
    }
    return i;
}
#endif

Kernel fastest_available() {
    if (available(Kernel::Avx2)) return Kernel::Avx2;
    if (available(Kernel::Sse2)) return Kernel::Sse2;
    return Kernel::Scalar;
}

// VCAL_ERFCX_KERNEL if set (Unavailable for a value that is not a kernel this CPU and build have), else the fastest.
Kernel initial_selection() {
    std::string value;
#if defined(_MSC_VER)
    char* raw = nullptr;
    std::size_t length = 0;
    if (_dupenv_s(&raw, &length, "VCAL_ERFCX_KERNEL") == 0 && raw != nullptr) {
        value = raw;
        std::free(raw);
    }
#else
    if (const char* raw = std::getenv("VCAL_ERFCX_KERNEL")) value = raw;
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
#if defined(VCAL_ERFCX_X86_64)
            return true;
#else
            return false;
#endif
        case Kernel::Avx2: {
#if defined(VCAL_ERFCX_AVX2)
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

bool nonnegative(const double* x, double* out, std::size_t n) {
    std::size_t done = 0;
    switch (selected()) {
        case Kernel::Scalar:
            break;
        case Kernel::Sse2:
#if defined(VCAL_ERFCX_X86_64)
            done = nonnegative_sse2(x, out, n);
#endif
            break;
        case Kernel::Avx2:
#if defined(VCAL_ERFCX_AVX2)
            done = nonnegative_avx2(x, out, n);
#endif
            break;
        case Kernel::Unavailable:
            return false;
    }
    for (std::size_t i = done; i < n; ++i) out[i] = generated::erfcx_nonnegative(x[i]);
    return true;
}

}  // namespace vcal::special::erfcx_batch
