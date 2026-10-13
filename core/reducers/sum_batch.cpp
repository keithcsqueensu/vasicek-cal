// SPDX-License-Identifier: Apache-2.0
//
// The grid reduction's sum kernels (P-6, D-404; sum_batch.hpp): the scalar kernel, the SSE2 lane type, the run-time
// selection and the dispatch. A call runs the selected kernel on as many whole blocks of its width as fit and
// finishes the remainder with the scalar code, which gives the same values: the results do not depend on the kernel
// or on where a block boundary falls.
#include "core/reducers/sum_batch.hpp"

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <string>

#include "common/cpu/avx2.hpp"
#include "core/reducers/weighted_sum.hpp"

#if defined(__x86_64__) || defined(_M_X64)
#define VCAL_REDUCE_X86_64 1
#include <emmintrin.h>

#include "core/reducers/sum_lanes.hpp"
#endif

namespace vcal::reducers::sum_batch {

#if defined(VCAL_REDUCE_AVX2)
// sum_batch_avx2.cpp: whole blocks of four only; each returns how many it did.
// Blocks of four whose lanes met a non-finite term are marked in dirty[i / 4] for the caller to redo.
long long weighted_sums_avx2(const double* L, long long T, long long K, const double* w, long long n, double* out,
                             unsigned char* dirty);
long long tile_bounds_avx2(const double* U, long long periods, long long tiles, const double* w, double margin,
                           double* bound);
#endif

namespace {

// The tile bound of one tile, as resample/bootstrap.hpp's scalar loop computes it.
double tile_bound_scalar(const double* U, std::int64_t periods, std::int64_t tiles, const double* w, double margin,
                         std::int64_t i) {
    double s = 0.0, mag = 0.0;
    for (std::int64_t j = 0; j < periods; ++j) {
        if (w[j] == 0.0) continue;
        const double x = w[j] * U[j * tiles + i];
        s += x;
        mag += std::fabs(x);
    }
    return s + margin * mag;
}

#if defined(VCAL_REDUCE_X86_64)
struct V128 {
    using T = __m128d;
    using M = __m128d;
    static constexpr int width = 2;
    static T set1(double v) { return _mm_set1_pd(v); }
    static T load(const double* p) { return _mm_loadu_pd(p); }
    static void store(double* p, T v) { _mm_storeu_pd(p, v); }
    static T add(T a, T b) { return _mm_add_pd(a, b); }
    static T sub(T a, T b) { return _mm_sub_pd(a, b); }
    static T mul(T a, T b) { return _mm_mul_pd(a, b); }
    static T abs(T a) { return _mm_andnot_pd(_mm_set1_pd(-0.0), a); }  // clears the sign bit, as fabs
    static M ge(T a, T b) { return _mm_cmpge_pd(a, b); }
    static M isfinite(T a) { return _mm_cmplt_pd(abs(a), _mm_set1_pd(HUGE_VAL)); }  // false for +-inf and NaN
    static M not_(M m) { return _mm_xor_pd(m, _mm_castsi128_pd(_mm_set1_epi64x(-1))); }
    static M or_(M a, M b) { return _mm_or_pd(a, b); }
    static M none_mask() { return _mm_setzero_pd(); }
    static bool none(M m) { return _mm_movemask_pd(m) == 0; }
    static T blend(M m, T a, T b) { return _mm_or_pd(_mm_and_pd(m, a), _mm_andnot_pd(m, b)); }  // m ? a : b
};

std::int64_t weighted_sums_sse2(const double* L, std::int64_t T, std::int64_t K, const double* w, std::int64_t k0,
                                std::int64_t n, double* out) {
    constexpr int kB = 4;  // four vectors of two interleaved
    std::int64_t i = 0;
    for (; i + 2 * kB <= n; i += 2 * kB) {
        if (!lanes_weighted_sums<V128, kB>(L + k0 + i, T, K, w, out + i)) {
            for (std::int64_t j = i; j < i + 2 * kB; ++j) out[j] = weighted_sum(L, T, K, w, k0 + j);
        }
    }
    for (; i + 2 <= n; i += 2) {
        if (!lanes_weighted_sums<V128, 1>(L + k0 + i, T, K, w, out + i)) {
            for (std::int64_t j = i; j < i + 2; ++j) out[j] = weighted_sum(L, T, K, w, k0 + j);
        }
    }
    return i;
}

std::int64_t tile_bounds_sse2(const double* U, std::int64_t periods, std::int64_t tiles, const double* w,
                              double margin, double* bound) {
    constexpr int kB = 4;
    std::int64_t i = 0;
    for (; i + 2 * kB <= tiles; i += 2 * kB) lanes_tile_bounds<V128, kB>(U + i, periods, tiles, w, margin, bound + i);
    for (; i + 2 <= tiles; i += 2) lanes_tile_bounds<V128, 1>(U + i, periods, tiles, w, margin, bound + i);
    return i;
}
#endif

Kernel fastest_available() {
    if (available(Kernel::Avx2)) return Kernel::Avx2;
    if (available(Kernel::Sse2)) return Kernel::Sse2;
    return Kernel::Scalar;
}

// VCAL_REDUCE_KERNEL if set (Unavailable for a value that is not a kernel this CPU and build have), else the fastest.
Kernel initial_selection() {
    std::string value;
#if defined(_MSC_VER)
    char* raw = nullptr;
    std::size_t length = 0;
    if (_dupenv_s(&raw, &length, "VCAL_REDUCE_KERNEL") == 0 && raw != nullptr) {
        value = raw;
        std::free(raw);
    }
#else
    if (const char* raw = std::getenv("VCAL_REDUCE_KERNEL")) value = raw;
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
#if defined(VCAL_REDUCE_X86_64)
            return true;
#else
            return false;
#endif
        case Kernel::Avx2: {
#if defined(VCAL_REDUCE_AVX2)
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

bool weighted_sums(const double* L, std::int64_t T, std::int64_t K, const double* w, std::int64_t k0, std::int64_t n,
                   double* out) {
    std::int64_t done = 0;
    switch (selected()) {
        case Kernel::Scalar:
            break;
        case Kernel::Sse2:
#if defined(VCAL_REDUCE_X86_64)
            done = weighted_sums_sse2(L, T, K, w, k0, n, out);
#endif
            break;
        case Kernel::Avx2:
#if defined(VCAL_REDUCE_AVX2)
            // In chunks of up to 256 points; blocks of four that met a non-finite term are redone with the scalar code.
            while (done < n) {
                constexpr std::int64_t kChunk = 256;
                unsigned char dirty[kChunk / 4] = {};
                const std::int64_t m = n - done < kChunk ? n - done : kChunk;
                const std::int64_t did = weighted_sums_avx2(L + k0 + done, T, K, w, m, out + done, dirty);
                for (std::int64_t blk = 0; blk < did / 4; ++blk) {
                    if (!dirty[blk]) continue;
                    for (std::int64_t j = done + 4 * blk; j < done + 4 * blk + 4; ++j) {
                        out[j] = weighted_sum(L, T, K, w, k0 + j);
                    }
                }
                done += did;
                if (did < m) break;  // the remainder (under four points) goes to the scalar code below
            }
#endif
            break;
        case Kernel::Unavailable:
            return false;
    }
    for (std::int64_t i = done; i < n; ++i) out[i] = weighted_sum(L, T, K, w, k0 + i);
    return true;
}

bool tile_bounds(const double* U, std::int64_t periods, std::int64_t tiles, const double* w, double margin,
                 double* bound) {
    std::int64_t done = 0;
    switch (selected()) {
        case Kernel::Scalar:
            break;
        case Kernel::Sse2:
#if defined(VCAL_REDUCE_X86_64)
            done = tile_bounds_sse2(U, periods, tiles, w, margin, bound);
#endif
            break;
        case Kernel::Avx2:
#if defined(VCAL_REDUCE_AVX2)
            done = tile_bounds_avx2(U, periods, tiles, w, margin, bound);
#endif
            break;
        case Kernel::Unavailable:
            return false;
    }
    for (std::int64_t i = done; i < tiles; ++i) bound[i] = tile_bound_scalar(U, periods, tiles, w, margin, i);
    return true;
}

}  // namespace vcal::reducers::sum_batch
