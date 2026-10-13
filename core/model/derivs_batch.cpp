// SPDX-License-Identifier: Apache-2.0
//
// The derivative node-function kernels (P-22, D-402; derivs_batch.hpp): the scalar kernel, the SSE2 lane type, the
// run-time selection and the dispatch. A call runs the selected kernel on as many whole blocks of its width as fit
// and finishes the remainder with the scalar kernel, which gives the same values: the results do not depend on the
// kernel or on where a block boundary falls.
#include "core/model/derivs_batch.hpp"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <string>

#include "common/cpu/avx2.hpp"
#include "core/model/derivs_lanes.hpp"

#if defined(__x86_64__) || defined(_M_X64)
#define VCAL_DERIVS_X86_64 1
#include <emmintrin.h>
#endif

namespace vcal::model::derivs_batch {

#if defined(VCAL_DERIVS_AVX2)
// derivs_batch_avx2.cpp: whole blocks of four only; returns how many it did.
std::size_t node_functions_avx2(const Params& p, const PeriodConstants& k, const double* z, const double* e,
                                const double* E, double* f, std::size_t n);
#endif

namespace {

PeriodConstants period_constants(const Params& p) {
    const double b2 = p.b * p.b;
    const double ab = p.a * p.b;
    PeriodConstants k;
    k.two_ab = 2.0 * p.a * p.b;
    k.two_b2 = 2.0 * b2;
    k.x1 = p.c1 / p.b;
    k.x11 = p.c2 / p.b;
    k.x12 = p.c1 / (2.0 * b2 * p.b);
    k.one_m_2r = 1.0 - 2.0 * p.rho;
    k.four_ab3 = 4.0 * ab * ab * ab;
    k.two_b2b2 = 2.0 * b2 * b2;
    return k;
}

// The scalar lane: one double through the same template, so the scalar kernel is the lanes' arithmetic too.
struct VScalar {
    using T = double;
    using M = bool;
    static T set1(double v) { return v; }
    static T add(T a, T b) { return a + b; }
    static T sub(T a, T b) { return a - b; }
    static T mul(T a, T b) { return a * b; }
    static T div(T a, T b) { return a / b; }
    static T neg(T a) { return -a; }
    static M lt(T a, T b) { return a < b; }
    static M gt(T a, T b) { return a > b; }
    static T blend(M m, T a, T b) { return m ? a : b; }
};

void node_functions_scalar(const Params& p, const PeriodConstants& k, const double* z, const double* e,
                           const double* E, double* f, std::size_t from, std::size_t n) {
    for (std::size_t i = from; i < n; ++i) {
        double out[5][1];
        lanes_node_functions<VScalar, 1>(p, k, z + i, e + i, E + i, out);
        for (int j = 0; j < 5; ++j) f[5 * i + j] = out[j][0];
    }
}

#if defined(VCAL_DERIVS_X86_64)
struct V128 {
    using T = __m128d;
    using M = __m128d;
    static T set1(double v) { return _mm_set1_pd(v); }
    static T add(T a, T b) { return _mm_add_pd(a, b); }
    static T sub(T a, T b) { return _mm_sub_pd(a, b); }
    static T mul(T a, T b) { return _mm_mul_pd(a, b); }
    static T div(T a, T b) { return _mm_div_pd(a, b); }
    static T neg(T a) { return _mm_xor_pd(a, _mm_set1_pd(-0.0)); }  // the sign bit, as unary minus
    static M lt(T a, T b) { return _mm_cmplt_pd(a, b); }
    static M gt(T a, T b) { return _mm_cmpgt_pd(a, b); }
    static T blend(M m, T a, T b) { return _mm_or_pd(_mm_and_pd(m, a), _mm_andnot_pd(m, b)); }  // m ? a : b
};

template <int kB>
void block_sse2(const Params& p, const PeriodConstants& k, const double* z, const double* e, const double* E,
                double* f, std::size_t i) {
    __m128d zv[kB], ev[kB], Ev[kB], out[5][kB];
    for (int j = 0; j < kB; ++j) {
        zv[j] = _mm_loadu_pd(z + i + 2 * j);
        ev[j] = _mm_loadu_pd(e + i + 2 * j);
        Ev[j] = _mm_loadu_pd(E + i + 2 * j);
    }
    lanes_node_functions<V128, kB>(p, k, zv, ev, Ev, out);
    for (int j = 0; j < kB; ++j) {  // node-major: each node's five values together
        double* q = f + 5 * (i + 2 * j);
        for (int m = 0; m < 5; ++m) {
            _mm_storel_pd(q + m, out[m][j]);
            _mm_storeh_pd(q + 5 + m, out[m][j]);
        }
    }
}

std::size_t node_functions_sse2(const Params& p, const PeriodConstants& k, const double* z, const double* e,
                                const double* E, double* f, std::size_t n) {
    constexpr int kB = 4;  // four vectors of two interleaved
    std::size_t i = 0;
    for (; i + 2 * kB <= n; i += 2 * kB) block_sse2<kB>(p, k, z, e, E, f, i);
    for (; i + 2 <= n; i += 2) block_sse2<1>(p, k, z, e, E, f, i);
    return i;
}
#endif

Kernel fastest_available() {
    if (available(Kernel::Avx2)) return Kernel::Avx2;
    if (available(Kernel::Sse2)) return Kernel::Sse2;
    return Kernel::Scalar;
}

// VCAL_DERIVS_KERNEL if set (Unavailable for a value that is not a kernel this CPU and build have), else the fastest.
Kernel initial_selection() {
    std::string value;
#if defined(_MSC_VER)
    char* raw = nullptr;
    std::size_t length = 0;
    if (_dupenv_s(&raw, &length, "VCAL_DERIVS_KERNEL") == 0 && raw != nullptr) {
        value = raw;
        std::free(raw);
    }
#else
    if (const char* raw = std::getenv("VCAL_DERIVS_KERNEL")) value = raw;
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
#if defined(VCAL_DERIVS_X86_64)
            return true;
#else
            return false;
#endif
        case Kernel::Avx2: {
#if defined(VCAL_DERIVS_AVX2)
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

bool node_functions(const Params& p, const double* z, const double* e, const double* E, double* f, std::size_t n) {
    const Kernel kern = selected();
    if (kern == Kernel::Unavailable) return false;
    const PeriodConstants k = period_constants(p);
    std::size_t done = 0;
    switch (kern) {
        case Kernel::Scalar:
            break;
        case Kernel::Sse2:
#if defined(VCAL_DERIVS_X86_64)
            done = node_functions_sse2(p, k, z, e, E, f, n);
#endif
            break;
        case Kernel::Avx2:
#if defined(VCAL_DERIVS_AVX2)
            done = node_functions_avx2(p, k, z, e, E, f, n);
#endif
            break;
        case Kernel::Unavailable:
            return false;
    }
    node_functions_scalar(p, k, z, e, E, f, done, n);
    return true;
}

}  // namespace vcal::model::derivs_batch
