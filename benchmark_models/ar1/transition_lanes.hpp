// SPDX-License-Identifier: Apache-2.0
//
// The AR(1) transition kernels' bodies (P-18, D-199; the contract is in transition.hpp), generic over a lane type V:
// one double (scalar), two (SSE2) or four (AVX2). Included only by transition.cpp (scalar, SSE2) and
// transition_avx2.cpp (AVX2), each with its own V.
//
// Everything here has internal linkage on purpose, as in dgp/philox_lanes.hpp: transition_avx2.cpp is compiled for
// AVX2 and runs only when the CPU has it, and an inline function with external linkage instantiated there could be
// merged by the linker with another file's copy and run on a CPU without AVX2. So this header includes only data
// (det_exp's generated constants) and <cstddef>; V supplies every operation.
//
// Bit for bit across V by construction: each lane performs exactly the scalar operations, in the same order, on
// the same inputs (no fused multiply-add: the project compiles with contraction off). The exp is det_exp
// (common/detmath/det_exp.hpp) with one change of method that does not change its value: 2^k is applied as
// (p 2^k1) 2^k2 with k1 = floor(k/2), k2 = k - k1, instead of ldexp(p, k). p is in [0.7, 1.42] and |k1| <= 578,
// so p 2^k1 is a normal number and exact, and the second product rounds once: the correctly rounded p 2^k, which
// is ldexp's result, subnormal and overflowing results included. tests/unit/ar1_transition_test.cpp checks it
// against det_exp at every kernel over the whole range, the edges and the subnormal band included.
#pragma once

#include <cstddef>

#include "common/detmath/generated/det_exp_constants.hpp"

namespace vcal::benchmark_models::ar1::kernel {
namespace {

constexpr double kZeroBelow = -750.0;  // transition.hpp's kExpZeroBelow

template <class V>
typename V::T lanes_exp(typename V::T x) {
    namespace gc = detmath::generated;
    using T = typename V::T;
    const auto nan = V::unord(x, x);
    const auto big = V::gt(x, V::set1(800.0));
    const auto small = V::lt(x, V::set1(-800.0));
    const T xc = V::vmax(V::vmin(x, V::set1(801.0)), V::set1(-801.0));  // keeps k in range on every lane
    const T kd = V::floor(V::add(V::mul(xc, V::set1(gc::kInvLn2)), V::set1(0.5)));
    const T r = V::sub(V::sub(xc, V::mul(kd, V::set1(gc::kLn2Hi))), V::mul(kd, V::set1(gc::kLn2Lo)));
    T p = V::set1(gc::kExpCoef[gc::kExpTerms]);
    for (int n = gc::kExpTerms - 1; n >= 0; --n) p = V::add(V::mul(p, r), V::set1(gc::kExpCoef[n]));
    const T k1 = V::floor(V::mul(kd, V::set1(0.5)));
    const T k2 = V::sub(kd, k1);
    T res = V::mul(V::mul(p, V::pow2(k1)), V::pow2(k2));
    res = V::blend(small, V::set1(0.0), res);
    res = V::blend(big, V::pinf(), res);
    return V::blend(nan, x, res);
}

// lanes_exp for arguments known to be <= 0 or -inf and never NaN: the exponent of every term of an online LSE
// (exp(-|d|)) and of the shifted pair weights. The same value as lanes_exp there, with the NaN and overflow cases
// dropped, and 2^k applied as one exact product p 2^k whenever every lane has k >= -1021 (then p 2^k >= 0.7 2^-1021
// is normal, so the product is exact and equals ldexp's result and the two-product form's).
template <class V>
typename V::T lanes_exp_nonpos(typename V::T x) {
    namespace gc = detmath::generated;
    using T = typename V::T;
    const auto small = V::lt(x, V::set1(-800.0));
    const T xc = V::vmax(x, V::set1(-801.0));
    const T kd = V::floor(V::add(V::mul(xc, V::set1(gc::kInvLn2)), V::set1(0.5)));
    const T r = V::sub(V::sub(xc, V::mul(kd, V::set1(gc::kLn2Hi))), V::mul(kd, V::set1(gc::kLn2Lo)));
    T p = V::set1(gc::kExpCoef[gc::kExpTerms]);
    for (int n = gc::kExpTerms - 1; n >= 0; --n) p = V::add(V::mul(p, r), V::set1(gc::kExpCoef[n]));
    T res;
    if (!V::any(V::lt(kd, V::set1(-1021.0)))) {
        res = V::mul(p, V::pow2(kd));
    } else {
        const T k1 = V::floor(V::mul(kd, V::set1(0.5)));
        res = V::mul(V::mul(p, V::pow2(k1)), V::pow2(V::sub(kd, k1)));
    }
    return V::blend(small, V::set1(0.0), res);
}

// Inputs within these bounds make every term finite: |z - pu| <= 2e100, so (z - pu)^2 c stays below 1e202 for
// c <= 1 / (2 (1 - 0.95^2)), and the sums of terms below 1e300 cannot overflow. Checked once per source or
// destination, so the per-term special cases drop out of the common path (P-18).
inline bool bounded(double x, double b) { return x <= b && x >= -b; }  // false for NaN and the infinities

// The online log-sum-exp of transition.hpp, one per lane.
template <class V>
struct LanesLse {
    using T = typename V::T;
    using M = typename V::M;
    T mx = V::ninf();
    T sum = V::set1(0.0);
    M nan = V::none();
    M inf = V::none();
    std::size_t terms = 0, below = 0;  // counted by the scalar kernel only (V::kCounts)

    void add(T x) {
        const M is_nan = V::unord(x, x);
        const M is_pinf = V::eq(x, V::pinf());
        const M is_ninf = V::eq(x, V::ninf());
        nan = V::bor(nan, is_nan);
        inf = V::bor(inf, is_pinf);
        const M ok = V::bnot(V::bor(V::bor(is_nan, is_pinf), is_ninf));
        const T d = V::sub(x, mx);
        const M gt = V::gt(x, mx);
        // exp(d) is exactly +0 and the sum is at least 1: adding it changes nothing (transition.hpp).
        const M skip = V::band(ok, V::band(V::bnot(gt), V::lt(d, V::set1(kZeroBelow))));
        if constexpr (V::kCounts) {
            ++terms;
            if (V::all(skip)) ++below;
        }
        if (V::all(V::bor(skip, V::bnot(ok)))) return;
        const T e = lanes_exp<V>(V::neg_abs(d));
        const T s_gt = V::add(V::mul(sum, e), V::set1(1.0));
        const T s_le = V::add(sum, e);
        sum = V::blend(ok, V::blend(gt, s_gt, s_le), sum);
        mx = V::blend(V::band(ok, gt), x, mx);
    }

    // add() for a term known finite on every lane: no special cases. The same result as add(): a skipped term's exp
    // is exactly +0, and sum + 0 = sum. The scalar kernel still skips it (and counts it); the vector kernels compute
    // it, since a skip only pays when every lane qualifies (0.1-0.6% of terms, P-18).
    void add_finite(T x) {
        const T d = V::sub(x, mx);
        const M gt = V::gt(x, mx);
        if constexpr (V::kCounts) {
            ++terms;
            if (V::all(V::band(V::bnot(gt), V::lt(d, V::set1(kZeroBelow))))) {
                ++below;
                return;
            }
        }
        const T e = lanes_exp_nonpos<V>(V::neg_abs(d));
        sum = V::blend(gt, V::add(V::mul(sum, e), V::set1(1.0)), V::add(sum, e));
        mx = V::blend(gt, x, mx);
    }

    // A term known to be -inf on every lane: no effect (counted as add() would count it).
    void add_neg_inf() {
        if constexpr (V::kCounts) ++terms;
    }

    void store(double* m, double* s, unsigned char* f) const {
        V::store(m, mx);
        V::store(s, sum);
        const unsigned bn = V::bits(nan), bi = V::bits(inf);
        for (int l = 0; l < V::kLanes; ++l) {
            f[l] = static_cast<unsigned char>(((bn >> l) & 1u) * 1u + ((bi >> l) & 1u) * 2u);
        }
    }
};

// Each output's sum is a chain: term j + 1 needs term j's max, and each exp is a chain of 17 multiply-adds. One
// block of lanes would leave the CPU waiting on that latency, so R independent blocks run interleaved (as the
// DGP's kernels interleave registers, P-17). Lanes and blocks never mix, so R changes no result.

// Destinations k0 .. k0 + R kLanes - 1 of the forward sum.
template <class V, int R>
void forward_block(const double* pu, const double* lhat, std::size_t np, const double* z, const double* a, double c,
                   double hlv, double* mx, double* sum, unsigned char* flags, std::size_t* terms, std::size_t* below) {
    using T = typename V::T;
    const T cc = V::set1(c), hv = V::set1(hlv);
    T zk[R], ak[R];
    LanesLse<V> acc[R];
    bool block_ok = bounded(c, 1e300) && bounded(hlv, 1e300);
    for (int q = 0; q < R; ++q) {
        zk[q] = V::load(z + q * V::kLanes);
        ak[q] = V::load(a + q * V::kLanes);
    }
    for (int i = 0; i < R * V::kLanes; ++i) block_ok = block_ok && bounded(z[i], 1e100) && bounded(a[i], 1e300);
    for (std::size_t j = 0; j < np; ++j) {
        const T puj = V::set1(pu[j]), lj = V::set1(lhat[j]);
        const bool pu_ok = block_ok && bounded(pu[j], 1e100);
        if (pu_ok && bounded(lhat[j], 1e300)) {
            for (int q = 0; q < R; ++q) {
                const T d = V::sub(zk[q], puj);
                acc[q].add_finite(V::add(lj, V::sub(ak[q], V::add(V::mul(V::mul(d, d), cc), hv))));
            }
        } else if (pu_ok && lhat[j] < -1e300 && lhat[j] - lhat[j] != 0.0) {  // lhat_j = -inf: every term -inf
            for (int q = 0; q < R; ++q) acc[q].add_neg_inf();
        } else {
            for (int q = 0; q < R; ++q) {
                const T d = V::sub(zk[q], puj);
                acc[q].add(V::add(lj, V::sub(ak[q], V::add(V::mul(V::mul(d, d), cc), hv))));
            }
        }
    }
    for (int q = 0; q < R; ++q) {
        acc[q].store(mx + q * V::kLanes, sum + q * V::kLanes, flags + q * V::kLanes);
        *terms += acc[q].terms;
        *below += acc[q].below;
    }
}

// Sources j0 .. j0 + R kLanes - 1 of the backward sum; e (if not null) points at row j0's first term, stride nn.
template <class V, int R>
void backward_block(const double* pu, const double* z, const double* a, const double* v, const double* b, std::size_t nn,
                    double c, double hlv, double* mx, double* sum, unsigned char* flags, double* e, std::size_t* terms,
                    std::size_t* below) {
    using T = typename V::T;
    const T cc = V::set1(c), hv = V::set1(hlv);
    T puj[R];
    LanesLse<V> acc[R];
    bool block_ok = bounded(c, 1e300) && bounded(hlv, 1e300);
    for (int q = 0; q < R; ++q) puj[q] = V::load(pu + q * V::kLanes);
    for (int i = 0; i < R * V::kLanes; ++i) block_ok = block_ok && bounded(pu[i], 1e100);
    double lane[V::kLanes];
    for (std::size_t k = 0; k < nn; ++k) {
        const T zk = V::set1(z[k]), ak = V::set1(a[k]), vk = V::set1(v[k]), bk = V::set1(b[k]);
        const bool r_ok = block_ok && bounded(z[k], 1e100) && bounded(a[k], 1e300);
        // 0: finite terms; 1: every term -inf (v or b is -inf, neither NaN nor +inf); 2: the general path.
        int cls = 2;
        if (r_ok && bounded(v[k], 1e300) && bounded(b[k], 1e300)) {
            cls = 0;
        } else if (r_ok && v[k] == v[k] && b[k] == b[k] && v[k] < 1e308 && b[k] < 1e308 &&
                   (v[k] < -1e308 || b[k] < -1e308)) {
            cls = 1;
        }
        for (int q = 0; q < R; ++q) {
            const T d = V::sub(zk, puj[q]);
            const T r = V::sub(ak, V::add(V::mul(V::mul(d, d), cc), hv));
            const T x = V::add(V::add(r, vk), bk);
            if (e != nullptr) {
                V::store(lane, x);
                for (int l = 0; l < V::kLanes; ++l) {
                    e[static_cast<std::size_t>(q * V::kLanes + l) * nn + k] = lane[l];
                }
            }
            if (cls == 0) {
                acc[q].add_finite(x);
            } else if (cls == 1) {
                acc[q].add_neg_inf();
            } else {
                acc[q].add(x);
            }
        }
    }
    for (int q = 0; q < R; ++q) {
        acc[q].store(mx + q * V::kLanes, sum + q * V::kLanes, flags + q * V::kLanes);
        *terms += acc[q].terms;
        *below += acc[q].below;
    }
}

// out[i] = det_exp(in[i] - shift) for kLanes values.
template <class V>
void exp_shifted_block(const double* in, double shift, double* out, std::size_t* terms, std::size_t* below) {
    const typename V::T x = V::sub(V::load(in), V::set1(shift));
    const auto zero = V::lt(x, V::set1(kZeroBelow));
    if constexpr (V::kCounts) ++*terms;
    if (V::all(zero)) {
        if constexpr (V::kCounts) ++*below;
        V::store(out, V::set1(0.0));
        return;
    }
    // Pair weights: every argument is at most 0 (the shift is their maximum) unless something is NaN.
    if (!V::any(V::bor(V::unord(x, x), V::gt(x, V::set1(0.0))))) {
        V::store(out, lanes_exp_nonpos<V>(x));
        return;
    }
    V::store(out, lanes_exp<V>(x));
}

// score_terms for destinations k0 .. k0 + kLanes - 1 of source row j.
template <class V>
void score_terms_block(const double* w, double uj, const double* z, double phi, double var, double* out) {
    using T = typename V::T;
    const T u = V::set1(uj), ph = V::set1(phi), vr = V::set1(var);
    const T a = V::sub(V::load(z), V::set1(phi * uj));
    const T t1 = V::div(V::mul(a, u), vr);
    const T t2 = V::div(V::mul(ph, V::mul(a, a)), V::set1(var * var));
    V::store(out, V::mul(V::load(w), V::add(V::sub(t1, t2), V::set1(phi / var))));
}

template <class V>
void exp_block(const double* in, double* out) {
    V::store(out, lanes_exp<V>(V::load(in)));
}

}  // namespace
}  // namespace vcal::benchmark_models::ar1::kernel
