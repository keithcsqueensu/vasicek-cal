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

    void store(double* m, double* s, unsigned char* f) const {
        V::store(m, mx);
        V::store(s, sum);
        const unsigned bn = V::bits(nan), bi = V::bits(inf);
        for (int l = 0; l < V::kLanes; ++l) {
            f[l] = static_cast<unsigned char>(((bn >> l) & 1u) * 1u + ((bi >> l) & 1u) * 2u);
        }
    }
};

// Destinations k0 .. k0 + kLanes - 1 of the forward sum.
template <class V>
void forward_block(const double* pu, const double* lhat, std::size_t np, const double* z, const double* a, double c,
                   double hlv, double* mx, double* sum, unsigned char* flags, std::size_t* terms, std::size_t* below) {
    using T = typename V::T;
    const T zk = V::load(z), ak = V::load(a), cc = V::set1(c), hv = V::set1(hlv);
    LanesLse<V> acc;
    for (std::size_t j = 0; j < np; ++j) {
        const T d = V::sub(zk, V::set1(pu[j]));
        const T r = V::sub(ak, V::add(V::mul(V::mul(d, d), cc), hv));
        acc.add(V::add(V::set1(lhat[j]), r));
    }
    acc.store(mx, sum, flags);
    *terms += acc.terms;
    *below += acc.below;
}

// Sources j0 .. j0 + kLanes - 1 of the backward sum; e (if not null) points at row j0's first term, stride nn.
template <class V>
void backward_block(const double* pu, const double* z, const double* a, const double* v, const double* b, std::size_t nn,
                    double c, double hlv, double* mx, double* sum, unsigned char* flags, double* e, std::size_t* terms,
                    std::size_t* below) {
    using T = typename V::T;
    const T puj = V::load(pu), cc = V::set1(c), hv = V::set1(hlv);
    LanesLse<V> acc;
    double lane[V::kLanes];
    for (std::size_t k = 0; k < nn; ++k) {
        const T d = V::sub(V::set1(z[k]), puj);
        const T r = V::sub(V::set1(a[k]), V::add(V::mul(V::mul(d, d), cc), hv));
        const T x = V::add(V::add(r, V::set1(v[k])), V::set1(b[k]));
        if (e != nullptr) {
            V::store(lane, x);
            for (int l = 0; l < V::kLanes; ++l) e[static_cast<std::size_t>(l) * nn + k] = lane[l];
        }
        acc.add(x);
    }
    acc.store(mx, sum, flags);
    *terms += acc.terms;
    *below += acc.below;
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
    V::store(out, lanes_exp<V>(x));
}

template <class V>
void exp_block(const double* in, double* out) {
    V::store(out, lanes_exp<V>(V::load(in)));
}

}  // namespace
}  // namespace vcal::benchmark_models::ar1::kernel
