// SPDX-License-Identifier: Apache-2.0
//
// The derivative node functions on kB vectors of lanes at once (P-22, D-402; derivs_batch.hpp), for a lane type V
// with set1, add, sub, mul, div, neg (sign flip), lt, gt and blend(m, a, b) = m ? a : b. Each lane performs
// BinomialLogIntegrandDerivs::operator()'s operations (binomial_mixture.hpp) in its order on the same operands; the
// per-period constants are its per-node subexpressions that do not depend on the node, computed once by the same
// expressions. kB vectors are interleaved so the divisions of one do not wait on another's.
//
// Included only by the kernel translation units (derivs_batch.cpp, derivs_batch_avx2.cpp). Internal linkage on
// purpose, as erfcx_lanes.hpp: derivs_batch_avx2.cpp is compiled for AVX2 and runs only when the CPU has it, so
// nothing instantiated there may be merged by the linker with another file's copy. This header includes only data.
#pragma once

#include "core/model/derivs_batch.hpp"
#include "core/special/generated/constants.hpp"

namespace vcal::model::derivs_batch {

// The operator()'s node-independent subexpressions, each by its expression there (period_constants() in
// derivs_batch.cpp).
struct PeriodConstants {
    double two_ab;     // 2.0 * a * b
    double two_b2;     // 2.0 * b2, b2 = b * b
    double x1;         // c1 / b
    double x11;        // c2 / b
    double x12;        // c1 / (2.0 * b2 * b)
    double one_m_2r;   // 1.0 - 2.0 * rho
    double four_ab3;   // 4.0 * ab * ab * ab, ab = a * b
    double two_b2b2;   // 2.0 * b2 * b2
};

namespace {

template <class V, int kB>
inline void lanes_node_functions(const Params& p, const PeriodConstants& k, const typename V::T* z,
                                 const typename V::T* e, const typename V::T* E, typename V::T (*f)[kB]) {
    using T = typename V::T;
    const T c = V::set1(p.c), a = V::set1(p.a), b = V::set1(p.b);
    const T zero = V::set1(0.0);
    const T k_neg = V::set1(special::constants::kSqrt2OverPi);  // inverse_mills_given_erfcx, x < 0
    const T k_pos = V::set1(special::constants::kInvSqrt2Pi);   // and x >= 0
    const T one = V::set1(1.0), half = V::set1(0.5);
    const T dflt = V::set1(p.defaults), surv = V::set1(p.survivors);
    const T two_ab = V::set1(k.two_ab), two_b2 = V::set1(k.two_b2), x1 = V::set1(k.x1), x11 = V::set1(k.x11);
    const T x12 = V::set1(k.x12), one_m_2r = V::set1(k.one_m_2r), four_ab3 = V::set1(k.four_ab3);
    const T two_b2b2 = V::set1(k.two_b2b2);
    const bool has_d = p.defaults > 0.0, has_s = p.survivors > 0.0;
    T x[kB], gx[kB], gxx[kB];
    for (int j = 0; j < kB; ++j) {
        x[j] = V::div(V::sub(c, V::mul(a, z[j])), b);  // (c - a * z) / b
        // inverse_mills_given_erfcx(t, e, E): kSqrt2OverPi / e if t < 0, else kInvSqrt2Pi * E / (1 - 0.5 * e * E);
        // the same two values serve t = x and t = -x.
        const T neg_branch = V::div(k_neg, e[j]);
        const T pos_branch = V::div(V::mul(k_pos, E[j]), V::sub(one, V::mul(V::mul(half, e[j]), E[j])));
        gx[j] = zero;
        gxx[j] = zero;
        if (has_d) {
            const T lam = V::blend(V::lt(x[j], zero), neg_branch, pos_branch);  // lambda(x)
            gx[j] = V::add(gx[j], V::mul(dflt, lam));
            gxx[j] = V::add(gxx[j], V::mul(dflt, V::mul(V::neg(lam), V::add(x[j], lam))));
        }
        if (has_s) {
            const T lam = V::blend(V::gt(x[j], zero), neg_branch, pos_branch);  // lambda(-x): -x < 0 iff x > 0
            gx[j] = V::sub(gx[j], V::mul(surv, lam));
            gxx[j] = V::add(gxx[j], V::mul(surv, V::mul(V::neg(lam), V::sub(lam, x[j]))));
        }
    }
    for (int j = 0; j < kB; ++j) {
        const T x2 = V::add(V::div(V::neg(z[j]), two_ab), V::div(x[j], two_b2));
        const T x22 = V::add(V::add(V::div(V::mul(z[j], one_m_2r), four_ab3), V::div(x2, two_b2)),
                             V::div(x[j], two_b2b2));
        const T g1 = V::mul(gx[j], x1), g2 = V::mul(gx[j], x2);
        f[0][j] = g1;
        f[1][j] = g2;
        f[2][j] = V::add(V::add(V::mul(V::mul(gxx[j], x1), x1), V::mul(gx[j], x11)), V::mul(g1, g1));
        f[3][j] = V::add(V::add(V::mul(V::mul(gxx[j], x2), x2), V::mul(gx[j], x22)), V::mul(g2, g2));
        f[4][j] = V::add(V::add(V::mul(V::mul(gxx[j], x1), x2), V::mul(gx[j], x12)), V::mul(g1, g2));
    }
}

}  // namespace

}  // namespace vcal::model::derivs_batch
