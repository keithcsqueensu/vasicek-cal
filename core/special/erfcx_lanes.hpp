// SPDX-License-Identifier: Apache-2.0
//
// The erfcx lane kernel's body (P-7, D-210; the contract is in erfcx_batch.hpp), generic over a lane type V: two
// doubles (SSE2) or four (AVX2). Included only by erfcx_batch.cpp (SSE2) and erfcx_batch_avx2.cpp (AVX2), each with
// its own V.
//
// Internal linkage on purpose, as in dgp/philox_lanes.hpp: erfcx_batch_avx2.cpp is compiled for AVX2 and runs only
// when the CPU has it, and an inline function with external linkage instantiated there could be merged by the linker
// with another file's copy and run on a CPU without AVX2. So this header includes only data (the generated table).
//
// Bit for bit with generated::erfcx_nonnegative (erfcx_chebyshev.hpp) by construction, lane by lane:
//   - the piece is the scalar code's: x <= 1, else x <= 3, else x <= 8, else the tail;
//   - t is its expression for that piece: (2x - centre2) / span, or 2 (8 / x) - 1 for the tail;
//   - Clenshaw's recurrence in chebyshev_eval's order, b0 = (c_k + (2t) b1) - b2 for k = width-1 .. 1, then
//     (c_0 + t b1) - b2, with the lane's own piece's coefficients selected exactly (a blend, not arithmetic). The
//     table pads each piece with +0 at its top, and a zero coefficient keeps b1 = b2 = +0 exactly while t is
//     finite, so a padded piece gives its unpadded series' value;
//   - the tail is divided by x.
// No fused multiply-add: these files compile with contraction off and strict FP (CMakeLists.txt).
#pragma once

#include "core/special/generated/erfcx_table.hpp"

namespace vcal::special::erfcx_batch {
namespace {

// B independent vectors at once: Clenshaw's recurrence is a chain of dependent operations, so one vector alone waits
// on latency; B chains interleaved keep the units busy. Each vector's operations are those of the scalar code, in
// its order, so B changes nothing in the values.
template <class V, int B>
void lanes_erfcx_nonnegative(const typename V::T (&x)[B], typename V::T (&out)[B]) {
    namespace g = vcal::special::generated;
    using T = typename V::T;
    using M = typename V::M;
    static_assert(g::kErfcxPieces == 4, "erfcx_lanes.hpp assumes three segments and a tail");
    M m1[B], m3[B], m8[B];
    T t[B], t2[B], b1[B], b2[B];
    const T two = V::set1(2.0);
    for (int j = 0; j < B; ++j) {
        m1[j] = V::le(x[j], V::set1(g::kErfcxUpper[0]));
        m3[j] = V::le(x[j], V::set1(g::kErfcxUpper[1]));
        m8[j] = V::le(x[j], V::set1(g::kErfcxUpper[2]));
        const T centre2 = V::blend(m1[j], V::set1(g::kErfcxCentre2[0]),
                                   V::blend(m3[j], V::set1(g::kErfcxCentre2[1]), V::set1(g::kErfcxCentre2[2])));
        const T span =
            V::blend(m1[j], V::set1(g::kErfcxSpan[0]), V::blend(m3[j], V::set1(g::kErfcxSpan[1]), V::set1(g::kErfcxSpan[2])));
        const T t_seg = V::div(V::sub(V::mul(two, x[j]), centre2), span);
        const T t_tail = V::sub(V::mul(two, V::div(V::set1(g::kErfcxTailX0), x[j])), V::set1(1.0));
        t[j] = V::blend(m8[j], t_seg, t_tail);
        t2[j] = V::mul(two, t[j]);
        b1[j] = V::set1(0.0);
        b2[j] = V::set1(0.0);
    }
    // When every lane of the block is in one piece (consecutive quadrature nodes usually are), its coefficients are
    // broadcast; otherwise each lane's own piece's coefficient is selected by blends. Both give each lane exactly
    // its piece's coefficient, so the two paths agree bit for bit.
    int piece = -1;  // the block's single piece, or -1 if mixed
    {
        bool all1 = true, all3 = true, all8 = true, none1 = true, none3 = true, none8 = true;
        for (int j = 0; j < B; ++j) {
            all1 = all1 && V::all(m1[j]);
            all3 = all3 && V::all(m3[j]);
            all8 = all8 && V::all(m8[j]);
            none1 = none1 && V::none(m1[j]);
            none3 = none3 && V::none(m3[j]);
            none8 = none8 && V::none(m8[j]);
        }
        if (all1) piece = 0;
        else if (none1 && all3) piece = 1;
        else if (none3 && all8) piece = 2;
        else if (none8) piece = 3;
    }
    if (piece >= 0) {
        const double* c = g::kErfcxCoef[piece];
        for (int k = g::kErfcxWidth - 1; k > 0; --k) {
            const T ck = V::set1(c[k]);
            for (int j = 0; j < B; ++j) {
                const T b0 = V::sub(V::add(ck, V::mul(t2[j], b1[j])), b2[j]);
                b2[j] = b1[j];
                b1[j] = b0;
            }
        }
        for (int j = 0; j < B; ++j) {
            const T r = V::sub(V::add(V::set1(c[0]), V::mul(t[j], b1[j])), b2[j]);
            out[j] = piece < 3 ? r : V::div(r, x[j]);
        }
        return;
    }
    auto coef = [&](int j, int k) {
        return V::blend(m1[j], V::set1(g::kErfcxCoef[0][k]),
                        V::blend(m3[j], V::set1(g::kErfcxCoef[1][k]),
                                 V::blend(m8[j], V::set1(g::kErfcxCoef[2][k]), V::set1(g::kErfcxCoef[3][k]))));
    };
    for (int k = g::kErfcxWidth - 1; k > 0; --k) {
        for (int j = 0; j < B; ++j) {
            const T b0 = V::sub(V::add(coef(j, k), V::mul(t2[j], b1[j])), b2[j]);
            b2[j] = b1[j];
            b1[j] = b0;
        }
    }
    for (int j = 0; j < B; ++j) {
        const T r = V::sub(V::add(coef(j, 0), V::mul(t[j], b1[j])), b2[j]);
        out[j] = V::blend(m8[j], r, V::div(r, x[j]));
    }
}

}  // namespace
}  // namespace vcal::special::erfcx_batch
