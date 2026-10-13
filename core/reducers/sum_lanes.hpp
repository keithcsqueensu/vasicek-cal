// SPDX-License-Identifier: Apache-2.0
//
// The grid reduction's sums on kB vectors of lanes at once (P-6, D-404; sum_batch.hpp), for a lane type V with
// set1, load, store, add, sub, mul, abs, ge, isfinite, not_, or_, none, none_mask and blend(m, a, b) = m ? a : b.
// One lane per grid point; each lane performs the scalar code's operations in their order on the same operands:
//   - weighted_sum (weighted_sum.hpp), on finite terms: x = w[t] * L[t][k], s = sum + x,
//     comp = comp + ((|sum| >= |x|) ? (sum - s) + x : (x - s) + sum), sum = s; the result is sum + comp. Both
//     branches are computed and the scalar's selected exactly. A block that meets a non-finite term is redone by the
//     caller with the scalar code (see below).
//   - the tile bounds: s = s + x and mag = mag + |x| with x = w[j] * U[j][tile], then s + margin * mag.
// A period with w[t] == 0 is skipped by every lane alike (w is shared). kB vectors are interleaved so the add
// latency of one chain hides behind the others.
//
// Included only by the kernel translation units (sum_batch.cpp, sum_batch_avx2.cpp). Internal linkage on purpose,
// as erfcx_lanes.hpp: sum_batch_avx2.cpp is compiled for AVX2 and runs only when the CPU has it, so nothing
// instantiated there may be merged by the linker with another file's copy. This header includes nothing.
#pragma once

namespace vcal::reducers::sum_batch {
namespace {

// L points at L[0 * K + k0]; the lanes' grid points are k0 .. k0 + kB * width - 1, contiguous. The finite path only:
// a lane's sum and compensation change exactly as the scalar's do while every term is finite. Returns false if any
// lane met a non-finite term, whose handling (a separate sum, no compensation) these lanes do not carry: the caller
// then recomputes the block with the scalar code, so the results are the scalar's either way. Non-finite terms are
// rare (a zero-default period at an extreme grid point), so the lanes stay small: three registers per vector.
template <class V, int kB>
inline bool lanes_weighted_sums(const double* L, long long T, long long K, const double* w, double* out) {
    using T_ = typename V::T;
    using M = typename V::M;
    constexpr int W = V::width;
    T_ sum[kB], comp[kB];
    M bad = V::none_mask();
    for (int j = 0; j < kB; ++j) {
        sum[j] = V::set1(0.0);
        comp[j] = V::set1(0.0);
    }
    for (long long t = 0; t < T; ++t) {
        if (w[t] == 0.0) continue;
        const T_ wt = V::set1(w[t]);
        const double* row = L + t * K;
        for (int j = 0; j < kB; ++j) {
            const T_ x = V::mul(wt, V::load(row + j * W));
            bad = V::or_(bad, V::not_(V::isfinite(x)));
            const T_ s = V::add(sum[j], x);
            const T_ br = V::blend(V::ge(V::abs(sum[j]), V::abs(x)), V::add(V::sub(sum[j], s), x),
                                   V::add(V::sub(x, s), sum[j]));
            comp[j] = V::add(comp[j], br);
            sum[j] = s;
        }
    }
    for (int j = 0; j < kB; ++j) V::store(out + j * W, V::add(sum[j], comp[j]));
    return V::none(bad);
}

// U points at U[0 * tiles + i0]; the lanes' tiles are i0 .. i0 + kB * width - 1, contiguous.
template <class V, int kB>
inline void lanes_tile_bounds(const double* U, long long periods, long long tiles, const double* w, double margin,
                              double* bound) {
    using T_ = typename V::T;
    constexpr int W = V::width;
    T_ s[kB], mag[kB];
    for (int j = 0; j < kB; ++j) {
        s[j] = V::set1(0.0);
        mag[j] = V::set1(0.0);
    }
    for (long long p = 0; p < periods; ++p) {
        if (w[p] == 0.0) continue;
        const T_ wp = V::set1(w[p]);
        const double* row = U + p * tiles;
        for (int j = 0; j < kB; ++j) {
            const T_ x = V::mul(wp, V::load(row + j * W));
            s[j] = V::add(s[j], x);
            mag[j] = V::add(mag[j], V::abs(x));
        }
    }
    const T_ m = V::set1(margin);
    for (int j = 0; j < kB; ++j) V::store(bound + j * W, V::add(s[j], V::mul(m, mag[j])));
}

}  // namespace
}  // namespace vcal::reducers::sum_batch
