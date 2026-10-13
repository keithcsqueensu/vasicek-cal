// SPDX-License-Identifier: Apache-2.0
//
// The grid reduction's sums in vector lanes (P-6, D-404): host-only kernels behind reducers::weighted_sums
// (weighted_sum.hpp) and the bounded argmax's tile bounds (resample/bootstrap.hpp).
//
// Three kernels, chosen at run time like the erfcx batch (P-7), the DGP's Bernoulli count (P-17) and the AR(1)
// transition (P-18): scalar, SSE2 (two lanes) and AVX2 (four, only on a CPU that has it). One lane per grid point;
// each lane performs the scalar code's operations in its order (sum_lanes.hpp), so every kernel gives the scalar
// values bit for bit (tests/unit/sum_batch_test.cpp). VCAL_REDUCE_KERNEL=scalar|sse2|avx2 forces one; a value this
// CPU or build lacks selects none, and the calls then refuse (return false), so a caller falls back to its scalar
// loop, and the tests skip visibly (exit code 77).
//
// Declarations only, and C++17: included by device-visible headers, but called only from host code.
#pragma once

#include <cstdint>

namespace vcal::reducers::sum_batch {

enum class Kernel { Scalar, Sse2, Avx2, Unavailable };

bool available(Kernel k);
const char* name(Kernel k);
Kernel selected();      // VCAL_REDUCE_KERNEL if set, else the fastest available; select() overrides
bool select(Kernel k);  // false if k is not available (Unavailable selects none: the calls then refuse)

// out[i] = reducers::weighted_sum(L, T, K, w, k0 + i) for i < n, bit for bit. Returns false, computing nothing, if no
// kernel is selected.
bool weighted_sums(const double* L, std::int64_t T, std::int64_t K, const double* w, std::int64_t k0, std::int64_t n,
                   double* out);

// For each tile i < tiles: s = sum over j < periods with w[j] != 0 of w[j] * U[j * tiles + i], and mag the same sum
// of |w[j] * U[j * tiles + i]|, each in ascending j with plain additions; bound[i] = s + margin * mag. The bounded
// argmax's tile bounds (resample/bootstrap.hpp), bit for bit. Returns false, computing nothing, if no kernel is
// selected.
bool tile_bounds(const double* U, std::int64_t periods, std::int64_t tiles, const double* w, double margin,
                 double* bound);

}  // namespace vcal::reducers::sum_batch
