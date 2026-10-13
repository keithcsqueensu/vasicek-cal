// SPDX-License-Identifier: Apache-2.0
//
// erfcx on many arguments at once, in vector lanes (P-7, D-212): host-only kernels behind the scalar erfcx of
// erfcx.hpp, for the quadrature's nodes. Every node of a Gauss-Hermite rule needs erfcx(|x| / sqrt 2) whichever
// tail log Phi takes it for (log_phi.hpp), so the binomial log-integrand computes its nodes' erfcx values first,
// here, and finishes each node with the scalar code (binomial_mixture.hpp).
//
// Three kernels, chosen at run time like the DGP's Bernoulli count (P-17) and the AR(1) transition (P-18): scalar,
// SSE2 (two lanes) and AVX2 (four, only on a CPU that has it). The scalar kernel is generated::erfcx_nonnegative
// itself; each lane of the others performs its operations in its order on the same coefficients (erfcx_lanes.hpp),
// so every kernel gives the same values bit for bit (tests/unit/erfcx_batch_test.cpp). VCAL_ERFCX_KERNEL=
// scalar|sse2|avx2 forces one; a value this CPU or build lacks selects none, and nonnegative() then refuses
// (returns false), so a caller falls back to its scalar path, and the tests skip visibly (exit code 77).
//
// Declarations only, and C++17: included by device-visible headers, but called only from host code.
#pragma once

#include <cstddef>

namespace vcal::special::erfcx_batch {

enum class Kernel { Scalar, Sse2, Avx2, Unavailable };

bool available(Kernel k);
const char* name(Kernel k);
Kernel selected();      // VCAL_ERFCX_KERNEL if set, else the fastest available; select() overrides
bool select(Kernel k);  // false if k is not available (Unavailable selects none: nonnegative() then refuses)

// out[i] = generated::erfcx_nonnegative(x[i]) bit for bit, for finite x[i] >= 0 (+0 and -0 included). A NaN or
// infinite x[i] gives an unspecified out[i]: callers keep those for the scalar path. Returns false, computing
// nothing, if no kernel is selected.
bool nonnegative(const double* x, double* out, std::size_t n);

}  // namespace vcal::special::erfcx_batch
