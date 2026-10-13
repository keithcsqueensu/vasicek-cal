// SPDX-License-Identifier: Apache-2.0
//
// The binomial log-integrand's five derivative node functions on many nodes at once, in vector lanes (P-22, D-402):
// host-only kernels behind BinomialLogIntegrandDerivs (binomial_mixture.hpp). After P-7 and P-21 the value's node
// batch leaves each node's e = erfcx(|x| / sqrt 2) and exp pair E in a per-call node cache; what the derivative
// functor still does per node is arithmetic on (z, e, E). The functor hands the whole cache to node_functions()
// once per quadrature call and reads each node's five values back in node order.
//
// Three kernels, chosen at run time like the erfcx batch (P-7) and the AR(1) transition (P-18): scalar, SSE2 (two
// lanes) and AVX2 (four, only on a CPU that has it). Each lane of every kernel performs the scalar functor's
// operations in its order (derivs_lanes.hpp): exact divisions, no reciprocals, no contraction, the inverse-Mills
// pair's two branches both computed and the one inverse_mills would take selected by x's sign. So every kernel
// gives the functor's values bit for bit (tests/unit/derivs_batch_test.cpp). VCAL_DERIVS_KERNEL=scalar|sse2|avx2
// forces one; a value this CPU or build lacks selects none, and node_functions() then refuses (returns false), so the
// functor keeps its per-node path, and the tests skip visibly (exit code 77).
//
// Declarations only, and C++17: included by device-visible headers, but called only from host code.
#pragma once

#include <cstddef>

namespace vcal::model::derivs_batch {

enum class Kernel { Scalar, Sse2, Avx2, Unavailable };

bool available(Kernel k);
const char* name(Kernel k);
Kernel selected();      // VCAL_DERIVS_KERNEL if set, else the fastest available; select() overrides
bool select(Kernel k);  // false if k is not available (Unavailable selects none: node_functions() then refuses)

// BinomialLogIntegrandDerivs's fields (binomial_mixture.hpp).
struct Params {
    double c, a, b, rho;
    double c1, c2;
    double defaults, survivors;
};

// For i < n: f[5 i + j] = the functor's j-th node function at z[i], given its cached e[i] and E[i], bit for bit, for
// nodes whose e[i] is not NaN (others get unspecified values: the functor keeps its own path for them). E[i] must be
// set wherever the functor's inverse_mills_given_erfcx would read it (its non-negative argument). Returns false,
// computing nothing, if no kernel is selected.
bool node_functions(const Params& p, const double* z, const double* e, const double* E, double* f, std::size_t n);

}  // namespace vcal::model::derivs_batch
