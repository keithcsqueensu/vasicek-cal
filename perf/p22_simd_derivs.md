# P-22: the derivative node functions in vector lanes (D-402)

## Expected, written before the build (2026-10-12)

- **Where the time is:** after P-21 (D-400), `log_contrib_derivs` costs the value plus about 2.5 µs per period
  (P-7's benchmark `p7_erfcx`, SSE2: 12.0 µs value, 14.5 µs with derivatives, 81 interior periods of the
  recovery grid, cloud Xeon, GCC 13, one thread). A microbenchmark of `BinomialLogIntegrandDerivs` alone on a
  filled node cache (128 nodes, n = 1,000, d = 10) gives 21.7 ns per node, 2.8 µs per period: the five node
  functions are all of the difference. Their own moment updates in the fold (`m_[j] += w·f[j]`) are about 4.9 ns
  per node and stay scalar.
- **The change:** once per quadrature call, the derivative functor computes every cached node's five node
  functions in vector lanes (two in SSE2, four in AVX2), each lane doing the scalar functor's operations in its
  order: x, the inverse-Mills pair from the cached e and E (both branches, blended by x's sign), g_x and g_xx,
  the x derivatives and the five functions. The existing fold then reads them in node order; the integrator and
  the accumulator do not change. Exact divisions throughout, no reciprocals, no contraction.
- **Expected:** the node functions 2.8 → 1.0-1.5 µs per period (seven divisions per node dominate, and vector
  division has only about twice the scalar throughput), so `log_contrib_derivs` **8-12% faster** at 0 < d < n,
  AVX2 a little ahead of SSE2; `log_contrib` unchanged. The composite rule's nodes (d ∈ {0, n}) are not cached
  (P-21) and do not benefit.
- **Stop rule:** if `log_contrib_derivs` gains under 3% on both kernels, P-22 is recorded as measured and not
  worth its code, and is not merged.
- **Bit for bit:** every lane value is the scalar expression on the same operands; scalar, SSE2 and AVX2 kernels
  are each checked against the scalar functor, with `VCAL_DERIVS_KERNEL` to force one, as P-7 and P-17.
