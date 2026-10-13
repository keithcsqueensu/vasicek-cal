# P-22: the derivative node functions in vector lanes (D-402): rejected

**Rejected (owner, 2026-10-13; D-402).** Built bit for bit and measured: 0-3% on `log_contrib_derivs` against the
8-12% expected, inside the noise. The code is kept on the branch `p22-simd-derivs` for reference and is not merged.

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

## What was built (branch `p22-simd-derivs`, kept for reference, not merged)

- `core/model/derivs_batch.{hpp,cpp}`, `derivs_batch_avx2.cpp`, `derivs_lanes.hpp`: the node functions on kB
  interleaved vectors, scalar, SSE2 and AVX2 kernels, `VCAL_DERIVS_KERNEL`, in `vcal_core_kernels`.
- The hook: `BinomialLogIntegrandDerivs` fills the node cache's functions at its first call after a batch (the E
  of one-sided periods by its own expression first) and reads node i's five values back; nodes without a cached e
  keep the per-node path.
- **Identity:** `unit_derivs_batch`, every kernel against the functor's own path under every erfcx kernel and
  none (69,681 cached nodes and 63,874 nodes through the functor per combination, 0 differences), and
  `log_contrib_derivs` on 400 periods with the parity rule and its doubled rule, 0 differences.

## Measured (2026-10-12): the expectation missed; the stop rule borderline

Cloud Xeon 2.8 GHz, GCC 13, one thread, Release.

- **The kernel alone** (the node functions of 128 cached nodes, n = 1,000, d = 10): scalar 16.5-16.9 ns per node,
  SSE2 8.7-8.8, AVX2 7.2-8.9. About 2× as expected.
- **`log_contrib_derivs`,** 36 interior periods (PD 0.3-10%, ρ 0.02-0.24, n 100-10,000), the parity rule, the three
  derivative paths interleaved in one process (A B C C B A, 30 rounds), median µs per period:

  | erfcx kernel | none (P-21's path) | derivs SSE2 | derivs AVX2 |
  |---|---|---|---|
  | AVX2 | 15.34 | 15.15 (−1.3%) | 14.84 (−3.2%) |
  | SSE2 | 14.91 | 14.63 (−1.9%) | 17.78 (+19%) |

  The interquartile ranges overlap in every row, and P-7's benchmark (`p7_erfcx`, three runs per kernel) shows
  no consistent difference either.
- **Why the kernel's 2× did not carry over:** in the per-node path the node functions' divisions overlap with the
  fold's `exp` and the cursor search (out-of-order execution hides much of their latency), so they cost less in
  place than alone. The batch saves about 1.2 µs of arithmetic per period but adds a separate pass that stores and
  reloads 640 values and loses that overlap. Storing the values node-major (one node's five together) instead of
  by function made no measurable difference. The AVX2 kernel next to the SSE2 erfcx kernel is slower, probably a
  transition or frequency effect between the two instruction sets.
- **Against the expectation:** missed. 8-12% was expected; the measured gain is 0-3%.
- **Against the stop rule (under 3% on both kernels):** borderline, not cleanly met. SSE2 is under 3% (−1.9%). AVX2
  is −3.2% in the interleaved benchmark, just over the line with overlapping interquartile ranges, and about 0 in
  `p7_erfcx`. And the AVX2 derivative kernel beside the SSE2 erfcx kernel is 19% slower, a combination the default
  selection never makes (both pick AVX2 where it exists) but `VCAL_*_KERNEL` overrides can.
- **Decision (owner, 2026-10-13): rejected.** A gain of at most 3%, inside the noise, does not pay for
  about 300 lines and a third kernel family to maintain. The derivative path's remaining cost is the fold and the
  value's own nodes; a real gain there needs the fold itself in lanes (its `exp` is libm, D-199).
