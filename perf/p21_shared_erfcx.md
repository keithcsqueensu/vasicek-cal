# P-21: one erfcx per node on the derivative path, and in the hint (D-400)

## Expected, written before the build (2026-10-12)

- **Where the time is:** after P-7 (D-212) the value's nodes take `erfcx` from one vector batch, but
  `log_contrib_derivs` was only 0-11% faster, because the derivative functor's `inverse_mills(x)` and
  `inverse_mills(-x)` recompute, per node, the same `erfcx(|x|/√2)` (twice) and the same exp pair the value just
  computed. The hint's Newton (`mixture_slopes`) calls `inverse_mills` four times per step at 0 < d < n (twice
  directly, twice inside `inverse_mills_slope`): four `erfcx` and two exp pairs for one pair of λ values. The SIMD
  sweep (D-213) put the hint at 10-30% of an interior period.
- **The change:** the value's batch leaves each node's `erfcx` and exp pair in a per-call cache, and the derivative
  functor finishes λ(±x) from them with `inverse_mills`'s own expression for each sign; the hint computes λ(x),
  λ(−x) as one pair from one `erfcx` and one exp pair, and λ' from those λ values. Bit for bit: every value is the
  expression it replaces, on the same operands.
- **Expected:** `log_contrib_derivs` 15-30% faster at 0 < d < n; `log_contrib` 5-15% faster (the hint's share).
  Measured with P-7's benchmark (`p7_erfcx`) against P-7's branch, AVX2 and SSE2.

## What was built

- **`inverse_mills_given_erfcx(x, e, E)` and `inverse_mills_pair(x)`** (`core/special/inverse_mills.hpp`):
  `inverse_mills`'s own expression for x's sign, given e = `erfcx(|x|/√2)` and E = `exp_neg_half_square(x)`; the
  pair takes one `erfcx` and one exp pair for λ(x) and λ(−x).
- **The node cache** (`model::NodeCache`, host only): with P-7's batch, `BinomialLogIntegrand::batch` leaves each
  node's z, e and (at 0 < d < n) E there; `log_phi_pair_given_erfcx` hands out the exp pair it computes.
  `log_contrib_derivs` makes one per call outside device code, and the derivative functor finds each node by a
  forward cursor (it is called in node order, skipping underflowed nodes). A node without a cached e (no kernel,
  a non-finite |x|, the composite rule) takes `inverse_mills_pair` or `inverse_mills`.
- **The hint:** `mixture_slopes` takes λ(±x) once, from `inverse_mills_pair` at 0 < d < n, and λ' from those values
  (`inverse_mills_slope_given`); before, four `inverse_mills` calls per Newton step.

## Identity (bit for bit)

- `unit_erfcx_batch`: `inverse_mills_pair` and the slope against two `inverse_mills` calls and the old slope on
  1,600,807 arguments (every erfcx piece, both signs, zeros, infinities, NaN); the derivative functor with and without
  the cache against a verbatim copy of the replaced functor on 57,512 nodes (d = 0, d = n, 0 < d < n, skipped
  nodes, an overflowing |x|), under every kernel; and P-7's `log_contrib` checks.
- **Against `main` before P-7, as built code:** a dump of 66,000 values (each period's hint mode, scale, centre and
  width, `log_contrib`, `log_contrib_derivs` with its score and Hessian; 6,000 random periods, the parity rule and its
  doubled rule) from `main` (26d8860... the P-7 base) and from this branch under each kernel: 0 differences.
- The committed goldens and the full non-slow suite (GCC 13, `-Werror`).

## Measured (2026-10-12)

Same machine and method as P-7 (cloud Xeon 2.8 GHz, GCC 13, one thread, `p7_erfcx --seconds 3`, three runs). "P-7" is
P-7's branch head built the same way.

| µs per period | P-7, SSE2 | P-7, AVX2 | P-21, SSE2 | P-21, AVX2 |
|---|---|---|---|---|
| `log_contrib` | 13.6-14.2 | 14.1-14.4 | 11.7-11.9 | 12.0-12.3 |
| `log_contrib_derivs` | 41.8-43.5 | 46.7-48.9 | 14.2-15.1 | 14.8-15.0 |

- **`log_contrib`: 15% faster** (the hint), at the top of the 5-15% expected.
- **`log_contrib_derivs`: 2.9-3.3× faster,** far above the 15-30% expected. The expectation took the SIMD sweep's
  27-30 ns per `erfcx` (D-213); on this machine a scalar `erfcx` costs about 50 ns, and each node's two
  `inverse_mills` calls also paid an exp pair and a division chain, so the derivative functor was most of the
  evaluation. It is now about the value plus the node functions' arithmetic. Checked bit for bit against `main` as
  above, because a gain this size could otherwise be a skipped computation.
- **Against `main`** (P-7 and P-21 together): `log_contrib` 21.6-22.8 → 11.7-12.3 µs (about 45% less time),
  `log_contrib_derivs` 48.1-49.6 → 14.2-15.1 µs (about 70% less). The derivative path is W* (40% of S-10), the
  analytic SE (13.5%) and the refinement; the next projection (D-187) measures the run.
- With the derivative functor no longer dominated by `inverse_mills`, AVX2 and SSE2 are now within a few percent on
  both measures.
