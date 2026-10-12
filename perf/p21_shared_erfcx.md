# P-21: one erfcx per node on the derivative path, and in the hint (D-214)

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
