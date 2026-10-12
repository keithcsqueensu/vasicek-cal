# P-7: erfcx in vector lanes, bit for bit (D-210)

## Expected, written before the build (2026-10-12)

- **Where the time is** (D-199's S-10 re-profile): `log_phi` is about 70% of a Gauss-Hermite evaluation of the
  binomial term, and `erfcx` about 62% of `log_phi`, so `erfcx` is about 43% of an evaluation. Quadrature
  evaluations of the binomial term dominate S-4b's W*, the analytic SE and the rest of calibrate.
- **The change:** every node of an adaptive Gauss-Hermite rule needs the same `erfcx(|x|/√2)` whichever tail it is
  for, so the nodes' `erfcx` values are computed first, in lanes (SSE2: 2, AVX2: 4), and the scalar code then
  finishes each node and folds the nodes in order as before. The Chebyshev pieces are padded with zero
  coefficients at the top, which leaves Clenshaw's recurrence exactly at +0 until a piece's own first coefficient,
  so lanes in different pieces share one pass. Bit for bit: the fits, surfaces and every downstream number are
  unchanged at every kernel.
- **Expected speed-up of one evaluation (AVX2):** `erfcx` 2.5-3.5× faster, so the evaluation 20-35% faster
  (43% × (1 − 1/2.5 … 1/3.5)). SSE2: smaller, 0-15%. A run dominated by these evaluations (S-10/S-4b) would then be
  about 15-25% shorter; the projector measures it.
- **Out of scope:** the score and Hessian's `inverse_mills` (a later item if this pays), the composite
  Gauss-Legendre rule for d ∈ {0, n}, and the device path (unchanged).
