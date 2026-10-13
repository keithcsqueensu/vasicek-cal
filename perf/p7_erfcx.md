# P-7: erfcx in vector lanes, bit for bit (D-212)

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

## What was built

- **`core/special/erfcx_batch.{hpp,cpp}`, `erfcx_batch_avx2.cpp`, `erfcx_lanes.hpp`:** `erfcx_batch::nonnegative(x,
  out, n)`, the scalar kernel (`generated::erfcx_nonnegative` itself) and lanes on SSE2 (2) and AVX2 (4), chosen at
  run time (`VCAL_ERFCX_KERNEL=scalar|sse2|avx2`, as P-17 and P-18). core/'s first compiled part, a host-only static
  library (`vcal_core_kernels`) linked through `vcal_core`, with strict FP and no contraction; only
  `erfcx_batch_avx2.cpp` is compiled for AVX2. Device code never calls it.
- **The table:** `tools/gen_special_tables.py` now also writes `core/special/generated/erfcx_table.hpp`, the same
  coefficients as `erfcx_chebyshev.hpp` in one array, each piece padded with +0 at its top. Every other generated
  file is byte-identical (`--check`).
- **The lanes,** per lane exactly the scalar code: the scalar piece test, its t, Clenshaw in `chebyshev_eval`'s order,
  the tail's division. Two speed measures that change no value: four independent vectors are interleaved (the
  recurrence is a dependent chain, so one vector waits on latency), and a block whose lanes all lie in one piece
  broadcasts that piece's coefficients instead of blending four.
- **The hook:** `BinomialLogIntegrand::batch(z, out, n)` computes the nodes' `erfcx(|x| / √2)` in one call and
  finishes each node with `log_phi_given_erfcx` / `log_phi_pair_given_erfcx`, which are `log_phi`'s and
  `log_phi_pair`'s own expressions (both now call them). `GaussHermiteFixed` and `GaussHermiteAdaptive` use it
  outside device code for rules of up to 256 nodes (the parity check's doubled rule included), with the node fold
  unchanged. A node whose |x| is not finite, and every node if no kernel is selected, takes the scalar path.

## Identity (bit for bit)

- **`unit_erfcx_batch`:** every kernel against the scalar erfcx on 800,401 arguments (a fine sweep of [0, 40],
  each piece boundary with 64 neighbouring doubles each side, the signed zeros, subnormals, huge values, 400,000
  log-uniform), at lengths 1-9 for the remainders; the integrand's batch against its `operator()` on 69,708 nodes
  (0 < d < n, d = 0, d = n, overflowing |x|); `log_contrib` and `log_contrib_derivs` on 400 periods through the
  parity rule and its doubled rule, every kernel against no kernel. 0 differences.
- **The committed goldens** under each kernel forced: `xref_core_regression_erfcx_*` and
  `unit_ar1_filter_erfcx_*` (the AR(1) filter's likelihood goldens and its φ = 0 identity). The full non-slow suite
  passes at the default kernel (AVX2 here) with GCC 13 and the new suites with Clang 18, `-Werror`.

## Measured (2026-10-12)

Machine: cloud VM, Intel Xeon @ 2.80 GHz (2 vCPU, AVX2 and AVX-512), GCC 13.3 `-O3`, contraction off, one
thread. `p7_erfcx --seconds 3`, three runs, kernels in alternating order; "main" is the same benchmark source built
on `main` (26d8860) with the same flags.

| Measure | main | scalar | SSE2 | AVX2 |
|---|---|---|---|---|
| `erfcx` alone, ns per argument (128 node-like arguments) | 48-52 | 51 | 21 | 14 |
| `log_contrib`, µs per period (81 recovery cells, 0 < d < n) | 21.6-22.8 | 18.9 | 14.3-15.1 | 14.5-15.0 |
| `log_contrib_derivs`, µs per period | 48.1-49.6 | 48.3-49.7 | 42.8-44.8 | 48.8-49.8 |

- **Against the expectation:** one evaluation (`log_contrib`) is **33% faster** at SSE2 and AVX2, at the top of the
  20-35% expected. The scalar kernel alone gives 15%: splitting the nodes' erfcx from their finishing lets the
  scalar code pipeline better. `erfcx` itself is 3.6× faster at AVX2, at the top of 2.5-3.5×.
- **Not as expected:** with the score and Hessian (`log_contrib_derivs`), SSE2 gains 11% and AVX2 nothing. That path
  is dominated by the derivative functor's `inverse_mills`, which P-7 does not touch, so its gain was always going
  to be smaller; that AVX2 gains less than SSE2 there on this VM is consistent with the AVX clock penalty on this
  Xeon, and is not explained further here. The default stays the fastest available kernel (AVX2); a machine where
  SSE2 is faster can force it with `VCAL_ERFCX_KERNEL=sse2`, with identical results.
- **For a run:** not measured on this 2-vCPU machine. The next full run's projection (D-187) measures it, on the
  machine that runs it, against this table.
