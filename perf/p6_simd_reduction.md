# P-6: the grid reduction in vector lanes (D-404)

## Expected, written before the build (2026-10-13)

- **Scope (D-199, with D-213's scope note):** `engine::weighted_sum`, the compensated sum over periods behind
  every grid reduction: the full-grid argmax of `reduce_weighted` (calibrate's grid search, the `Full` reference
  path), the bounded argmax's visited tiles and its tile bounds (`resample::detail::argmax_bounded`, the bootstrap
  default since D-172), and the unit-weight surface sums in `posterior.hpp` and `profile.hpp`. Lanes over grid
  points: the recovery grid is 61 (PD) × 41 (ρ) with ρ contiguous, so a 4 × 4 tile row is four adjacent k, one
  AVX2 vector or two SSE2 vectors.
- **Where the time is:** the scalar `weighted_sum` costs 2.4 ns per term at T = 13 and 40, and 3.4 ns at T = 100.
  That is a full-grid sum of 77, 248 and 842 µs (microbenchmark on random L, one thread). About 7 cycles per term:
  latency-bound, because each Neumaier step waits on the last. In the runs:
  - **S-10 (parametric bootstrap):** the grid reduction was 5-7% of the run (D-199's profile).
  - **W × L bootstrap tools:** in `surface_cache_spike --profile` (subset scenarios 37, 4, 49), the ArgMax
    reduction is 0.2 s of the 0.4 s of per-replicate work.
- **The change:** a lane per grid point, each doing `weighted_sum`'s operations in its order: the same products
  w[t]·L[t][k], the same Neumaier branches (both computed, the one the scalar takes selected exactly), and the
  same non-finite handling (a lane's non-finite terms go to its own `nonfinite` sum and skip compensation). Several
  vectors are interleaved to hide the add latency. Scalar, SSE2 and AVX2 kernels, `VCAL_REDUCE_KERNEL`, in
  `vcal_core_kernels`, as P-7. The callers compute a run of adjacent k at once and push them to their reducer in
  the same ascending order.
- **Expected:** the reduction **3-5× faster** in the microbenchmark (four lanes, with the latency hidden). At run
  level: **S-10 3-5% faster**; the W × L bootstrap's per-replicate work **20-35% faster**; calibrate's grid search
  faster by its reduction share. Results are bit for bit.
- **Stop rule:** under 1.5× on the reduction microbenchmark, or under 2% on S-10's tool and under 10% on the
  W × L tool's bootstrap phase: recorded as measured, not merged.
- **Bit for bit:** each lane's arithmetic is the scalar's on the same operands, with no contraction and no
  reassociation. The checks: every kernel against `weighted_sum` on random, adversarial (huge, tiny, mixed-sign,
  cancelling) and non-finite surfaces at every length and offset; `replicate_estimates` (Full and Bounded) and
  calibrate bit for bit under every kernel; and the committed goldens under each forced kernel.

## What was built

- **`core/reducers/weighted_sum.hpp`:** `weighted_sum`, moved unchanged from `engine/surface.hpp` (`engine::` keeps
  the name through a using-declaration), and `weighted_sums(L, T, K, w, k0, n, out)`, which outside device code
  runs the kernels and otherwise the scalar loop.
- **`core/reducers/sum_batch.{hpp,cpp}`, `sum_batch_avx2.cpp`, `sum_lanes.hpp`:** scalar, SSE2 and AVX2 kernels
  for `weighted_sums` and the bounded argmax's tile bounds, selected by `VCAL_REDUCE_KERNEL`, in
  `vcal_core_kernels`. Four vectors are interleaved.
- **A lean finite path:** the lanes carry only the sum and the compensation, plus one mask of non-finite terms.
  A block whose lanes met a non-finite term is redone with the scalar code. The first version carried the
  scalar's non-finite sum in every lane as well; it reached only 1.5-2.3× (AVX2), because four values per vector
  spilled registers.
- **The callers:** `reduce_weighted` (chunks of 64 adjacent k, pushed in ascending k), the bounded argmax (each
  tile row of up to four adjacent k, and the tile bounds), `grid_posterior` and `summed_surface`.

## Identity (bit for bit)

- **`unit_sum_batch`:**
  - every kernel against `weighted_sum` on 56,754 sums: random, adversarial (huge, tiny, subnormal, signed zeros,
    cancelling) and non-finite surfaces, with zero, integer and fractional weights, at run lengths 1-40 and
    offsets that leave remainders;
  - 522 tile bounds against the scalar loop;
  - calibrate, and 64 bootstrap replicates with both searches, under every kernel against none.
  All 0 differences. `unit_reduce_kernel_*` covers each forced kernel and an unavailable one (exit 77).
- **The committed regression estimates** under each forced kernel (`xref_core_regression_reduce_*`), the full
  non-slow suite (79 of 79, GCC 13, `-Werror`), and the slow full-against-bounded argmax check
  (`study_argmax_check`).

## Measured (2026-10-13)

Cloud Xeon 2.8 GHz, GCC 13, one thread, Release.

| The reduction, recovery grid (2,501 points), ns per term | scalar | SSE2 | AVX2 |
|---|---|---|---|
| T = 13, full grid | 1.68 | 0.91 | 0.50 (3.4×) |
| T = 40, full grid | 1.77 | 0.95 | 0.48 (3.7×) |
| T = 100, full grid (2 MB surface, out of cache) | 2.25 | 1.65 | 1.02 (2.2×) |
| runs of 4 (tile rows), T = 13 / 40 / 100 | 1.81 / 1.79 / 3.13 | | 0.76 / 0.63 / 1.71 |

- **The W × L bootstrap tool** (`surface_cache_spike --profile`, subset scenarios 37, 4 and 49, 20 replicates):
  the ArgMax reduction 0.9 → 0.6 s (−33%), the cached arm's total 4.8 → 4.4 s.
- **S-10's tool** (`study_parametric_bootstrap`, scenarios 37 and 4, one replicate, A B B A): 51.8, 53.0 s on
  `main` against 52.6, 52.4 s. No measurable change.
- **Against the expectation:** the reduction met it (3.3-3.7× with AVX2 at T ≤ 40, 2.2× at T = 100). The W × L
  bootstrap phase met it (−33%). S-10's tool missed (about 0% against 3-5%). D-199's 5-7% share of S-10 was
  probably mostly reduction work the re-profile counted elsewhere: calibrate's one full-grid argmax per refit is
  50-250 µs, against milliseconds of integrals per refit.
- **Against the stop rule:** the W × L bootstrap phase is over 10%, so it is not triggered. Adopted.
- **Who gains:** the W × L bootstrap studies (replicate estimates through the bounded argmax), the posterior's grid
  sums (the Bayesian studies) and the profile's summed surface. The parametric bootstrap and the AR(1) MLE barely
  do.
