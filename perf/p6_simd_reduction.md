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
