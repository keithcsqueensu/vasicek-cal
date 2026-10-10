# P-18: the AR(1) filter's transition sums on kernels (D-199)

The AR(1) filter (`benchmark_models/ar1/filter.hpp`, M6) spends most of an evaluation and its score in log-sum-exps
over node pairs: about 128 × 128 terms per pair of adjacent periods, five passes' worth per evaluation plus score.
P-18 moves those sums onto kernels (`benchmark_models/ar1/transition*.{hpp,cpp}`) that run one lane per output node,
with the in-house exp shared with the DGP (`common/detmath/det_exp.hpp`).

## Why the exp changed

A bit-identical vector kernel cannot vectorise libm's exp. Measured on the filter, `std::exp` was 59-63% of an
evaluation plus score on every (T, n) cell, so a kernel that kept it would have been capped well below 2×. The
filter had no pinned results yet, so the owner chose to move its transition to `det_exp` (D-199): 2 ulp from mpmath,
built from correctly rounded operations, and identical on every platform.

| (n, T) | evaluation + score | transition terms | that many `std::exp` |
|---|---|---|---|
| 100, 20 | 15.6 ms | 3.36M | 9.1 ms (59%) |
| 100, 100 | 83.2 ms | 18.4M | 51.6 ms (62%) |
| 1,000, 100 | 43.0 ms | 9.18M | 26.3 ms (61%) |
| 10,000, 100 | 36.8 ms | 8.11M | 23.3 ms (63%) |

## Speed

One evaluation (three passes) plus one score, best of three, φ = 0.6, on one DGP panel per cell (Core Ultra 9
275HX, MSVC 19.51). "libm filter" is the filter before P-18 on today's core (the P-7 extension included in both
columns), so the comparison isolates the transition:

| (n, T) | libm filter | scalar | SSE2 | AVX2 |
|---|---|---|---|---|
| 100, 20 | 14.2 ms | 26.9 | 19.1 | 9.5 |
| 100, 100 | 77.7 | 147.2 | 104.1 | 51.7 |
| 1,000, 20 | 7.4 | 13.8 | 9.8 | 5.0 |
| 1,000, 100 | 39.0 | 73.4 | 52.4 | 26.4 |
| 10,000, 20 | 6.5 | 12.4 | 8.9 | 4.4 |
| 10,000, 100 | 33.6 | 64.6 | 46.1 | 23.0 |

AVX2 is 1.46-1.50× the libm filter. The scalar fallback is about 1.9× slower and SSE2 about 1.35× slower, accepted by
the owner: correctness and identical results across paths come first, and every x86-64 machine has at least SSE2.

**The forward kernel alone,** a 128 × 128 block: libm's online LSE 3.38 ns per term; the kernels 7.25 (scalar),
5.20 (SSE2) and 2.18 (AVX2). `det_exp`'s 17-term Horner chain defines its bits and is about half of the AVX2
kernel's operations; it is the cap.

**How the AVX2 kernel got there:**

| step | ns per term |
|---|---|
| per-term special cases, two-product scaling, the all-lanes skip | 2.99 |
| interleaving four blocks | 2.99: throughput-bound, not latency-bound |
| exact fast paths: inputs classified once per source or destination; one exact product for 2^k when every lane is normal | 2.18 |

The score's pair terms (two divisions per pair) were then vectorised with exact divisions and summed in order, taking
the n = 100, T = 20 evaluation plus score from 10.2 to 9.5 ms.

**Evaluations and scores (owner, in place of changing the Hessian refresh):** a fit makes about 86 evaluations
and 71 scores (12 fits, range 71-103 and 63-85). At AVX2 cost an evaluation is 5.2 ms and a score 4.3 ms at
n = 100, T = 20: evaluations are about 59% of a fit's filter time, scores 41%.

**Skipped exps (owner's item 1):** the share of transition exps whose argument is below -750, where exp is exactly
+0, is 0.6% at φ = 0, 0.4% at 0.6 and 0.1% at 0.9, over about 327M terms at each φ (six (T, n) cells, four panels
each, an evaluation and a score at two points). The skip is exact and kept in the scalar kernel; the vector kernels
compute those terms, because a skip pays only when every lane qualifies.

## The S-35 M6 arm

`perf/project_runtime.py --run` on `study_lrpd_serial_m6 --boot-sample 8` (the comparator sampled at 8 refits per
job and scaled to B = 199; no rows written), 24 threads, efficiency 0.997:

| | projection | 90% range | with the ±5% allowance | accepted at |
|---|---|---|---|---|
| before P-18 (d68c364) | 46.1 h | 42.4-50.0 h | 40.3-52.5 h | 4 replicates per scenario |
| after P-18 | 34.9 h | 32.9-37.0 h | 31.2-38.9 h | 8 replicates per scenario |

24% shorter, and still about twice the registration's estimate of 18 h. The arm's scale is the owner's decision,
before the validation run (D-197).

**Where the 34.9 h goes** (the 4-replicate sample kept, `phase_split`): the PD and ρ profile intervals on every fit
50.9% (17.8 h), the point estimate 26.0% (9.1 h), the bootstrap comparator 20.5% (7.2 h; 900 replicates × 199
refits), φ's profile on the subset 2.6% (0.9 h); the binomial start, the simulation and the fills under 0.1%. By
cell, T = 100 is 62% of the run (n = 100: 27.5%, n = 1,000: 22.9%, n = 10,000: 11.4%).

## Identity

- Every kernel against a plain reference, bit for bit: the forward and backward LSEs, the shifted exps and the score
  terms at sizes that leave remainders, with -inf, NaN and +inf terms; the vector `det_exp` on 462,050 arguments and
  its non-positive path on 381,000 more (`unit_ar1_transition`; each kernel forced: `unit_ar1_kernel_*`).
- The filter under every kernel: likelihood, placement check, per-period terms and score, 360 numbers
  (`ar1_filter_is_the_same_under_every_kernel`); the filter suite with each kernel forced
  (`unit_ar1_filter_kernel_*`); 12 full fits and 4 profile sets, every field, the same under all three kernels.
- φ = 0 is still the per-period likelihood bit for bit (48 points). Against the libm filter: goldens within 7.2e-15
  (10 of 66 rows moved), ref 1.1e-13, smoothness 4.4e-13, the score 5.6e-10, scipy 2.9e-14 / 7.5e-15 / 1.8e-14.
- The owner's items 1-6, against the `det_exp` filter: the snapshot reproduced bit for bit under every kernel.
- `det_exp` shared: every DGP golden, both reference-panel hashes and the mirror test unchanged.
