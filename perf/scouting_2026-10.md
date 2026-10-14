# Performance scouting, October 2026: the re-profile on `main` and six candidates (D-500)

**Proposals only.** Nothing here is built, benchmarked as a candidate, or merged. Each candidate still needs its own
decision entry, benchmark and identity or tolerance test under the [P-item rules](README.md#rules-for-every-p-item-d-177)
before it merges, and the owner assigns its P-number and its place in the order. The numbers below are shares and
estimates from a cloud machine with GCC; anything adopted is measured on the owner's laptop, interleaved A B C C B A.

Contents: [Summary](#summary) · [Method](#method) · [The re-profile](#the-re-profile) ·
[Candidates](#candidates) · [Checked, not proposed](#checked-not-proposed) · [Reproducing the re-profile](#reproducing-the-re-profile)

## Summary

The S-10 re-profile of D-199 (W* 39.8%, argmax + refine + SE 27.7%, analytic SE 13.5%, simulate 12.5%) is out of date.
On current `main` (c70f9b4, with P-6, P-7, P-15, P-17 and P-21), the S-10/S-4 run's projected thread-seconds split as
**W* (S-4b) 31.6%, argmax + refine + SE 26.7%, simulate 21.2%, analytic SE 11.5%, surface 4.5%, Jeffreys 3.2%,** the rest
under 1%. The T = 100, n = 10,000 cell is 36% of the run. Inside the integration phases the derivative evaluation
(`log_contrib_derivs`, one per distinct observation per evaluation point) is the unit of cost: a bootstrap fit takes one
value pass, one Hessian pass (calibrate), one more Hessian pass (the analytic SE) and, on 199 of 999 panels, about 13
Hessian passes for W*. libm's `exp`, `log` and `log1p` are a quarter of the instructions in that cell (the in-house
`log`/`exp` decision of D-199, not a new item).

The S-35 M6 arm's tool splits as **the PD and ρ profile intervals ~60-66% of a job, the point estimate ~30-34%,** the
comparator and φ's profile on the subset only (D-199 put them at 20.5% and 2.6% of the arm). The AR(1) transition
kernels (P-18) are about 72% of the instructions at T = 100, n = 100, the costliest cell (27% of the arm), where a
one-sided period's 256 composite nodes make each transition pair 2-4× a 128 × 128 one. In the n = 10,000 cells, **about
half of all filter evaluations are rejected line-search points:** after BFGS updates the Hessian goes bad, the full step
falls downhill, and the search halves it five to seven times.

| Rank | Candidate | Acts on | Changes results | Expected run-level gain |
|---|---|---|---|---|
| 1 | [C1: calibrate's value pass taken from its Hessian pass](#c1-calibrates-value-pass-taken-from-its-hessian-pass) | S-10 argmax + refine + SE (26.7%); every tool that calibrates | no | S-10 about 8% |
| 2 | [C2: bootstrap fits read cached rows in place, argmax bounded by cached tile maxima](#c2-bootstrap-fits-read-cached-rows-in-place-argmax-bounded-by-cached-tile-maxima) | S-10 surface (4.5%) and the K × T mechanics inside calibrate | no | S-10 about 6-9% |
| 3 | [C3: W* without repeated evaluations, and P-20's concrete form](#c3-w-without-repeated-evaluations-and-p-20s-concrete-form) | S-10 W* (31.6%) | tier A no; tier B within solver tolerance | S-10 2% + 2% (after P-19) + about 5% |
| 4 | [C4: the M6 search's line search](#c4-the-m6-searchs-line-search) | M6 fit and profiles (~95% of the arm) | within solver tolerance | M6 arm about 7% (range 4-11%) |
| 5 | [C5: one `erfcx` batch for the composite rule, and its nodes cached for the derivatives](#c5-one-erfcx-batch-for-the-composite-rule-and-its-nodes-cached-for-the-derivatives) | one-sided periods (d ∈ {0, n}): S-10's n = 100 cells, M6's binomial terms | no | S-10 about 1%; M6 under 1% |
| 6 | [C6: fewer composite panels in the AR(1) filter](#c6-fewer-composite-panels-in-the-ar1-filter) | M6 transitions into and out of one-sided periods (n = 100 cells, 43% of the arm) | within a stated tolerance (observed ≤ 3.4e-10 in l); breaks the φ = 0 bit-for-bit identity | M6 arm 10-15%; **not recommended as is** |

The interior rule's node count is **not** a candidate: 64 or 32 Gauss-Hermite nodes reproduce the parity filter to rounding
near the data, but the filter itself is unconverged in the box's corners at n = 10,000 (its own doubled rule disagrees by
up to 1.7 nats there), and a reduced rule disagrees by 1-100 nats; see C6.

## Method

- **Machine:** a cloud VM, Intel Xeon @ 2.1 GHz (2 vCPU, AVX2), Ubuntu 24.04, GCC 13.3, `cpu-release` preset (`-O3`,
  contraction off) with `-g -fno-omit-frame-pointer` added for the profiler. Every kernel selected at run time was AVX2
  (`dgp_kernel=avx2`; the erfcx, reduce and AR(1) kernels likewise). Timings here are noisy at about 3% and the owner's
  runs use MSVC; shares are the result, not seconds.
- **Phase shares:** `study_parametric_bootstrap --replicates 2 --scenarios 36,...,44 --threads 1 --timing-out` and
  `study_lrpd_serial_m6 --replicates 2 --scenarios 36,...,44 --threads 1 --boot-sample 8 --timing-out`: one scenario per
  (T, n) cell (PD 1%, ρ 0.12; ids 36-44 cover the nine cells), two replicates each, one thread, 724 s and 686 s of wall
  time. Each cell's cost was scaled by its nine scenarios and the phase schedule of D-187 (Jeffreys 0.05, the S-4 widths
  0.20, fills once; the M6 comparator 0.1 on the subset's scenarios only). `perf/project_runtime.py` on the same samples
  projects 5.84 h × 9 and 117 h × 9 on one thread, since it reads one scenario per cell as a nine-scenario run.
- **Inside the phases:** valgrind's callgrind (`--read-inline-info=yes`, so inlined header code is attributed to its own
  function), on a throwaway driver that repeats one replicate's per-bootstrap-panel work (simulate, surface from the
  cache, calibrate with the sampled check, the analytic SE, W* on 199 of 999) with the cache warm, at (100, 10,000) and
  (100, 100); and on the M6 tool itself at (100, 100), φ = 0.6, one replicate. Callgrind counts instructions, not time:
  an AVX2 instruction counts once for four lanes and a memory instruction counts once whatever it waits on, so the
  memory-bound shares below understate time and the kernel shares overstate it a little. `perf` is not available on
  this machine; gprof would have needed `-pg`, which changes the code.
- **Counts:** throwaway drivers counted panel evaluations in W* (polish against the two inner solves, and exact repeats
  of an evaluation point), the M6 search's evaluations, scores, Hessian refreshes and line-search rejections, and the
  AR(1) filter's value under reduced rules. The drivers are described in the last section and were not committed.

## The re-profile

### S-10/S-4: `study_parametric_bootstrap`, by phase and cell

Shares of the projected full run (81 × 1,000 replicates, B = 999), with the phase schedule applied:

| Phase | Share | D-199 (3876889) |
|---|---|---|
| W* (S-4b) | **31.6%** | 39.8% |
| argmax + refine + SE | **26.7%** | 27.7% |
| simulate | **21.2%** | 12.5% |
| analytic SE | 11.5% | 13.5% |
| surface (row lookups and copies) | 4.5% | — |
| Jeffreys (replicates 0-49) | 3.2% | — |
| profile, W0, S-4 widths, cache fill | 1.4% together | 6.5% (the rest) |

P-7 and P-21 took the integration phases down and P-17 did not move, so simulate's share rose; on this VM's AVX2 the
Philox kernel may also run relatively slower than on the owner's laptop, where D-199 measured 12.5%.

| (T, n) | share of run | simulate | surface | argmax + refine + SE | analytic SE | W* | Jeffreys |
|---|---|---|---|---|---|---|---|
| 20, 100 | 3.5% | 4% | 3% | 21% | 13% | 56% | 2% |
| 20, 1,000 | 5.3% | 6% | 3% | 31% | 14% | 40% | 4% |
| 20, 10,000 | 7.9% | 28% | 4% | 23% | 11% | 31% | 3% |
| 40, 100 | 4.3% | 6% | 5% | 29% | 13% | 44% | 2% |
| 40, 1,000 | 8.0% | 8% | 4% | 32% | 14% | 38% | 3% |
| 40, 10,000 | 14.9% | 29% | 3% | 23% | 11% | 29% | 3% |
| 100, 100 | 6.3% | 10% | 11% | 32% | 9% | 34% | 3% |
| 100, 1,000 | 13.5% | 11% | 6% | 32% | 13% | 32% | 4% |
| 100, 10,000 | **36.3%** | 31% | 4% | 24% | 11% | 26% | 3% |

### Inside the S-10 phases (callgrind, instructions)

**(100, 10,000), one replicate's bootstrap work, cache warm.** Of the instructions: simulate 31% (the AVX2 Philox count,
`count_blocks_avx2`, 30%); `panel_derivs` 51% in all, through three callers: calibrate's Hessian 11.4%, the analytic SE
11.4% and W* 28%; calibrate's value pass (`log_contrib` per distinct observation at the estimate, plus the doubled rule on
5% of fits) 9.3%; calibrate's full-grid argmax (`reduce_weighted`, P-6 lanes) 1.8%; the surface phase 6.3% (the row
copies 3.2% and the zero fill of `L.assign` 3.1%). Inside the derivative evaluations: the node batch (P-7's erfcx lanes
plus `log_phi_pair` with its libm `log`, `log1p` and `exp`) 35%, the derivative functor 10%, and the fold's `exp` per
node. libm in total: `exp` 15.9%, `log` 5.0%, `log1p` 4.3% of all instructions.

**(100, 100), the same.** The integrals are cheap (6.3 distinct observations per panel) and the K × T mechanics show:
the row copies 17.4% and the zero fill 17.3% of instructions (memory-bound, so less in time: the timing table's 11%),
the full-grid argmax 10%, `panel_derivs` 41% (calibrate 11%, the SE 6%, W* 24%), the value pass 6.7%, one-sided periods'
composite rule 13.4% (of which the nodes' scalar `g(z)` 9.3%, `sinh` 1.6%). 50% of this cell's periods are one-sided.

**W*: evaluations per bootstrap panel** (`profile_lr_statistic`; 40 panels per scenario, replicates 0-1):

| scenario (T, n, PD) | polish | the two inner solves at θ̂ | exact repeats of a point |
|---|---|---|---|
| 36 (20, 100, 1%) | 7.0 | 11.3 | 8.2% |
| 38 (20, 10,000, 1%) | 6.0 | 7.8 | — |
| 42 (100, 100, 1%) | 6.6 | 7.8 | 6.9% |
| 44 (100, 10,000, 1%) | 5.9 | 7.2 | 7.6% |
| 9 (20, 100, 0.1%) | 50.9 | 9.2 | — |
| 17 (100, 10,000, 0.1%) | 5.9 | 7.0 | — |

So W* is about 13 Hessian passes per panel, 45% of them in the polish, against calibrate's one value pass and one
Hessian pass; and about one evaluation per panel is an exact repeat of a point already evaluated (`start_axis1`
re-evaluates the polished maximum, whose derivatives the last inner solve computed). The PD 0.1%, n = 100 cell's polish
is an outlier (50 evaluations: the surface is nearly flat along ρ there) in a cell that is 1% of the run.

### S-35's M6 arm: `study_lrpd_serial_m6`, by phase and cell

Per job (a replicate at φ ∈ {0, 0.3, 0.6}), one thread; the comparator (φ = 0.6, replicates 0-99, 199 refits each) and
φ's profile run on the subset's scenarios only, of which the sample met two (37 and 43), so their shares are taken from
D-199 (20.5% and 2.6% of the arm) rather than from this sample.

| (T, n) | share of the arm (fit + PD/ρ profiles) | M6 fit | PD and ρ profiles | one-sided periods |
|---|---|---|---|---|
| 20, 100 | 5.0% | 33% | 67% | 42% |
| 20, 1,000 | 3.8% | 30% | 70% | 11% |
| 20, 10,000 | 3.3% | 30% | 70% | 0% |
| 40, 100 | 10.9% | 30% | 70% | 54% |
| 40, 1,000 | 5.8% | 35% | 65% | — |
| 40, 10,000 | 6.5% | 31% | 69% | 0% |
| 100, 100 | **26.6%** | 34% | 66% | 47% |
| 100, 1,000 | **24.3%** | 39% | 61% | 2% |
| 100, 10,000 | 13.9% | 34% | 66% | 0% |

The cell weights agree with D-199's (T = 100: n = 100 27.5%, n = 1,000 22.9%, n = 10,000 11.4%). The n = 100 cells are
the dearest per job (80 s at (100, 100) against 42 s at (100, 10,000)), because a one-sided period's composite rule has
256 nodes and the transition into and out of it has 256 × 128 or 256 × 256 terms.

**Inside (callgrind at (100, 100), φ = 0.6):** `log_likelihood` (three passes) 55%, `score` 44%; the transition kernels
72% (forward 30%, backward 30%, `exp_shifted` 12%), the binomial terms about 5% (the composite rule's scalar nodes 3%,
the Gauss-Hermite batch 1.4%), `weighted_moments` and the kernels' libm `fmax` 2.9%, libm `exp` 2.5%, vector resizes and
zero fills 2.5%, the central-difference Hessian (`Search::hessian`) 3.5% of its own.

**The search's counts** (evaluations and scores per fit and per PD + ρ profile, and the line search):

| scenario (T, n) | per fit: evaluations / scores | per profile pair | line-search evaluations rejected | iterations with a downhill full step |
|---|---|---|---|---|
| 38 (20, 10,000) | 100 / 70 | 314 / 144 | **540 of 1,045 (52%)** | 93 of 505 (77 with a BFGS-updated Hessian) |
| 41 (40, 10,000) | 90 / 68 | 271 / 126 | 300 of 604 (50%) | 57 of 304 |
| 44 (100, 10,000) | 86 / 71 | 292 / 127 | 163 of 307 (53%) | 33 of 144 |
| 40 (40, 1,000) | 87 / 69 | 183 / 122 | 122 of 411 (30%) | — |
| 43 (100, 1,000) | 76 / 63 | 135 / 96 | 34 of 157 (22%) | — |
| 39 (40, 100) | 59 / 50 | 154 / 111 | 67 of 308 (22%) | — |
| 42 (100, 100) | 64 / 57 | 116 / 98 | 6 of 130 (5%) | — |

A rejected line-search point costs a full three-pass evaluation and no score. In the n = 10,000 cells the full Newton
step goes downhill in about one iteration in five, almost always with a BFGS-updated Hessian (77 of 93 at scenario 38,
against 4 with a freshly differenced one), and each such iteration then halves the step 5-7 times (27 of 505 steps were
accepted at α < 1/64). Weighted by cell cost, about a quarter of the arm's line-search evaluations are rejected points,
roughly 10% of the arm's time.

## Candidates

Ranked by expected run-level gain, with the cost of a results change weighed as the rules ask. Portability throughout:
GCC 11, MSVC 14.39 and CUDA 11.8; AVX2 chosen at run time with SSE2 and scalar fallbacks; device-visible code stays C++17.

### C1: calibrate's value pass taken from its Hessian pass

- **What:** `engine/calibrate.hpp` only. `calibrate_from_surface` integrates each distinct observation at the estimate
  twice: once for `l_primary` (the log-likelihood, the quadrature check's primary side, and the check's rounding scale)
  through `log_contrib`, and once more inside `panel_derivs` for the analytic Hessian (P-15), whose per-observation value
  `log_contrib_derivs` returns is `log_contrib`'s bit for bit (the objective's contract, D-170: the same nodes and fold;
  checked again here on 28,800 (observation, θ, rule) cases under the scalar, SSE2 and AVX2 erfcx kernels, 0
  differences). For a fit that reaches the Hessian (no grid-edge or flat flag) the Hessian pass runs first and `l_primary`
  is taken from its per-observation values, in period order through `first[]` as now; the doubled rule's `l_check` pass
  is unchanged; flagged fits and objectives without derivatives keep the value pass. `panel_derivs_natural` would return
  (or write) the per-observation values, a signature change inside `engine/derivs.hpp`.
- **Where it acts:** the argmax + refine + SE phase, 26.7% of S-10. At (100, 10,000) the value pass is 9.3% of the
  instructions against the Hessian's 11.4% and the phase's 22.6%, so about 40% of the phase; at (100, 100) 26% of it.
  The doubled rule on 5% of fits is about a tenth of the pass and stays. Arithmetic: 0.267 × (0.3 to 0.4) × 0.9 ≈ **7-10%
  of the run, 8% expected.** The recovery harness, the M6 tool's binomial starts and comparator refits, and every
  `calibrate` call take the same saving in their calibrate phase.
- **Results:** no (bit for bit): every number is the one it replaces, from the same evaluation point and the same
  summation order; the check still compares `l_primary` with `l_check` term by term.
- **Stop rule:** under 4% on `study_parametric_bootstrap` at scenarios 37 and 4 (one thread, A B B A), recorded and not
  merged. (P-6's 3-5% expectation missed on this tool because its reduction share was small; this item's share is
  measured above, so a miss would mean the measurement was wrong.)
- **Identity test:** `calibrate_takes_its_value_pass_from_the_hessian_pass_bit_for_bit`: on the recovery subset's panels
  and 64 bootstrap panels per cell, every `Estimate2` field against a verbatim copy of the old path, including the check's
  max, total and flagged count and the quadrature-check-off path (P-16); the committed goldens and the S-10 tool's rows
  (scenarios 37 and 4, as P-14) byte for byte.
- **Risks:** none for portability or parity (host-only engine code; the device path does not run calibrate). The order of
  `calibrate_from_surface`'s steps changes (Hessian before the check); `quad_check_flagged` must still be set from the
  same comparison.
- **Prior art:** P-14 (each distinct observation integrated once per evaluation point) and P-19 (one computation shared
  by two users at the same point). P-14 removed repeats across periods within a pass; this removes the repeat between the
  two passes. P-19 (in flight) shares calibrate's Hessian with the studentised SE; with both, an unflagged bootstrap fit
  integrates each distinct observation once with derivatives instead of three times (value, Hessian, SE), about
  3.5 → 1.25 in value-pass units on the integration part of calibrate plus the SE phase.

### C2: bootstrap fits read cached rows in place, argmax bounded by cached tile maxima

- **What:** a surface *view* for fits whose rows come from `SurfaceRowCache` (P-1): T row pointers instead of a
  materialised T × K `L`. `evaluate_surface_cached` would return the pointers; `calibrate_from_surface` would take either
  a contiguous `L` or a view through one accessor; `reducers::weighted_sums` (P-6) gains a row-pointer entry whose lanes
  do the same operations in the same order (`core/reducers/weighted_sum.hpp`, `sum_batch*.cpp`, `sum_lanes.hpp`);
  `refine_2d`'s nine stencil sums and `summed_surface` (W*'s and the profile's brackets) read through the view. With the
  copy gone, calibrate's full-grid argmax (`reduce_weighted` over K × T) becomes the bounded argmax of D-172
  (`resample::detail::argmax_bounded`, unit weights, one replicate), with each row's tile maxima `U` (K / 16 doubles) and
  its finiteness flag computed once at fill time and kept in the cache entry (`engine/surface_cache.hpp`), so the per-fit
  cost is the tile bounds (T × 176 terms on the recovery grid) plus the visited tiles, not K × T. Files: `engine/
  surface_cache.hpp`, `engine/calibrate.hpp`, `engine/refine.hpp`, `engine/profile.hpp`, `engine/surface.hpp`,
  `core/reducers/*`, `resample/bootstrap.hpp`, and the S-10 tool's `fit_panel`.
- **Where it acts:** the surface phase (4.5% of S-10: a 2 MB zero fill plus a 2 MB copy per bootstrap fit at T = 100) and
  the argmax inside calibrate (1.8% of instructions at (100, 10,000), 10% at (100, 100)). By cell: about 20% of the
  (T, 100) cells (14% of the run), about 9% of the (T, 1,000) cells (27%), about 5% of the (T, 10,000) cells (59%):
  2.8 + 2.4 + 3.0 ≈ **8%, with 6-9% the range** (the memory-bound shares are instruction counts, and could be more or less
  in time).
- **Results:** no (bit for bit). The bounded argmax returns the full grid's argmax exactly (index, value and tie-break,
  with the NaN guard falling back to the full grid; D-172, `study_argmax_check`); the stencil and summed sums read the same
  values in the same order.
- **Stop rule:** under 3% on the S-10 tool at scenarios 51 (T = 100, n = 100) and 44 (T = 100, n = 10,000), one thread,
  A B B A: recorded and not merged.
- **Identity test:** `calibrate_on_a_row_view_is_calibrate_bit_for_bit` on the recovery subset and bootstrap panels
  (every field of `Estimate2`, the profile intervals and W* through the view against the contiguous path);
  `study_argmax_check` extended to unit weights through the view; the committed goldens; the S-10 tool's rows byte for
  byte on scenarios 51 and 29 (P-14's pair).
- **Risks:** the largest diff of the six (an accessor through the engine's hot loops); the row-pointer kernel must stay
  bit for bit with the contiguous one on every kernel path, tested as P-6's were. The cache entry grows by K / 16 doubles
  per row (1.3 KB on the recovery grid; the row is 20 KB). The GPU backend (M4) keeps the contiguous layout, which the
  view does not touch. The per-lookup mutex and key string in `SurfaceRowCache::row` are unchanged (measured small).
- **Prior art:** P-1 (rows cached), P-5 (bounded argmax, so far only in resampling), P-6 (lanes over adjacent k). It
  differs from P-23 (rejected: the row fill is once per scenario) in acting on the per-fit copy and reduction, which run
  999 times per replicate.

### C3: W* without repeated evaluations, and P-20's concrete form

W* (`profile_lr_statistic`, 31.6% of S-10; also the engine's S-4b option in `resample/bartlett.hpp`) is about 13
Hessian passes per panel: a nested polish (Newton on the PD profile, each value an inner Newton over ρ, about 6
evaluations), then the inner solve at θ̂ for each parameter (about 7 together). Three tiers, the first two bit for bit.

- **Tier A, memoised evaluation points (`engine/profile.hpp`, `NewtonProfile`):** keep the full `PanelDerivs` of the last
  few evaluated points and return the memo when `derivs(u0, u1)` is asked for exactly the same pair. The one repeat per
  panel measured above is `start_axis1`'s re-evaluation at the polished maximum, which the last inner solve already
  computed (7-8% of W*'s evaluations). `profile_intervals_newton` makes the same call and gains the same. Bit for bit:
  the memo returns the numbers the recomputation would produce (a pure function of the point, the panel and the rule).
  **Expected: about 2% of S-10.** Stop rule: under 1.5% on the S-10 tool's W* phase, not merged (it is a few lines, so the
  bar is low). Test: every W* and profile interval on the recovery subset and 64 bootstrap panels per cell, bit for bit
  against the memo disabled; `evaluations` falls by the repeat count. Risk: none for portability; the ulp round trip of
  D-199 means a point reached through `to_scaled(value)` is not the same pair as `r.scaled`, so the memo only catches
  exact repeats, as intended.
- **Tier A′, with P-19:** once `Estimate2` carries the calibrate's scaled coordinates (P-19's pass-through), the polish's
  first evaluation, `inner(0, û₀)` started at `w₀ = û₁`, is calibrate's Hessian evaluation at the same pair. Passing that
  `PanelDerivs` into `NewtonProfile` as a seeded memo entry saves one more pass per panel, bit for bit. **About 2% of
  S-10**, conditional on P-19's form; not to be built against P-19's branch (in flight).
- **Tier B, P-20 made concrete (within solver tolerance):** replace the nested polish in `profile_lr_statistic` with a
  two-dimensional safeguarded Newton on (u₀, u₁) from the estimate, using the analytic Hessian `panel_derivs` already
  returns (the refined vertex is within O(h²) ≈ 1e-3 of the maximum in u, so two steps reach `kProfileMaxTol`, three
  evaluations against about six), and seed the inner solves' `last[a]` from the maximum's Hessian as `start_axis1` does.
  l_max moves by at most the tolerance's effect (≤ tol²|H|, far below the residuals the tests bound), so W*, k_b and
  the S-4b intervals solved from k_b change within solver tolerance: a golden refresh for `bartlett_test` and a flip
  check of every recovery fit's S-4b interval (hours on the owner's laptop) are the cost. `profile_intervals_newton`
  keeps its polish, so every profile interval and golden is untouched. **Expected: polish 6 → 3 evaluations, about 20%
  of W*, about 5% of S-10** (less at PD 0.1%, n = 100, where the polish's 50 evaluations are the flat ρ direction).
  Stop rule: under 10% on the W* phase, not merged (the flip check does not pay for less). Test: the W* values against
  the nested polish within `TOL_PROFILE_SOLVER_AGREEMENT_U`-style tolerances in u on the recovery subset, the flip counts
  reported even when zero.
- **Prior art:** P-2 (profile solves by Newton, within solver tolerance, 0 flips) and P-19 (one evaluation shared). The
  owner's P-20 note names a warm start; within a replicate the 199 panels share only θ̂, and the inner solves already
  start from the polished maximum along the envelope slope, so the warm start worth having is the one above: the
  maximum from the Hessian the fit already has, not from a neighbouring panel.
- **Risks:** tier B changes a `native` engine option's numbers (S-4b, ABI 0.4) and `resample/bartlett.hpp` must take the
  same path; `kProfileMaxTol` is in u, so the test is in u.

### C4: the M6 search's line search

- **What:** `benchmark_models/ar1/search.hpp`, `Search::newton`. When the full Newton step is rejected, the backtracking
  halves α up to 40 times; measured above, a downhill full step is followed by 5-7 halvings, and it happens in about one
  iteration in five at n = 10,000, almost always after BFGS updates. Two changes, each within solver tolerance: (i) on a
  rejected step, the next trial α comes from the quadratic through l(u), the slope and l(u + s) (then the cubic, as a
  textbook backtracking line search), clamped to [0.1α, 0.5α], so one or two trials replace five to seven; (ii) a
  downhill full step (ln < l, not merely short of the Armijo bound) counts as two backtracks, so the next iteration
  refreshes the differenced Hessian rather than stepping again from the BFGS matrix that produced it. Neither touches the
  Hessian refresh policy the owner kept in D-199 beyond (ii)'s trigger, and (ii) can be left out if the owner prefers.
- **Where it acts:** the M6 fit and the PD/ρ profiles, 95% of the arm (the comparator's refits are the same search). The
  rejected evaluations are 50-53% of line-search evaluations at n = 10,000, 22-30% at n = 1,000, 5-22% at n = 100; weighted
  by cell cost about a quarter of the arm's line-search evaluations, and evaluations are 55% of the time: about 10% of the
  arm. (i) alone would keep roughly one trial per rejection: **about 7% of the arm, 4-11% the range;** (ii) is for the
  measurement to decide (a refresh costs 2m evaluations and scores, so it pays only where the rejections cluster).
- **Results:** within solver tolerance. The maximiser is still found to `kStepTolerance` (1e-10 in u) and the profile ends
  to `kProfileRootTol`, but the iterates differ, so every M6 estimate, SE, profile end and the comparator move at the
  tolerance's level: a golden refresh of the AR(1) fit goldens and a flip check of the M6 validation run's W1-W9
  (`studies/m6-validation`, 3.5 h on the owner's laptop; the S-35 M6 arm's registration is not amended, D-197). That is
  the cost weighed against 4-11% of a 35 h arm.
- **Stop rule:** under 5% on `study_lrpd_serial_m6` at scenarios 44 and 43, φ = 0.6, one thread, A B B A: recorded and
  not merged.
- **Tolerance test:** the twelve fits and four profile sets of the AR(1) snapshot against the current search within
  `kStepTolerance` and `kProfileRootTol` in u (a `TOL_AR1_SEARCH_AGREEMENT_U` entry), the log-likelihoods within
  rounding, and the evaluation counts reported; the flip counts in the decision entry even when zero.
- **Risks:** none for portability (scalar host code); a line search that is too greedy could accept a step the halving
  would have rejected, so the Armijo test stays and the clamp keeps α inside the halving's range. The comparator's refits
  (the same search) move too.
- **Prior art:** P-2 (a better solver for the same equations, within solver tolerance). The owner's D-199 note on the
  refresh policy covered the Hessian, not the line search, and no measurement of rejected evaluations existed; the
  counts above are the new argument.

### C5: one `erfcx` batch for the composite rule, and its nodes cached for the derivatives

- **What:** `core/quadrature/composite_legendre.hpp`: the 256 node abscissae `z_c + w sinh(u)` are known before any
  `g(z)` is evaluated, so the rule can call the integrand's `batch` (P-7's erfcx lanes) for all nodes, then fold them in
  the same order, as `GaussHermiteAdaptive` does; and the batch fills the `NodeCache` (P-21), so the derivative functor on
  a one-sided period finishes λ(∓x) from the value's erfcx instead of calling `inverse_mills` per node. The extent
  bisection's 113-116 scalar `g` calls stay as they are (each depends on the last). P-7's "not done" item.
- **Where it acts:** one-sided periods (d ∈ {0, n}): 33-92% of periods in S-10's n = 100 cells (14% of the run), 2-12% at
  n = 1,000, under 1% at n = 10,000; 42-54% of M6's periods at n = 100 (43% of the arm). At S-10's (100, 100) the composite
  rule is 13.4% of instructions, its nodes' scalar `g` 9.3%; P-7 took a third off a batched node, and the derivative path's
  per-node `inverse_mills` is about as much again, so about 30% of the composite rule's cost: 4% of that cell, **about 1%
  of S-10** (2-3% of the n = 100 cells). In M6 the binomial terms are about 5% of a job, so **under 1% of the arm.**
- **Results:** no (bit for bit): `batch` is `operator()` node by node (P-7's identity), and the fold is unchanged.
- **Stop rule:** under 10% on a one-sided period's `log_contrib` and `log_contrib_derivs` in `p7_erfcx` (extended to
  d ∈ {0, n}), not merged.
- **Identity test:** `p7_erfcx`'s identity extended to one-sided periods under every kernel; the composite goldens of the
  convergence table; the S-10 rows at scenarios 36 and 9 byte for byte.
- **Risks:** the rule is device-visible C++17: the batch path stays behind `#if !defined(__CUDA_ARCH__)` as in
  `gauss_hermite.hpp`; the derivative functor's cursor assumes node order, which the composite fold keeps.
- **Prior art:** P-7 and P-21, extended to the rule they left out; D-213's measurement of one-sided periods.

### C6: fewer composite panels in the AR(1) filter

- **What:** the AR(1) filter (`benchmark_models/ar1/filter.hpp`) integrates a one-sided period with the parity split
  rule's composite rule, 16 panels × 16 points = 256 nodes, so a transition into or out of such a period has 256 × 128 or
  256 × 256 terms against 128 × 128. With 8 panels (128 nodes) the transitions are uniform. The interior rule would stay
  at 128 (below).
- **Where it acts:** M6's n = 100 cells (43% of the arm), where 42-54% of periods are one-sided: the expected pair size
  falls from about 188² to 128² (×0.46), the transition kernels are 72% of the cell, so the cell is about 0.72 × 0.46 +
  0.28 = 0.61 of its time: **about 15% of the arm, 10-15% the range.**
- **Results:** within a stated tolerance. Measured here with the filter's value under 8 panels against 16 on 4 panels ×
  4 φ × 9 (PD, ρ) points per scenario, box corners included: ≤ 2.3e-13 at n = 100 (where it matters), ≤ 3.4e-10 at the box
  corners with n = 10,000 (the same level as the doubled rule's own disagreement at n = 100-1,000). **But it breaks the
  φ = 0 bit-for-bit identity with the binomial likelihood** (D-195, `ar1_phi_zero_is_the_per_period_likelihood_bit_for_bit`, whose panels include one-sided periods, and the strongest
  check of the M6 validation, W2-W3 in D-200), which the parity rule's 16 panels make exact; the identity would become
  a within-tolerance one. A reduced rule during the search only, with the estimate's final evaluation at the parity
  rule, keeps the identity at the estimate but not along the search.
- **Not proposed, interior nodes:** 64 and 32 Gauss-Hermite nodes reproduce the parity filter within 1e-14 to 1e-12
  near the data on every cell measured (the doubled rule's own level), but at the box's corners (ρ = 0.5 with PD = 1e-4 or
  0.2, n = 10,000, data at φ = 0.9) the filter itself is unconverged, its doubled rule 0.17-1.7 nats off, and 64 nodes
  are 1-140 nats off. The search and the profiles evaluate such points only when they wander, but a rule that is wrong
  there cannot carry a stated tolerance over the box.
- **Stop rule:** under 10% on the (100, 100) cell, not merged; and **not recommended** unless the owner accepts the loss of
  the φ = 0 identity for the n = 100 gain.
- **Tolerance test:** a `TOL_AR1_COMPOSITE_PANELS_LL` entry (1e-9 proposed) checked by the filter's doubled-rule check at
  every estimate, plus the AR(1) goldens refreshed and the M6 validation's flip check.
- **Risks:** the one above; nothing for portability.
- **Prior art:** none among the P-items (all kept the parity rule); D-118's choice of 16 panels was for the surface's own
  rounding floor up to ρ = 0.99, which the filter in the recovery box does not need.

## Checked, not proposed

- **The in-house `log`/`exp` (D-199, a decision of its own):** at (100, 10,000) libm's `exp` is 15.9%, `log` 5.0% and `log1p`
  4.3% of S-10's per-panel instructions, nearly all in the node batch (`log_phi_pair`'s log, log1p and exp pair) and the
  fold's exp per node; in M6 libm `exp` is 2.5% (`weighted_moments`). This is the number that decision was waiting for;
  it is not re-proposed as an item.
- **Simulate (21% of S-10):** P-17's AVX2 Philox is already the whole phase (30% of instructions at (100, 10,000));
  the only lever left is wider lanes (AVX-512), which the owner's CPU lacks, or a different DGP, which D-109/D-149 fix.
- **The M6 comparator's refits** (20.5% of the arm): a refit is a full scan over six φ values plus a polish; a warm start
  from the parent's φ̂ would skip the scan but could land on another local maximum, a results change beyond tolerance.
- **M6's `fmax` and scratch resizes:** libm's `fmax` is 2.9% of M6's instructions (`weighted_moments` and the kernels) and
  vector resizes with zero fills 2.5%. An exact inline maximum (NaN and signed-zero handling as `fmax`'s) and sized-once
  scratch would be bit for bit and worth 2-4% of the arm; a tidy-up for an M6 PR rather than an item.
- **The M6 polish's first evaluation** repeats the scan's last accepted point (one evaluation and one score per polished
  maximum, about 1.5% of a fit): bit for bit, small; fold into C4's PR.
- **P-9, P-12, P-13, P-22, P-23:** not re-proposed; nothing above covers an argument their measurements missed.
- **Threads and scheduling:** both tools hand out work dynamically and ran at efficiency 0.997-1.00 (D-186, D-208); the
  Jeffreys tables (built once per n at start-up with all threads; about 10 minutes single-threaded here) and the row
  cache's mutex per lookup are not measurable at run level.

## Reproducing the re-profile

Build `cpu-release` with `-DCMAKE_CXX_FLAGS="-g -fno-omit-frame-pointer"`. Then:

- `study_parametric_bootstrap --replicates 2 --scenarios 36,37,38,39,40,41,42,43,44 --threads 1 --out rows.csv
  --timing-out s10.csv` and `study_lrpd_serial_m6 --replicates 2 --scenarios 36,...,44 --threads 1 --boot-sample 8 --out
  m6.csv --timing-out m6.csv`; shares per phase and cell from the samples with D-187's schedule (`project_runtime.py`
  reads one scenario per cell as a nine-scenario run; multiply by nine, and weight the M6 subset-only phases by the
  subset's cells).
- Inside the phases: `valgrind --tool=callgrind --read-inline-info=yes` on a driver that, for one scenario and replicate
  0, fits the panel through `evaluate_surface_cached` and `calibrate_from_surface`, then for b = 0..B−1 simulates the
  bootstrap panel at θ̂ (`kBootSeed`, replicate b), fits it the same way with the check on b % 20 == 0, takes
  `analytic_se_scaled`, and `profile_lr_statistic` on the first 199 B/999; run the loop twice and instrument only the
  second (`CALLGRIND_START_INSTRUMENTATION`) so the fills are excluded. For M6, callgrind on the tool itself
  (`--phis 0.6 --replicates 1 --scenarios 42 --threads 1`).
- Counts: a driver that calls `NewtonProfile::polish` and the two `inner` solves itself and reads `evaluations`, with a
  map of evaluated (u₀, u₁) pairs for the repeats; the M6 search copied with counters at its line search, Hessian
  refreshes and convergence exits; and the filter's `log_likelihood` under `SplitRule`s built from
  `gauss_hermite_rule(64 | 32)` and `CompositeLegendre{gauss_legendre_rule(16), 8}` against the parity rule.
