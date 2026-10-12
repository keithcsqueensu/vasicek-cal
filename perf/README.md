# perf/: performance measurement

Every performance item (P-item) of the performance track (D-177), its benchmark and its write-up, plus
the runtime-projection method used before every full run. P-items change how results are computed,
never what; they measure speed and never decide a result. The studies' own performance questions
(S-16 FP32 search, S-17 GPU scaling) stay studies, registered in [`studies/`](../studies/README.md).

Contents: [P-items](#p-items) · [Proposed SIMD items](#proposed-simd-items-2026-10-12) ·
[Rules for every P-item](#rules-for-every-p-item-d-177) · [Projecting a full run](#projecting-a-full-run-d-187) ·
[Benchmarks](#benchmarks)

## P-items

| ID | Item | Status | Decision | Changes results? | Write-up |
|---|---|---|---|---|---|
| P-1 | Surface rows cached by (n, d) | done | D-167 | no | `study_surface_cache_spike` (benchmarks below) |
| P-2 | Profile and q solves by safeguarded Newton | done | D-170 | within solver tolerance; 0 flips | `study_surface_cache_spike` |
| P-3 | Study fits compute only scored arms | done | D-171 | no | `study_surface_cache_spike` |
| P-4 | Bootstrap reduction over distinct rows | done | D-172 | within rounding; 0 flips | `study_surface_cache_spike` |
| P-5 | Bounded argmax, full-grid fallback | done | D-172 | no | `study_surface_cache_spike`, `study_argmax_check` |
| P-6 | SIMD reduction across grid points | planned, fourth (5–7% of the run) | D-199 | no | — |
| P-7 | SIMD `erfcx` (narrowed; `log_phi_pair` extension done) | next | D-199 | no | — |
| P-8 | GPU backend | M4 | D-168 | within tolerance | — |
| P-9 | Memoise identical panels | measured, not built (~1%) | D-189 | no | [p9_repeat_rate.md](p9_repeat_rate.md) |
| P-10 | Dynamic scheduling across threads | planned | D-186 | no | — |
| P-11 | Checkpoint and resume | part 1 done; part 2 open | D-179 | no | [p11_checkpoint.md](p11_checkpoint.md) |
| P-12 | Warm-started bootstrap search | **rejected** (competes with P-6; extra fallback path) | D-199 | — | — |
| P-13 | LTO/PGO for local study builds | planned, third | D-199 | no | — |
| P-14 | Calibrate's checks per distinct observation | done | D-178 | no | [p14_calibrate_distinct.md](p14_calibrate_distinct.md) |
| P-15 | Calibrate's SEs from the analytic Hessian | done | D-185 | within tolerance (a correction) | [p15_analytic_se.md](p15_analytic_se.md) |
| P-16 | Sampled, reported quadrature check on bootstrap fits | done | D-180 | no | [p16_sampled_check.md](p16_sampled_check.md) |
| P-17 | SIMD Philox for the DGP's Bernoulli draws | done | D-194 | no | [p17_simd_philox.md](p17_simd_philox.md) |
| P-18 | AR(1) transition sums on kernels with shared `det_exp` | done | D-199 | within rounding (≤ 7.2e-15) | [p18_ar1_transition.md](p18_ar1_transition.md) |
| P-19 | Analytic SE shared with calibrate (~13% of S-10) | future, with next golden refresh | D-199 | within rounding; flip check | — |
| P-20 | Warm-start S-4b's W* refits (40% of S-10) | future | D-199 | yes; own flip check | — |
| P-21 | One `erfcx` per node on the derivative path | **proposed** | — | no | [below](#p-21-one-erfcx-per-node-on-the-derivative-path) |
| P-22 | SIMD derivative node functions | **proposed** (after P-7, P-21) | — | no | [below](#p-22-simd-derivative-node-functions) |
| P-23 | Hint and extent search across grid points | **proposed** | — | no | [below](#p-23-hint-and-extent-search-across-grid-points) |

- **Order of the open items (D-199):** P-7 next, then P-13 (third) and P-6 (fourth); P-19 with the next golden
  refresh; P-20 later, with its own flip check.
- **Also future (D-199):** the in-house `log`/`exp` for `log_phi`, a decision of its own (within rounding;
  would make results platform-independent).
- **Rejected:** P-12 (above), and sequential stopping of Monte Carlo runs (a design change, not a
  performance item).
- **The runtime projection** ([below](#projecting-a-full-run-d-187)) is **validated** on the M6 validation run
  (3.499 h against 3.3–3.7 h; D-202), and in range again on the S-35 M6 arm (34.94 h against 32.9–37.0 h;
  D-208).

## Proposed SIMD items (2026-10-12)

A sweep of the CPU paths for vectorisation not yet covered by P-6, P-7, P-17 or P-18. These are
proposals for the owner's order, not decisions: none has a `DECISIONS.md` entry yet, and each would need
the [rules below](#rules-for-every-p-item-d-177) before it merges. P-7 is in progress and its scope is not
touched here; P-22 and P-23 reuse its vector `erfcx` as a building block. All three keep results bit for
bit: one lane per independent evaluation, each lane performing the scalar operations in the scalar order,
FMA contraction off, every kernel path tested as P-17.

**Measurements** (`simd_scout`; parity rule, PD 1%, ρ 0.12; one thread on a cloud VM, Intel Xeon 2.1 GHz,
GCC 13.3, Release; indicative only, to be repeated on the owner's machine before any decision):

| (n, d) | `log_contrib` | `log_contrib_derivs` | derivs / value | hint alone | hint / value | integrand evaluations (nodes + search) |
|---|---|---|---|---|---|---|
| 100, 1 | 14.8 µs | 24.9 µs | 1.68 | 1.7 µs | 12% | 128 + 0 |
| 1,000, 5 | 14.6 | 25.7 | 1.76 | 1.5 | 10% | 128 + 0 |
| 10,000, 20 | 17.0 | 26.7 | 1.57 | 5.1 | 30% | 128 + 0 |
| 100, 0 | 50.2 | 61.7 | 1.23 | 1.3 | 2.5% | 256 + 113 |
| 1,000, 0 | 52.0 | 71.6 | 1.38 | 1.2 | 2.3% | 256 + 115 |
| 10,000, 0 | 53.8 | 64.6 | 1.20 | 1.6 | 2.9% | 256 + 116 |

Per call: `log_phi_pair` 68–99 ns (two runs), `inverse_mills(x)` + `inverse_mills(-x)` 90 ns, `erfcx` 27–30 ns,
the composite rule's node transform (`sinh`, `cosh`, `log`) 28–32 ns. Two things stand out. The derivative
path pays about as much for its inverse-Mills pair as for the value itself. And a one-sided period
(d ∈ {0, n}, common in low-PD panels) costs about 3.4× an interior one, a third of it in the extent search.

### P-21: one `erfcx` per node on the derivative path

**What.** At each node `log_contrib_derivs` calls `log_phi_pair(x)` (one `erfcx(|x|/√2)`, one exp pair
`exp(−s.hi/2)·exp(−s.lo/2)`), then `BinomialLogIntegrandDerivs` calls `inverse_mills(x)` and
`inverse_mills(−x)`, which recompute the same `erfcx(|x|/√2)` twice and the same exp pair once. The node's x
is the same expression in both (`(c − √ρ z)/√(1 − ρ)`). So λ(±x) can be finished from the pair's own
e = `erfcx(|x|/√2)` and E = `exp(−s.hi/2)·exp(−s.lo/2)`, each by `inverse_mills`'s own expression for its
sign: λ(−|x|) = √(2/π)/e and λ(|x|) = (1/√(2π))·E/(1 − 0.5·e·E), as `log_phi_pair` does for log Φ. Three
`erfcx` and four `exp` per node become one and two. `detail::mixture_slopes` (the hint's Newton) calls the same
inverse-Mills pair at 0 < d < n and shares the change.

**Changes results:** no. Each value is the expression it replaces, on the same operands; the test is
`inverse_mills_pair_is_two_inverse_mills_calls_bit_for_bit` on the `log_phi_pair` test's arguments, plus
`log_contrib_derivs` (value, score, Hessian) bit for bit on the recovery grid's distinct rows.

**Cost of the change.** The accumulator calls the derivative functor only for nodes whose weight does not
underflow, after the value. The tails must be handed from the value to the derivatives: either the integrand
caches them per node, or one fused integrand returns the value and the five node functions together (it
then computes them for underflowed nodes too, discarded as before). Scalar code; no kernel paths.

**Expected:** about a third of a derivative evaluation at 0 < d < n (about 65 of 195 ns per node, from the
per-call costs above), less at d ∈ {0, n} where only one inverse-Mills call is made. That path is the
profile solves: W* of S-4b (39.8% of the S-10 run, D-199), the analytic SE (13.5%) and the refinement.
The run-level saving should be projected with the D-187 method before it is built.

### P-22: SIMD derivative node functions

**What.** After P-7 and P-21, everything `BinomialLogIntegrandDerivs` does per node is arithmetic on
(z, x, λ(x), λ(−x)): x_ρ, x_ρρ, g_x, g_xx and the five node functions, then `m_[j] += w·f[j]` and the
rescale `m_[j] *= r`. Lanes over nodes (four in AVX2, two in SSE2), the node functions written to scratch
and folded in node order by the existing accumulator, as P-18 does with its pair terms. Divisions stay
exact divisions, no reciprocals. The remaining libm calls (the exp pair, the fold's exp) stay scalar.

**Changes results:** no; scalar, SSE2 and AVX2 kernels, each bit for bit against the scalar functor, the
`VCAL_*_KERNEL` override and build-info record as P-17 and P-18.

**Depends on** P-7's node batching (the same block of nodes, the same scratch) and on P-21 (without it the
lanes wait on three scalar `erfcx`).

### P-23: hint and extent search across grid points

**What.** Two serial searches sit in front of every integral: the hint's bracket and Newton
(`binomial_mixture_hint`, 10–30% of an interior period's cost, 30% at n = 10,000) and the composite rule's
extent bisections (113–116 of about 370 integrand evaluations at d ∈ {0, n}). Neither vectorises within
one integral: each step depends on the last. Across grid points they do. A surface row is one observation
at every grid point, so neighbouring grid points share n and d, take the same branches (one-sided or not)
and differ only in (PD, ρ). The row fill (`evaluate_surface`, `SurfaceRowCache::row`) would evaluate four
grid points per call, each lane running its own search in the scalar order, finished lanes masked until
the slowest lane stops (Newton's and the bisection's iteration counts differ per lane). The predicate
`g(z) − z²/2 < target` is one `log_phi` per lane, which P-7's vector `erfcx` serves.

**Changes results:** no: each lane's iterates are the scalar iterates. Test: every row of the recovery
subset's distinct observations, bit for bit against the scalar row fill, under every kernel, including
lanes that stop at the iteration cap and the masked tail when K is not a multiple of the lane count.

**Considered and not proposed:** a four-way section search inside one integral (4 probes per step, about
2.3× fewer steps). It returns the same end point only if `below` is monotone along the ray in floating
point, which rounding near the target does not guarantee; it would be a within-rounding change with its own
flip check, while lanes over grid points give the same gain bit for bit.

**Cost of the change.** The row fill becomes a blocked loop over grid points with a scalar tail; the
`parallel_for` stays over blocks, so results still do not depend on the thread count. The profile solves
evaluate off-grid points one at a time and do not benefit.

### Scope note for P-6

P-6's entry describes a SIMD reduction across grid points, as in `reduce_weighted`, but since D-172 the engine's default is the bounded argmax
(`ArgmaxSearch::Bounded`); `reduce_weighted` runs only when its guard fails. The 5–7% of the run that P-6
targets is therefore mostly `weighted_sum` inside visited tiles. Proposed for P-6's own decision entry:

- **Lanes = one tile row.** Tiles are 4 × 4 grid points (`kArgmaxTile`), and k is contiguous along the ρ
  axis, so a tile row is four adjacent k: one AVX2 vector, two SSE2. L is read at stride K per period with
  no gather; each lane keeps its own Neumaier sum and non-finite tally in the scalar order. On the recovery
  grid (41 ρ points) the last tile column is one point wide and is masked.
- **The tile bounds,** `sum_j w_j U[j][tile]` and its magnitude for every tile of every replicate, with
  lanes over tiles (`U` is contiguous in the tile index).
- **The `U` table,** a per-tile `fmax`, exact in any order.
- **The summed-surface calls** with unit weights in `posterior.hpp`, `profile.hpp` and
  `conditional_pd.hpp`, which call the same `weighted_sum`.

Not SIMD, noted while reading: each replicate sorts all tile bounds (176 on the recovery grid) and usually
visits a few; a heap pops them in the same order (bound descending, then index) for less work. Measure
before deciding.

### Unlocked only by the in-house `log`/`exp` decision (D-199)

These remain scalar under bit-for-bit rules because they call libm; they belong in that decision's
cost-benefit, not in new P-items: `log` and `log1p` in `log_phi`, its exp pair, the online log-sum-exp's
`exp` per node, the composite rule's `sinh`, `cosh` and `log` per node (28–32 ns, about a third of a
one-sided node), the posterior's mass exponentials and the AR(1) filter's `weighted_moments`. With an
in-house `log`/`exp`, a whole `log_contrib` could run with lanes over grid points and the online fold kept
per lane, with no scratch buffer: P-23's layout extended to the integral itself.

### Checked and not proposed

- **Bootstrap index draws** (`resample/weights.hpp`): scalar Philox, one block per two draws. P-17's
  lanes would apply bit for bit, but S-10 draws about 5 × 10⁴ blocks per panel against thousands of
  quadratures per distinct row: well under 1%.
- **Grid-point setup** (`make_vasicek1f` per grid point: `probit(PD)` and two square roots, recomputed for
  every ρ): tabulating per axis value is exact, but it is tens of ns against 15–50 µs per evaluation.
- **The DGP's factor draws** (one polar Box–Muller pair per period): T per panel, negligible beside the
  Bernoulli count P-17 already vectorised.

## Rules for every P-item (D-177)

1. **A decision entry** in `DECISIONS.md` before it merges, stating whether it changes results:
   **no** (bit for bit), **within rounding** (summation order only), **within solver tolerance** (a
   solver's stopping rule), or **within a stated tolerance** (a `TOL_` entry in the register).
2. **A benchmark here,** reproducible from one command, with the machine, compiler and thread count
   recorded beside the numbers, before and after.
3. **An identity or tolerance test** in CTest: bit-for-bit equality with the path it replaces, or
   agreement to the stated tolerance on panels chosen to be hard for it.
4. **Anything that changes results** gets a golden refresh with a coverage-flip check over every
   recovery fit (PD, ρ and q coverage by each interval, flags, edge counts), one commit per P-item, the
   flip counts reported in the decision entry even when zero.

## Projecting a full run (D-187)

Before any full study run, project its runtime from a timing sample, not from the study subset. The
subset is chosen for statistical coverage (T × ρ); runtime follows T × n. S-10/S-4's projection from the
subset missed by 1.375× (5 h projected, 6.4 h actual, D-186), because the costliest cell (T = 100,
n = 10⁴) had no subset scenario. The study subset is unchanged: it stays the basis for predictions and
paired verdicts. Only the runtime projection changes.

1. **Sample every cell.** Run the tool on a few replicates of every scenario of the full matrix with
   `--timing-out FILE` (`study_parametric_bootstrap`, `study_bayes_coverage`, `study_estimator_pass`,
   `study_severe_period`, `study_jackknife_run` and `study_lrpd_serial` take it). Each (scenario, replicate) job's thread-seconds
   are recorded by phase (`tests/harness/timing_sample.hpp`).
2. **Follow the phase schedule.** Each phase carries its share of a full run's replicates: `1` for every
   replicate, a fraction for an arm run on the first replicates only (S-10's Jeffreys check on 0–49 of
   1,000: 0.05; the jackknife run's exact refits on 0–4: 0.005), or `once` for a one-time cost (filling the
   surface-row cache). A sample runs every arm on every sampled replicate, so an unweighted projection is
   biased whenever arms are sampled. A fractional phase is averaged over the sampled replicates that run it
   in a full run.
3. **Weight by cell.** A cell's cost per replicate times its scenario count times the full run's replicates,
   summed over cells. Wall time is that over (threads × efficiency), where the efficiency is measured from
   real runs, not assumed: thread-seconds / (threads × wall) of the last full run, which covers scenario
   overlap, cache fills and the tail where threads idle. S-10/S-4: 555,004 / (24 × 23,190) = **0.997**,
   the default (`--efficiency`). Each validation re-measures it, and a run that falls below it says so.
4. **Report a range, and sample until it is tight.** The range is the 5th–95th percentile of a bootstrap over
   jobs within each cell. `perf/project_runtime.py --run "<tool>"` starts at 2 replicates per scenario and
   doubles until the range is within ±10% of the projection (`--target`) **and** the projection agrees with
   the previous sample size's within the same ±10% (D-192), up to `--max`. It prints the subset-based
   projection alongside for comparison. The agreement check is a safety net for biases the range cannot
   see.
5. **Count fixed costs once (D-193).** A surface-row cache fill is paid once per distinct observation for
   the whole run. Every tool with `--timing-out` measures each job's fills (`FillMeter`, from the cache's
   per-thread tally) and records them in a `cache fill` phase with share `once`, subtracted from the phase
   they happened in. The projection is then the fixed part once plus the per-replicate part × the count.
   The fixed part is the sample's fills, which cover only the rows the sample met. In S-35 that was about
   half of the full run's 205 s of fills: 1% of the run.
6. **Sample at sustained clocks (D-193).** A short sample runs at boost clocks and projects low. S-35's
   32-replicate sample cost 0.215 s per job cold and 0.279 s straight after sustained load, against 0.275 s
   in the full run. So `--run` first discards `--warmup` seconds of runs (default 120). It then accepts a
   sample only if the sample itself ran at least `--min-wall` seconds (default 120), as well as meeting the
   range and agreement checks (step 4).
7. **State the machine-state allowance.** The 90% range covers sampling noise only. Two 128-replicate S-35
   samples, taken on different occasions, projected 24,059 and 25,990 thread-seconds, about ±4% apart. The
   report therefore also prints the range widened by `--allowance` (default ±5%). On S-35, run from cold,
   the method stopped at 128 replicates per scenario (a 134 s sample): 0.28 h, and 0.27–0.29 h with the
   allowance, against 0.293 h actual. That is a check on the run the rules were built from, not a
   validation.

**Calibrated on S-10, not yet validated** (`perf/projection_s10.md`). With 2 replicates per scenario
(18 jobs per cell) the projection was 6.13 h (90% range 5.70–6.56 h) against the actual 6.44 h. With 1 per
scenario (9 per cell) it was 5.34 h (5.00–5.68 h): a 17% miss, outside its range, because single
replicates' costs are skewed (0.44–2.53× their scenario's mean). The bootstrap range is a sampling range
only, and it under-covered at 9 jobs per cell. The method was tuned against this one run, so its next
full run is a test: the projection is recorded before the run and compared with the actual wall time after
it. If the actual lands in range again, the method is called validated.

**Second test, S-35 (D-192): a miss.** Its projection, recorded before the run from a 2-replicate sample,
was 0.39 h (0.35–0.43 h); the run took 0.293 h at efficiency 1.00. The sample charged the surface-row
cache's first fills, which `study_lrpd_serial` does not time as a `once` phase, to too few jobs: 8 and 32
replicates project 0.28 h. The agreement rule above was added for it, and stops at 16 replicates with
0.29 h (0.28–0.30 h). The method stays calibrated, not validated.

**Third test, M6 validation (D-202): in range.** With D-193's rules (warm-up, doubling with the agreement
check, a minimum sample wall time, fixed costs added once), the projection recorded before the run was 3.5 h
(3.3–3.7 h; 3.2–3.9 h with the machine-state allowance), from a sample on another seed. The run took 3.499 h,
projection / actual 1.005, on an otherwise idle machine. The method is now validated, for tools whose phases
it times.

**Fourth test, the S-35 M6 arm (D-208): in range.** The projection recorded before the run was 34.9 h (90%
range 32.9–37.0 h; D-199). The run took 125,767 s (34.94 h) in three checkpointed segments (54,728, 45,468 and
25,571 s; the gaps between them, 35 s, excluded), projection / actual 0.999. That makes two full runs in range
after D-193, one about 3.5 h and one about 35 h long. During the first segment the machine also ran single-threaded
work: the scorer's editing, its seconds-long self-test, and one CI-config check.

## Benchmarks

| P-item | Benchmark | What it reports |
|---|---|---|
| P-1 to P-5 | `study_surface_cache_spike --profile` (`tests/studies/surface_cache_spike.cpp`; kept there, where it was built, D-167) | per-phase CPU seconds on the recovery subset (surface, calibrate, profile, q interval, bootstrap reduction vs refinement), evaluations per fit, cached and uncached, and that every fit equals `recovery::fit` bit for bit |
| P-5 | `study_argmax_check` (`tests/studies/argmax_check.cpp`; slow CTest at R = 40) | bounded against full-grid argmax, every bootstrap and jackknife row of the subset |
| P-9 | `p9_repeat_rate` (`perf/p9_repeat_rate.cpp`), with a timing sample of the run it would speed up (`perf/p9_repeat_rate.md`) | repeated panels per scenario, replicate and (T, n) cell, under the ordered and the sorted key; the saving they project |
| P-15 | `study_parametric_bootstrap` before and after, in ABBA order; `recovery_harness --save-fits` before and after for the flip check (`perf/p15_analytic_se.md`) | calibrate-phase CPU seconds; every recovery fit's fields and coverage flips |
| P-18 | `perf/p18_ar1_transition.md`: one evaluation plus one score per (T, n) cell, the libm filter against each transition kernel; the transition kernel's ns per term; `study_lrpd_serial_m6` projected before and after | ms per evaluation plus score by kernel; the S-35 M6 arm's projection; identity: `unit_ar1_transition`, `unit_ar1_kernel_*`, `unit_ar1_filter_kernel_*` |
| P-21 to P-23 | `simd_scout` (`perf/simd_scout.cpp`; built, run by hand) | per period, single-threaded: `log_contrib`, `log_contrib_derivs` and the hint alone in µs; integrand evaluations split into the rule's nodes and the composite rule's extent search; per-call ns of `log_phi_pair`, the inverse-Mills pair, `erfcx` and the composite node transform |
| P-17 | `p17_simd_philox` (`perf/p17_simd_philox.cpp`; CTest `p17_simd_philox_identity` runs it briefly); `study_parametric_bootstrap` before and after, in ABBA order (`perf/p17_simd_philox.md`) | ms per panel by Bernoulli kernel (scalar, SSE2, AVX2) on the recovery grid's (T, n) cells, every panel checked against `simulate_panel`; the S-10 tool's simulate phase and rows |

New benchmarks for P-6 onward go in this directory. P-9 started with a measurement, not a build: the
repeat rate of identical panels and the saving it projects (`p9_repeat_rate`, `perf/p9_repeat_rate.md`,
D-189: about 1%, not built).
