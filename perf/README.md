# perf/: performance measurement

Every performance item (P-item) of the performance track (D-177), its benchmark and its write-up, plus
the runtime-projection method used before every full run. P-items change how results are computed,
never what; they measure speed and never decide a result. The studies' own performance questions
(S-16 FP32 search, S-17 GPU scaling) stay studies, registered in [`studies/`](../studies/README.md).

Contents: [P-items](#p-items) · [Rules for every P-item](#rules-for-every-p-item-d-177) ·
[Projecting a full run](#projecting-a-full-run-d-187) · [Benchmarks](#benchmarks)

## P-items

| ID | Item | Status | Decision | Changes results? | Write-up |
|---|---|---|---|---|---|
| P-1 | Surface rows cached by (n, d) | done | D-167 | no | `study_surface_cache_spike` (benchmarks below) |
| P-2 | Profile and q solves by safeguarded Newton | done | D-170 | within solver tolerance; 0 flips | `study_surface_cache_spike` |
| P-3 | Study fits compute only scored arms | done | D-171 | no | `study_surface_cache_spike` |
| P-4 | Bootstrap reduction over distinct rows | done | D-172 | within rounding; 0 flips | `study_surface_cache_spike` |
| P-5 | Bounded argmax, full-grid fallback | done | D-172 | no | `study_surface_cache_spike`, `study_argmax_check` |
| P-6 | SIMD reduction across grid points | done: the reduction 3.3–3.7× (AVX2), the W × L bootstrap phase −33%; S-10 ~0% | D-199, D-404 | no | [p6_simd_reduction.md](p6_simd_reduction.md) |
| P-7 | SIMD `erfcx` (narrowed; `log_phi_pair` extension done) | done: `log_contrib` 33% faster (SSE2, AVX2) | D-199, D-212 | no | [p7_erfcx.md](p7_erfcx.md) |
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
| P-21 | One `erfcx` per node on the derivative path and in the hint | done: `log_contrib_derivs` ~3× faster, `log_contrib` 15% (after P-7) | D-213, D-400 | no | [p21_shared_erfcx.md](p21_shared_erfcx.md) |
| P-22 | SIMD derivative node functions | **rejected** (measured 0–3%, expected 8–12%; ~300 lines and a third kernel family) | D-213, D-402 | — | [p22_simd_derivs.md](p22_simd_derivs.md) |
| P-23 | Hint and extent search across grid points | **proposed** | D-213 | no | [p21_p23_simd_sweep.md](p21_p23_simd_sweep.md#p-23-hint-and-extent-search-across-grid-points) |

- **Proposed, not yet ordered (D-213):** P-23 (P-21 done, D-400; P-22 rejected, D-402), from a SIMD sweep, with a scope
  note for P-6 ([p21_p23_simd_sweep.md](p21_p23_simd_sweep.md)).
- **Order of the open items (D-199):** P-7 done (D-212); next P-13 (third) and P-6 (fourth); P-19 with the next golden
  refresh; P-20 later, with its own flip check.
- **Also future (D-199):** the in-house `log`/`exp` for `log_phi`, a decision of its own (within rounding;
  would make results platform-independent).
- **Rejected:** P-12 and P-22 (above), and sequential stopping of Monte Carlo runs (a design change, not a
  performance item).
- **The runtime projection** ([below](#projecting-a-full-run-d-187)) is **validated** on the M6 validation run
  (3.499 h against 3.3–3.7 h; D-202), and in range again on the S-35 M6 arm (34.94 h against 32.9–37.0 h;
  D-208).

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
| P-7 | `p7_erfcx` (`perf/p7_erfcx.cpp`; CTest `p7_erfcx_identity` runs it briefly), the same source built on `main` for "before" (`perf/p7_erfcx.md`) | µs per `log_contrib` and `log_contrib_derivs` by erfcx kernel on the recovery grid, every kernel checked against the no-kernel path first; identity: `unit_erfcx_batch`, `unit_erfcx_kernel_*`, `xref_core_regression_erfcx_*`, `unit_ar1_filter_erfcx_*` |
| P-21 | `p7_erfcx` against P-7's branch built the same way; a cross-build dump of hints, values and derivatives against `main` (`perf/p21_shared_erfcx.md`) | µs per `log_contrib` and `log_contrib_derivs`; identity: `unit_erfcx_batch` (`inverse_mills_pair_*`, `derivative_functor_*`) |
| P-18 | `perf/p18_ar1_transition.md`: one evaluation plus one score per (T, n) cell, the libm filter against each transition kernel; the transition kernel's ns per term; `study_lrpd_serial_m6` projected before and after | ms per evaluation plus score by kernel; the S-35 M6 arm's projection; identity: `unit_ar1_transition`, `unit_ar1_kernel_*`, `unit_ar1_filter_kernel_*` |
| P-17 | `p17_simd_philox` (`perf/p17_simd_philox.cpp`; CTest `p17_simd_philox_identity` runs it briefly); `study_parametric_bootstrap` before and after, in ABBA order (`perf/p17_simd_philox.md`) | ms per panel by Bernoulli kernel (scalar, SSE2, AVX2) on the recovery grid's (T, n) cells, every panel checked against `simulate_panel`; the S-10 tool's simulate phase and rows |
| P-21 to P-23 | `simd_scout` (`perf/simd_scout.cpp`; built, run by hand; `perf/p21_p23_simd_sweep.md`) | per period, single-threaded: `log_contrib`, `log_contrib_derivs` and the hint alone in µs; integrand evaluations split into the rule's nodes and the composite rule's extent search; per-call ns of `log_phi_pair`, the inverse-Mills pair, `erfcx` and the composite node transform |

New benchmarks for P-6 onward go in this directory. P-9 started with a measurement, not a build: the
repeat rate of identical panels and the saving it projects (`p9_repeat_rate`, `perf/p9_repeat_rate.md`,
D-189: about 1%, not built).
