# perf/: performance measurement

Performance benchmarks for the P-items of the performance track (`STATE.md`, "Performance"; D-177).
They measure speed; they are not studies, and they never decide a result. The studies' performance
questions (S-16 FP32 search, S-17 GPU scaling) stay studies, registered in `studies/`.

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
   projection alongside for comparison. The agreement check catches a bias the range cannot: a one-time
   cost that a tool does not time as a `once` phase is charged to too few jobs in a small sample.

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

## Benchmarks

| P-item | Benchmark | What it reports |
|---|---|---|
| P-1 to P-5 | `study_surface_cache_spike --profile` (`tests/studies/surface_cache_spike.cpp`; kept there, where it was built, D-167) | per-phase CPU seconds on the recovery subset (surface, calibrate, profile, q interval, bootstrap reduction vs refinement), evaluations per fit, cached and uncached, and that every fit equals `recovery::fit` bit for bit |
| P-5 | `study_argmax_check` (`tests/studies/argmax_check.cpp`; slow CTest at R = 40) | bounded against full-grid argmax, every bootstrap and jackknife row of the subset |
| P-15 | `study_parametric_bootstrap` before and after, in ABBA order; `recovery_harness --save-fits` before and after for the flip check (`perf/p15_analytic_se.md`) | calibrate-phase CPU seconds; every recovery fit's fields and coverage flips |

New benchmarks for P-6 onward go in this directory. P-9 starts with a measurement, not a build: the
repeat rate of identical panels per scenario on the recovery subset and the saving it projects.
