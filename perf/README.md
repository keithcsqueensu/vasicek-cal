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

## Benchmarks

| P-item | Benchmark | What it reports |
|---|---|---|
| P-1 to P-5 | `study_surface_cache_spike --profile` (`tests/studies/surface_cache_spike.cpp`; kept there, where it was built, D-167) | per-phase CPU seconds on the recovery subset (surface, calibrate, profile, q interval, bootstrap reduction vs refinement), evaluations per fit, cached and uncached, and that every fit equals `recovery::fit` bit for bit |
| P-5 | `study_argmax_check` (`tests/studies/argmax_check.cpp`; slow CTest at R = 40) | bounded against full-grid argmax, every bootstrap and jackknife row of the subset |

New benchmarks for P-6 onward go in this directory. P-9 starts with a measurement, not a build: the
repeat rate of identical panels per scenario on the recovery subset and the saving it projects.
