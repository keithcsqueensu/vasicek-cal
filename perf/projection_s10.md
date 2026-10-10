# The runtime projection, calibrated on S-10 (D-187)

The cell-weighted timing-sample projection (`perf/README.md`, "Projecting a full run") checked against the
one full run with recorded timings: S-10/S-4 (D-181). That run was 81 scenarios × 1,000 replicates, and
took 23,190 s wall and 555,004 thread-seconds on 24 threads, with an Intel Core Ultra 9 275HX and MSVC
19.51 Release.

## The check

- **The binary:** the run's own, d800d34 (SHA-256 5d80ae13…f7dd), from before P-15. It was run with
  `--replicates 1` and then `--replicates 2` over all 81 scenarios.
- **The data:** that binary predates `--timing-out`, so the projection used its per-scenario phase table,
  the same thread-seconds summed over a scenario's sampled replicates.
- **The arithmetic:** each cell's mean cost per replicate, with the phase schedule:
  - Jeffreys × 0.05 and the S-4 widths × 0.20;
  - the surface phase counted once, as the cache fill;
  - the other phases × 1;
  - then times its 9 scenarios × 1,000 replicates.
- **The range:** a bootstrap over the cell's sampled scenarios. The tool's `--timing-out` now records jobs,
  and `perf/project_runtime.py` resamples those.

| Sample | Jobs per cell | Sample wall | Projection | 90% range | Actual | Actual / projection |
|---|---|---|---|---|---|---|
| 1 replicate per scenario | 9 | 107 s | 5.33 h | 4.98–5.66 h | 6.44 h, outside | 1.21 |
| 2 replicates per scenario | 18 | 135 s | 6.11 h | 5.68–6.54 h | 6.44 h, inside | 1.05 |
| subset-based (81/9), 1 per scenario | — | — | 3.91 h | — | 6.44 h | 1.64 |
| subset-based (81/9), 2 per scenario | — | — | 4.60 h | — | 6.44 h | 1.40 |

## Why one replicate per scenario missed

The miss is replicate noise, not a short-run effect:
- **Single replicates disagree:** replicate 0 alone projects 460,195 thread-seconds (excluding the surface
  phase), and replicate 1 alone 595,950. The full run, on the same accounting, was 538,207.
- **Single-replicate costs are skewed:** across scenarios, one replicate costs 0.44–2.53× its scenario's
  full-run mean.

With 9 jobs per cell the sample mean is noisy, and the bootstrap range under-states it for costs this
skewed. Hence the default of 2 replicates per scenario (≥ 18 jobs per cell in the recovery matrix), and the
rule of sampling until the range is within ±10%.

## Status

**Calibrated on S-10.** This is one run, and the method's choices (the phase schedule, the default sample
size) were made against it. The next full study run is its test: the projection and its range are recorded
before the run starts and compared with the actual wall time after it ends. If the actual lands in range
again, the method is called validated (D-187).

For reference, the current code (P-15, D-185) projects the same S-10/S-4 run at 4.02 h (90% range
3.85–4.18 h, ±4%) from a 114 s sample at 2 replicates per scenario: `python perf/project_runtime.py --run
"build/cpu-release/tests/study_parametric_bootstrap.exe"`.
