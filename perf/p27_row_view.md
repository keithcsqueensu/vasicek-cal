# P-27: bootstrap fits on cached rows in place, calibrate's argmax bounded by cached tile maxima (D-410)

The rest of the scouting note's C2 ([scouting_2026-10.md](scouting_2026-10.md#c2-bootstrap-fits-read-cached-rows-in-place-argmax-bounded-by-cached-tile-maxima);
D-500), after P-26 removed the zero fill (D-409). No change in results.

## The ceiling, measured before the build (2026-11-06)

`main` at 99c923a (P-24 to P-26 merged), cloud Xeon, 2 vCPU (AVX2), GCC 13, Release;
`study_parametric_bootstrap --replicates 10 --scenarios 42,44,51 --threads 2`, twice, with a throwaway timer around
`calibrate_from_surface`'s full-grid argmax (`reduce_weighted`; not committed). Per-replicate work is weighted as the
full run weights its phases (Jeffreys ×0.05, the S-4 widths ×0.20). All three scenarios have T = 100 on the recovery
grid (K = 61 × 41 = 2,501), so the argmax costs the same per fit in each: 5.49 and 5.42 s over 30,030 fits, about
0.18 ms per fit, a third to each scenario.

| Scenario (T, n) | per-replicate work | surface (the copy) | calibrate's argmax | ceiling |
|---|---|---|---|---|
| 42 (100, 100) | 9.46, 9.64 s | 0.82, 0.84 (8.7%) | about 1.82 (19%) | about 28% |
| 51 (100, 100) | 10.04, 10.14 s | 0.74, 0.74 (7.3%) | about 1.82 (18%) | about 25% |
| 44 (100, 10,000) | 55.77, 54.89 s | 1.93, 1.88 (3.4%) | about 1.82 (3.3%) | about 6.7% |

The argmax is over half the fit phase in the n = 100 cells (1.82 of 3.2 s), more than the scouting note's instruction
counts suggested (10% of the cell's instructions). Over the run, with the scouting note's cell weights (n = 100 cells
14%, n = 1,000 27%, n = 10,000 59%; the n = 1,000 cells taken halfway, about 12%), the ceiling is **about 10%**,
well above the stop rule's 4% for building at all.

## Expected, written before the build (2026-11-06)

- **What goes:** for every bootstrap fit, the T × K copy of cached rows into `L` (2 MB at T = 100), and calibrate's
  full-grid argmax over K × T terms. What stays: the T row lookups in the cache; the tile bounds, T × 176 terms (176
  tiles of 4 × 4 on the recovery grid, D-172's tiling) with the T × 176 tile maxima gathered into one buffer; the
  visited tiles, 16 points × T each; and the copy for the 199 of 999 panels that W* fits from a contiguous surface.
- **Expected:** the surface phase falls by about 80% (only W*'s panels copy) and calibrate's argmax by about 85-90%.
  **S-10 per-replicate work: about −20% at (100, 100), about −5% at (100, 10,000); about 8-9% of the run.**
- **Stop rule:** under 3% of per-replicate work on scenarios 51 and 44 together (A B B A, the same machine state):
  recorded, not merged.
- **Bit for bit:** the bounded search returns the full grid's argmax exactly (index, value and tie-break), the
  property D-172 proved and its slow tests check; with unit weights and one replicate it is calibrate's
  `reduce_weighted`. Its guard keeps exactness: a panel with any non-finite value in a row takes the full grid, here
  by copying the rows and calling the contiguous path unchanged. The stencil sums through row pointers repeat
  `weighted_sum`'s operations in its order.
- **Not touched:** the contiguous `calibrate_from_surface`, `evaluate_surface`, the GPU backend, the C ABI, every
  caller that holds a contiguous surface. The view is for fits whose rows come from `SurfaceRowCache`.
- **Checks:** a unit test of the row-pointer fit against the contiguous one on the recovery subset's panels,
  bootstrap panels and panels with non-finite rows (the fallback), every `Estimate2` field; the non-slow suite; the
  S-10 tool's rows against `main` byte for byte; the recovery harness's fits at R = 50, every field.

## What was built

- **`engine/argmax_tiles.hpp`:** D-172's tiling constants, moved unchanged from `resample/bootstrap.hpp` (which now
  takes them from here), and `tile_maxima`, a row's maximum over each 4 × 4 tile, computed as the bootstrap's search
  computes it.
- **`SurfaceRowCache`** (`engine/surface_cache.hpp`): each entry also keeps its row's tile maxima and whether every
  value is finite, computed once at fill. `entry()` returns them; `row()` is unchanged. An entry is 1.4 KB larger on
  the recovery grid (176 tiles; the row is 20 KB).
- **`SurfaceRows`, `surface_rows_cached`, `copy_rows`:** a panel's surface as T pointers into the cache, with the
  lookups of `evaluate_surface_cached` (so the cache's fills, hits and misses are unchanged), and the copy it stands
  for.
- **`weighted_sum_rows`** (`core/reducers/weighted_sum.hpp`): `weighted_sum` through row pointers, the same operations
  in the same order.
- **`calibrate_from_rows`** (`engine/calibrate.hpp`): D-172's bounded search on unit weights (the tile maxima gathered
  into one T × 176 buffer, the tile bounds through P-6's kernel, the tiles visited in descending order of their
  bound, each point summed through the rows), then the rest of calibrate. That rest is now `detail::
  calibrate_from_argmax`, shared verbatim with `calibrate_from_surface`, whose argmax and stencil are unchanged. A
  panel with a non-finite value in any row copies its rows and runs `calibrate_from_surface` (D-172's guard).
- **`study_parametric_bootstrap`:** every bootstrap fit uses `calibrate_from_rows`; only W*'s 199 of 999 panels copy
  their rows, since `profile_lr_statistic` takes a contiguous surface. The panel's own fit, once per replicate, is
  unchanged.

## Identity

- **`calibrate_from_rows_is_calibrate_from_surface_bit_for_bit`** (`unit_surface_cache`): every `Estimate2` field the
  same, and every cached row's tile maxima equal to `tile_maxima` on it.
  - 324 cached panels: two recovery replicates of all 81 scenarios, and a bootstrap-like panel drawn at each estimate.
  - 60 hand-made surfaces: smooth, flat with ties, huge and tiny values, and 15 with an infinite or NaN value (the
    fallback). 0 differ.
- **The non-slow suite** passes with the committed goldens (81 of 81, GCC 13, warnings as errors), and so do the slow
  `study_argmax_check`, `study_parametric_bootstrap_smoke` and `study_parametric_bootstrap_resume`.
- **The S-10 tool's rows** of all four timing runs below: byte-identical to `main`'s.
- **The recovery harness** (which fits through `calibrate_cached`, the contiguous path, and so checks the split of
  calibrate and the moved constants): `--replicates 50 --save-fits` on `main` and this branch, every field of all
  4,050 fits identical (the two 4-byte padding slots of the raw `rc::Fit` records skipped), and the printed summaries
  identical.

## Measured (2026-11-06)

Cloud Xeon, 2 vCPU (AVX2), GCC 13, Release; `study_parametric_bootstrap --replicates 10 --scenarios 51,44
--threads 2`, `main` (99c923a) and this branch A B B A; thread-seconds summed over the 10 replicates; per-replicate
work weighted as the full run weights its phases.

| | (100, 100), scenario 51: `main` | P-27 | change | (100, 10,000), scenario 44: `main` | P-27 | change |
|---|---|---|---|---|---|---|
| per-replicate work | 10.45, 10.56 | 8.70, 8.67 | **−17.3%** | 56.81, 55.44 | 53.73, 53.07 | **−4.9%** |
| surface | 0.90, 0.91 | 0.33, 0.33 | −63% | 2.21, 1.87 | 0.90, 0.77 | −59% |
| argmax + refine + SE | 3.41, 3.50 | 2.31, 2.31 | −33% | 10.24, 9.85 | 8.56, 8.48 | −15% |
| W* (S-4b) | 4.34, 4.36 | 4.28, 4.24 | −2% | 17.10, 17.04 | 17.62, 17.24 | +2% |
| simulate (unchanged code) | 1.32, 1.32 | 1.30, 1.32 | −1% | 24.12, 23.59 | 23.57, 23.48 | −1% |

- **Against the stop rule** (under 3% on scenarios 51 and 44 together): −6.8%, not triggered.
- **Against the expectation:** per-replicate work held: −17% against about −20% at (100, 100), −5% against about −5%
  at (100, 10,000). Inside it, both phases fell less than expected:
  - **The surface phase:** −60% against about −80%. The row lookups (a mutex and a key string per distinct
    observation) are a larger part of what stays than the copy of W*'s panels.
  - **Calibrate's argmax:** from about 1.8 s to about 0.7 s at (100, 100), about 60% of it removed against 85-90%
    expected. The bounded search's own cost (gathering the T × 176 tile maxima, the tile bounds, the sort, and the
    visited tiles summed point by point through the rows) is about 40% of the full grid's.
- **Run level:** with the scouting note's cell weights (the n = 1,000 cells taken between the two measured), about
  7%, a little under the 8-9% expected. The owner's laptop measures it on the next full run.
- **Left:** the row lookups could go too (the panel's distinct observations resolved to cache entries once, and the
  bootstrap's panels looked up by (n, d) in a flat table), and the visited tiles could be summed in lanes through the
  rows. Neither is proposed until the laptop's full-run phase split shows what they are worth.
