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
