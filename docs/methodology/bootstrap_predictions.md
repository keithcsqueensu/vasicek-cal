# Bootstrap coverage: predictions made before the run (M2b, D-136)

Written and committed **before** the recovery harness computed any bootstrap interval, at the
owner's request, so that the review compares results with expectations instead of rationalising
whatever comes out. The git history shows this file predates the results. The results, and how
they compare, will be recorded in [recovery.md](recovery.md) and D-137; this file is not edited
after the run.

## What is run

- **Scenarios:** every scenario of the recovery matrix (81 × 1,000 simulated panels; see
  [recovery.md](recovery.md)).
- **Bootstrap:** the iid bootstrap of periods, B = 999 replicates per panel. Each panel's
  bootstrap stream uses seed = recovery seed XOR (scenario << 32 | replicate).
- **Interval:** 95% percentile interval (type-7 quantiles), with grid-edge replicates kept at the
  bound.
- **Coverage:** over all 1,000 panels. Nothing is deferred.
- **Verdict:** the same band as for the profile interval, 0.927–0.973.
- **Not in this study:** the moving-block bootstrap, whose coverage study belongs with the AR(1)
  data in M6 (D-133).

## Scenario groups

Groups are fixed by the M1.8 results (Wald verdict per scenario):

| Group | Scenarios | Count |
|---|---|---|
| A: near-uninformative or near a bound (Wald DEFERRED) | 0, 1, 2, 3, 4, 6, 7, 9, 10, 12, 15, 18, 19, 20, 21, 22, 23, 24, 27, 28, 30, 33, 36, 39, 45, 46, 48, 54, 57, 60, 72 | 31 |
| B: informative, T = 20 | 11, 29, 37, 38, 47, 55, 56, 63, 64, 65, 73, 74 | 12 |
| C: informative, T = 40 | 5, 13, 14, 31, 32, 40, 41, 49, 50, 58, 59, 66, 67, 68, 75, 76, 77 | 17 |
| D: informative, T = 100 | 8, 16, 17, 25, 26, 34, 35, 42, 43, 44, 51, 52, 53, 61, 62, 69, 70, 71, 78, 79, 80 | 21 |

## Predictions

1. **Group A, ρ: mostly undercoverage, often severe; not conservative.**
   - **Prediction:** at least two thirds of these 31 ρ verdicts fall below the band, many below
     0.90.
   - **Reason:** the bootstrap is inconsistent when the parameter is at or near the boundary of
     the parameter space (Andrews 2000). When ρ̂ sits near its lower bound, resampled estimates
     pile up at the bound, and the percentile interval shrinks towards it.
   - **The contrast:** this is the opposite of the profile interval, which was *conservative*
     in 20 such settings because truncation widened it.
2. **Group A, PD: mixed.** Undercoverage where PD̂ itself often sits at its lower bound (n = 100,
   PD 0.1%: scenarios 0, 3, 6, 9, 12, 15, 18, 21, 24). Mostly in the band elsewhere.
3. **Group B (T = 20), both parameters: undercoverage, more than profile.**
   - **Prediction:** coverage of about 0.88–0.93, so most of these 24 verdicts fall below the
     band.
   - **Reason:** percentile intervals from 20 periods are too short, like Wald's.
   - **For ρ:** the misses fall mostly below the truth (the interval lies entirely under it),
     because the bootstrap distribution is centred on a downward-biased ρ̂ and the percentile
     method does not correct for that bias.
4. **Group C (T = 40):** PD mostly in the band. ρ mixed, with some verdicts below it
   (about 0.91–0.93).
5. **Group D (T = 100):** both parameters in the band.
6. **CONSERVATIVE verdicts (above the band):** none or very few (at most 5).

**Totals implied, out of 162 verdicts:**

- below the band: roughly 50–90;
- in the band: roughly 70–110;
- above the band: at most 5.
