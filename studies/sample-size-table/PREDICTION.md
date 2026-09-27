# S-2 Sample-size planning table: predictions made before the run

Pre-registration for study S-2 (`studies/README.md`, D-148). This file is committed **before**
the table is computed. No RMSE value in the recovery summary was read in writing it: the
predictions come from large-n theory, stated below, and from facts already recorded in earlier
entries (RMSE falls with T in every cell, D-124; which scenarios are near-uninformative, D-136).
The results will be compared with it prediction by prediction, and misses will be reported, not
explained away.

**This file is not edited after the run.** A change of method goes in a new, dated section
appended at the end, committed before the changed run.

## The question

How many years (periods) of data are needed to estimate ρ within ±0.05 (absolute), and PD within
±25% of itself (relative), both at 95%, as a function of PD, ρ and the number of obligors n?

## What is run

No new fits. A stdlib-only Python script, `studies/sample-size-table/sample_size_table.py`, reads
the committed recovery summary (`tests/golden/recovery/summary.csv`: 81 scenarios, R = 1,000,
bias and RMSE over all replicates, grid-edge estimates at their grid value, D-121). For each of the
27 cells (PD, ρ, n) and each parameter:

1. **Points used.** The three scenarios T ∈ {20, 40, 100}. A scenario is **mostly at a bound** if
   more than half its replicates are on the grid edge (the summary's `edge` count > 500). Such
   scenarios are left out of the fit. If all three are, the cell is **not estimable at any T
   studied**.
2. **Scaling check.** Ordinary least squares of log RMSE on log T over the points used gives a
   slope b and the largest absolute residual. **The scaling holds** if three points are used,
   b ∈ [−0.75, −0.30], and every residual is at most 0.15 (15% in RMSE). With two points it is
   **unchecked**.
3. **The implied T.** With the slope fixed at −½, C = exp(mean over the points used of
   (log RMSE + ½ log T)), so RMSE(T) ≈ C/√T, and T* = (1.96·C/target)², with target 0.05 for ρ
   and 0.25·PD for PD. T* is reported rounded up to a whole year. It is also reported under the
   fitted slope, where the scaling is checked and holds: T* = T₁·(1.96·RMSE(T₁)/target)^(−1/b) from
   the fitted line.
4. **Labels.** T* outside [20, 100] is marked an **extrapolation**. A cell where the scaling fails
   carries that label beside its T*, since T* then rests on a law the data do not follow.

RMSE includes bias. Bias falls like 1/T, faster than the spread, so at T = 20 it slightly inflates
C, which makes T* a little conservative. That is acceptable for a planning aid, and the scaling
check shows where it is not small.

The optional check in the original entry (refitting the subset at one implied T) is **not
planned**. It would be added under an appended prediction if the owner asks.

## Mechanism: the large-n benchmark

As n → ∞ the default rate reveals x_t = Φ⁻¹(DR_t) = (c − √ρ·Z_t)/√(1 − ρ), with c = Φ⁻¹(PD). So
x_t is normal with variance σ² = ρ/(1 − ρ), and:

- **ρ:** Var(σ̂²) ≈ 2σ⁴/T, and dρ/dσ² = (1 − ρ)², so SD(ρ̂) ≈ √2·ρ(1 − ρ)/√T;
- **PD:** ĉ = √(1 − ρ)·mean(x) has SD √ρ/√T, and dPD/PD = λ(c)·dc with λ(c) = φ(c)/Φ(c), so the
  relative SD of PD̂ is λ(c)·√ρ/√T.

This gives the benchmark T* (years), with no binomial noise and no bounds:

| PD | ρ | T* for ρ ± 0.05 | T* for PD ± 25% |
|---|---|---|---|
| 0.1% | 0.02 / 0.12 / 0.24 | 2 / 35 / 103 | 14 / 84 / 168 |
| 1% | 0.02 / 0.12 / 0.24 | 2 / 35 / 103 | 9 / 53 / 105 |
| 5% | 0.02 / 0.12 / 0.24 | 2 / 35 / 103 | 6 / 32 / 63 |

**At finite n** the binomial noise adds to the systematic variance and raises T*, most where the
expected number of defaults per period, n·PD, is small. At n·PD ≲ 1 (PD 0.1% with n ≤ 1,000, and
PD 1% with n = 100) it dominates. **Near a bound** the RMSE is capped by the distance to the
bound, so it falls more slowly than T^(−½) and the scaling check can fail. At ρ = 0.02, the
±0.05 target is lenient: even an estimate stuck at the floor (ρ = 0.001) is 0.019 from the truth.

## Predictions

"Informative" cells are those whose three scenarios are all outside group A of D-136 (not
near-uninformative, not near a bound).

| # | Prediction |
|---|---|
| Q1 | **Scaling holds** (step 2) for PD in at least 19 of 27 cells and for ρ in at least 16 of 27. Failures concentrate in cells with ρ = 0.02 or n·PD ≤ 1: at least two-thirds of the failures of each parameter are in such cells |
| Q2 | **Not estimable at any T studied:** at most 3 cells, all with PD = 0.1% and n = 100 |
| Q3 | **Monotone in the design:** fixed-slope T* for ρ falls as n rises (at fixed PD, ρ) and rises with ρ (at fixed PD, n); T* for PD falls as n rises and as PD rises. At least 90% of the adjacent pairs, in each direction, go the predicted way |
| Q4 | **ρ at 0.02:** T* ≤ 20 (at or below the studied range) in at least 7 of the 9 cells |
| Q5 | **ρ at 0.24:** T* > 60 in all 9 cells, and T* > 100 (beyond the studied range) in at least 6 |
| Q6 | **ρ at 0.12, n = 10,000, PD ≥ 1%:** T* between 25 and 70 in both cells; with n = 100, T* > 100 in all 3 cells |
| Q7 | **Large-n benchmark,** in the 6 cells with n = 10,000 and PD ≥ 1%. For PD, the fitted C/PD is within ±35% of λ(c)·√ρ in at least 5 of the 6. For ρ, the fitted C is within ±35% of √2·ρ(1 − ρ) in at least 3 of the 4 cells with ρ ≥ 0.12. The two ρ = 0.02 cells are reported but not counted, since the floor at 0.001 caps their RMSE |
| Q8 | **PD ± 25%, n = 10,000:** T* in [22, 50] for (5%, 0.12); [36, 85] for (1%, 0.12); [75, 170] for (1%, 0.24); > 100 for (0.1%, 0.24); ≤ 20 for (5%, 0.02) and (1%, 0.02) |
| Q9 | **PD at n·PD ≤ 1** (PD 0.1% with n ∈ {100, 1,000}; PD 1% with n = 100): T* > 100 for PD ± 25% in every such cell that is estimable |
| Q10 | **Headline:** no cell with ρ ≥ 0.12 reaches both targets within T ≤ 40, apart from at most two cells with PD = 5% and n ≥ 1,000 |

**What would be a surprise worth its own diagnosis:** a fitted slope steeper than −0.75 in an
informative cell. That would mean the T = 20 RMSE is inflated by more than bias, for example by a
heavy tail of estimates, and T* read off the −½ law would then be too pessimistic.

## Verdicts

S-2 is a planning aid, not a verdict family. It changes no pinned verdict. Each cell's result is
a T* with its labels (extrapolation, scaling fails or unchecked, not estimable).

## Cost

Seconds: one pass over 81 rows.

## How the comparison will be reported

Each prediction Q1–Q10 gets a row: prediction, result, held or not. A prediction that holds in
direction but misses in size is recorded as not held, with the size of the miss. The table, the
comparison and the monitoring implication go in S-2's entry in `studies/README.md` and in a
D-entry.
