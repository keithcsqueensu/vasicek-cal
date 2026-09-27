# S-2 Sample-size planning table: results

Run on 2026-09-27, after the prediction ([`PREDICTION.md`](PREDICTION.md)) was merged, with
`python studies/sample-size-table/sample_size_table.py` (stdlib only) on the committed recovery
summary (`tests/golden/recovery/summary.csv`, 81 scenarios × R = 1,000). No new fits; it runs in
under a second.

**How to read the grids:**
- **T\*** is the number of years (periods, n obligors each) at which 1.96 × RMSE meets the target,
  from RMSE ≈ C/√T fitted over T ∈ {20, 40, 100}.
- **Only unmarked figures are inside the verified range.** The T^(−½) law was verified at T = 20, 40
  and 100, and held in every cell where it could be checked (24 of 27, slope −0.40 to −0.59).
  - **†** marks a figure beyond 100 years;
  - **‡** marks a figure below 20 years.

  Both are extrapolations of the law, not observations. They are probably of the right size, but
  they were not observed.
- **§** marks a cell where the law is unchecked: only T = 100 was usable, since at T = 20 and 40 more
  than half the estimates are on the grid edge.
- **The relative ρ target (±25%) is a post-hoc view,** added at the owner's request after the run
  (D-305). With the slope fixed at −½, T\* scales as 1/target², so it follows from the same fitted C
  by arithmetic: T\*(relative) = T\*(absolute) × (0.05/(0.25·ρ))². It was therefore not
  pre-registered, since there was nothing left to predict, and it involves no new fits.
- **"not estimable"** marks a cell where no T was usable.

### Years for ρ within ±0.05 absolute (95%)

| PD | n | n PD | rho = 0.02 | rho = 0.12 | rho = 0.24 |
|---|---|---|---|---|---|
| 0.001 | 100 | 0.1 | not estimable | 1,803 † § | 3,371 † § |
| 0.001 | 1,000 | 1 | 24 | 181 † | 537 † |
| 0.001 | 10,000 | 10 | 3 ‡ | 56 | 206 † |
| 0.01 | 100 | 1 | 47 | 225 † | 507 † |
| 0.01 | 1,000 | 10 | 4 ‡ | 57 | 184 † |
| 0.01 | 10,000 | 100 | 2 ‡ | 36 | 119 † |
| 0.05 | 100 | 5 | 11 ‡ | 77 | 215 † |
| 0.05 | 1,000 | 50 | 2 ‡ | 40 | 114 † |
| 0.05 | 10,000 | 500 | 2 ‡ | 35 | 110 † |

- † Beyond 100 years: an extrapolation of the T^(-1/2) law, verified only at T = 20, 40 and 100.
- ‡ Below 20 years: shorter than any history studied, also an extrapolation.
- § Unchecked: only T = 100 was usable (at T = 20 and 40 most estimates are on the grid edge), so
  the law could not be checked in this cell.

### Years for ρ within ±25% relative (95%; post hoc, D-305)

| PD | n | n PD | rho = 0.02 | rho = 0.12 | rho = 0.24 |
|---|---|---|---|---|---|
| 0.001 | 100 | 0.1 | not estimable | 5,007 † § | 2,341 † § |
| 0.001 | 1,000 | 1 | 2,368 † | 502 † | 373 † |
| 0.001 | 10,000 | 10 | 272 † | 155 † | 144 † |
| 0.01 | 100 | 1 | 4,660 † | 624 † | 352 † |
| 0.01 | 1,000 | 10 | 354 † | 157 † | 128 † |
| 0.01 | 10,000 | 100 | 140 † | 100 | 83 |
| 0.05 | 100 | 5 | 1,065 † | 214 † | 150 † |
| 0.05 | 1,000 | 50 | 183 † | 111 † | 79 |
| 0.05 | 10,000 | 500 | 124 † | 95 | 76 |

- † Beyond 100 years: an extrapolation of the T^(-1/2) law, verified only at T = 20, 40 and 100.
- ‡ Below 20 years: shorter than any history studied, also an extrapolation.
- § Unchecked: only T = 100 was usable (at T = 20 and 40 most estimates are on the grid edge), so
  the law could not be checked in this cell.

### Years for PD within ±25% relative (95%)

| PD | n | n PD | rho = 0.02 | rho = 0.12 | rho = 0.24 |
|---|---|---|---|---|---|
| 0.001 | 100 | 0.1 | not estimable | 804 † § | 1,187 † § |
| 0.001 | 1,000 | 1 | 81 | 224 † | 591 † |
| 0.001 | 10,000 | 10 | 21 | 149 † | 506 † |
| 0.01 | 100 | 1 | 71 | 134 † | 250 † |
| 0.01 | 1,000 | 10 | 16 ‡ | 81 | 201 † |
| 0.01 | 10,000 | 100 | 11 ‡ | 71 | 188 † |
| 0.05 | 100 | 5 | 19 ‡ | 48 | 97 |
| 0.05 | 1,000 | 50 | 7 ‡ | 38 | 86 |
| 0.05 | 10,000 | 500 | 6 ‡ | 39 | 82 |

- † Beyond 100 years: an extrapolation of the T^(-1/2) law, verified only at T = 20, 40 and 100.
- ‡ Below 20 years: shorter than any history studied, also an extrapolation.
- § Unchecked: only T = 100 was usable (at T = 20 and 40 most estimates are on the grid edge), so
  the law could not be checked in this cell.

## What the table says

- **ρ within ±0.05 absolute is a long-history target** once ρ is moderate:
  - at ρ = 0.12 it takes about 35–40 years with 50 or more expected defaults a year (n·PD ≥ 50),
    about 56 with 10, 77 with 5, and 180–225† with 1;
  - at ρ = 0.24 it needs more than 100 years everywhere studied, 110–120† even for large
    portfolios. Those figures are extrapolated; what the data show directly is that the target
    is not met at T = 100.

  So a typical 20–30 year history does not pin ρ to ±0.05 unless ρ is small. Where n·PD ≥ 50, the
  large-n formula T\* ≈ (1.96·√2·ρ(1 − ρ)/0.05)² is close: C within 8%, so T\* within about 16%.
- **At ρ = 0.02 the ±0.05 target is lenient,** because the estimate cannot go below 0.001, so an
  estimate is never far from the truth on the low side. The small T\* there (2–4‡ years with
  n·PD ≥ 10) says the error is small in absolute terms, not that ρ is known relatively: the interval
  still reaches down to the floor.
- **Much of the high-ρ result is the target's shape (the relative view, D-305).** ±0.05 is about
  ±20% of ρ at 0.24 but about ±40% at 0.12, and ±250% at 0.02. Held to a common ±25% relative
  target, the order reverses. For large portfolios (n·PD ≥ 50):
  - ρ = 0.24 needs 76–83 years, inside the verified range, against 110–119† for ±0.05;
  - ρ = 0.12 needs 95–111 (the top end†), against 35–40;
  - ρ = 0.02 needs 124–183†.

  This is what large-n theory predicts: the relative SD of ρ̂ is √2·(1 − ρ)/√T, which falls as ρ
  rises, giving 118, 95 and 71 years at ρ = 0.02, 0.12 and 0.24. Relative to its size, ρ is
  estimated best when it is large. With few defaults the finite-n penalty dominates either way:
  at n·PD = 1 every ρ cell needs more than 350† years.
- **PD within ±25%:**
  - PD depends on ρ as strongly as ρ itself does: at 5% PD it takes about 6–7 years at ρ = 0.02,
    38–39 at 0.12 and 82–86 at 0.24 (large portfolios);
  - at ρ = 0.12, 1% PD takes 71–81 years and 0.1% PD about 150†.

  The factor cycle, not the number of obligors, limits PD's precision once n·PD ≥ 10.
- **The finite-n penalty.** The fitted C against the large-n benchmark shows how much binomial noise
  costs. For ρ at 0.12 and 0.24, C is:
  - 1.00–1.08× the benchmark at n·PD ≥ 50;
  - 1.27–1.42× at n·PD = 10, so 1.6–2.0× the years;
  - 2.2–2.6× at n·PD = 1, so 5–7× the years (all these cells are extrapolated†).

  For PD, C is 1.02–1.34× the benchmark at n·PD ≥ 50, 1.22–1.74× at n·PD = 10, and 1.54–2.85× at
  n·PD = 1.

## Prediction against result

| # | Prediction | Result | Held |
|---|---|---|---|
| Q1 | Scaling holds for PD in ≥ 19 of 27 cells, for ρ in ≥ 16; failures mostly at ρ = 0.02 or n·PD ≤ 1 | holds in 24 of 27 for both; no cell fails; the other 3 are unchecked (2) or not estimable (1), all at PD 0.1%, n = 100 | yes |
| Q2 | Not estimable: at most 3 cells, all PD 0.1% with n = 100 | 1 cell: (0.1%, 0.02, 100) | yes |
| Q3 | Monotone in n, ρ and PD in ≥ 90% of adjacent pairs, each direction | ρ vs n 17/17, ρ vs ρ 17/17, PD vs PD 17/17, PD vs n 16/17 (94%; (5%, 0.12): 37.7 at n = 1,000 against 38.0 at 10,000, a tie within noise) | yes |
| Q4 | ρ = 0.02: T\* ≤ 20 in ≥ 7 of 9 cells | 6 of 9. Misses: (0.1%, 1,000) at 24, (1%, 100) at 47, and the not-estimable cell | **no**, by 1 cell |
| Q5 | ρ = 0.24: T\* > 60 in all 9, > 100 in ≥ 6 | > 100 in all 9 (110†–3,371†; that the target is not met by T = 100 is observed, the size beyond it extrapolated) | yes |
| Q6 | ρ = 0.12, n = 10,000, PD ≥ 1%: T\* in [25, 70]; at n = 100, T\* > 100 in all 3 | 36 and 35: held. n = 100: 1,803† §, 225† and **77** at PD 5% | **no**: (5%, 0.12, 100) needs 77 years, not > 100 |
| Q7 | Large-n benchmark within ±35%: PD in ≥ 5 of 6, ρ in ≥ 3 of 4 | PD 6 of 6 (ratios 1.02–1.34); ρ 4 of 4 (1.00–1.08); the two ρ = 0.02 cells, reported only, at 1.09 and 1.02 | yes |
| Q8 | PD ± 25%, n = 10,000: (5%, 0.12) in [22, 50]; (1%, 0.12) in [36, 85]; (1%, 0.24) in [75, 170]; (0.1%, 0.24) > 100; (5%, 0.02) and (1%, 0.02) ≤ 20 | 39, 71, **188†**, 506†, 6‡, 11‡ | **no**: (1%, 0.24) needs 188, above the range by 18 (the fitted slope, −0.56, gives 163) |
| Q9 | PD at n·PD ≤ 1: T\* > 100 in every estimable cell | > 100 (all †) in 6 of 8; **81** at (0.1%, 0.02, 1,000) and **71** at (1%, 0.02, 100) | **no** |
| Q10 | No cell with ρ ≥ 0.12 reaches both targets within 40 years, apart from at most two cells with PD 5% and n ≥ 1,000 | exactly those two: (5%, 0.12, 1,000) at 40/38 and (5%, 0.12, 10,000) at 35/39 | yes |

**Six held, four missed.** The misses, not explained away:

- **Q4, Q6 and Q9 all underestimate how well the estimators do with few defaults when ρ is low or
  PD is 5%.** At ρ = 0.02 there is almost no factor variance, so the PD estimate is close to the
  pooled binomial one: relative SD ≈ 1/√(n·PD·T), which gives T\* ≈ 61/(n·PD) years. That is 61 at
  n·PD = 1, against the observed 71–81 and a prediction of > 100. The prediction treated n·PD ≤ 1 as
  uniformly hopeless, which holds at ρ ≥ 0.12 but not at ρ = 0.02. The same effect puts
  (1%, 0.02, 100)'s ρ at 47 years rather than ≤ 20, and (5%, 0.12, 100) at 77 rather than > 100.
- **Q8 misses in size at (1%, 0.24, 10,000):** 188† against an upper limit of 170. Both figures
  lie beyond the verified range: the miss is in the extrapolated law, and the data show only that
  the target is not met at T = 100. The prediction
  allowed up to 1.6× the benchmark (105), and the observed factor is 1.8×. The slope there
  is −0.56: RMSE at T = 20 carries extra skew, which the fixed −½ law turns into a conservative T\*.
- **No surprise of the kind the prediction named:** no slope is steeper than −0.59, let alone −0.75.

## Full output

`sample_size_table.py`, verbatim: the per-cell fits and the large-n benchmark. C for PD is shown
relative to PD, and the `rho_rel` rows are the post-hoc relative target (D-305), which share the
`rho` rows' C.

| PD | rho | n | n PD | param | points used | slope | max resid | scaling | C | large-n C | T* (slope -1/2) | T* (fitted slope) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.001 | 0.02 | 100 | 0.1 | rho | 0 | | | | | | not estimable at any T studied | |
| 0.001 | 0.02 | 100 | 0.1 | rho_rel | 0 | | | | | | not estimable at any T studied | |
| 0.001 | 0.02 | 100 | 0.1 | pd | 0 | | | | | | not estimable at any T studied | |
| 0.001 | 0.02 | 1,000 | 1 | rho | 3 | -0.396 | 0.019 | holds | 0.1241 | 0.02772 | 24 | 21 |
| 0.001 | 0.02 | 1,000 | 1 | rho_rel | 3 | -0.396 | 0.019 | holds | 0.1241 | 0.02772 | 2,368 † | 6,753 † |
| 0.001 | 0.02 | 1,000 | 1 | pd | 3 | -0.491 | 0.018 | holds | 1.143 x PD | 0.4762 x PD | 81 | 82 |
| 0.001 | 0.02 | 10,000 | 10 | rho | 3 | -0.498 | 0.018 | holds | 0.04204 | 0.02772 | 3 ‡ | 3 ‡ |
| 0.001 | 0.02 | 10,000 | 10 | rho_rel | 3 | -0.498 | 0.018 | holds | 0.04204 | 0.02772 | 272 † | 274 † |
| 0.001 | 0.02 | 10,000 | 10 | pd | 3 | -0.485 | 0.014 | holds | 0.5828 x PD | 0.4762 x PD | 21 | 21 |
| 0.001 | 0.12 | 100 | 0.1 | rho | 1 |  |  | unchecked | 1.083 | 0.1493 | 1,803 † | - |
| 0.001 | 0.12 | 100 | 0.1 | rho_rel | 1 |  |  | unchecked | 1.083 | 0.1493 | 5,007 † | - |
| 0.001 | 0.12 | 100 | 0.1 | pd | 1 |  |  | unchecked | 3.615 x PD | 1.166 x PD | 804 † | - |
| 0.001 | 0.12 | 1,000 | 1 | rho | 3 | -0.504 | 0.015 | holds | 0.3429 | 0.1493 | 181 † | 179 † |
| 0.001 | 0.12 | 1,000 | 1 | rho_rel | 3 | -0.504 | 0.015 | holds | 0.3429 | 0.1493 | 502 † | 493 † |
| 0.001 | 0.12 | 1,000 | 1 | pd | 3 | -0.577 | 0.036 | holds | 1.905 x PD | 1.166 x PD | 224 † | 180 † |
| 0.001 | 0.12 | 10,000 | 10 | rho | 3 | -0.511 | 0.017 | holds | 0.1903 | 0.1493 | 56 | 56 |
| 0.001 | 0.12 | 10,000 | 10 | rho_rel | 3 | -0.511 | 0.017 | holds | 0.1903 | 0.1493 | 155 † | 151 † |
| 0.001 | 0.12 | 10,000 | 10 | pd | 3 | -0.508 | 0.054 | holds | 1.556 x PD | 1.166 x PD | 149 † | 146 † |
| 0.001 | 0.24 | 100 | 0.1 | rho | 1 |  |  | unchecked | 1.481 | 0.258 | 3,371 † | - |
| 0.001 | 0.24 | 100 | 0.1 | rho_rel | 1 |  |  | unchecked | 1.481 | 0.258 | 2,341 † | - |
| 0.001 | 0.24 | 100 | 0.1 | pd | 1 |  |  | unchecked | 4.393 x PD | 1.65 x PD | 1,187 † | - |
| 0.001 | 0.24 | 1,000 | 1 | rho | 3 | -0.508 | 0.032 | holds | 0.5908 | 0.258 | 537 † | 517 † |
| 0.001 | 0.24 | 1,000 | 1 | rho_rel | 3 | -0.508 | 0.032 | holds | 0.5908 | 0.258 | 373 † | 361 † |
| 0.001 | 0.24 | 1,000 | 1 | pd | 3 | -0.507 | 0.001 | holds | 3.099 x PD | 1.65 x PD | 591 † | 571 † |
| 0.001 | 0.24 | 10,000 | 10 | rho | 3 | -0.507 | 0.004 | holds | 0.3661 | 0.258 | 206 † | 202 † |
| 0.001 | 0.24 | 10,000 | 10 | rho_rel | 3 | -0.507 | 0.004 | holds | 0.3661 | 0.258 | 144 † | 141 † |
| 0.001 | 0.24 | 10,000 | 10 | pd | 3 | -0.586 | 0.054 | holds | 2.867 x PD | 1.65 x PD | 506 † | 352 † |
| 0.01 | 0.02 | 100 | 1 | rho | 3 | -0.448 | 0.009 | holds | 0.1741 | 0.02772 | 47 | 48 |
| 0.01 | 0.02 | 100 | 1 | rho_rel | 3 | -0.448 | 0.009 | holds | 0.1741 | 0.02772 | 4,660 † | 8,028 † |
| 0.01 | 0.02 | 100 | 1 | pd | 3 | -0.492 | 0.007 | holds | 1.074 x PD | 0.3769 x PD | 71 | 72 |
| 0.01 | 0.02 | 1,000 | 10 | rho | 3 | -0.501 | 0.023 | holds | 0.04793 | 0.02772 | 4 ‡ | 4 ‡ |
| 0.01 | 0.02 | 1,000 | 10 | rho_rel | 3 | -0.501 | 0.023 | holds | 0.04793 | 0.02772 | 354 † | 352 † |
| 0.01 | 0.02 | 1,000 | 10 | pd | 3 | -0.457 | 0.016 | holds | 0.5048 x PD | 0.3769 x PD | 16 ‡ | 15 ‡ |
| 0.01 | 0.02 | 10,000 | 100 | rho | 3 | -0.500 | 0.004 | holds | 0.03018 | 0.02772 | 2 ‡ | 2 ‡ |
| 0.01 | 0.02 | 10,000 | 100 | rho_rel | 3 | -0.500 | 0.004 | holds | 0.03018 | 0.02772 | 140 † | 140 † |
| 0.01 | 0.02 | 10,000 | 100 | pd | 3 | -0.528 | 0.013 | holds | 0.4038 x PD | 0.3769 x PD | 11 ‡ | 11 ‡ |
| 0.01 | 0.12 | 100 | 1 | rho | 3 | -0.504 | 0.012 | holds | 0.3821 | 0.1493 | 225 † | 222 † |
| 0.01 | 0.12 | 100 | 1 | rho_rel | 3 | -0.504 | 0.012 | holds | 0.3821 | 0.1493 | 624 † | 610 † |
| 0.01 | 0.12 | 100 | 1 | pd | 3 | -0.484 | 0.004 | holds | 1.473 x PD | 0.9233 x PD | 134 † | 139 † |
| 0.01 | 0.12 | 1,000 | 10 | rho | 3 | -0.528 | 0.014 | holds | 0.1912 | 0.1493 | 57 | 56 |
| 0.01 | 0.12 | 1,000 | 10 | rho_rel | 3 | -0.528 | 0.014 | holds | 0.1912 | 0.1493 | 157 † | 146 † |
| 0.01 | 0.12 | 1,000 | 10 | pd | 3 | -0.498 | 0.002 | holds | 1.142 x PD | 0.9233 x PD | 81 | 81 |
| 0.01 | 0.12 | 10,000 | 100 | rho | 3 | -0.502 | 0.004 | holds | 0.1531 | 0.1493 | 36 | 37 |
| 0.01 | 0.12 | 10,000 | 100 | rho_rel | 3 | -0.502 | 0.004 | holds | 0.1531 | 0.1493 | 100 | 100 |
| 0.01 | 0.12 | 10,000 | 100 | pd | 3 | -0.507 | 0.008 | holds | 1.069 x PD | 0.9233 x PD | 71 | 70 |
| 0.01 | 0.24 | 100 | 1 | rho | 3 | -0.505 | 0.019 | holds | 0.5742 | 0.258 | 507 † | 496 † |
| 0.01 | 0.24 | 100 | 1 | rho_rel | 3 | -0.505 | 0.019 | holds | 0.5742 | 0.258 | 352 † | 346 † |
| 0.01 | 0.24 | 100 | 1 | pd | 3 | -0.505 | 0.012 | holds | 2.014 x PD | 1.306 x PD | 250 † | 245 † |
| 0.01 | 0.24 | 1,000 | 10 | rho | 3 | -0.491 | 0.026 | holds | 0.3457 | 0.258 | 184 † | 189 † |
| 0.01 | 0.24 | 1,000 | 10 | rho_rel | 3 | -0.491 | 0.026 | holds | 0.3457 | 0.258 | 128 † | 131 † |
| 0.01 | 0.24 | 1,000 | 10 | pd | 3 | -0.530 | 0.030 | holds | 1.805 x PD | 1.306 x PD | 201 † | 184 † |
| 0.01 | 0.24 | 10,000 | 100 | rho | 3 | -0.476 | 0.021 | holds | 0.2774 | 0.258 | 119 † | 125 † |
| 0.01 | 0.24 | 10,000 | 100 | rho_rel | 3 | -0.476 | 0.021 | holds | 0.2774 | 0.258 | 83 | 85 |
| 0.01 | 0.24 | 10,000 | 100 | pd | 3 | -0.555 | 0.013 | holds | 1.748 x PD | 1.306 x PD | 188 † | 163 † |
| 0.05 | 0.02 | 100 | 5 | rho | 3 | -0.453 | 0.003 | holds | 0.08323 | 0.02772 | 11 ‡ | 10 ‡ |
| 0.05 | 0.02 | 100 | 5 | rho_rel | 3 | -0.453 | 0.003 | holds | 0.08323 | 0.02772 | 1,065 † | 1,480 † |
| 0.05 | 0.02 | 100 | 5 | pd | 3 | -0.466 | 0.016 | holds | 0.5516 x PD | 0.2917 x PD | 19 ‡ | 18 ‡ |
| 0.05 | 0.02 | 1,000 | 50 | rho | 3 | -0.509 | 0.002 | holds | 0.03442 | 0.02772 | 2 ‡ | 2 ‡ |
| 0.05 | 0.02 | 1,000 | 50 | rho_rel | 3 | -0.509 | 0.002 | holds | 0.03442 | 0.02772 | 183 † | 178 † |
| 0.05 | 0.02 | 1,000 | 50 | pd | 3 | -0.501 | 0.022 | holds | 0.3343 x PD | 0.2917 x PD | 7 ‡ | 7 ‡ |
| 0.05 | 0.02 | 10,000 | 500 | rho | 3 | -0.465 | 0.025 | holds | 0.02837 | 0.02772 | 2 ‡ | 1 ‡ |
| 0.05 | 0.02 | 10,000 | 500 | rho_rel | 3 | -0.465 | 0.025 | holds | 0.02837 | 0.02772 | 124 † | 135 † |
| 0.05 | 0.02 | 10,000 | 500 | pd | 3 | -0.484 | 0.030 | holds | 0.2978 x PD | 0.2917 x PD | 6 ‡ | 6 ‡ |
| 0.05 | 0.12 | 100 | 5 | rho | 3 | -0.505 | 0.010 | holds | 0.2237 | 0.1493 | 77 | 77 |
| 0.05 | 0.12 | 100 | 5 | rho_rel | 3 | -0.505 | 0.010 | holds | 0.2237 | 0.1493 | 214 † | 210 † |
| 0.05 | 0.12 | 100 | 5 | pd | 3 | -0.483 | 0.009 | holds | 0.8817 x PD | 0.7145 x PD | 48 | 48 |
| 0.05 | 0.12 | 1,000 | 50 | rho | 3 | -0.495 | 0.011 | holds | 0.1607 | 0.1493 | 40 | 40 |
| 0.05 | 0.12 | 1,000 | 50 | rho_rel | 3 | -0.495 | 0.011 | holds | 0.1607 | 0.1493 | 111 † | 112 † |
| 0.05 | 0.12 | 1,000 | 50 | pd | 3 | -0.499 | 0.024 | holds | 0.7836 x PD | 0.7145 x PD | 38 | 38 |
| 0.05 | 0.12 | 10,000 | 500 | rho | 3 | -0.488 | 0.018 | holds | 0.149 | 0.1493 | 35 | 34 |
| 0.05 | 0.12 | 10,000 | 500 | rho_rel | 3 | -0.488 | 0.018 | holds | 0.149 | 0.1493 | 95 | 97 |
| 0.05 | 0.12 | 10,000 | 500 | pd | 3 | -0.464 | 0.015 | holds | 0.7864 x PD | 0.7145 x PD | 39 | 38 |
| 0.05 | 0.24 | 100 | 5 | rho | 3 | -0.509 | 0.022 | holds | 0.3738 | 0.258 | 215 † | 209 † |
| 0.05 | 0.24 | 100 | 5 | rho_rel | 3 | -0.509 | 0.022 | holds | 0.3738 | 0.258 | 150 † | 146 † |
| 0.05 | 0.24 | 100 | 5 | pd | 3 | -0.542 | 0.013 | holds | 1.25 x PD | 1.011 x PD | 97 | 91 |
| 0.05 | 0.24 | 1,000 | 50 | rho | 3 | -0.503 | 0.025 | holds | 0.2716 | 0.258 | 114 † | 113 † |
| 0.05 | 0.24 | 1,000 | 50 | rho_rel | 3 | -0.503 | 0.025 | holds | 0.2716 | 0.258 | 79 | 79 |
| 0.05 | 0.24 | 1,000 | 50 | pd | 3 | -0.486 | 0.003 | holds | 1.18 x PD | 1.011 x PD | 86 | 88 |
| 0.05 | 0.24 | 10,000 | 500 | rho | 3 | -0.504 | 0.025 | holds | 0.2668 | 0.258 | 110 † | 109 † |
| 0.05 | 0.24 | 10,000 | 500 | rho_rel | 3 | -0.504 | 0.025 | holds | 0.2668 | 0.258 | 76 | 76 |
| 0.05 | 0.24 | 10,000 | 500 | pd | 3 | -0.509 | 0.027 | holds | 1.152 x PD | 1.011 x PD | 82 | 81 |
