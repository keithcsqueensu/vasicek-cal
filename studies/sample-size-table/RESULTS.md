# S-2 Sample-size planning table: results

Run on 2026-09-27, after the prediction ([`PREDICTION.md`](PREDICTION.md)) was merged, with
`python studies/sample-size-table/sample_size_table.py` (stdlib only) on the committed recovery
summary (`tests/golden/recovery/summary.csv`, 81 scenarios × R = 1,000). No new fits; it runs in
under a second.

**How to read the grids:**
- **T\*** is the number of years (periods, n obligors each) at which 1.96 × RMSE meets the target,
  from RMSE ≈ C/√T fitted over T ∈ {20, 40, 100}.
- **"(extrap.)"** marks T\* outside the studied 20–100 years, so it rests on the T^(−½) law beyond
  the data. The law held in every cell where it could be checked (24 of 27, slope −0.40 to −0.59).
- **"(unchecked)"** marks a cell where only T = 100 was usable: at T = 20 and 40, more than half the
  estimates are on the grid edge.
- **"not estimable"** marks a cell where no T was usable.

### Years for ρ within ±0.05 (95%)

| PD | n | n PD | rho = 0.02 | rho = 0.12 | rho = 0.24 |
|---|---|---|---|---|---|
| 0.001 | 100 | 0.1 | not estimable | 1,803 (extrap.) (unchecked) | 3,371 (extrap.) (unchecked) |
| 0.001 | 1,000 | 1 | 24 | 181 (extrap.) | 537 (extrap.) |
| 0.001 | 10,000 | 10 | 3 (extrap.) | 56 | 206 (extrap.) |
| 0.01 | 100 | 1 | 47 | 225 (extrap.) | 507 (extrap.) |
| 0.01 | 1,000 | 10 | 4 (extrap.) | 57 | 184 (extrap.) |
| 0.01 | 10,000 | 100 | 2 (extrap.) | 36 | 119 (extrap.) |
| 0.05 | 100 | 5 | 11 (extrap.) | 77 | 215 (extrap.) |
| 0.05 | 1,000 | 50 | 2 (extrap.) | 40 | 114 (extrap.) |
| 0.05 | 10,000 | 500 | 2 (extrap.) | 35 | 110 (extrap.) |

### Years for PD within ±25% (95%)

| PD | n | n PD | rho = 0.02 | rho = 0.12 | rho = 0.24 |
|---|---|---|---|---|---|
| 0.001 | 100 | 0.1 | not estimable | 804 (extrap.) (unchecked) | 1,187 (extrap.) (unchecked) |
| 0.001 | 1,000 | 1 | 81 | 224 (extrap.) | 591 (extrap.) |
| 0.001 | 10,000 | 10 | 21 | 149 (extrap.) | 506 (extrap.) |
| 0.01 | 100 | 1 | 71 | 134 (extrap.) | 250 (extrap.) |
| 0.01 | 1,000 | 10 | 16 (extrap.) | 81 | 201 (extrap.) |
| 0.01 | 10,000 | 100 | 11 (extrap.) | 71 | 188 (extrap.) |
| 0.05 | 100 | 5 | 19 (extrap.) | 48 | 97 |
| 0.05 | 1,000 | 50 | 7 (extrap.) | 38 | 86 |
| 0.05 | 10,000 | 500 | 6 (extrap.) | 39 | 82 |

## What the table says

- **ρ within ±0.05 is a long-history target** once ρ is moderate:
  - at ρ = 0.12 it takes about 35–40 years with 50 or more expected defaults a year (n·PD ≥ 50),
    about 56 with 10, 77 with 5, and 180–225 with 1;
  - at ρ = 0.24 it is beyond 100 years everywhere studied (110–120 even for large portfolios).

  So a typical 20–30 year history does not pin ρ to ±0.05 unless ρ is small. Where n·PD ≥ 50, the
  large-n formula T\* ≈ (1.96·√2·ρ(1 − ρ)/0.05)² is close: C within 8%, so T\* within about 16%.
- **At ρ = 0.02 the ±0.05 target is lenient,** because the estimate cannot go below 0.001, so an
  estimate is never far from the truth on the low side. The small T\* there (2–4 years with
  n·PD ≥ 10) says the error is small in absolute terms, not that ρ is known relatively: the interval
  still reaches down to the floor.
- **PD within ±25%:**
  - PD depends on ρ as strongly as ρ itself does: at 5% PD it takes about 6–7 years at ρ = 0.02,
    38–39 at 0.12 and 82–86 at 0.24 (large portfolios);
  - at ρ = 0.12, 1% PD takes 71–81 years and 0.1% PD about 150.

  The factor cycle, not the number of obligors, limits PD's precision once n·PD ≥ 10.
- **The finite-n penalty.** The fitted C against the large-n benchmark shows how much binomial noise
  costs. For ρ at 0.12 and 0.24, C is:
  - 1.00–1.08× the benchmark at n·PD ≥ 50;
  - 1.27–1.42× at n·PD = 10, so 1.6–2.0× the years;
  - 2.2–2.6× at n·PD = 1, so 5–7× the years.

  For PD, C is 1.02–1.34× the benchmark at n·PD ≥ 50, 1.22–1.74× at n·PD = 10, and 1.54–2.85× at
  n·PD = 1.

## Prediction against result

| # | Prediction | Result | Held |
|---|---|---|---|
| Q1 | Scaling holds for PD in ≥ 19 of 27 cells, for ρ in ≥ 16; failures mostly at ρ = 0.02 or n·PD ≤ 1 | holds in 24 of 27 for both; no cell fails; the other 3 are unchecked (2) or not estimable (1), all at PD 0.1%, n = 100 | yes |
| Q2 | Not estimable: at most 3 cells, all PD 0.1% with n = 100 | 1 cell: (0.1%, 0.02, 100) | yes |
| Q3 | Monotone in n, ρ and PD in ≥ 90% of adjacent pairs, each direction | ρ vs n 17/17, ρ vs ρ 17/17, PD vs PD 17/17, PD vs n 16/17 (94%; (5%, 0.12): 37.7 at n = 1,000 against 38.0 at 10,000, a tie within noise) | yes |
| Q4 | ρ = 0.02: T\* ≤ 20 in ≥ 7 of 9 cells | 6 of 9. Misses: (0.1%, 1,000) at 24, (1%, 100) at 47, and the not-estimable cell | **no**, by 1 cell |
| Q5 | ρ = 0.24: T\* > 60 in all 9, > 100 in ≥ 6 | > 100 in all 9 (110–3,371) | yes |
| Q6 | ρ = 0.12, n = 10,000, PD ≥ 1%: T\* in [25, 70]; at n = 100, T\* > 100 in all 3 | 36 and 35: held. n = 100: 1,803, 225 and **77** at PD 5% | **no**: (5%, 0.12, 100) needs 77 years, not > 100 |
| Q7 | Large-n benchmark within ±35%: PD in ≥ 5 of 6, ρ in ≥ 3 of 4 | PD 6 of 6 (ratios 1.02–1.34); ρ 4 of 4 (1.00–1.08); the two ρ = 0.02 cells, reported only, at 1.09 and 1.02 | yes |
| Q8 | PD ± 25%, n = 10,000: (5%, 0.12) in [22, 50]; (1%, 0.12) in [36, 85]; (1%, 0.24) in [75, 170]; (0.1%, 0.24) > 100; (5%, 0.02) and (1%, 0.02) ≤ 20 | 39, 71, **188**, 506, 6, 11 | **no**: (1%, 0.24) needs 188, above the range by 18 (the fitted slope, −0.56, gives 163) |
| Q9 | PD at n·PD ≤ 1: T\* > 100 in every estimable cell | > 100 in 6 of 8; **81** at (0.1%, 0.02, 1,000) and **71** at (1%, 0.02, 100) | **no** |
| Q10 | No cell with ρ ≥ 0.12 reaches both targets within 40 years, apart from at most two cells with PD 5% and n ≥ 1,000 | exactly those two: (5%, 0.12, 1,000) at 40/38 and (5%, 0.12, 10,000) at 35/39 | yes |

**Six held, four missed.** The misses, not explained away:

- **Q4, Q6 and Q9 all underestimate how well the estimators do with few defaults when ρ is low or
  PD is 5%.** At ρ = 0.02 there is almost no factor variance, so the PD estimate is close to the
  pooled binomial one: relative SD ≈ 1/√(n·PD·T), which gives T\* ≈ 61/(n·PD) years. That is 61 at
  n·PD = 1, against the observed 71–81 and a prediction of > 100. The prediction treated n·PD ≤ 1 as
  uniformly hopeless, which holds at ρ ≥ 0.12 but not at ρ = 0.02. The same effect puts
  (1%, 0.02, 100)'s ρ at 47 years rather than ≤ 20, and (5%, 0.12, 100) at 77 rather than > 100.
- **Q8 misses in size at (1%, 0.24, 10,000):** 188 against an upper limit of 170. The prediction
  allowed up to 1.6× the benchmark (105), and the observed factor is 1.8×. The slope there
  is −0.56: RMSE at T = 20 carries extra skew, which the fixed −½ law turns into a conservative T\*.
- **No surprise of the kind the prediction named:** no slope is steeper than −0.59, let alone −0.75.

## Full output

`sample_size_table.py`, verbatim: the per-cell fits (C for PD is shown relative to PD) and the
large-n benchmark.

| PD | rho | n | n PD | param | points used | slope | max resid | scaling | C | large-n C | T* (slope -1/2) | T* (fitted slope) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.001 | 0.02 | 100 | 0.1 | rho | 0 | | | | | | not estimable at any T studied | |
| 0.001 | 0.02 | 100 | 0.1 | pd | 0 | | | | | | not estimable at any T studied | |
| 0.001 | 0.02 | 1,000 | 1 | rho | 3 | -0.396 | 0.019 | holds | 0.1241 | 0.02772 | 24 | 21 |
| 0.001 | 0.02 | 1,000 | 1 | pd | 3 | -0.491 | 0.018 | holds | 1.143 x PD | 0.4762 x PD | 81 | 82 |
| 0.001 | 0.02 | 10,000 | 10 | rho | 3 | -0.498 | 0.018 | holds | 0.04204 | 0.02772 | 3 (extrap.) | 3 (extrap.) |
| 0.001 | 0.02 | 10,000 | 10 | pd | 3 | -0.485 | 0.014 | holds | 0.5828 x PD | 0.4762 x PD | 21 | 21 |
| 0.001 | 0.12 | 100 | 0.1 | rho | 1 |  |  | unchecked | 1.083 | 0.1493 | 1,803 (extrap.) | - |
| 0.001 | 0.12 | 100 | 0.1 | pd | 1 |  |  | unchecked | 3.615 x PD | 1.166 x PD | 804 (extrap.) | - |
| 0.001 | 0.12 | 1,000 | 1 | rho | 3 | -0.504 | 0.015 | holds | 0.3429 | 0.1493 | 181 (extrap.) | 179 (extrap.) |
| 0.001 | 0.12 | 1,000 | 1 | pd | 3 | -0.577 | 0.036 | holds | 1.905 x PD | 1.166 x PD | 224 (extrap.) | 180 (extrap.) |
| 0.001 | 0.12 | 10,000 | 10 | rho | 3 | -0.511 | 0.017 | holds | 0.1903 | 0.1493 | 56 | 56 |
| 0.001 | 0.12 | 10,000 | 10 | pd | 3 | -0.508 | 0.054 | holds | 1.556 x PD | 1.166 x PD | 149 (extrap.) | 146 (extrap.) |
| 0.001 | 0.24 | 100 | 0.1 | rho | 1 |  |  | unchecked | 1.481 | 0.258 | 3,371 (extrap.) | - |
| 0.001 | 0.24 | 100 | 0.1 | pd | 1 |  |  | unchecked | 4.393 x PD | 1.65 x PD | 1,187 (extrap.) | - |
| 0.001 | 0.24 | 1,000 | 1 | rho | 3 | -0.508 | 0.032 | holds | 0.5908 | 0.258 | 537 (extrap.) | 517 (extrap.) |
| 0.001 | 0.24 | 1,000 | 1 | pd | 3 | -0.507 | 0.001 | holds | 3.099 x PD | 1.65 x PD | 591 (extrap.) | 571 (extrap.) |
| 0.001 | 0.24 | 10,000 | 10 | rho | 3 | -0.507 | 0.004 | holds | 0.3661 | 0.258 | 206 (extrap.) | 202 (extrap.) |
| 0.001 | 0.24 | 10,000 | 10 | pd | 3 | -0.586 | 0.054 | holds | 2.867 x PD | 1.65 x PD | 506 (extrap.) | 352 (extrap.) |
| 0.01 | 0.02 | 100 | 1 | rho | 3 | -0.448 | 0.009 | holds | 0.1741 | 0.02772 | 47 | 48 |
| 0.01 | 0.02 | 100 | 1 | pd | 3 | -0.492 | 0.007 | holds | 1.074 x PD | 0.3769 x PD | 71 | 72 |
| 0.01 | 0.02 | 1,000 | 10 | rho | 3 | -0.501 | 0.023 | holds | 0.04793 | 0.02772 | 4 (extrap.) | 4 (extrap.) |
| 0.01 | 0.02 | 1,000 | 10 | pd | 3 | -0.457 | 0.016 | holds | 0.5048 x PD | 0.3769 x PD | 16 (extrap.) | 15 (extrap.) |
| 0.01 | 0.02 | 10,000 | 100 | rho | 3 | -0.500 | 0.004 | holds | 0.03018 | 0.02772 | 2 (extrap.) | 2 (extrap.) |
| 0.01 | 0.02 | 10,000 | 100 | pd | 3 | -0.528 | 0.013 | holds | 0.4038 x PD | 0.3769 x PD | 11 (extrap.) | 11 (extrap.) |
| 0.01 | 0.12 | 100 | 1 | rho | 3 | -0.504 | 0.012 | holds | 0.3821 | 0.1493 | 225 (extrap.) | 222 (extrap.) |
| 0.01 | 0.12 | 100 | 1 | pd | 3 | -0.484 | 0.004 | holds | 1.473 x PD | 0.9233 x PD | 134 (extrap.) | 139 (extrap.) |
| 0.01 | 0.12 | 1,000 | 10 | rho | 3 | -0.528 | 0.014 | holds | 0.1912 | 0.1493 | 57 | 56 |
| 0.01 | 0.12 | 1,000 | 10 | pd | 3 | -0.498 | 0.002 | holds | 1.142 x PD | 0.9233 x PD | 81 | 81 |
| 0.01 | 0.12 | 10,000 | 100 | rho | 3 | -0.502 | 0.004 | holds | 0.1531 | 0.1493 | 36 | 37 |
| 0.01 | 0.12 | 10,000 | 100 | pd | 3 | -0.507 | 0.008 | holds | 1.069 x PD | 0.9233 x PD | 71 | 70 |
| 0.01 | 0.24 | 100 | 1 | rho | 3 | -0.505 | 0.019 | holds | 0.5742 | 0.258 | 507 (extrap.) | 496 (extrap.) |
| 0.01 | 0.24 | 100 | 1 | pd | 3 | -0.505 | 0.012 | holds | 2.014 x PD | 1.306 x PD | 250 (extrap.) | 245 (extrap.) |
| 0.01 | 0.24 | 1,000 | 10 | rho | 3 | -0.491 | 0.026 | holds | 0.3457 | 0.258 | 184 (extrap.) | 189 (extrap.) |
| 0.01 | 0.24 | 1,000 | 10 | pd | 3 | -0.530 | 0.030 | holds | 1.805 x PD | 1.306 x PD | 201 (extrap.) | 184 (extrap.) |
| 0.01 | 0.24 | 10,000 | 100 | rho | 3 | -0.476 | 0.021 | holds | 0.2774 | 0.258 | 119 (extrap.) | 125 (extrap.) |
| 0.01 | 0.24 | 10,000 | 100 | pd | 3 | -0.555 | 0.013 | holds | 1.748 x PD | 1.306 x PD | 188 (extrap.) | 163 (extrap.) |
| 0.05 | 0.02 | 100 | 5 | rho | 3 | -0.453 | 0.003 | holds | 0.08323 | 0.02772 | 11 (extrap.) | 10 (extrap.) |
| 0.05 | 0.02 | 100 | 5 | pd | 3 | -0.466 | 0.016 | holds | 0.5516 x PD | 0.2917 x PD | 19 (extrap.) | 18 (extrap.) |
| 0.05 | 0.02 | 1,000 | 50 | rho | 3 | -0.509 | 0.002 | holds | 0.03442 | 0.02772 | 2 (extrap.) | 2 (extrap.) |
| 0.05 | 0.02 | 1,000 | 50 | pd | 3 | -0.501 | 0.022 | holds | 0.3343 x PD | 0.2917 x PD | 7 (extrap.) | 7 (extrap.) |
| 0.05 | 0.02 | 10,000 | 500 | rho | 3 | -0.465 | 0.025 | holds | 0.02837 | 0.02772 | 2 (extrap.) | 1 (extrap.) |
| 0.05 | 0.02 | 10,000 | 500 | pd | 3 | -0.484 | 0.030 | holds | 0.2978 x PD | 0.2917 x PD | 6 (extrap.) | 6 (extrap.) |
| 0.05 | 0.12 | 100 | 5 | rho | 3 | -0.505 | 0.010 | holds | 0.2237 | 0.1493 | 77 | 77 |
| 0.05 | 0.12 | 100 | 5 | pd | 3 | -0.483 | 0.009 | holds | 0.8817 x PD | 0.7145 x PD | 48 | 48 |
| 0.05 | 0.12 | 1,000 | 50 | rho | 3 | -0.495 | 0.011 | holds | 0.1607 | 0.1493 | 40 | 40 |
| 0.05 | 0.12 | 1,000 | 50 | pd | 3 | -0.499 | 0.024 | holds | 0.7836 x PD | 0.7145 x PD | 38 | 38 |
| 0.05 | 0.12 | 10,000 | 500 | rho | 3 | -0.488 | 0.018 | holds | 0.149 | 0.1493 | 35 | 34 |
| 0.05 | 0.12 | 10,000 | 500 | pd | 3 | -0.464 | 0.015 | holds | 0.7864 x PD | 0.7145 x PD | 39 | 38 |
| 0.05 | 0.24 | 100 | 5 | rho | 3 | -0.509 | 0.022 | holds | 0.3738 | 0.258 | 215 (extrap.) | 209 (extrap.) |
| 0.05 | 0.24 | 100 | 5 | pd | 3 | -0.542 | 0.013 | holds | 1.25 x PD | 1.011 x PD | 97 | 91 |
| 0.05 | 0.24 | 1,000 | 50 | rho | 3 | -0.503 | 0.025 | holds | 0.2716 | 0.258 | 114 (extrap.) | 113 (extrap.) |
| 0.05 | 0.24 | 1,000 | 50 | pd | 3 | -0.486 | 0.003 | holds | 1.18 x PD | 1.011 x PD | 86 | 88 |
| 0.05 | 0.24 | 10,000 | 500 | rho | 3 | -0.504 | 0.025 | holds | 0.2668 | 0.258 | 110 (extrap.) | 109 (extrap.) |
| 0.05 | 0.24 | 10,000 | 500 | pd | 3 | -0.509 | 0.027 | holds | 1.152 x PD | 1.011 x PD | 82 | 81 |
