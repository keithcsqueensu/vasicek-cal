# The shared jackknife run (S-3, S-5, S-21): predictions made before the run

Pre-registration for studies S-3 (`jackknife-bias-rho`), S-5 (`bca-intervals`) and S-21
(`period-influence`), which share one jackknife pass (D-150). This file is committed **before**
any jackknife estimate, BCa interval or influence has been computed on a recovery panel. The
results will be compared with it prediction by prediction, and misses will be reported, not
explained away. It lives in `studies/jackknife-bias-rho/`; the S-5 and S-21 entries point here.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a
new, dated section appended at the end, committed before the run it concerns.

**Evidence used in writing it:** only the pinned recovery summary (M1.8, M2a, M2b; D-121, D-131,
D-136, D-137).

## What is run

- **Panels:** the recovery panels and fits, unchanged (seeds, replicates, box, grid, parity
  integrator and check). The study subset first (scenarios 29, 37, 72, 4, 68, 49, 7, 43, 51,
  R = 1,000), then the full matrix (81 × 1,000) to pin verdicts. The model is correctly specified:
  the targets are the true PD and ρ.
- **The jackknife pass:** for each panel, the T delete-one-period estimates θ̂₍₋ₜ₎ from W × L
  (M2b's jackknife weights, `resample::replicate_estimates`): the 3 × 3 refinement on the
  reweighted surface, exactly what `calibrate` gives on the reduced panel. Grid-edge replicates are
  kept at their grid value and counted. The leave-two-out pass uses the T(T − 1)/2 pairs the same
  way (190, 780 and 4,950 rows at T = 20, 40, 100).
- **One fitting run pins the results:** the full-matrix run writes its outputs directly, and the
  reviewed verdict labels for the new verdict families are kept in a separate file, so that a
  review never requires a refit.
- **Exact refits (the resolution check):** for replicates 0–4 of every scenario, each
  leave-one-out panel is also maximised off the grid (the polished maximum of the profile code),
  and the difference from the refined estimate is reported in SE units.

### S-3: the bias-corrected ρ̂ (a `native` option)

- **Estimate:** ρ̃ = T·ρ̂ − (T − 1)·mean_t ρ̂₍₋ₜ₎, on the natural scale as stated in the index
  entry, with ρ̂ the refined estimate (the same method as the jackknife estimates, so the
  correction measures bias rather than a method difference). ρ̃ outside the box is set to the
  bound and counted.
- **Interval:** the parity profile interval for ρ, shifted in logit ρ by logit ρ̃ − logit ρ̂. An
  end truncated at a bound of the box stays at the bound.
- **Reported:** the bias and RMSE of ρ̃ against ρ̂'s, over all replicates; coverage of the shifted
  interval over all replicates, with its verdict under the D-131 band and policy; pairwise, the
  replicates it gains and loses against the parity interval.

### S-5: BCa intervals (a `native` option)

- **In the logit coordinate**, for PD and for ρ, from the existing B = 999 iid bootstrap
  replicates (D-136) and the jackknife above:
  - z₀ = Φ⁻¹((#{u*_b < û} + ½·#{u*_b = û}) / B′), B′ the finite replicates;
  - a = Σ(ū₍·₎ − u₍₋ₜ₎)³ / (6·(Σ(ū₍·₎ − u₍₋ₜ₎)²)^{3/2}), over the jackknife estimates;
  - α₁,₂ = Φ(z₀ + (z₀ ± z₀.₉₇₅)/(1 − a(z₀ ± z₀.₉₇₅))), the ends being type-7 quantiles of the
    replicate estimates at α₁ and α₂.
- **Not computed** when z₀ is infinite (every replicate on one side of û), the jackknife has no
  spread, or a denominator is not positive. Such an interval does not cover (D-131).
- **Reported:** coverage over all replicates and its verdict, pairwise against the percentile
  interval, with the number not computed.

### S-21: period influence (descriptive; no verdicts)

- **Leave-one-out influence:** Δ_t = (ρ̂₍₋ₜ₎ − ρ̂)/SE(ρ̂), with the Hessian SE (D-119), on the
  replicates without a Wald flag; the same for PD. Per replicate: the largest |Δ_t|, its period,
  that period's realised factor Z_t (from the DGP) and its default rate.
- **Leave-two-out:** the largest |Δ_{s,t}| over pairs, and its ratio to the largest |Δ_t|.
- **Reported per scenario:** medians and 90th percentiles of those maxima, and the share of
  replicates whose most influential period is the one with the most extreme Z_t, adverse or
  benign.

## Mechanism behind the predictions

- **ρ̂'s bias is small against its spread** (pinned): bias/SD is −0.10 to −0.25 in group B,
  −0.07 to −0.18 in C, and −0.04 to −0.13 in D. Removing a bias of 0.2 SD moves the coverage of a
  well-calibrated interval by about half a point, so **S-3 cannot move coverage much**. The
  jackknife removes the O(1/T) term of the bias, and costs variance: its RMSE should be no
  better at T = 20.
- **The percentile interval's shortfall for ρ is mostly not bias** (pinned): group B ρ coverage is
  0.803–0.887, far more than a 0.2 SD bias explains. D-137's diagnosis is skewness plus a
  resampling variance of (T − 1)/T. BCa corrects median bias (z₀) and skewness (a), not the
  variance shortfall, and it cannot help where the bootstrap is inconsistent (the boundary
  breakdown of group A).
- **Influence:** ρ is a variance-like parameter, so a period's influence grows roughly with
  (Z_t² − 1): the most influential period should usually be the one with the most extreme factor.
  In SE units, the largest of T influences is about max(Z²−1)/√(2T): about 0.6 at T = 20, 0.55 at
  T = 40, 0.5 at T = 100.

## Predictions, full matrix

Groups as in D-136: A (31), B (12), C (17), D (21).

### S-3

| # | Prediction |
|---|---|
| J1 | **Bias:** in groups B–D (50 scenarios), |bias(ρ̃)| < |bias(ρ̂)| in at least 40 |
| J2 | **RMSE:** RMSE(ρ̃) > RMSE(ρ̂) in at least 9 of the 12 group B scenarios; within ±3% of it in at least 15 of the 21 group D scenarios |
| J3 | **Coverage:** in groups B–D, the shifted interval's coverage differs from the parity profile's by at most 0.015 in every scenario, and is higher in at least 30 of 50 |
| J4 | **The ρ small-T findings (29, 55, 74 in B; 72 in A):** at most 2 of the 4 move into the band |
| J5 | **Broken PASS verdicts:** at most 2 ρ profile PASS verdicts in groups B–D fall out of the band with the shift |
| J6 | **Clamping:** ρ̃ is set to a bound in at least 5% of replicates in at least 10 of the 31 group A scenarios, and in under 1% of replicates in every group C and D scenario |

### S-5

| # | Prediction |
|---|---|
| J7 | **ρ, groups B–D:** BCa coverage exceeds percentile coverage in at least 45 of 50; BCa is still below the band in at least 8 of the 12 group B scenarios; in group D, at least 14 of 21 in the band |
| J8 | **PD, groups B–D:** BCa is in the band in at least 35 of 50 (percentile today: 26) |
| J9 | **Group A (boundary breakdown, 57 verdicts):** BCa does not fix it; at least 45 of the 57 stay below the band, and BCa is not computed in at least 2% of replicates in at least 10 scenarios |
| J10 | **Totals over 162:** BCa below the band 60–100 (percentile: 125); above the band at most 5 |

### S-21

| # | Prediction |
|---|---|
| J11 | **Size:** in groups B–D, the median over replicates of the largest |Δ_t| for ρ is 0.45–0.85 SE at T = 20, 0.4–0.75 at T = 40 and 0.3–0.65 at T = 100, and falls with T for every (PD, ρ, n) cell outside group A |
| J12 | **Which period:** in groups B–D, the most influential period for ρ is the one with the most extreme Z_t in at least 60% of replicates in every scenario |
| J13 | **Pairs:** the median ratio of the largest leave-two-out |Δ| to the largest leave-one-out is 1.3–1.9 in every group B–D scenario |
| J14 | **Resolution:** the refined leave-one-out estimates differ from the exact off-grid maxima by less than 0.05 SE in every checked replicate outside group A |

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not. A prediction that holds in direction
but misses in size is recorded as not held, with the size of the miss. Results and the comparison
go in `studies/README.md` (the S-3, S-5 and S-21 entries) and in a D-entry. The pinned parity
verdicts are unchanged: S-3 and S-5 are `native` options reported beside them.
