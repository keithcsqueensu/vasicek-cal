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

## Addendum, 2026-09-27: a bias-corrected Wald interval for q (S-23's mitigation 3a)

Appended after S-23 was pinned (D-156) and before any jackknife run, under this file's rule for
additions. Nothing above changes. **Evidence used:** S-23's pinned results
(`tests/golden/recovery/summary.csv`, the q columns) and its diagnosis of the delta-method Wald
findings for q.

### What S-23 found, and why this arm

S-23's delta-method Wald interval for the 99.9% conditional PD q is below the band in 17 of the 50
assessed scenarios. The diagnosis is **low estimates with narrow intervals**:

- q̂ is biased low, by −0.13 to −0.36 SE on average, carrying ρ̂'s downward bias;
- the delta-method SE moves with the estimate, so the lowest estimates get the narrowest intervals:
  across replicates, the error in logit(q̂) and its SE correlate at +0.64 to +0.88;
- removing the mean bias exactly (an oracle shift of the pinned estimates) lifts coverage only to
  0.916–0.947. That would bring 15 of the 17 into the band, all but 29 (0.916) and 64 (0.924).

So both the centre and the SE matter. The arm below tests each with the jackknife this run
computes anyway, at no extra cost.

### Definitions

Both sub-arms use s = logit(q), in the coordinate of S-23's delta-method interval.

- **Centre:** s̃ = logit q(PD̂, ρ̃), with ρ̃ S-3's corrected estimate (set to the bound where it
  leaves the box) and PD̂ the parity estimate.
- **(a) delta-method SE:** s̃ ± 1.96·SE_Δ, with SE_Δ S-23's delta-method SE of logit(q̂) at the
  estimate.
- **(b) jackknife SE:** s̃ ± 1.96·SE_J, with SE_J² = ((T − 1)/T)·Σₜ(s₍₋ₜ₎ − s̄)² and
  s₍₋ₜ₎ = logit q(PD̂₍₋ₜ₎, ρ̂₍₋ₜ₎), from the same delete-one estimates.
- **Coverage:** on the replicates without a Wald flag, the set S-23's Wald interval uses, so that the
  three are compared pairwise. The verdict follows S-23's Wald rule. (b) is also reported over all
  replicates, since it needs no Hessian.
- **Reported:** coverage and verdict per scenario, and misses below and above. The mechanism is
  also reported directly: the correlation across replicates between the error s̃ − s_true and each
  SE.

### Predictions

| # | Prediction |
|---|---|
| J15 | **(a) helps, partially:** among S-23's 17 below-band Wald scenarios, (a)'s coverage exceeds S-23's Wald coverage in at least 14, and (a) is in the band in 5–12 of them. That is fewer than the oracle's 15: the jackknife removes only part of the bias, adds variance, and keeps the delta-method SE's link to the estimate |
| J16 | **(b) covers at least as well as (a):** in at least 14 of the 17, and (b) is in the band in at least as many of the 17 as (a). If (a) and (b) come out the same, the simpler (a) is preferred |
| J17 | **The mechanism, tested directly:** the error–SE correlation is lower for SE_J than for SE_Δ in at least 12 of the 17 |
| J18 | **What it costs elsewhere:** among S-23's 32 Wald PASS scenarios, at most 3 leave the band under (a) and at most 3 under (b). (b) is above the band in at most 5 of the 50 assessed scenarios (the jackknife variance tends to be too large) |

## Addendum, 2026-09-28: a polished arm, after exploring the study subset

Appended after the subset exploration of the registered arms (R = 1,000 on scenarios 29, 37, 72,
4, 68, 49, 7, 43, 51) and before the polished arm below has run on it. Nothing above changes: the
registered arms run on the full matrix exactly as defined, and J1–J18 are scored as written.

### What the subset showed

- **Delete-one estimates are not exact.** The jackknife estimates come from the grid refinement,
  and the resolution check (J14's refits) found gaps of up to 0.52 SE where a delete-one fit's
  refinement is rejected (the estimate is then a grid point), and 0.12 SE where it is accepted, at
  scenario 29 (ρ = 0.02, n = 10⁴). Elsewhere on the subset the gap is at most 0.05 SE.
- **S-3 amplifies them.** ρ̃ multiplies the delete-one mean by T − 1. The shifted interval's
  coverage fell well below the parity profile's at 29 (0.873 against 0.927), where the resolution
  error is largest, and at 4 and 7 (group A), where ρ̃ was set to the box's floor in 36% and 21% of
  replicates.
- So the registered S-3 arm measures the grid as well as the jackknife. The arm below separates
  them.

**Seen before writing this:** the subset summary of the registered arms, and a 4-replicate smoke
test of the polished arm on scenario 29, run to check the code (shifted interval covering 4 of 4,
ρ̃'s RMSE 0.0051 on those 4). Nothing else of the polished arm.

### The polished arm (subset only)

- **Estimates:** ρ̂_p and PD̂_p are the full panel's exact off-grid maximum (the profile code's
  polished maximum). ρ̂_p₍₋ₜ₎ and PD̂_p₍₋ₜ₎ are each delete-one panel's, from the same code on
  that panel. ρ̃_p = T·ρ̂_p − (T − 1)·mean ρ̂_p₍₋ₜ₎, set to a bound of the box where it leaves it.
- **S-3's shifted interval:** the parity profile interval shifted in logit ρ by
  logit ρ̃_p − logit ρ̂_p, truncated ends staying at the box, as for the registered arm.
- **J15–J18's sub-arms:** (a) and (b) centred at logit q(PD̂_p, ρ̃_p); (b)'s jackknife SE from the
  polished delete-one estimates.
- **Cost:** one profile per delete-one panel, about 3 h on the subset. The full matrix is not run.

### Predictions

The six informative subset scenarios are 29, 37, 68, 49 (T = 20, 40) and 43, 51 (T = 100).

| # | Prediction |
|---|---|
| J19 | **Resolution explains 29:** the polished arm's shifted-interval coverage at scenario 29 is at least 0.91 (registered arm 0.873), and ρ̃_p's RMSE there is at least 20% below the registered ρ̃'s (0.0099) |
| J20 | **Where resolution is small, polishing changes little:** in 37, 68, 49, 43 and 51, the polished and registered arms' shifted-interval coverages differ by at most 0.01 |
| J21 | **Clamping is not a resolution effect:** in 4 and 7 the polished arm's shifted-interval coverage is still at least 0.03 below the parity profile's (0.986 and 0.964) |
| J22 | **With exact estimates, S-3's shift is harmless in the informative scenarios:** in at least 5 of the 6, the polished arm's shifted-interval coverage is within 0.015 of the parity profile's |
| J23 | **q, sub-arm (a):** in each of the six, the polished (a)'s coverage is at least the registered (a)'s minus 0.005, and it is higher at 29 |
