# S-34 Sensitivity to severe new periods: predictions made before the run

Pre-registration for study S-34 (`studies/README.md`, D-152). This file is committed **before**
any S-34 quantity has been computed: no recovery panel has been extended or refitted, on the
subset or anywhere else. The results will be compared with it prediction by prediction, and misses
will be reported, not explained away.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a
new, dated section appended at the end, committed before the run it concerns.

**Evidence used in writing it:**

- the model alone, in the large-n limit: [`large_n_reference.py`](large_n_reference.py) (standard
  library only; it reads no recovery panel), which gives every reference number below;
- the scenario definitions (PD, ρ, T, n) and the arithmetic of each added period's default count;
- results already pinned: the groups of D-136, the profile coverage of PD, ρ and q (D-121, D-156),
  and D-159's lesson that grid-refined estimates are not exact enough to difference.

## What is run

- **Panels:** the recovery panels, unchanged (seeds, replicates, box, grid, parity integrator and
  its check). The study subset first (scenarios 29, 37, 72, 4, 68, 49, 7, 43, 51, R = 1,000), then
  the full matrix (81 × 1,000). The model is correctly specified and the added periods are possible
  draws from it, so the truth does not change.
- **The added periods:** each panel is extended by k = 1 or 2 periods, each with the scenario's n
  obligors, at a fixed adverse factor level z_a:
  - **1-in-100:** z_a = −Φ⁻¹(0.99) = −2.3263; **1-in-1,000:** z_a = −Φ⁻¹(0.999) = −3.0902.
  - In the engine's convention (higher Z = better conditions; `core/model/vasicek.hpp`), the
    period's conditional PD is p_a = Φ((Φ⁻¹(PD) − √ρ·z_a)/√(1 − ρ)), at the scenario's true PD and ρ.
  - Its default count is the median of Binomial(n, p_a): the smallest d with P(D ≤ d) ≥ ½.
  - So the added periods are the same for every replicate of a scenario, and four variants are run
    per replicate: 1 × 1-in-100, 2 × 1-in-100, 1 × 1-in-1,000, 2 × 1-in-1,000.
- **Estimates:** PD̂ and ρ̂ are the **exact off-grid maxima** (the profile code's polished maximum)
  of the original and of each extended panel, and q̂ = q(PD̂, ρ̂), the 99.9% conditional PD of
  S-23. Not the grid refinement: D-159 found its error large enough to matter in differences.
- **Intervals:** the parity profile intervals for PD, ρ and q (S-23's) of each extended panel.
- **Reported per scenario and variant:**
  - **Shift in SE units:** (θ̂_ext − θ̂)/SE(θ̂) for PD and ρ, with the original panel's Hessian SE
    (D-119); for q, (logit q̂_ext − logit q̂)/SE, with S-23's delta-method SE of logit q̂. Medians over
    the replicates without a Wald flag on the original panel.
  - **Relative shift:** medians of θ̂_ext/θ̂ − 1, over all replicates.
  - **Interval end points:** medians of each end's ratio, extended over original.
  - **Exceedance:** the share of replicates whose q̂_ext lies above the original panel's q profile
    interval (its upper end), and the share whose original q̂ lies below the extended interval.
  - **Coverage of the unchanged truth** by the extended panel's profile intervals, over all
    replicates, beside the pinned coverage of the original panel.
  - Fits flagged on an extended panel are counted and reported.
- **One fitting run pins the results** (D-156's process change). S-34 is descriptive: it adds no
  verdict family, so there are no labels to review.

## Mechanism behind the predictions

- **The large-n limit.** As n → ∞ each period's probit default rate x_t = (Φ⁻¹(PD) − √ρ·z_t)/√(1 − ρ)
  is observed exactly: x_t ~ N(m, s²) with s² = ρ/(1 − ρ). The MLE is the mean μ̂ and the 1/T
  variance v̂ of the x_t, with ρ̂ = v̂/(1 + v̂) and Φ⁻¹(q̂) = μ̂ + Φ⁻¹(0.999)·√v̂. A severe period
  adds x_a = m + s·u (u = 2.3263 or 3.0902): it raises μ̂ by about s·u/(T + 1) and v̂ by about
  s²(u² − 1)/(T + 1). In SE units both are ∝ 1/√T.
- **q̂'s shift in SE units does not depend on PD or ρ** in that limit, because q̂ is a quantile of
  the fitted normal for x. ρ̂'s barely does. The reference (medians, 20,000 draws):

  | T | Added | q̂ shift (SE) | ρ̂ shift (SE) | P(q̂_ext above the old upper end) | q coverage, before → after (Wald) |
  |---|---|---|---|---|---|
  | 20 | 1 × 1-in-100 | 0.81 | 0.68 | 0.011 | 0.908 → 0.987 |
  | 20 | 2 × 1-in-100 | 1.46 | 1.20 | 0.206 | 0.901 → 0.953 |
  | 20 | 1 × 1-in-1,000 | 1.40 | 1.28 | 0.154 | 0.908 → 0.964 |
  | 20 | 2 × 1-in-1,000 | 2.45 | 2.24 | 0.759 | 0.901 → 0.718 |
  | 40 | 1 × 1-in-100 | 0.59 | 0.48 | 0.000 | 0.926 → 0.972 |
  | 40 | 2 × 1-in-100 | 1.11 | 0.91 | 0.016 | 0.927 → 0.951 |
  | 40 | 1 × 1-in-1,000 | 1.02 | 0.92 | 0.005 | 0.926 → 0.955 |
  | 40 | 2 × 1-in-1,000 | 1.91 | 1.73 | 0.450 | 0.927 → 0.797 |
  | 100 | 1 × 1-in-100 | 0.38 | 0.31 | 0.000 | 0.939 → 0.961 |
  | 100 | 2 × 1-in-100 | 0.73 | 0.60 | 0.000 | 0.941 → 0.948 |
  | 100 | 1 × 1-in-1,000 | 0.67 | 0.60 | 0.000 | 0.939 → 0.953 |
  | 100 | 2 × 1-in-1,000 | 1.29 | 1.16 | 0.004 | 0.941 → 0.874 |

  As a rule of thumb, one 1-in-100 period moves q̂ by about (3.6–3.8)/√T SE and one 1-in-1,000
  period by about (6.3–6.7)/√T; two periods move it 1.75–1.94 times as far as one.
- **Finite n dilutes it.** With few defaults a period's count is a noisy reading of its factor
  (S-21, D-159: the most extreme-factor period is the most influential in about a third of
  histories at n·PD ≤ 5). The mixture likelihood then pulls the added period's factor towards the
  centre, and the original SE is larger. Both make the shift smaller in SE units where n·PD is
  small; where n·PD ≥ 100 the limit should be close.
- **Why coverage of the unchanged truth can rise.** The added periods are adverse, so they push ρ̂
  and q̂ up, against their downward small-T bias (S-23), and they widen the interval. One severe
  period therefore moves a low estimate towards the truth more often than it pushes a good one
  past it. Only when the push is large against the SE (two 1-in-1,000 periods, small T) does
  coverage fall.
- **Where the "severe" period is not severe.** In scenarios 0, 3 and 6 (group A: PD 0.1%, ρ 0.02,
  n = 100) the median count of the added period is 0 at both severities, no more than an ordinary
  period's. Adding zero-default periods lowers the default rate. In groups B–D the added period
  has at least 5 defaults (1-in-100) and 9 (1-in-1,000).

## Predictions, full matrix

Groups as in D-136: A (31), B (12), C (17), D (21). A scenario-variant is one scenario under one of
the four variants. SE-unit shifts are for q unless stated.

| # | Prediction |
|---|---|
| K1 | **Direction:** in groups B–D, the median shifts of PD̂, ρ̂ and q̂ are positive in every scenario-variant (200 of 200, for each of the three) |
| K2 | **T:** in groups B–D, q̂'s median shift in SE units falls from T = 20 to 40 to 100 in every (PD, ρ, n) cell with all three T outside group A, for each of the four variants |
| K3 | **Two periods against one:** in groups B–D, q̂'s median shift for two periods is 1.55–2.05 times that for one period of the same severity in at least 90 of the 100 scenario-severity pairs |
| K4 | **Size where n·PD ≥ 100** (the 18 B–D scenarios with n = 10⁴ and PD ≥ 1%): q̂'s median shift is within ±25% of the reference table in at least 58 of their 72 scenario-variants (80%) |
| K5 | **Dilution where n·PD ≤ 10** (the 23 B–D scenarios with n·PD ∈ {1, 5, 10}): q̂'s median shift is below the reference in at least 69 of their 92 scenario-variants (75%) |
| K6 | **ρ̂ where n·PD ≥ 100:** ρ̂'s median shift is within ±25% of the reference's ρ̂ column in at least 58 of the 72 scenario-variants |
| K7 | **Exceedance of the original q interval:** q̂_ext lies above the original profile interval's upper end in at most 5% of replicates for 1 × 1-in-100 in every B–D scenario; in at least 40% of replicates for 2 × 1-in-1,000 in at least 8 of the 12 group B scenarios; and in at most 10% of replicates for every variant in every group D scenario |
| K8 | **Coverage of the unchanged truth (q, profile):** after 1 × 1-in-100 it is at least the pinned coverage in at least 40 of the 50 B–D scenarios; after 2 × 1-in-1,000 it is at least 0.05 below the pinned coverage in at least 8 of the 12 group B scenarios, and below it in at least 15 of the 21 group D scenarios |
| K9 | **Relative size at T = 20:** 1 × 1-in-1,000 raises q̂ by a median of at least 15% in at least 10 of the 12 group B scenarios |
| K10 | **The non-severe variants (scenarios 0, 3, 6):** q̂'s median relative shift is at most 0 under all four variants in all three scenarios |

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not. A prediction that holds in direction
but misses in size is recorded as not held, with the size of the miss. The rows are computed by a
script committed before the full-matrix run's results exist. Results, the comparison, a mitigation
and a monitoring implication go in `studies/README.md` (the S-34 entry) and in a D-entry, with a
`MANIFEST.json` giving commits, commands, times and hashes. The pinned verdicts are unchanged.
