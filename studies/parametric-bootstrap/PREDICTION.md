# S-10 Parametric bootstrap intervals, against profile and Jeffreys: predictions made before the run

Pre-registration for study S-10 (`studies/README.md`, D-148, D-168). This file is committed **before**
any parametric-bootstrap panel has been simulated. The results will be compared with it prediction by
prediction, and misses will be reported, not explained away. S-4 (the Bartlett-corrected profile
threshold, [`../bartlett-profile/PREDICTION.md`](../bartlett-profile/PREDICTION.md)) shares this run:
its usable factor (S-4b) is computed from the same bootstrap panels.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a new,
dated section appended at the end, committed before the run it concerns.

**Evidence used in writing it:** pinned results only, through [`reference.py`](reference.py): the
recovery summary (profile and iid percentile coverage, the bias and spread of logit ρ̂), S-13's
R = 10,000 profile coverage at T = 20 (D-158), and S-9's per-replicate rows (Jeffreys equal-tailed
coverage on the same panels, D-175). No parametric-bootstrap panel has been simulated.

## The question, and why it is asked this way

The parametric bootstrap is the one remaining candidate interval method (S-23's mitigation 3) and the
usual remedy for small-sample profile intervals. It is not run as a standalone study: the question is
whether it beats the two intervals that already work. So every panel gets all three methods in the
same run, and the predictions are about coverage **relative to** the profile interval and the Jeffreys
equal-tailed interval on the same replicates, and about width, since an interval can cover better by
being wider, which costs something when the estimate feeds a stress figure.

The comparators on the recovery panels (`reference.py`, pooled coverage per group):

| group | PD: profile, Jeffreys | ρ: profile, Jeffreys |
|---|---|---|
| A (31, near-uninformative) | 0.9548, 0.9517 | 0.9645, 0.9722 |
| B (12, T = 20) | 0.9401, 0.9422 | 0.9371, 0.9430 |
| C (17, T = 40) | 0.9468, 0.9489 | 0.9448, 0.9482 |
| D (21, T = 100) | 0.9478, 0.9489 | 0.9479, 0.9495 |

In group B the Jeffreys interval covers more than the profile interval on the same replicates (ρ: 109
replicates covered by Jeffreys only, 38 by the profile only; paired z = 5.9) at essentially the same
width (median width ratio 1.003), but it still undercovers ρ by about 0.7 points.

## What is run

- **Panels:** the recovery panels, unchanged (81 scenarios × replicates 0–999). The study subset first,
  to measure the cost; then the full matrix. One fitting run pins the results (D-156).
- **Per panel:** the parity fit (θ̂, its Hessian SEs); the pinned profile interval; the Jeffreys
  equal-tailed interval with the resolution rule (see "Comparators" below); and the parametric
  bootstrap:
  - **B = 999 bootstrap panels** simulated at θ̂ with the recovery DGP (D-109) and the panel's own n and
    T, under a separate seed (`S10PBOOT`, 0x53313050424F4F54), scenario id as the panel's, replicate =
    1,000·r + b for bootstrap panel b of replicate r. Each is fitted as the recovery panels are (the
    parity grid through the row cache, D-167; refinement; Hessian SEs).
  - **Percentile interval,** per parameter in logit: the 2.5% and 97.5% points (type 7) of the 999
    bootstrap estimates; grid-edge estimates kept at their grid value, non-finite ones left out, as
    for the iid bootstrap (D-136).
  - **Studentised (bootstrap-t) interval,** per parameter on the logit scale (u = logit PD, logit ρ):
    t*_b = (u*_b − û)/se*_b, and the interval is [û − t*₀.₉₇₅·se, û − t*₀.₀₂₅·se] mapped back to the
    natural scale. **The SE, for the panel (se) and for every bootstrap panel (se*_b) alike, is the
    analytic one:** the observed information from the analytic score and Hessian of the panel
    log-likelihood (D-170's posterior-moment derivatives, `panel_derivs`) at the refined estimate, in
    the logit coordinates, inverted as a 2 × 2 matrix, the square root of its diagonal entry for the
    parameter. This is not the SE `calibrate` reports (D-119: central differences of the objective);
    the two agree closely away from the bounds, and the run reports their largest relative difference
    on the panels. Bootstrap fits where the information is not positive definite are left out; if more
    than 10% are, the interval is not computed (and does not cover, D-131).
- **Comparators in the same rows:** the profile interval is recomputed (it reproduces the pinned ends,
  checked). The Jeffreys interval (with the resolution rule, as S-9 ran it) costs about 5 CPU seconds
  per panel, so it is **not** recomputed in full: the run recomputes it on replicates 0–49 of every
  scenario (4,050 panels) and requires them to equal S-9's committed rows bit for bit; the other
  replicates are joined from S-9's rows (`studies/bayes-coverage/fits.parquet`, the same panels and
  engine). A single mismatch stops the comparison until it is explained.
- **Reported, per scenario, parameter and method:** coverage and its class against the D-131 band;
  misses below and above the truth; median logit width, and its ratio to the profile interval's on
  the same replicates; paired counts against the profile and Jeffreys intervals (covered by one and
  not the other).

## Mechanism behind the predictions

- **The percentile interval inherits the bias twice.** The bootstrap distribution of θ̂* is centred
  near θ̂ + b, so its percentile interval sits about 2b from the truth. At T = 20 logit ρ̂ is biased by
  −0.07 to −0.23 of its spread (`reference.py`), which alone predicts 0.925–0.948 (mean about 0.935).
  The profile interval's own shortfall at T = 20 (skew and a variable curvature) comes on top, so the
  percentile interval should cover ρ **less** than the profile interval at T = 20. It will cover far
  more than the iid percentile interval (0.80–0.89 in group B, D-137), whose shortfall came from
  resampling T periods, which the parametric bootstrap does not do.
- **The studentised interval corrects location and scale together.** Pivoting on t* removes the bias
  and the dependence of the spread on the parameter to second order, so it should reach about 95% at
  T = 20, at a cost in width: t*'s quantiles are wider than ±1.96 when the SE is noisy, as at T = 20.
- **At T ≥ 40 and away from the bounds all three methods converge,** and the comparison is about width.
- **Near the bounds (group A)** percentile and studentised intervals break down as the iid bootstrap
  did (bootstrap estimates pile on the box's bound; SEs are not finite): reported, not predicted.

## Predictions

Groups as in D-136. "Pooled" coverage is over all replicates of the group's scenarios.

| # | Prediction |
|---|---|
| K1 | **Percentile below profile for ρ at T = 20:** in group B, pooled percentile ρ coverage is 0.910–0.935 and below the profile's (0.9371); below the profile's in at least 9 of the 12 scenarios |
| K2 | **Studentised at least as good as Jeffreys for ρ at T = 20, but wider:** in group B, pooled studentised ρ coverage (logit ρ, analytic Hessian SE, as defined above) is at least Jeffreys' (0.9430) and at most 0.960; its median width ratio to the profile interval is 1.03–1.20 in at least 10 of 12 scenarios |
| K3 | **Convergence at T ≥ 40:** in groups C and D (38 scenarios), for PD and ρ, the percentile and studentised pooled coverages are each within 0.010 of the profile's, and their median width ratios to the profile interval lie in 0.95–1.08 in at least 32 of 38 scenarios per method and parameter |
| K4 | **Far better than the iid bootstrap:** in groups B–D, the parametric percentile ρ coverage exceeds the pinned iid percentile coverage in at least 45 of 50 scenarios |
| K5 | **Neither bootstrap interval beats Jeffreys at equal width:** for each of the percentile and studentised intervals, it is not the case that it *beats* Jeffreys and has *equal width* (terms defined below) |

**Terms in K5, fixed before the run:**

- **Beats:** on ρ in group B (12 scenarios × 1,000 replicates, paired on the same panels), the paired
  z-statistic of the bootstrap interval against Jeffreys is above 2: z = (n₁₀ − n₀₁)/√(n₁₀ + n₀₁), with
  n₁₀ the replicates the bootstrap interval covers and Jeffreys does not, and n₀₁ the reverse.
- **Equal width:** the median over the 12 group B scenarios of the per-scenario ratio of median logit ρ
  widths (bootstrap over Jeffreys) lies within 0.98–1.02.

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not; a prediction that holds in direction but
misses in size is not held. The rows are computed by a script committed before the full run's results
exist. The results update the interval comparison in `studies/README.md` (profile, Jeffreys,
BCa, and the two parametric intervals, with coverage and width side by side) and give the recommendation
for which interval to use when; with a D-entry and a `MANIFEST.json`.

## Addendum, 2026-09-29: how the run and its scoring are computed, fixed before the subset run

Appended before the subset run; no prediction, threshold or definition above changes. It fixes what
the text leaves open, and is implemented in `study_parametric_bootstrap`
(`tests/studies/parametric_bootstrap.cpp`) and [`compare.py`](compare.py), committed with it.

- **Engine additions for the run** (tested before it): `engine::profile_lr_statistic` returns W at
  given values from the same polished maximum and inner maximisation as the profile intervals (W = 2c
  at the solved ends to 2.3e-9); `engine::analytic_se_scaled` is the analytic-observed-information SE
  (it agrees with calibrate's D-119 SE to 6.5e-5 relative, `TOL_ANALYTIC_SE_VS_CENTRAL_REL`).
- **Checks before scoring:** the Jeffreys interval on replicates 0–49 must equal S-9's rows bit for
  bit, and the recomputed profile ends S-9's rows; either failing stops the scoring. W0 ≤ c against the
  profile interval's coverage is counted and reported.
- **Pooling:** "pooled" coverage is over all replicates of the group's scenarios; K3 and H4 pool groups
  C and D together (38 scenarios), per method and parameter.
- **S-4a** covers when W0 ≤ c·k_a, k_a the mean of the scenario's finite W0; **S-4b** when W0 ≤ c·k_b.
  Their widths use replicates 0–199, where the corrected intervals are solved; a width ratio is the
  median of the method's logit widths over the median of the profile's on the same replicates.
- **The subset run** (the study subset of D-150, 9 scenarios × 1,000) is scored with the same script;
  every figure is printed with its denominators, and "held" reads "subset" when a prediction's
  registered denominator is incomplete. It is reported as such, not as the study's result.
- **Disclosure:** the tool was smoke-run on 2 replicates of scenarios 29 and 43 to check the code and the
  scoring pipeline (the checks passed); those rows are not results.
