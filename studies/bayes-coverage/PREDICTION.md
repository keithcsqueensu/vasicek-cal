# S-9 Frequentist coverage of grid-Bayesian credible intervals: predictions made before the run

Pre-registration for study S-9 (`studies/README.md`, D-148). This file is committed **before** the
grid-Bayesian estimator exists (M3) and before any posterior has been computed on a recovery panel.
The results will be compared with it prediction by prediction, and misses will be reported, not
explained away.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a
new, dated section appended at the end, committed before the run it concerns. M3 implements the
estimator as defined below; if M3 has to define it differently, that change is appended here,
dated, before S-9 runs.

**Evidence used in writing it:** results already pinned: the groups of D-136, the profile coverage
of PD and ρ (D-121, D-131), S-13's finding that T = 20 profile intervals undercover by about one
point with most misses below the truth (D-158), and each scenario's spread of logit(estimate)
across replicates against the parity grid's spacing ([`grid_spread.py`](grid_spread.py)).

## What is run

- **Panels:** the recovery panels, unchanged. The study subset first, then the full matrix
  (81 × 1,000). The model is correctly specified, so the targets are the true PD and ρ.
- **The posterior,** on the grid, from the per-period surfaces L the MLE already computes:
  - the log-posterior at grid point k is log π(θ_k) + Σ_t L[t, k] + log J_k, where J_k is the cell's
    area on the natural scale (the grid is logit-spaced, so a density on (PD, ρ) is weighted by
    ∏ v(1 − v)·Δu), normalised by `LogSumExpPosterior`;
  - **priors:** (a) **flat** on the natural scale over the box; (b) **Jeffreys**, √det I(PD, ρ),
    with I the Fisher information of one binomial-mixture period with the scenario's n, computed
    once per n at every grid point by the parity rule (T periods scale I by T, which cancels);
  - **marginals:** each axis's posterior mass at its grid points, spread uniformly in the logit
    coordinate over each point's cell (half-cells at the box's ends), so the marginal CDF is
    piecewise linear and continuous.
- **The resolution rule (part of the estimator, for M3 to implement):** the grid must have **at
  least k = 4 points per posterior SD along each axis**, in the logit coordinate. The parity grid
  locates the posterior; the estimator then checks the rule and, where it fails, refines:
  1. a local grid in (logit PD, logit ρ), spacing at most SD/4 per axis, spanning ±8 SD around the
     posterior mode (clipped to the box, and at least ±2 parity spacings), where SD starts from the
     MLE's Hessian SE in logit units (D-119) and the parity-grid posterior's marginal SD, whichever
     is smaller;
  2. the posterior is recomputed on the local grid (the per-period contributions are evaluated at
     its points), its marginal SDs are measured again, and the step repeats while the rule fails;
  3. after three refinements that still fail the rule, or if the parity-grid posterior puts more
     than 10⁻⁶ of its mass outside the local grid, the estimator refuses and flags the fit.

  Why 4: for a normal posterior, cell-uniform marginals with 4 points per SD put the 95% interval's
  ends within 0.01 SD of the exact ones and its true mass within 0.001 of 0.95 at any grid offset
  (2 points per SD: 0.04 SD and 0.005; 1 point: 0.16 SD and 0.016). The rule scales the grid with
  the data instead of fixing one resolution.
- **Intervals, 95%:** equal-tailed (the 2.5% and 97.5% points of the marginal CDF) and HPD (the
  shortest set of cells, with the partial cell at each end, on the logit scale). Coverage is over
  all replicates; a fit the estimator refuses (the resolution rule) counts as not covering, as in
  D-131, and refusals are reported.
- **Two arms, both on the full matrix:** the estimator with the resolution rule; and, as a
  diagnostic, the parity grid (61 × 41) with the rule switched off.
- **Reported:** coverage and verdict (D-131 band) per scenario, parameter, prior and interval type,
  pairwise against the pinned profile interval; the side of misses; median widths in logit.
- **One fitting run pins the results** (D-156's process change); out-of-band verdicts are reviewed
  in a separate labels file, as in D-155.

## Mechanism behind the predictions

- **Bernstein–von Mises:** where the data are informative the posterior is close to normal around
  the MLE with the inverse-information covariance, so equal-tailed credible intervals behave like
  the profile interval. Where the data are thin, the prior matters.
- **The grid is too coarse in places.** In groups B–D the estimator's spread across replicates is
  0.25–2.99 grid spacings for PD and 0.82–3.12 for ρ. Where it is under half a spacing (PD in 8,
  32, 34, 35, 58, 59, 61, 62: large n·PD at T = 40 and 100), almost all posterior mass sits on one
  or two grid points; spreading it uniformly over their cells makes the interval about a cell wide,
  several posterior SDs, so it **overcovers**. The resolution rule removes this. Since the spread is
  under 4 spacings everywhere in groups B–D (at most 3.12), the rule refines almost every fit there:
  the parity grid only locates the posterior.
- **The priors pull differently at small T.** Near ρ = 0 Jeffreys behaves like 1/ρ and pulls ρ
  down; the flat prior puts more mass on larger ρ than the logit grid's own spacing does and pulls
  ρ up. ρ̂ is biased down at T = 20 and the profile interval's misses there are mostly below the
  truth (S-13), so the flat prior should cover ρ better at T = 20 than Jeffreys.

## Predictions

Groups as in D-136: A (31), B (12), C (17), D (21). A verdict is one scenario and one parameter.

| # | Prediction |
|---|---|
| F1 | **Without the rule the grid overcovers where it is coarse:** in the diagnostic arm (parity grid, rule off), the equal-tailed PD interval is above the band under both priors in at least 6 of the 8 coarse scenarios (8, 32, 34, 35, 58, 59, 61, 62) |
| F2 | **The rule refines almost everywhere and rarely refuses:** in groups B–D the estimator refines in at least 95% of replicates, and refuses in at most 0.1% of replicates in every scenario |
| F3 | **With the rule, credible intervals match the profile:** in groups B–D, the Jeffreys equal-tailed coverage is within 0.015 of the profile coverage in at least 80 of the 100 verdicts (50 scenarios × PD, ρ) |
| F4 | **The flat prior helps ρ at T = 20:** with the rule, in group B, the flat prior's equal-tailed ρ coverage is at least Jeffreys' in at least 9 of the 12 scenarios |
| F5 | **HPD against equal-tailed:** with the rule, in groups B–D, the HPD interval's median logit width is below the equal-tailed interval's for ρ in at least 45 of the 50 scenarios, and their ρ coverages differ by at most 0.02 in at least 40 of the 50 |

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not. A prediction that holds in direction
but misses in size is recorded as not held, with the size of the miss. The rows are computed by a
script committed before the full-matrix run's results exist. Results, the comparison, a mitigation
and a monitoring implication go in `studies/README.md` (the S-9 entry) and in a D-entry, with a
`MANIFEST.json`. The pinned verdicts are unchanged.

## Addendum, 2026-09-28: the resolution rule as implemented (M3.3), before S-9 runs

Appended under this file's rule for changes M3 has to make, before any S-9 run; nothing above is
edited, and F1–F5 are unchanged. **Evidence:** the estimator's unit tests and its reference panels,
which use a validation seed of their own (M3BAYVAL), not the recovery panels.

As registered, the local grid's **extent** was ±8 SD with SD the smaller of the Hessian SE and the
posterior SD. At small n and T the ρ posterior is wide and skewed towards the box's floor, so the
Hessian SE at the mode is far smaller than the posterior SD, and a ±8 SE grid dropped real posterior
mass: 3 of 16 reference fits (T = 20, n = 100) were refused. That would have made F2 fail for a
reason unrelated to the question. As implemented (`engine/posterior.hpp`, D-166):

- **Extent:** the hull of the posterior mean ± 8 posterior SDs and the edges of the cells holding the
  current marginal's 5·10⁻⁸ and 1 − 5·10⁻⁸ quantiles, at least ±2 current spacings, clipped to the
  box. Centred on the posterior **mean**, not the mode.
- **Spacing:** at most s/5, s the smaller of the posterior SD and (on the first refinement) the
  Hessian SE: 5 points per SD rather than 4, because the cell-uniform SD shrinks slightly with the
  cells, and a grid built at exactly 4 would fail the rule it was built for.
- **Refusal:** after three local grids that still fail the rule, if a local grid would need more than
  2,001 points on an axis, or if the current grid's mass outside the next grid exceeds 10⁻⁶.
- **HPD:** the shortest interval holding 95% of the cell-uniform marginal (a partial cell at each end),
  as registered; it can lie up to half a cell from the smooth HPD (`TOL_POSTERIOR_HPD_END_SD`).

The rule itself (at least 4 grid points per posterior SD on each axis) is unchanged.

## Addendum, 2026-09-28: how the run and its scoring are computed, fixed before the run

Appended before any S-9 run on the recovery panels. No prediction, threshold or estimator definition
changes; this fixes what the text leaves open, so the scoring cannot be chosen after the results are
seen. It is implemented in [`compare.py`](compare.py), committed with this addendum.

- **The tool:** `study_bayes_coverage` (`tests/studies/bayes_coverage.cpp`) fits each panel as the
  recovery harness does (calibrate on the parity grid through the row cache; the pinned profile
  intervals) and computes `engine::grid_posterior` in four arms: flat or Jeffreys, with the rule or
  without it. The Jeffreys table is computed once per n on the parity grid. One row per (scenario,
  replicate); [`export_fits.py`](export_fits.py) writes `fits.parquet`, and `compare.py` scores F1–F5
  from that file and the D-136 groups.
- **The diagnostic arm:** `grid_posterior` gained `resolution_rule = false`, which returns the
  posterior on the parity grid as it stands, never refined or refused (D-173). With the rule on, the
  estimator is unchanged: its reference results reproduce with a worst difference of 0.
- **Coverage and verdicts:** over all replicates; a refused fit or one without a finite posterior
  does not cover (D-131). The band is 0.95 ± 3.29·√(0.95·0.05/1,000). The profile coverage in F3 is
  the pinned one.
- **F2** is scored for each prior: refinement pooled over the B–D replicates, refusals per B–D
  scenario (at most 1 of 1,000). Held if both priors hold. **F5** is scored for each prior, with the
  median logit width over the replicates that have an interval; held if both priors hold.
- **Disclosure:** the tool was run on 4 replicates of scenarios 1, 29, 43 and 68 and 2 of 1, 8 and 29,
  to check the code and time it (about 5 CPU seconds per panel with the rule, both priors); those
  rows were not scored.
