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
- **Intervals, 95%:** equal-tailed (the 2.5% and 97.5% points of the marginal CDF) and HPD (the
  shortest set of cells, with the partial cell at each end, on the logit scale). Coverage is over
  all replicates; every posterior interval is computed, so none is DEFERRED.
- **Two grids:** the parity grid (61 × 41) on the full matrix; a **fine grid** (241 × 161, the same
  box and scale, 4 × finer per axis) on the study subset only.
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
  several posterior SDs, so it **overcovers**. The fine grid removes this.
- **The priors pull differently at small T.** Near ρ = 0 Jeffreys behaves like 1/ρ and pulls ρ
  down; the flat prior puts more mass on larger ρ than the logit grid's own spacing does and pulls
  ρ up. ρ̂ is biased down at T = 20 and the profile interval's misses there are mostly below the
  truth (S-13), so the flat prior should cover ρ better at T = 20 than Jeffreys.

## Predictions

Groups as in D-136: A (31), B (12), C (17), D (21). A verdict is one scenario and one parameter.

| # | Prediction |
|---|---|
| F1 | **The grid overcovers where it is coarse:** on the parity grid, the equal-tailed PD interval is above the band under both priors in at least 6 of the 8 coarse scenarios (8, 32, 34, 35, 58, 59, 61, 62) |
| F2 | **Where the grid is adequate it matches the profile:** on the parity grid, in the verdicts whose spread is at least 1.5 spacings (16 for PD, 25 for ρ; `grid_spread.py`), the Jeffreys equal-tailed coverage is within 0.015 of the profile coverage in at least 80% of the 41 |
| F3 | **The fine grid matches the profile:** on the fine grid, in the six informative subset scenarios (29, 37, 68, 49, 43, 51), the Jeffreys equal-tailed coverage for PD and for ρ is within 0.02 of the profile coverage in at least 10 of the 12 verdicts |
| F4 | **The flat prior helps ρ at T = 20:** on the parity grid, in group B, the flat prior's equal-tailed ρ coverage is at least Jeffreys' in at least 9 of the 12 scenarios |
| F5 | **HPD against equal-tailed:** in groups B–D on the parity grid, the HPD interval's median logit width is below the equal-tailed interval's for ρ in at least 45 of the 50 scenarios, and their ρ coverages differ by at most 0.02 in at least 40 of the 50 |

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not. A prediction that holds in direction
but misses in size is recorded as not held, with the size of the miss. The rows are computed by a
script committed before the full-matrix run's results exist. Results, the comparison, a mitigation
and a monitoring implication go in `studies/README.md` (the S-9 entry) and in a D-entry, with a
`MANIFEST.json`. The pinned verdicts are unchanged.
