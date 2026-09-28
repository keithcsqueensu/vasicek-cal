# S-15 Simulation-based calibration of the grid-Bayesian estimator: predictions made before the run

Pre-registration for study S-15 (`studies/README.md`, D-148, D-151). This file is committed
**before** the grid-Bayesian estimator exists (M3) and before any posterior has been computed. The
results will be compared with it prediction by prediction, and misses will be reported, not
explained away.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a
new, dated section appended at the end, committed before the run it concerns. The estimator is the
one S-9 defines ([`../bayes-coverage/PREDICTION.md`](../bayes-coverage/PREDICTION.md)); a change
there made before S-9 runs applies here too.

**Evidence used in writing it:** S-9's definitions and its grid-spread reference; the model's
standard errors in logit units (below). No posterior, SBC draw or recovery panel.

## What is run

Simulation-based calibration (Talts et al. 2018): if the posterior is computed exactly under the
prior the parameters are drawn from, the posterior CDF at the true value is uniform on [0, 1].

- **Settings (n, T):** (1,000, 20), (10⁴, 40) and (1,000, 100). N = 1,000 draws each.
- **Draws:** (PD, ρ) from the prior on the parity grid's box: a grid cell is drawn with the prior's
  mass, then a point uniformly within it in the logit coordinate (the same cell-uniform convention
  as S-9's marginals). A panel is simulated with the recovery DGP (D-109) under a separate stream:
  scenario id 1000 + setting, replicate = draw.
- **Priors:** flat on the natural scale (every setting); Jeffreys (setting (1,000, 20) only).
- **Posterior:** as in S-9, on the parity grid; and on S-9's fine grid (241 × 161) at (1,000, 100).
- **Rank statistic:** the marginal posterior CDF at the true value, for PD and for ρ. It is
  continuous under the cell-uniform convention, so there are no ties.
- **Reported:** a 20-bin histogram per parameter, its χ² p-value against uniform, and the **tail
  share**, the fraction of ranks below 0.025 or above 0.975 (0.05 when calibrated; the Monte Carlo
  band at N = 1,000 is 0.05 ± 0.023).

## Mechanism behind the predictions

- **The grid is the only approximation.** Draws come from the grid's own cell-uniform prior, so any
  departure from uniformity comes from evaluating the likelihood at grid points only and spreading
  each point's mass over its cell.
- **Where the posterior is narrower than a cell, SBC fails in a known shape.** Spreading mass over a
  cell makes the posterior too wide, so the true value's rank piles up in the middle: a hump, with a
  tail share below 0.05. In logit units the posterior SD of ρ̂ is about √(2/T), 0.32 at T = 20,
  0.22 at 40 and 0.14 at 100, against a spacing of 0.173. For PD it depends on (PD, ρ); under the
  flat prior most draws have PD above 5%, where the logit SD at T = 100 is about 0.1, against a
  spacing of 0.130. So (1,000, 20) should be calibrated and (1,000, 100) should not, on the parity
  grid; the fine grid, with spacings a quarter as wide, should restore it.

## Predictions

| # | Prediction |
|---|---|
| G1 | **(1,000, 20), parity grid, flat prior:** for PD and for ρ, the χ² p-value is at least 0.01 and the tail share is within 0.027–0.073 |
| G2 | **(1,000, 20), Jeffreys prior:** the same holds |
| G3 | **(1,000, 100), parity grid:** for ρ and for PD the tail share is below 0.027, a hump |
| G4 | **(1,000, 100), fine grid:** for PD and for ρ, the tail share is within 0.027–0.073 |
| G5 | **(10⁴, 40), parity grid:** the tail share for ρ lies between those of (1,000, 20) and (1,000, 100) on the parity grid |

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not. The rows are computed by a script
committed before the run's results exist. Results, the comparison, a mitigation and a monitoring
implication go in `studies/README.md` (the S-15 entry) and in a D-entry, with a `MANIFEST.json`.
