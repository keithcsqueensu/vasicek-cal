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
- **Posterior:** S-9's estimator, with its resolution rule (at least 4 grid points per posterior SD
  on each axis, refining locally); and, as a diagnostic, the parity grid with the rule switched
  off. Under the rule, a draw's rank is taken on the grid the estimator ends on, with the prior's
  mass on that grid.
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
  spacing of 0.130. So on the parity grid without the rule (1,000, 20) should be calibrated and
  (1,000, 100) should not; the resolution rule should restore calibration everywhere.

## Predictions

| # | Prediction |
|---|---|
| G1 | **(1,000, 20), with the rule, flat prior:** for PD and for ρ, the χ² p-value is at least 0.01 and the tail share is within 0.027–0.073 |
| G2 | **(1,000, 20), with the rule, Jeffreys prior:** the same holds |
| G3 | **(1,000, 100), parity grid without the rule:** for ρ and for PD the tail share is below 0.027, a hump |
| G4 | **(1,000, 100), with the rule:** for PD and for ρ, the χ² p-value is at least 0.01 and the tail share is within 0.027–0.073 |
| G5 | **(10⁴, 40), parity grid without the rule:** the tail share for ρ lies between those of (1,000, 20) and (1,000, 100) without the rule |

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not. The rows are computed by a script
committed before the run's results exist. Results, the comparison, a mitigation and a monitoring
implication go in `studies/README.md` (the S-15 entry) and in a D-entry, with a `MANIFEST.json`.

## Addendum, 2026-09-28: how the run and its scoring are computed, fixed before the run

Appended before the run. No prediction or threshold changes; this fixes what the text leaves open.
It is implemented in `study_bayes_sbc` (`tests/studies/bayes_sbc.cpp`) and [`compare.py`](compare.py),
committed with this addendum.

- **The draw's random stream** (the text fixes the panel's stream, not the parameters'): Philox4x32-10
  keyed by the seed "S15SBCDR" (0x5331355342434452), counter (setting, draw, prior, block), words to
  doubles as the DGP does. Block 0's words 0–1 pick the cell by inverting the prior's cumulative cell
  masses in grid order; block 0's words 2–3 and block 1's words 0–1 place the point uniformly in the
  cell on the PD and ρ axes. The prior's cell masses are the posterior code's own (`masses` with a zero
  log-likelihood), so the draws follow the grid's cell-uniform prior exactly. The flat and Jeffreys
  draws at (1,000, 20) share the panel stream (scenario 1000 + setting, replicate = draw), with
  different parameters.
- **The rank:** `engine::marginal_cdf` on the final marginal (new, D-173), at the true value in logit.
  A draw the estimator refuses has no rank: it is left out of that arm's histogram and tail share,
  and refusals are reported. The diagnostic arm is `grid_posterior` with `resolution_rule = false`.
- **Scoring:** the 20-bin χ² test has 19 degrees of freedom. The tail-share range 0.027–0.073 is
  inclusive; "below 0.027" is strict. G5 uses the flat prior (the only prior run at (10⁴, 40) and
  (1,000, 100)) and holds when (10⁴, 40) lies strictly between the other two, in either order.
- **Disclosure:** the tool was run on 4 draws per setting and prior to check the code; not scored.
