# Grid-Bayesian estimation (M3)

This note describes, in plain prose, the grid-Bayesian estimator of PD and ρ (`engine/posterior.hpp`):
the posterior on the estimation grid, its priors, the resolution rule that keeps the grid fine enough,
and its credible intervals. It was defined in study S-9's registration (D-161) and amended there before
S-9 ran (addendum of 2026-09-28; D-166). The script `validation/scipy/grid_posterior.py` replicates it
(section 6).

## 1. The posterior on the grid

Everything is in the grid's logit coordinates u = (logit PD, logit ρ), where grid points are equally
spaced. With the per-period log-likelihoods l_t the engine already computes (the surface L), the
log posterior density at grid point k is

    log π_u(u_k) + Σ_t l_t(θ_k),

π_u being the prior's density in u. Each point's mass is its density times its cell's area (cells are
centred on the points, half-cells at a grid's ends), and the masses are normalised. Marginals are the
masses summed over the other axis, each point's mass spread uniformly over its cell (the cell-uniform
convention), so a marginal's CDF is piecewise linear and continuous: quantiles need no tie-breaking.
Sums run in a fixed order, so results do not depend on the thread count.

## 2. Priors

- **Flat** on the natural scale over the box: π_u = PD(1 − PD)·ρ(1 − ρ), the Jacobian of the logit.
- **Jeffreys:** √det I_u(PD, ρ), I_u the Fisher information, in u, of one binomial-mixture period
  with the panel's n (T periods multiply it by T, which cancels). I_u = Σ_d P(d)·g_d g_dᵀ, g_d the
  score of log P(d) by central differences (step 10⁻⁴ in u), over the counts with P(d) ≥ 10⁻¹⁴ × the
  largest; log P(d) is the binomial-mixture l_t, by the parity integrator. It is tabulated once per n
  on the parity grid and interpolated bilinearly in log at other points. Against scipy's adaptive quad
  the table agrees within 1.45·10⁻⁹ in log.

## 3. The resolution rule

A grid posterior is only as fine as its grid: where the posterior is narrower than a few cells, its
quantiles are set by the grid, not the data. On the parity grid (61 × 41) the estimator's spread in
groups B–D is 0.25–2.99 spacings for PD (`studies/bayes-coverage/grid_spread.py`), so a fixed grid
cannot serve. **The rule: at least 4 grid points per posterior SD along each axis.** For a normal
posterior, cell-uniform marginals at 4 points per SD put the 95% interval's ends within 0.010 SD of the
exact ones and its true mass within 0.001 of 0.95 at any grid offset (2 points: 0.04 SD and 0.005).

The parity grid locates the posterior. While the rule fails, the posterior is recomputed on a local
grid, the per-period log-likelihoods evaluated at its points:

- **extent,** per axis: the hull of the posterior mean ± 8 posterior SDs and the edges of the cells
  holding the marginal's 5·10⁻⁸ and 1 − 5·10⁻⁸ quantiles, at least ±2 current spacings, within the
  box. The tail term matters for skewed posteriors, such as ρ's towards the box's floor at small n;
- **spacing:** at most s/5, s the smaller of the posterior SD and, on the first refinement, the MLE's
  Hessian SE in u. Five rather than four points per SD, because the cell-uniform SD shrinks slightly
  with the cells;
- **refusal,** flagged, never an approximate answer: after three local grids that still fail, if a
  local grid would need more than 2,001 points on an axis, or if more than 10⁻⁶ of the current grid's
  mass lies outside the next grid.

## 4. Intervals

95% intervals per parameter, on the natural scale: **equal-tailed**, the marginal CDF's 2.5% and 97.5%
points; and **HPD**, the shortest interval holding 95% of the cell-uniform marginal, with a partial cell
at each end (found exactly: its length is piecewise linear in its lower end, so the minimum is where
an end crosses a cell edge). Near its minimum that length is almost flat, so candidates up to half a
cell apart can differ by rounding alone; the choice is made canonical, so that a platform's libm
cannot move it: candidates within 10⁻⁹ (relative) of the shortest count as shortest, and among them
the one centred closest to the marginal's median is taken. On an exactly normal posterior the equal-tailed ends are within 0.009 SD
and the HPD ends within 0.071 SD of the exact ones (the HPD of a step density can sit up to half a cell
from the smooth HPD while its width is right).

## 5. Validation in this step

- A posterior known exactly (a test objective whose likelihood is normal in u), with SDs from 0.02 to
  0.8 in u, refined and unrefined: ends, mean and SD as in section 4 (`TOL_POSTERIOR_*`).
- A bimodal posterior with a narrow Hessian SE at one mode is kept whole; a posterior piled against the
  box's floor is not refused; a posterior too narrow to resolve in three refinements is refused.
- Identical results for any thread count.
- Reference results (`tests/golden/posterior/`, `posterior_reference`): a Jeffreys table (n = 100) and
  16 fits, both priors, on panels of a validation seed of its own (M3BAYVAL), not the recovery panels,
  so S-9's and S-15's results are first computed in their registered runs. Recomputed on every
  platform by `unit_posterior`.

## 6. Replicating with scipy

`validation/scipy/grid_posterior.py` recomputes the Jeffreys table with scipy's adaptive quad (1.45·10⁻⁹
in log) and re-runs every reference fit with a fixed 200-node Gauss–Hermite rule in numpy: the parity
posterior and both priors, the rule's decision on the parity grid (exactly the engine's), the rule on
the engine's final grid, and both intervals (within 3.0·10⁻⁵ in logit). It runs in the CI job
"validation (scipy)".

## 7. What is not here yet

- Frequentist coverage of these intervals (S-9) and simulation-based calibration (S-15): their
  registered runs.
- An interval for the 99.9% conditional PD from the posterior.
- The C ABI (v0 exposes the binomial MLE only).
