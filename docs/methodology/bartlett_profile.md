# The Bartlett-corrected profile-likelihood interval (native option)

This note describes, in plain prose, the Bartlett-corrected profile-likelihood interval
(`resample/bartlett.hpp`; C ABI `vcal_calibrate_bartlett`, 0.4). It is study S-4's S-4b (D-176, D-181),
made a `native` option by D-182 and D-184. It is not the preferred interval. The preferred interval
remains the Jeffreys equal-tailed credible interval with the resolution rule (`grid_bayesian.md`, D-182).

## 1. What it corrects

The 95% profile-likelihood interval for a parameter θ_a (PD or ρ) is the set of values whose profile
log-likelihood P_a(θ_a) is within c/2 of the maximum l_max, where c = 3.841458820694124 is the 0.95
quantile of χ²₁. The threshold rests on the likelihood-ratio statistic

    W_a = 2 (l_max − P_a(θ_a))

being χ²₁ at the true θ_a, so that E[W_a] = 1. With short histories it is not. In the recovery
scenarios at T = 20 its mean is about 1.10 for both parameters (the median over the T = 20 scenarios,
S-4's H2), and the profile interval
covers about a point short of 95% (S-13).

A Bartlett correction rescales the statistic by its mean, so that W_a / k_a is closer to χ²₁. The
interval becomes

    { θ_a : P_a(θ_a) ≥ l_max − k_a c / 2 },

which is the profile interval solved at the threshold k_a c / 2 instead of c / 2. Its ends are found
exactly as the profile interval's (`binomial_mixture_mle.md`), and are truncated at the grid's bounds
in the same way.

## 2. Estimating k_a by parametric bootstrap

The true mean E[W_a] depends on the unknown parameters, so it is estimated at the estimate
θ̂ = (PD̂, ρ̂):

1. **Simulate B panels.** Panel b (b = 0, …, B − 1) is drawn from the one-factor DGP (`dgp.md`) at
   (PD̂, ρ̂), with the original panel's obligor counts n_t period by period. The DGP's key is the caller's
   seed, the caller's stream as its scenario field, and first_panel + b as its replicate field. Panel b
   is therefore a pure function of (seed, stream, first_panel + b, PD̂, ρ̂, n_1, …, n_T).
2. **Refit each panel** exactly as the original: the same grid, quadrature rules and calibration.
3. **Compute W\*_{a,b} = 2 (l\*_max − P\*_a(θ̂_a))**, panel b's profile statistic at the value θ̂_a
   that generated it.
4. **Average:** k_a is the mean of the finite W\*_{a,b}, summed in panel order.

k_a falls back to 1, the plain profile interval, in two cases:

- **On the bound:** θ̂_a is on the grid's bound, where the bootstrap distribution is truncated.
- **Too few finite:** more than 10% of the W\*_{a,b} are not finite.

If the estimate itself is unusable (flat or numeric flags), no bootstrap is run. The factor is then
NaN, and the intervals are the profile intervals, flagged not computed as they would be anyway.

B defaults to 199, the count S-4 used. Each estimate therefore costs about B fits. The bootstrap
panels share the n_t, so their surface rows repeat, and they are computed through a row cache (D-167),
whose results are identical.

## 3. Results and reproducibility

- **Reproducible:** the result is deterministic in (panel, grid, seed, stream, first_panel, B), and
  bitwise identical for any thread count.
- **Tied to the study run:** with seed 0x53313050424F4F54, stream = the recovery scenario,
  first_panel = 1000 × the replicate and B = 199, it is S-4's S-4b for that replicate. The unit test
  `bartlett_reproduces_the_s4_run` checks five replicates of the committed run bit for bit
  (`tests/golden/bartlett/s4_run.csv`).
- **What S-4 found** (D-181), on the 81 recovery scenarios × 1,000 replicates:
  - **at T = 20:** coverage of ρ 0.947 and of PD 0.945, against the profile interval's 0.937 and 0.940
    and Jeffreys' 0.943 and 0.942;
  - **in the band:** no verdict below the band in the 50 informative scenarios;
  - **width:** about 4–5% wider than the profile interval at T = 20, and within 1% of it at T = 100.
- **Why it is not the default** (D-182): S-4's registered rule judged it against Jeffreys, and found it
  a match bought with width. Changing the default needs a new registered study with an explicit
  width–coverage rule.

## 4. Replication in a short script

Given a fitter that returns l_max and the profile log-likelihood (as `validation/scipy/conditional_pd_profile.py`
does for a profile interval) and the DGP's line-for-line Python implementation (`tools/gen_dgp_tables.py`;
`dgp.md`), the factor is:

```python
W = []
for b in range(B):
    d = simulate_panel(seed, stream, first_panel + b, pd_hat, rho_hat, n)
    fit_b = fit(n, d)
    W.append(2 * (fit_b.loglik_max - fit_b.profile(a, theta_hat[a])))
finite = [w for w in W if math.isfinite(w)]
k = 1.0 if on_bound or B - len(finite) > 0.10 * B else sum(finite) / len(finite)
```

The interval is then the profile interval at threshold k · c / 2. The engine's ends are bitwise only
against the engine. An independent script agrees to the profile-endpoint tolerances
(`tolerances.md`).
