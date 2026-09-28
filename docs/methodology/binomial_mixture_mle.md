# Binomial-mixture maximum likelihood (the parity profile)

This note describes, in plain prose, how vasicek-cal estimates the probability of default (PD)
and the asset correlation ρ of the one-factor Vasicek model from a panel of default counts. It
covers the numerical method, the ranges over which it has been validated, what each run checks
and reports about itself, and what is known about its statistical behaviour. Anyone replicating
the fit should not need anything else. Decision numbers (D-nnn) point to `DECISIONS.md` for the reasoning, but
are not required reading.

The accompanying script `validation/scipy/binomial_mixture_mle.py` reproduces the method with
ordinary scipy tools (section 8).

## 1. Model and likelihood

A panel has periods t = 1, …, T. In period t there are n_t obligors, of whom d_t default.
Conditional on a systematic factor Z_t ~ N(0, 1), independent across periods, each obligor
defaults independently with probability

    p(z) = Φ( (Φ⁻¹(PD) − √ρ · z) / √(1 − ρ) ).

So d_t, given Z_t = z, is Binomial(n_t, p(z)), and the log-likelihood of period t is

    l_t(PD, ρ) = log C(n_t, d_t) + log ∫ p(z)^d_t (1 − p(z))^(n_t − d_t) φ(z) dz.

The panel log-likelihood is ℓ(PD, ρ) = Σ_t l_t. The binomial coefficient does not depend on the
parameters, but it is kept, so that log-likelihoods, and anything built on them such as AIC,
BIC and likelihood ratios, are the true values (D-070). Everything is computed in log space:
log Φ comes from a dedicated function accurate in both tails, never from log(Φ(x)).

## 2. The integral

The integrand exp(h(z)), with h(z) = d log p(z) + (n − d) log(1 − p(z)) − z²/2, has two shapes,
and parity uses a different rule for each (D-118).

**0 < d < n: adaptive Gauss–Hermite, N = 128 nodes.** h has a single peak, which for large n is
very narrow. The rule is centred on the peak and scaled to its width (Liu and Pierce, 1994):

- **Mode z\*:** the root of h′, found by Newton's method with bisection safeguards, starting
  from where p(z) = (d + ½)/(n + 1).
- **Scale:** σ = 1/√(−h″(z\*)).
- **Nodes:** z_i = z\* + σ·u_i, where u_i and w_i are the standard 128-point Gauss–Hermite
  nodes and weights (tabulated from mpmath at 50 digits, D-080).
- **Integral:** log ∫ … = log σ + log Σ w_i exp(h(z_i) + u_i²/2) − log √(2π), with the sum
  taken as a log-sum-exp.

**d = 0 or d = n: composite Gauss–Legendre, 16 panels × 16 points.** Now h has no interior peak.
The integrand is the normal prior cut off by a steep "survival" factor, (1 − p)ⁿ for d = 0 or
pⁿ for d = n. A polynomial rule centred on a mode resolves that cliff slowly. With 128 nodes,
adaptive Gauss–Hermite is off by up to 1.9·10⁻⁵ in a period's log-likelihood at ρ = 0.5, and by
6·10⁻³ at ρ = 0.9 (D-116). The rule:

- **Domain:** where h is within e⁻⁵⁰ of its maximum. Each end is found by doubling a step away
  from the peak, then exactly 64 bisections.
- **Transition point and width:**
  - z_c is where the survival factor equals ½. For d = 0 that is where p = 1 − 2^(−1/n), about
    ln 2 / n; for d = n it is the mirror image.
  - w is the factor's width there, 1/(n·β·λ), with β = √ρ/√(1 − ρ) and λ the inverse Mills ratio.
- **Nodes:** substitute z = z_c + w·sinh(u), so points crowd at the cliff and spread
  geometrically into the bulk. Divide the u-range into 16 equal panels with 16 Gauss–Legendre
  points each.

The panel count and the placement are the same at every parameter value. The log-likelihood
surface is therefore a smooth function of the parameters, with no count to jump. That is tested:
second differences agree with an independent implementation to 4e-14 of their size.

**Why not ordinary (fixed-node) Gauss–Hermite (D-037).** A fixed rule places its nodes where the
prior φ(z) has its mass, not where the integrand has it. For large n the integrand is a spike of
width about 1/√n, far narrower than the node spacing. Take n = 10⁶, d = 1,000, PD = 0.1%,
ρ = 0.12:

- the true period log-likelihood is −7911.0799;
- fixed-node Gauss–Hermite with 128 nodes gives −7923.29, wrong by 12.2;
- the adaptive rule with only 32 nodes is exact to all printed digits.

An error of 12 in one period's log-likelihood is several times the 1.92 that decides a 95%
interval. The error does not show up as noise: it is a smooth bias, so the fit converges
confidently to the wrong answer. That is why parity uses the adaptive rule, and why a replication
using `numpy.polynomial.hermite.hermgauss` without recentring will not reproduce the engine for
large obligor counts.

## 3. Validated ranges and accuracy

Precision is measured in units of machine epsilon ε ≈ 2.2e-16. For log-likelihoods the unit is
scaled to the size of the terms, ε·max(1, |log C| + |log I|). The two terms can nearly cancel
(n = 10⁶, d = 350,000: each about 6.5·10⁵, their sum −15), and no finite-precision method can do
better than ε times the size of the terms.

- **Against mpmath at 50 digits** (golden set): every case is within 11.7 ε for zero- and
  all-default periods and 2.2 ε otherwise. The set covers PD 10⁻⁶ to 0.5, ρ 10⁻⁴ to 0.99, n from
  1 to 10⁶, d = 0, d = n and interior d.
- **Against the independent reference implementation** (`ref/`, which shares no code or
  algorithm): 4,296 periods, all within 35 ε (term-scaled). They span PD 10⁻⁶ to 0.5, ρ 10⁻⁴ to
  0.9, n from 1 to 10⁶, with d ∈ {0, n, ≈ n·PD, ≈ 3n·PD}. The largest absolute error anywhere is
  2.3e-10, on a period whose terms are 6.5·10⁵ in size, i.e. rounding.
- **Estimates against `ref/`:** within 0.026 SE on six simulated panels. SEs agree within 1.6%
  on the four panels not flagged near a bound.

Outside these ranges (for example ρ > 0.99, or n > 10⁶) nothing is claimed. The per-run check
below still applies there.

## 4. The per-run check (D-092, D-120)

Every fit re-evaluates each period at the estimate with the rule doubled: Gauss–Hermite with 256
nodes, and the composite rule with 32 panels.

- **Flagged period:** one whose two values differ by more than max(10⁻¹⁰, 64·ε·(|log C| + |log I|)).
  The second term is the rounding floor for large obligor counts.
- **Run flag:** any flagged period sets `quadrature unconverged` on the fit.
- **Always reported:** the per-period maximum and the total difference.

**Does the check see real errors?** It is tested on the golden set and against `ref/`: no period
whose true error exceeds the threshold goes unflagged. It is also tested on a deliberately
coarse version of the same rules, which is inaccurate on 42 golden periods. There the check
never reported less than 0.99998 of the true error. The tolerance assumes it reports at least
half, so twice the reported figures bounds the true error. That factor is empirical, not proven.

## 5. Materiality

A likelihood-ratio 95% interval moves by a change of 1.92 in log-likelihood.

- **Per period:** the largest numerical error found anywhere in the validated range is 2.3e-10,
  about 10¹⁰ times smaller.
- **Per panel:** even a panel of 10,000 such periods would carry at most 2.3e-6.
- **Conclusion:** within the validated ranges the numerical error of parity is immaterial to
  any estimate, interval or test.

Outside them, the per-run check's reported maximum and total are the measure. Doubled
(section 4), they bound the error to compare with 1.92.

## 6. Estimation

1. **Grid.** The log-likelihood surface is evaluated on a rectangular grid in logit coordinates,
   u = log(v/(1 − v)), for both PD and ρ. The box and resolution are set by the caller. The
   validation runs use PD ∈ [10⁻⁴, 0.2] with 61 points and ρ ∈ [10⁻³, 0.5] with 41 points. ρ is
   capped at 0.5 by default (D-089).
   - Periods with identical (n, d) are evaluated once (D-122), with no effect on the results.
   - Sums over periods are compensated and taken in a fixed order, so results are bitwise
     identical for any number of threads.
2. **Grid maximum.** The grid point with the largest log-likelihood, ties going to the lower
   index.
3. **Refinement (D-095).** The exact quadratic, including the PD–ρ cross term, is fitted through
   the 3 × 3 block of grid points around the maximum, and its vertex is taken. If the vertex lies
   outside that block, the grid point is kept and the fit is flagged `refinement rejected`. The
   refinement is accurate to about 0.01 standard errors.
4. **Standard errors (D-119).** The Hessian of ℓ is taken at the refined estimate by central
   differences in logit coordinates. The step is 0.15 times the standard error implied by the
   grid stencil; results move by less than 5·10⁻⁵ between steps of 0.10 and 0.20. The covariance
   is the inverse of the negative Hessian, converted to PD and ρ by the delta method.

**Flags.** Each is reported on the fit, and none is resolved silently:

| Flag | Meaning | What to do |
|---|---|---|
| grid edge | the maximum is on the edge of the box (e.g. no defaults at all) | no SE; the estimate is the bound. Widen the box if the edge is artificial. |
| flat surface | the curvature is not negative definite | no SE; the data do not identify both parameters |
| refinement rejected | the quadratic's vertex is outside the stencil | estimate = grid point; SEs still reported |
| quadrature unconverged | a period failed the per-run check (section 4) | the reported maximum and total bound the error |
| numeric | a non-finite value appeared on the surface | investigate the inputs |
| **near bound** | the estimate is within 2 SEs of a bound of the box, in logit coordinates | the Wald interval is unreliable here; use the profile-likelihood interval (section 6a) |
| **ρ not identified** | no period has n ≥ 2 (D-302). A single obligor's likelihood is E[p(Z)] = PD or 1 − PD, with no ρ in it | ρ̂ is still reported but is meaningless: rounding picks it (S-1). Its profile interval is the whole box. PD is unaffected |

The near-bound flag matters in practice. Low-default panels and small correlations routinely put
ρ̂ within two standard errors of zero. There the likelihood is far from quadratic, and a
symmetric Wald interval misstates the uncertainty in either direction.

## 6a. Profile-likelihood intervals (M2; D-128–D-130)

The 95% profile-likelihood interval for PD is the set of PD values whose profile log-likelihood,
maximised over ρ, is within c = 1.92073 of the maximum. The same holds for ρ, maximised over PD.
c is half the 95% point of χ² with one degree of freedom. The interval needs no standard error,
follows the likelihood's actual shape, and is the right interval near a bound.

1. **The maximum** ℓ_max is found off the grid by Brent's method, as a nested maximisation. Each
   value of the PD profile is itself a maximisation over ρ.
2. **Each end:**
   - Walk outwards from the maximum along the grid. Grid values of the profile can only
     underestimate it, so a grid point above the threshold is inside the interval.
   - At the first grid point below the threshold, evaluate the profile exactly.
   - Solve profile = ℓ_max − c by Brent's root finder, to 10⁻⁹ in logit coordinates. The residual
     in log-likelihood is reported: at most 2·10⁻⁹ in the tests, against a threshold of 1.92.
3. **A bound of the box:** if the profile is still above the threshold at the bound, the end is
   the bound itself, flagged truncated. It is never extrapolated. At small ρ, or with very few
   defaults, the lower ρ end is routinely the bound.

An independent implementation (`ref/`: golden-section over the whole range and bisection) agrees
with the engine's endpoints to 1.5·10⁻¹⁰ in logit coordinates, and on every truncation.

**Boundary theory (D-129).** When the true ρ is 0, on the edge of the parameter space, the
likelihood-ratio statistic for *testing* ρ = 0 is not χ²₁. It is an equal mixture of χ²₀ and χ²₁
(Chernoff 1954; Self and Liang 1987), with 5% critical value 2.71 rather than 3.84. For
*confidence intervals* the engine deliberately uses χ²₁ with truncation at the bound, the
standard practice. A test of "no correlation" should use the mixture, and is not offered by the
parity profile. The box's lower ρ bound is 10⁻³, so intervals are truncated there.

## 6b. Which interval to use

**Profile-likelihood intervals (section 6a) are the recommended method for inference in this
model.** In the recovery study they cover as advertised in 136 of 162 settings. Where they do not,
they are either conservative (near-uninformative data) or within a point of the target at 20 to
40 periods. Wald intervals undercover at small T and are unreliable near a bound.

Bootstrap percentile intervals are provided for comparison, and **are not recommended for ρ**.
They undercover for ρ at every sample size studied, and for PD with 40 or fewer periods, for two
reasons:

- **Boundary breakdown.** Near a bound the bootstrap is inconsistent (Andrews 2000): resampled
  estimates pile onto the bound and intervals collapse.
- **No bias or skewness correction.** Elsewhere ρ̂ is biased downwards and skewed, and the
  percentile interval corrects neither, so its intervals sit below the truth.

Coverage by group (1,000 panels per setting; the band is 0.927–0.973; per-setting verdicts are
pinned in the recovery goldens, D-137):

| Group | Parameter | Scenarios | PASS | Below the band | Bootstrap coverage | Zero-width intervals | Profile coverage |
|---|---|---|---|---|---|---|---|
| A: near a bound / near-uninformative | PD | 31 | 5 | 26 | 0.762–0.944 | 622 | 0.930–0.986 |
| A | ρ | 31 | 0 | 31 | 0.431–0.922 | 2,878 | 0.924–0.999 |
| B: T = 20 | PD | 12 | 1 | 11 | 0.875–0.936 | 0 | 0.920–0.953 |
| B | ρ | 12 | 0 | 12 | 0.803–0.887 | 0 | 0.921–0.948 |
| C: T = 40 | PD | 17 | 7 | 10 | 0.899–0.955 | 0 | 0.924–0.964 |
| C | ρ | 17 | 0 | 17 | 0.874–0.909 | 0 | 0.934–0.956 |
| D: T = 100 | PD | 21 | 18 | 3 | 0.909–0.948 | 0 | 0.939–0.955 |
| D | ρ | 21 | 6 | 15 | 0.908–0.938 | 0 | 0.938–0.961 |

## 6c. Intervals for the 99.9% conditional PD (S-23)

Capital and stress calculations use the PD conditional on an adverse factor level, not PD or ρ
themselves. At the 0.1% adverse quantile of the factor that quantity is

    q(PD, ρ) = Φ( (Φ⁻¹(PD) + √ρ · Φ⁻¹(0.999)) / √(1 − ρ) ),   Φ⁻¹(0.999) = 3.0902.

The engine's convention is that a higher Z means better conditions, so the adverse level is
Z = −3.09. `engine/conditional_pd.hpp` gives q, its point estimate q̂ = q(PD̂, ρ̂), and three
95% intervals.

**Profile likelihood (recommended).** The interval is the set of values c whose profile
log-likelihood P_q(c), the maximum of ℓ(PD, ρ) over the box subject to q(PD, ρ) = c, is within
1.92073 of ℓ_max: the same threshold as section 6a. On the curve q = c the PD is fixed by ρ:

    PD_c(ρ) = Φ( √(1 − ρ) · Φ⁻¹(c) − √ρ · Φ⁻¹(0.999) ),

so P_q(c) is a maximisation over ρ alone, over the values of ρ for which PD_c(ρ) lies in the box.

1. **Inner maximum.** A guide comes first: the grid's surface, interpolated at points of the curve
   (eight per grid step of ρ). It chooses a bracket, and Brent's method then maximises the
   likelihood itself over logit ρ. The bracket widens until the maximum is interior or on a bound.
   A bound is an end of the ρ axis, or a point where PD_c(ρ) reaches an end of the PD axis.
2. **Each end.** q is not a grid axis, so the end is bracketed by walking outwards from the
   maximum in logit(q), one PD grid step at a time. Points where the guide is above the
   threshold are passed over. The bracket's own ends are always evaluated exactly, since the
   guide is not a bound. Brent's root finder then solves P_q(c) = ℓ_max − 1.92073 to 10⁻⁹ in
   logit(q), and the residual is reported: at most 8.8·10⁻¹¹ in the unit test, and asserted
   below 10⁻⁷ on every recovery fit.
3. **The box.** q attainable in the box runs from 1.455·10⁻⁴ (PD = 10⁻⁴, ρ = 10⁻³) to 0.9713
   (PD = 0.2, ρ = 0.5). An interval that is still above the threshold there ends at that limit,
   flagged *truncated*, and is never extrapolated. An end whose inner maximiser lies on a bound
   of the box is flagged *box-limited*. Such an end can be held by a bound through the nuisance,
   ρ usually, without q itself reaching its range. Truncation implies box-limited.

**Delta-method Wald.** q̂'s standard error in logit(q) comes from the covariance of the
estimates at the Hessian (section 6) and the gradient of logit(q) with respect to logit PD and
logit ρ. The interval is logit(q̂) ± 1.96·SE. It is unreliable wherever the parameters' Wald
intervals are: on the grid edge, near a bound, or on a flat surface.

**Bootstrap percentile.** q evaluated at each iid bootstrap replicate's estimates (B = 999),
with type-7 quantiles. Percentile intervals are invariant under monotone transformations, so
this is the same as taking them in logit(q).

**Recovery results (S-23, D-156).** Over the 81 recovery scenarios × 1,000 replicates:

- **Profile interval for q:** as reliable as those for PD and ρ. In every informative scenario its
  coverage is within 0.015 of theirs. Two scenarios at T = 20 are below the band. 15 are above it,
  where the data are nearly uninformative and the ends are held by the box.
- **Delta-method Wald interval:** undercovers in 17 of the 50 assessed scenarios, and not only at
  T = 20. Its SE is right, but q̂ is biased low (median q̂/q − 1 between −1% and −11%), because
  ρ̂'s small-sample bias carries into q.
- **Bootstrap percentile interval:** undercovers in 73 of 81 scenarios, for the reasons found for ρ.

Use the profile interval for q. The study's entry in `studies/README.md` lists the mitigations and
the monitoring implication.

A second scipy script (section 8) solves the profile end points on its own and checks the
engine's.

## 7. Known statistical behaviour (M1.8, M2a, M2b)

The recovery harness fits 1,000 simulated panels for each of 81 settings: PD 0.1–5%, ρ 0.02–0.24,
T from 20 to 100 periods, and n from 100 to 10,000 obligors ([recovery.md](recovery.md),
[recovery_results.md](recovery_results.md)).

- **The engine is correct.** Error falls as the number of periods grows, in every setting and
  for both parameters. No fit of the 81,000 needed its quadrature flagged.
- **ρ̂ is biased downwards at small T.** Where the data identify ρ (coverage not deferred),
  the bias is 0.0002–0.014, largest with 20 periods. In low-information settings, where most
  fits sit at a bound, it is larger in either direction (−0.08 to +0.07). This is a property of
  the maximum-likelihood estimator, not of the engine, and it is reported, not corrected.
- **Wald intervals undercover a little at small T.** With 20 to 40 periods, 95% intervals for PD
  cover 91–93% in some settings. Here the standard error, taken from the curvature, is 1–9%
  smaller than the actual spread of the estimates. A t quantile on T − 1 degrees of freedom
  recovers part of it: 0.924–0.943 in those settings, reported as a `native` comparison (D-131).
- **Wald intervals for ρ can overcover.** 97–100% in five settings, because ρ̂ is skewed
  towards zero in logit coordinates.
- **Near a bound, Wald coverage is not assessed.** Where at least 5% of fits are flagged near a
  bound (62 of 162 setting × parameter combinations, mostly low-default or low-correlation
  settings), the Wald verdict is deferred.
- **Profile-likelihood intervals (section 6a) cover as advertised in most settings.** Their
  coverage counts every fit, with nothing deferred: 136 of 162 verdicts are inside the band.
  - They fix 15 of the 19 Wald shortfalls and 41 of the 62 deferred verdicts.
  - In 20 near-uninformative settings they are **conservative** (0.977–0.999), because truncation
    at a bound widens them.
  - In 6 settings with 20 to 40 periods they still **undercover** slightly: 0.920–0.927 against a
    band from 0.927. For ρ, that follows ρ̂'s downward bias.
  - Two remedies are on the backlog as options, not changes: a Bartlett-type threshold correction
    and a bias-corrected ρ̂.

The 19 out-of-band verdicts (12 for PD, 7 for ρ) are listed individually, each with its
diagnosis and its standard-error ratio, in `recovery_results.md`. The parity sign-off
([parity_signoff.md](parity_signoff.md)) lists them together with the 62 deferred verdicts.

## 8. Replicating with scipy (D-126)

`validation/scipy/binomial_mixture_mle.py`, 125 lines, deliberately does *not* copy the
engine's numerics. It uses:

- `scipy.stats.binom.logpmf` for the binomial term, in log space;
- `scipy.integrate.quad` for the integral. The integration range is centred on the integrand's
  mode (found by Brent's method) and extends until the integrand has dropped by e⁻⁵⁰; the
  requested relative accuracy is 10⁻¹².
- `scipy.optimize.minimize` (Nelder–Mead) for the fit, in logit coordinates over the same box.

It checks three things against the engine's committed results on the six simulated panels
of M1.7:

1. the panel log-likelihood at the engine's 702 surface points. Worst relative difference
   4.1e-13; tolerance 10⁻¹¹.
2. the log-likelihood at scipy's own optimum minus that at the engine's estimate. Worst 3.9e-4,
   always positive; tolerance 10⁻³. That is the cost of the engine's grid refinement, which is
   accurate to about 0.025 SE.
3. the estimates, in units of the engine's standard error. Worst 0.026; tolerance 0.06.

The tolerances are the `TOL_SCIPY_*` entries of `tests/tolerances.toml`, which the script reads
directly. They describe what a straightforward independent implementation achieves, not the
ε-level agreement between the engine and its reference implementation. To run:

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/binomial_mixture_mle.py

It takes about five minutes. CI runs it on every commit with the same pinned versions.

**The conditional PD's profile interval (S-23).** `validation/scipy/conditional_pd_profile.py`
solves the end points of section 6c on its own, on the 162 recovery replay panels
(`tests/golden/recovery/replay_panels.csv`), and compares them with the engine's
(`replay.csv`). It uses:

- `scipy.special` (`log_ndtr`, `gammaln`) for the binomial term;
- the trapezoid rule for each period's integral, on a window around the integrand's mode. The
  window is placed from a Gaussian approximation and widened until both ends are e⁻⁶⁰ below
  the maximum. It is not quadrature code shared with the engine;
- nested bounded Brent over the box for ℓ_max, bounded Brent over the feasible stretch of logit ρ
  for the inner maximum, and `brentq` for each end point.

**Two pitfalls for anyone writing their own script.** Both showed up while this one was being
written, and both make an independent check disagree with a correct engine:

- **Nelder–Mead stalls on a bound of the box.** With bounds, `scipy.optimize.minimize(...,
  method="Nelder-Mead")` clips the simplex. When the maximum lies just inside the box (ρ̂ a little
  above 10⁻³ with few defaults), the simplex collapses onto the bound and stops there. On one
  replay panel ℓ_max came out 1.5·10⁻⁶ too low, which moved both interval ends by about 10⁻⁶ in
  logit(q). The script therefore maximises by nested one-dimensional searches over the whole box
  instead.
- **Bounded Brent stops short of a bound.** `minimize_scalar(..., method="bounded")` never
  evaluates at the bounds, and its tolerance includes a relative term √ε·|x| (about 10⁻⁷ in logit
  ρ). When the maximum lies *on* a bound, which is common for box-limited ends and for ρ̂ = 10⁻³,
  it returns a point up to that distance inside. With a steep likelihood that costs up to 10⁻⁷ in
  log-likelihood. The script evaluates the bound itself whenever the search ends within 10⁻⁶ of
  it, and keeps the better value.

Both errors are small in log-likelihood but systematic, and the second is exactly where the
interesting intervals are, at the box. A replication that skips these fixes should expect
disagreements of 10⁻⁷–10⁻⁶ in logit(q) at box-limited ends, and should not read them as engine
errors.

The end points agree within 1.2·10⁻⁹ in logit(q) (`TOL_SCIPY_Q_PROFILE_ENDPOINT_S` = 3·10⁻⁹).
Ends the engine truncated at the limit of q are inside the interval here too, and the
box-limited flags agree on every end. It takes about three minutes.

    uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/conditional_pd_profile.py
