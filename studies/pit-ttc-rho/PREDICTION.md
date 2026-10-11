# S-36 Rating philosophy and fitted correlation: design and predictions

Status: **DRAFT, not registered.** This is a proposal for the owner. It is registered only when the owner
approves it and it is merged. No generator code exists yet; nothing has been simulated or fitted, and nothing
may be until this file is registered (D-148). The reference values below come from the model alone
([`reference.py`](reference.py), [`reference.csv`](reference.csv)): quadrature, no panels.

The design and the predictions are written together, before the generator, so that the predictions are blind
by construction. Implementation readings that the code later forces (as with D-203) go in an appended,
dated addendum; no prediction, threshold or scope changes after registration.

**Owner's choices so far (2026-10-11):** the rating rule has momentum (a rating moves only part of the way
towards its target each period); the seed, the grade labels and the grid were left to the proposer and are
set below.

## Question

How does the correlation ρ fitted to a grade's default history change as the rating system moves from
through-the-cycle (TTC) to point-in-time (PIT), and how does rating momentum change that? Miu and Ozdemir
(2005) argue that the correlation a grade's history implies is higher under TTC. That direction is the
candidate prediction; this study gives it a size, separates the rating's philosophy from its speed, and asks
what a single-grade fit recovers at realistic T and n.

## The DGP

**Obligors.** Each period has n = 10,000 obligors, drawn afresh from a fixed population distribution.
Obligor i in period t falls in TTC grade k of the master scale with the grade scenario matrix's share
(D-151), and its TTC PD is log-uniform within the grade, on [PD_k/√1.5, PD_k·√1.5). Its TTC probit is
c_i = Φ⁻¹(TTC PD_i).

**The factor and defaults.** One factor Z_t per period, independent N(0, 1), as in the base DGP. Obligor i
defaults when U_i < p_i(Z_t), with the base model's conditional PD p_i(z) = Φ(x_i(z)),
x_i(z) = (c_i − √ρ·z)/√(1 − ρ) (a higher z is better). Every obligor has the scale's asset correlation ρ
(L 0.24, M 0.12, H 0.02).

**The rating rule: philosophy λ and speed θ.** The rating's target is a convex combination of the obligor's
TTC and PIT probits, y_i = (1 − λ)·c_i + λ·x_i; λ = 0 is TTC, λ = 1 is PIT. With momentum, the rating moves
only part of the way towards the target each period, at speed θ ∈ (0, 1]. Because the target's cyclical part
is common to all obligors, the rating uses the factor smoothed at that speed,

    M_t = (1 − θ)·M_{t−1} + θ·Z_t,    s_i,t = (1 − λ)·c_i + λ·(c_i − √ρ·M_t)/√(1 − ρ),

so θ = 1 is no momentum and smaller θ is a slower rating. M starts stationary: M_{−1} ~ N(0, θ/(2 − θ)),
drawn once per replicate, independent of the Z_t.

*The simplification, stated:* this is momentum through the systematic part only. Obligors are redrawn each
period, so there are no obligor-level migrations, survivorship or rating histories (S-33 covers those). For
the question asked, the grade's correlation, the systematic part is what matters: in the narrow-grade limit
below, an obligor-level partial adjustment of surviving obligors gives the same grade default rates.

**Grades.** Grade g holds the obligors with s in [b_{g−1}, b_g), where the b are the probits of the geometric
midpoints between adjacent master-scale PDs; grades 1 and 17 are open-ended. A grade keeps its master-scale
label under every rule. Under a hybrid or slow rule the label is no longer the grade's long-run PD, and that
gap is reported (Q9).

**Common random numbers.** The factor draws, M_{−1}, the obligors' TTC PDs and their default uniforms do not
depend on (λ, θ): the rule only relabels obligors into grades. Every rule is fitted on the same defaults, so
the comparison across rules is paired, and the portfolio-level panel (all grades together) is identical
under every rule, which Q1 checks.

**Recorded data.** Per grade and period, n_gt (the obligors rated g) and d_gt (their defaults). A period with
n_gt = 0 carries no count and drops out of that grade's fit.

**Draws.** On Philox, addressed like the base DGP's, with the seed **`0x5049545454435248` ("PITTTCRH")**,
used by no other study (checked against `DECISIONS.md` and the code). The runtime projection's timing sample
uses **`0x5049545454435453` ("PITTTCTS")**. Per (scenario, replicate, period): the factor first, then per
obligor its TTC grade, its position within the grade and its default uniform, each from its own block range;
M_{−1} from a period index of its own. Exact arithmetic as in `docs/methodology/dgp.md` (the in-house log, exp
and Φ⁻¹). As every new DGP variant, it ships with a prose specification, a line-for-line Python mirror and a
reference-panel hash (D-149, D-110 to D-112), before this study runs.

## The runs

- **Scenarios:** the grade scenario matrix (D-151): scales L, M, H × T ∈ {20, 40, 100}, 9 scenarios,
  n = 10,000 per period, R = 1,000 replicates.
- **Rules (13):** TTC (λ = 0, where θ has no effect), and λ ∈ {0.25, 0.5, 0.75, 1} × θ ∈ {1, 0.5, 0.25},
  all on every replicate (common random numbers).
- **Fits:** the parity binomial-mixture MLE, per grade (17 per rule), with profile intervals for PD and ρ;
  and once per replicate, the portfolio-level fit to all grades together.
- **A HAC interval beside the profile (owner, 2026-10-11).** For every grade fit, a Wald interval in logit
  units from a heteroskedasticity- and autocorrelation-consistent (HAC) sandwich SE: V = H⁻¹ S H⁻¹, with H
  the observed information of the fit and S the Newey–West long-run variance of the per-period scores
  (each period's gradient of its log-likelihood at the estimate), Bartlett kernel, lag
  L = ⌊4·(T/100)^(2/9)⌋ (2 at T = 20, 3 at T = 40, 4 at T = 100), fixed now. The Bartlett kernel keeps S positive
  semi-definite, so a negative autocorrelation shrinks the SE without making it undefined. Fits at the box's
  bound or flagged near-bound get no HAC interval (counted). This is the candidate `native` option S-35
  named (D-192): it needs its own decision, a methodology note, a scipy replication and a unit test against
  a hand-computed panel before this study runs, and it changes no parity result.
- **Informative grades (fixed now, from the model):** grades 2–16 whose expected defaults per period under
  TTC, n·share_k·PD_k, are at least 5: grades 9–14 of L, 8–16 of M and 9–16 of H (23 grades). The edge
  grades 1 and 17, and the thin grades, are reported, not predicted.
- **Cost:** about 9 × 1,000 × (13 × 17 + 1) ≈ 2.0 million single-grade fits, plus the obligor simulation
  (10⁴ draws per period, once per replicate, shared by every rule). Projected from a timing sample before the
  run (D-187, D-193). If the projection is too long for the laptop, the θ = 0.25 rules at T = 20 are dropped
  first, by an addendum before the run.

## Targets (D-150: two targets)

This is a misspecification study: a single-grade Vasicek model with independent periods is fitted to a grade
whose membership moves with the factor, and, with momentum, with its past. So:

- **The primary target is the pseudo-true value** ρ_g(λ, θ), the value the fit converges to on a very long
  history. Each period's composition and default rate have a stationary distribution, so the very-long-panel
  value is well defined, and an iid fit's target depends only on the period's marginal distribution (as in
  S-35). `reference.py` computes it at n → ∞ by quadrature over (Z_t, M_{t−1}): the first two moments of the
  grade's probit default rate over the populated periods. The binomial MLE at n = 10,000 targets about the
  same value (S-27), which the tolerances allow for.
- **The distance from the generating ρ** (the obligors' asset correlation) is reported beside it: it is what
  a practitioner reading a grade's ρ̂ as an asset correlation gets wrong.

**The closed form (narrow grades).** With a = 1/√(1 − ρ), k = 1 − λ + λa and β = λa/k, inside a grade the PIT
probit is a·s_g/k + a√ρ·(β·M_t − Z_t). So the grade's probit default rate is Gaussian, with

    v = ρ/(1 − ρ) · (1 − 2βθ + β²θ/(2 − θ)),    ρ_eff = v/(1 + v),

and lag-j autocorrelation (1 − θ)^(j−1) · γ₁, with γ₁ = (1 − θ)·βθ·(β/(2 − θ) − 1)·ρ/((1 − ρ)·v), which is
negative: a rating that lags the cycle leaves better obligors in a grade after a bad year, and worse ones
after a good year. Without momentum (θ = 1), v = ρ(1 − λ)²/((1 − ρ)k²): ρ at λ = 0, 0 at λ = 1. The
finite-grade values in `reference.csv` agree with the closed form to within 0.002 on every informative grade
under every rule.

**ρ_eff by rule** (closed form):

| Scale (ρ) | TTC | λ = 0.25: θ = 1, 0.5, 0.25 | λ = 0.5 | λ = 0.75 | λ = 1 (PIT) |
|---|---|---|---|---|---|
| L (0.24) | 0.240 | 0.142, 0.191, 0.216 | 0.064, 0.151, 0.196 | 0.016, 0.118, 0.181 | 0, 0.095, 0.169 |
| M (0.12) | 0.120 | 0.069, 0.094, 0.107 | 0.031, 0.072, 0.096 | 0.008, 0.056, 0.087 | 0, 0.044, 0.081 |
| H (0.02) | 0.020 | 0.011, 0.016, 0.018 | 0.005, 0.012, 0.016 | 0.001, 0.009, 0.014 | 0, 0.007, 0.013 |

Lag-1 autocorrelation of the grade's probit default rate: 0 without momentum; at θ = 0.5, −0.07 (λ = 0.25)
to −0.25 (λ = 1); at θ = 0.25, −0.05 to −0.13. (The box's lower bound for ρ is 10⁻³, so ρ_eff at λ = 1,
θ = 1, and at λ = 0.75, θ = 1 on scale H, is at or below it.)

**What the table says before any run.** Momentum undoes most of what PIT-ness does to the fitted correlation:
a slow PIT rating (θ = 0.25) gives a ρ̂ of 65–70% of the asset correlation on all three scales, while a
fast one gives about none. So "PIT or TTC" is not enough to say what ρ a grade's history implies; the
rating's speed matters as much.

**The grade's long-run PD against its label** (`pd_pseudo` / master-scale PD, informative grades): 1.01 under
TTC; without momentum, below the label for a hybrid rule (scale L at λ = 0.5: 0.80–0.89) and back at it under
PIT (0.98–1.00); with momentum, above it, most at slow PIT (scale L at λ = 1, θ = 0.25: 1.30–1.76). A lagging
rating keeps obligors in their grade after they have worsened.

## Predictions

Strict pass/fail unless marked. "Median" is over replicates; coverage is against the D-131 band
(0.927–0.973 at R = 1,000). "Informative grades" are the 23 fixed above.

| # | Prediction |
|---|---|
| Q1 | **Check (stops scoring if it fails):** in every replicate, the portfolio-level fit (all grades) is bit-for-bit identical under the 13 rules. Defaults do not depend on the rule, only grade membership does |
| Q2 | **TTC recovers the asset correlation:** at T = 100, in at least 20 of the 23 informative grades the ρ profile interval's coverage of the generating ρ is in the band (under TTC the pseudo-true value is the generating ρ to within 0.002) |
| Q3 | **Philosophy (Miu and Ozdemir), without momentum:** in every scenario and every informative grade, at θ = 1, median ρ̂ falls from λ = 0 to λ = 0.5, and is non-increasing across the five λ, ties at the box's floor allowed |
| Q4 | **Speed:** in every scenario at T = 100 and every informative grade, for each λ > 0, median ρ̂ rises as θ falls from 1 to 0.5 to 0.25, ties at the floor allowed |
| Q5 | **Size, at T = 100:** in every informative grade under every rule whose pseudo-true ρ is at least 0.01, median ρ̂ / ρ_g(λ, θ) is between 0.85 and 1.05. Below 1 allows for ρ̂'s known small-sample downward bias (M1; S-3) |
| Q6 | **Fast PIT sits at the floor:** at T = 100, at λ = 1 and θ = 1, pooled over the informative grades, the share of fits with ρ̂ at the box's lower bound (flagged) is between 0.35 and 0.65 (a pseudo-true value below the bound puts about half the estimates on it) |
| Q7 | **The interval no longer covers the asset correlation:** at T = 100, θ = 1 and λ ≥ 0.5, in every informative grade, the ρ profile interval covers the generating ρ in under half of the replicates |
| Q8 | **Coverage of the pseudo-true value, at T = 100 (expected ranges, not pass/fail):** pooled over the informative grades, (a) without momentum, at λ = 0.25 and 0.5, the ρ profile interval covers ρ_g with frequency 0.90–0.97; (b) with momentum, at θ = 0.5 and λ ≥ 0.5, the PD profile interval covers PD_g with frequency above 0.95: the negative autocorrelation makes the mean default rate less variable than an iid fit assumes (the factor part of its variance shrinks to about 0.4 of the iid value; the binomial part, independent across periods, does not shrink, so the overcoverage is smaller where defaults are few) |
| Q10 | **The HAC interval, at T = 100 (expected ranges, not pass/fail):** pooled over the informative grades, (a) with momentum, at θ = 0.5 and λ ≥ 0.5, the HAC interval for PD covers PD_g with frequency 0.90–0.97, closer to 0.95 than the profile interval's (Q8b): it sees the negative autocorrelation the iid likelihood misses; (b) under TTC, where periods are independent, its PD and ρ coverage of the generating values is within 0.03 of the profile interval's, at most 0.03 below it (the cost of estimating S with T = 100 periods) |
| Q9 | **The grade's long-run PD, at T = 100:** in every informative grade under every rule, median PD̂ / PD_g(λ, θ) (the pseudo-true PD in `reference.csv`) is between 0.95 and 1.05 |

**Reported without a prediction:** everything at T = 20 and 40 (where ρ̂'s small-T bias and the floor
dominate); the edge and thin grades; the share of empty grade-periods by rule; the generating-ρ distance per
grade; the realised autocorrelation of each grade's default rates against the closed form; the HAC interval
at T = 20 and 40, its width against the profile interval's, and its ρ coverage under every rule; and, for
capital,
each grade's 99.9% conditional PD from its fitted (PD̂, ρ̂) beside the one from the generating ρ and the
master-scale PD (S-23's interval), the practical consequence of reading a grade's ρ̂ as an asset correlation.

## Verdicts it could change

None of the pinned verdicts: it is a new DGP with its own reference values.

## Monitoring implication

To be filled in when the study finishes (D-152). The candidates: a check that the correlation a calibration
uses matches both the rating philosophy and the rating speed of the grades it is applied to, and a check of a
grade's realised default rate against its label, which a lagging rating moves.

## Still for the owner before registration

1. The momentum form: the smoothed-factor rule above (momentum through the systematic part, obligors
   redrawn), or obligor-level partial adjustment with persistent obligors (adds migrations and survivorship,
   overlaps S-33, and costs a persistent-population generator).
2. The informative-grade rule and the bands in Q2 to Q10.
3. The HAC interval's settings (Bartlett kernel, the lag rule, Wald in logit units).
