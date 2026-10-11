# S-36 Rating philosophy and fitted correlation: design and predictions

Status: **DRAFT, not registered.** This is a proposal for the owner. It is registered only when the owner
approves it and it is merged. No generator code exists yet; nothing has been simulated or fitted, and nothing
may be until this file is registered (D-148). The reference values below come from the model alone
([`reference.py`](reference.py), [`reference.csv`](reference.csv)): quadrature, no panels.

The design and the predictions are written together, before the generator, so that the predictions are blind
by construction. Implementation readings that the code later forces (as with D-203) go in an appended,
dated addendum; no prediction, threshold or scope changes after registration.

## Question

How does the correlation ρ fitted to a grade's default history change as the rating system moves from
through-the-cycle (TTC) to point-in-time (PIT)? Miu and Ozdemir (2005) argue that the correlation a grade's
history implies is higher under TTC. That direction is the candidate prediction; this study gives it a
size, and asks what a single-grade fit recovers at realistic T and n.

## The DGP (proposed; owner to approve each choice marked ◆)

**Obligors.** Each period has n = 10,000 obligors, drawn afresh from a fixed population distribution
(◆ no persistence or migration between periods; S-33 covers rating dynamics). Obligor i in period t has a
TTC PD on the master scale of the grade scenario matrix (D-151): it falls in TTC grade k with the scale's
share, and its PD is log-uniform within the grade, on [PD_k/√1.5, PD_k·√1.5) (◆). Its TTC probit is
c_i = Φ⁻¹(TTC PD_i).

**The factor and defaults.** One factor Z_t per period, independent N(0, 1), as in the base DGP. Obligor i
defaults when U_i < p_i(Z_t), with the base model's conditional PD,
p_i(z) = Φ(x_i(z)), x_i(z) = (c_i − √ρ·z)/√(1 − ρ) (a higher z is better). Every obligor shares the same
asset correlation ρ, the scale's (L 0.24, M 0.12, H 0.02).

**The rating rule.** Obligor i is rated on a hybrid score, a convex combination of its TTC and PIT probits
(◆):

    s_i(Z_t) = (1 − λ)·c_i + λ·x_i(Z_t),   λ ∈ [0, 1].

At λ = 0 the rating is the obligor's TTC grade, fixed through the cycle, so a grade's default rate moves
with Z. At λ = 1 it is the PIT probit itself, so a grade's conditional PD stays inside the grade's band and
its default rate barely moves. Grade g holds the obligors with s in [b_{g−1}, b_g), where the b are the
probits of the geometric midpoints between adjacent master-scale PDs; grades 1 and 17 are open-ended. A
grade keeps its master-scale label at every λ (◆; under a hybrid rule the label is no longer the grade's
long-run PD, which is part of what is reported).

**Common random numbers.** The factor draws, the obligors' TTC PDs and their default uniforms do not depend
on λ: λ only relabels obligors into grades. So every λ is fitted on the same defaults, the comparison across
λ is paired, and the portfolio-level panel (all grades together) is identical at every λ, which Q1 checks.

**Recorded data.** Per grade and period, n_gt (the obligors rated g) and d_gt (their defaults). A period with
n_gt = 0 carries no count and drops out of that grade's fit.

**Draws.** On Philox, addressed like the base DGP's, with a new reserved seed (◆ to be reserved in the
decision that registers this file) so that no draw is shared with any existing study: per (scenario,
replicate, period), the factor first, then per obligor its TTC grade, its position within the grade and its
default uniform, each from its own block range. Exact arithmetic as in `docs/methodology/dgp.md` (the
in-house log, exp and Φ⁻¹). As every new DGP variant, it ships with a prose specification, a line-for-line
Python mirror and a reference-panel hash (D-149, D-110 to D-112), before this study runs.

## The runs

- **Scenarios:** the grade scenario matrix (D-151): scales L, M, H × T ∈ {20, 40, 100}, 9 scenarios,
  n = 10,000 per period, R = 1,000 replicates (◆).
- **λ ∈ {0, 0.25, 0.5, 0.75, 1}** (◆), all five on every replicate (common random numbers).
- **Fits:** the parity binomial-mixture MLE, per grade (17 per λ), with profile intervals for PD and ρ; and
  once per replicate, the portfolio-level fit to all grades together.
- **Informative grades (fixed now, from the model):** grades 2–16 whose expected defaults per period at
  λ = 0, n·share_k·PD_k, are at least 5. That is grades 9–14 of L, 8–16 of M and 9–16 of H (23 grades).
  The edge grades 1 and 17, and the thin grades, are reported, not predicted.
- **Cost:** about 9 × 1,000 × (5 × 17 + 1) ≈ 774,000 single-grade fits, plus the obligor simulation (10⁴
  draws per period, once per replicate, shared by the five λ). Projected from a timing sample before the
  run (D-187, D-193), on another seed.

## Targets (D-150: two targets)

This is a misspecification study: a single-grade Vasicek model is fitted to a grade whose membership moves
with the factor. So:

- **The primary target is the pseudo-true value** ρ_g(λ), the value the fit converges to on a very long
  history. With the period-to-period composition stationary, the very-long-panel value is well defined.
  `reference.py` computes it at n → ∞ (the Vasicek-rate pseudo-true values, the first two moments of the
  grade's probit default rate over the populated periods); the binomial MLE at n = 10,000 targets about
  the same value (S-27), which the tolerances below allow for.
- **The distance from the generating ρ** (the obligors' asset correlation) is reported beside it: it is
  what a practitioner reading a PIT grade's ρ̂ as an asset correlation gets wrong.

**The closed form.** Inside a narrow grade, the grade's probit default rate is linear in z, and its
Vasicek-implied correlation is

    ρ_eff(λ) = v/(1 + v),   v = ρ(1 − λ)² / ((1 − ρ)·k²),   k = 1 − λ + λ/√(1 − ρ),

so ρ_eff = ρ at λ = 0 and 0 at λ = 1. The finite-grade values in `reference.csv` agree with it to within
0.002 on every informative grade at every λ. They depart from it on the edge grades (1 and 17,
open-ended) and on grades near the scale's ends at λ ≥ 0.5, where a grade empties in good or bad years (in
`reference.csv`, `p_empty` is the factor probability of a grade-period with fewer than one expected obligor,
on a 64-node Gauss–Hermite grid, so it is coarse). Under the hybrid rule the long-run PD of
a grade also falls below its label (at λ = 0.5 on scale L, 0.38–0.96 of the master-scale PD on grades 1–16),
because a hybrid score is not the probit of any PD; at λ = 1 it returns to the label (0.97–1.02 on grades
2–16).

| Scale (ρ) | λ = 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| L (0.24) | 0.240 | 0.142 | 0.064 | 0.016 | 0 |
| M (0.12) | 0.120 | 0.069 | 0.031 | 0.008 | 0 |
| H (0.02) | 0.020 | 0.011 | 0.005 | 0.001 | 0 |

(ρ_eff, closed form. The box's lower bound for ρ is 10⁻³, so every value from λ = 0.75 on H and λ = 1
everywhere is at or below it.)

## Predictions

Strict pass/fail unless marked. "Median" is over replicates; coverage is against the D-131 band
(0.927–0.973 at R = 1,000).

| # | Prediction |
|---|---|
| Q1 | **Check (stops scoring if it fails):** in every replicate, the portfolio-level fit (all grades) is bit-for-bit identical at the five λ. Defaults do not depend on λ, only grade membership does |
| Q2 | **λ = 0 recovers the asset correlation:** at T = 100, in at least 20 of the 23 informative grades the ρ profile interval's coverage of the generating ρ is in the band (at λ = 0 the pseudo-true value is the generating ρ to within 0.002) |
| Q3 | **The direction (Miu and Ozdemir):** in every scenario and every informative grade, median ρ̂ falls from λ = 0 to λ = 0.5, and is non-increasing across the five λ, ties at the box's floor allowed |
| Q4 | **The size, at T = 100:** in every informative grade at each (λ, scale) whose pseudo-true ρ is at least 0.01 (L at λ ≤ 0.75, M at λ ≤ 0.5, H at λ ≤ 0.25), median ρ̂ / ρ_g(λ) is between 0.85 and 1.05. Below 1 allows for ρ̂'s known small-sample downward bias (M1; S-3) |
| Q5 | **At λ = 1 the fit sits at the floor:** at T = 100, pooled over the informative grades, the share of fits with ρ̂ at the box's lower bound (flagged) is between 0.35 and 0.65 (a pseudo-true value below the bound puts about half the estimates on it) |
| Q6 | **The interval no longer covers the asset correlation:** at T = 100 and λ ≥ 0.5, in every informative grade, the ρ profile interval covers the generating ρ in under half of the replicates |
| Q7 | **Coverage of the pseudo-true value, at T = 100 (expected range, not pass/fail):** at λ = 0.25 and 0.5, pooled over the informative grades, the ρ profile interval covers ρ_g(λ) with frequency 0.90–0.97. A shortfall would mean the composition change adds variance the single-grade likelihood does not see |
| Q8 | **The grade's long-run PD, at T = 100:** in every informative grade and every λ, median PD̂ / PD_g(λ) (the pseudo-true PD in `reference.csv`) is between 0.95 and 1.05 |

**Reported without a prediction:** everything at T = 20 and 40 (where ρ̂'s small-T bias and the floor
dominate); the edge and thin grades; the share of empty grade-periods by λ; the generating-ρ distance per
grade; and, for capital, each grade's 99.9% conditional PD from its fitted (PD̂, ρ̂) beside the one from the
generating ρ and the master-scale PD (S-23's interval), which is the practical consequence of reading a PIT
grade's ρ̂ as an asset correlation.

## Verdicts it could change

None of the pinned verdicts: it is a new DGP with its own reference values (ρ̂ against λ).

## Monitoring implication

To be filled in when the study finishes (D-152). The candidate: a check that the correlation a calibration
uses matches the rating philosophy of the grades it is applied to.

## Owner decisions needed before registration

1. The rating rule: a convex combination of TTC and PIT probits (above), or another form (for example a
   partial adjustment of the previous period's grade, which adds rating momentum).
2. The population: obligors redrawn every period, and log-uniform PDs within a TTC grade.
3. Grade labels kept at the master scale under every λ.
4. The λ grid and the scenario set (the grade scenario matrix at R = 1,000).
5. The informative-grade rule (interior, at least 5 expected defaults per period at λ = 0), and the bands in
   Q2 to Q8.
6. A seed to reserve, in the decision that registers this file.
