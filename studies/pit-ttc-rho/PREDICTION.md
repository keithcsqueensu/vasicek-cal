# S-36 Rating philosophy and fitted correlation: design and predictions

Status: **REGISTERED (D-401, owner's approval 2026-10-12).** No generator code exists yet; nothing has been
simulated or fitted, and nothing may be before this file is merged (D-148). The reference values below come
from the model alone ([`reference.py`](reference.py), [`reference.csv`](reference.csv),
[`reference_idio.csv`](reference_idio.csv)): quadrature, no panels, with a Monte Carlo check of arm I
([`reference_idio_mc.csv`](reference_idio_mc.csv)).

The design and the predictions are written together, before the generator, so that the predictions are blind
by construction. Implementation readings that the code later forces (as with D-203) go in an appended,
dated addendum; no prediction, threshold or scope changes after registration.

**Owner's choices so far:**

- 2026-10-11: the rating rule has momentum (a rating moves only part of the way towards its target each
  period); the seed, the grade labels and the grid were left to the proposer and are set below.
- 2026-10-12: momentum is both systematic and idiosyncratic, so the study has two arms: the main arm
  (systematic momentum only) and arm I (systematic and idiosyncratic momentum). The informative-grade rule
  and the bands were left to the proposer and are fixed below. The interval tests are core: the coverage
  predictions (Q8, Q10, Q11) are strict.

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

**Arm I: idiosyncratic momentum.** In the main arm the rating sees nothing obligor-specific beyond c_i, so
only the systematic part can lag. Arm I adds obligor-specific news the rating can see and lag. A share
κ = 0.5 of each obligor's idiosyncratic variance is an observable driver η_i,t, stationary N(0, 1) and AR(1)
in time with persistence ψ = 0.8; the rest stays unobservable noise:

    x_i,t(z) = (c_i − √ρ·z − √((1 − ρ)κ)·η_i,t) / √((1 − ρ)(1 − κ)).

The asset correlation is still ρ, and integrating η out gives the base model's conditional PD, so a TTC
grade's default rate is exactly as in the main arm. The rating tracks the driver with a lag of its own,
N_i,t = (1 − θ_I)·N_i,t−1 + θ_I·η_i,t, as M_t tracks Z_t at speed θ_S (the main arm's θ):

    s_i,t = (1 − λ)·c_i + λ·(c_i − √ρ·M_t − √((1 − ρ)κ)·N_i,t) / √((1 − ρ)(1 − κ)).

θ_I = 1 is a rating that sees today's news; smaller θ_I lags it. As in the main arm, obligors are drawn
afresh each period, now with (η, N) from their stationary joint law: N ~ N(0, V_N) and η | N ~ N(μN, σ²),
with V_N = θ_I²/(1 − (1 − θ_I)²) · (1 + (1 − θ_I)ψ)/(1 − (1 − θ_I)ψ), Cov(N, η) = θ_I/(1 − (1 − θ_I)ψ),
μ = Cov/V_N and σ² = 1 − Cov²/V_N. A population of persistent obligors without survivorship would give the
same cross-section in every period; this draw gives it without simulating each obligor's path (S-33 covers
survivorship). κ and ψ are fixed now, not varied: κ = 0.5 makes the observable news as large as the noise
the rating cannot see, and ψ = 0.8 makes it persistent enough that a lag matters.

**Grades.** Grade g holds the obligors with s in [b_{g−1}, b_g), where the b are the probits of the geometric
midpoints between adjacent master-scale PDs; grades 1 and 17 are open-ended. A grade keeps its master-scale
label under every rule. Under a hybrid or slow rule the label is no longer the grade's long-run PD, and that
gap is reported (Q9).

**Common random numbers.** The factor draws, M_{−1}, the obligors' TTC PDs and their default uniforms do not
depend on the rule: the rule only relabels obligors into grades. Within an arm, every rule is fitted on the
same defaults, so the comparison across rules is paired, and the portfolio-level panel (all grades
together) is identical under every rule of the arm, which Q1 checks. The two arms share the factor, the TTC
PDs and the default uniforms too; their defaults differ only through η in x.

**Recorded data.** Per grade and period, n_gt (the obligors rated g) and d_gt (their defaults). A period with
n_gt = 0 carries no count and drops out of that grade's fit.

**Draws.** On Philox, addressed like the base DGP's, with the seed **`0x5049545454435248` ("PITTTCRH")**,
used by no other study (checked against `DECISIONS.md` and the code). The runtime projection's timing sample
uses **`0x5049545454435453` ("PITTTCTS")**. Per (scenario, replicate, period): the factor first, then per
obligor its TTC grade, its position within the grade and its default uniform, each from its own block range;
M_{−1} from a period index of its own. Arm I adds two standard normals per obligor, G₁ and G₂, from two
further block ranges: N = √V_N·G₁ and η = μN + σG₂. Exact arithmetic as in `docs/methodology/dgp.md` (the in-house log, exp
and Φ⁻¹). As every new DGP variant, it ships with a prose specification, a line-for-line Python mirror and a
reference-panel hash (D-149, D-110 to D-112), before this study runs.

## The runs

- **Scenarios:** the grade scenario matrix (D-151): scales L, M, H × T ∈ {20, 40, 100}, 9 scenarios,
  n = 10,000 per period, R = 1,000 replicates.
- **Rules, main arm (13):** TTC (λ = 0, where θ has no effect), and λ ∈ {0.25, 0.5, 0.75, 1} × θ ∈ {1, 0.5,
  0.25}, all on every replicate (common random numbers).
- **Rules, arm I (13):** TTC (a control: the rating ignores η, so the grade's default rates are the main
  arm's TTC law), and λ ∈ {0.5, 1} × θ_S ∈ {1, 0.25} × θ_I ∈ {1, 0.5, 0.25}.
- **Fits:** the parity binomial-mixture MLE, per grade (17 per rule), with profile intervals for PD and ρ;
  and once per replicate and arm, the portfolio-level fit to all grades together.
- **A HAC interval beside the profile (owner, 2026-10-11).** For every grade fit, a Wald interval in logit
  units from a heteroskedasticity- and autocorrelation-consistent (HAC) sandwich SE: V = H⁻¹ S H⁻¹, with H
  the observed information of the fit and S the Newey–West long-run variance of the per-period scores
  (each period's gradient of its log-likelihood at the estimate), Bartlett kernel, lag
  L = ⌊4·(T/100)^(2/9)⌋ (2 at T = 20, 3 at T = 40, 4 at T = 100), fixed now. The Bartlett kernel keeps S positive
  semi-definite, so a negative autocorrelation shrinks the SE without making it undefined. Fits at the box's
  bound or flagged near-bound get no HAC interval (counted). This is the candidate `native` option S-35
  named (D-192): it needs its own decision, a methodology note, a scipy replication and a unit test against
  a hand-computed panel before this study runs, and it changes no parity result.
- **Informative grades (fixed now, from the model; proposer's rule, 2026-10-12):**
  - **The grades:** grades 2–16 whose expected defaults per period under TTC, n·share_k·PD_k, are at least
    5: grades 9–14 of L, 8–16 of M and 9–16 of H (23 grades), the same in both arms. The edge grades 1 and
    17 (open-ended, so their PDs are not the label's) and the thin grades are reported, not predicted.
  - **Why 5:** below it, a grade has a period with no defaults often enough (e^−5 ≈ 0.7% at the mean, far
    more in good years) that its ρ̂ is driven by zero counts, and the binomial noise in the period's default
    rate is as large as the factor signal on scale H. The rule uses the TTC population so that it does not
    depend on the rule being scored.
  - **Grade–rule pairs:** a pair whose reference probability of an empty period (`p_empty`, a period in
    which the grade's expected count is below 1) exceeds 0.01 is reported, not predicted: an empty period
    drops out of the fit, so its target is no longer the reference's. That excludes 11 of the main arm's
    299 pairs (fast hybrid and PIT rules, λ ≥ 0.75 and θ = 1, on the top grades of L and M, and one at
    λ = 1, θ = 0.5) and 4 of arm I's 299 (λ = 1, θ_S = 1, θ_I < 1, on L). The lists are fixed by the
    committed CSVs.
- **Cost:** about 9 × 1,000 × (26 × 17 + 2) ≈ 4.0 million single-grade fits, plus the obligor simulation
  (10⁴ draws per period, once per replicate and arm, shared by the arm's rules; arm I draws two normals more
  per obligor). Projected from a timing sample before the run (D-187, D-193). If the projection is too long
  for the laptop, arm I at T = 20 is dropped first, then the main arm's θ = 0.25 rules at T = 20, by an
  addendum before the run.

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

**Arm I.** Inside a narrow grade the PIT probit is a·s_g/k + a√ρ·(β·M_t − Z_t) + a√((1 − ρ)κ)·(β·N_i − η_i),
now with a = 1/√((1 − ρ)(1 − κ)), k = 1 − λ + λa and β = λa/k. The last term varies across the grade's
obligors and diversifies, which divides the factor variance by D = 1 + a²(1 − ρ)κ·V_I, with
V_I = Var(βN − η) = β²V_N − 2β·Cov(N, η) + 1:

    v = a²ρ·(1 − 2βθ_S + β²θ_S/(2 − θ_S)) / D,    ρ_eff = v/(1 + v).

With the observable news a larger share of the score, the grade is wide in the default driver, and the
narrow-grade limit is a poor guide for the fast hybrid rules (λ = 0.5, θ_S = 1: 0.074 against a finite-grade
0.046–0.055 on L). The targets are therefore the finite-grade values in `reference_idio.csv`, a quadrature
over (Z_t, M_{t−1}) and the population of (c, N), with η integrated out exactly given N. A Monte Carlo of
the grade construction (150,000 obligors, `--mc-check`) agrees with them on every rule and scale
within its own noise (about ±0.02 on L, judged by the TTC control, whose answer is known).

**Arm I, ρ_g by rule** (finite-grade, range over the informative grades):

| Rule (λ, θ_S, θ_I) | L (ρ 0.24) | M (ρ 0.12) | H (ρ 0.02) |
|---|---|---|---|
| TTC | 0.239–0.240 | 0.119–0.120 | 0.020 |
| 0.5, 1, 1 | 0.046–0.055 | 0.019–0.024 | 0.003–0.004 |
| 0.5, 1, 0.5 | 0.044–0.051 | 0.017–0.023 | 0.002–0.003 |
| 0.5, 1, 0.25 | 0.044–0.049 | 0.015–0.022 | 0.002–0.003 |
| 0.5, 0.25, 1 | 0.289 | 0.148–0.149 | 0.025 |
| 0.5, 0.25, 0.5 | 0.268–0.269 | 0.135–0.137 | 0.023 |
| 0.5, 0.25, 0.25 | 0.243–0.244 | 0.120–0.122 | 0.020 |
| 1, 1, any θ_I | below 0.001 | below 0.001 | below 0.001 |
| 1, 0.25, 1 | 0.288 | 0.149 | 0.026 |
| 1, 0.25, 0.5 | 0.266 | 0.135 | 0.023 |
| 1, 0.25, 0.25 | 0.234 | 0.116 | 0.019 |

**What the tables say before any run.** Momentum undoes most of what PIT-ness does to the fitted correlation:
a slow PIT rating (θ = 0.25) gives a ρ̂ of 65–70% of the asset correlation on all three scales, while a
fast one gives about none. So "PIT or TTC" is not enough to say what ρ a grade's history implies; the
rating's speed matters as much. Arm I adds two things:

- **Idiosyncratic momentum lowers the fitted correlation:** under slow systematic momentum, each step of
  θ_I (1 to 0.5 to 0.25) lowers ρ_g by 7–16%. A rating that lags obligor news leaves obligors with stale news in a grade;
  their spread in default risk is idiosyncratic and diversifies.
- **The fitted correlation can exceed the asset correlation:** a rating that sees obligor news promptly
  (θ_I = 1) but lags the cycle (θ_S = 0.25) gives a ρ_g about 1.2 times the generating ρ on every scale,
  whether the rule is hybrid or PIT. Sorting on the news makes the grade homogeneous in it, so the factor is
  a larger share of what is left.

**The grade's long-run PD against its label** (`pd_pseudo` / master-scale PD, informative grades): 1.01 under
TTC; without momentum, below the label for a hybrid rule (scale L at λ = 0.5: 0.80–0.89) and back at it under
PIT (0.98–1.00); with momentum, above it, most at slow PIT (scale L at λ = 1, θ = 0.25: 1.30–1.76). A lagging
rating keeps obligors in their grade after they have worsened. In arm I the gap is larger: 1.03–4.85 under
slow systematic momentum (`reference_idio.csv`).

## Predictions

All strict pass/fail: a miss is recorded as a miss. "Median" is over replicates. "The band" is the D-131
coverage band, 0.927–0.973 at R = 1,000. "Informative grades" are the 23 fixed above, less the excluded
grade–rule pairs; "pooled" is over a scale's informative grades. ρ_g and PD_g are the pseudo-true values of
the arm's reference CSV. Q1–Q9 cover both arms unless they say otherwise.

| # | Prediction |
|---|---|
| Q1 | **Check (stops scoring if it fails):** in every replicate, the portfolio-level fit (all grades) is bit-for-bit identical under the 13 rules of each arm. Defaults do not depend on the rule, only grade membership does |
| Q2 | **TTC recovers the asset correlation:** at T = 100, in each arm, in at least 20 of the 23 informative grades the ρ profile interval's coverage of the generating ρ is in the band (under TTC the pseudo-true value is the generating ρ to within 0.002) |
| Q3 | **Philosophy (Miu and Ozdemir), without momentum (main arm):** in every scenario and every informative grade, at θ = 1, median ρ̂ falls from λ = 0 to λ = 0.5, and is non-increasing across the five λ, ties at the box's floor allowed |
| Q4 | **Speed (main arm):** in every scenario at T = 100 and every informative grade, for each λ > 0, median ρ̂ rises as θ falls from 1 to 0.5 to 0.25, ties at the floor allowed |
| Q5 | **Size, at T = 100:** in every informative grade under every rule whose ρ_g is at least 0.01, median ρ̂ / ρ_g is between 0.85 and 1.05. Below 1 allows for ρ̂'s known small-sample downward bias (M1; S-3) |
| Q6 | **Fast PIT sits at the floor:** at T = 100, at λ = 1 and θ = 1 (main arm) and at λ = 1, θ_S = θ_I = 1 (arm I), pooled, the share of fits with ρ̂ at the box's lower bound (flagged) is between 0.35 and 0.65 (a pseudo-true value below the bound puts about half the estimates on it) |
| Q7 | **The interval no longer covers the asset correlation (main arm):** at T = 100, θ = 1 and λ ≥ 0.5, in every informative grade, the ρ profile interval covers the generating ρ in under half of the replicates |
| Q8 | **ρ interval, coverage of the pseudo-true value:** at T = 100, for every (scale, rule) of either arm whose ρ_g is at least 0.01 in every informative grade, the pooled coverage of ρ_g by the ρ profile interval is between 0.90 and 0.97. Momentum's negative autocorrelation inflates the variance of the sample variance by 1 + 2Σγ_j² (at most about 1.17 here), so it moves this coverage by at most about 0.02 |
| Q9 | **The grade's long-run PD, at T = 100:** in every informative grade under every rule, median PD̂ / PD_g is between 0.95 and 1.05 |
| Q10 | **PD interval, coverage of the pseudo-true value, at T = 100:** (a) without systematic momentum (θ = 1, or θ_S = 1 in arm I, where periods are independent), for every (scale, rule) whose ρ_g is at least 0.01 in every informative grade, the pooled coverage of PD_g by the PD profile interval is in the band; (b) with systematic momentum (θ ≤ 0.5 in the main arm, θ_S = 0.25 in arm I) and λ ≥ 0.5, on scales L and M, it is above 0.96 for every (scale, rule): the negative autocorrelation makes the mean default rate less variable than an iid fit assumes. Scale H is reported: there the binomial part of the variance, which does not shrink, dominates |
| Q11 | **The HAC interval, PD only, at T = 100:** (a) on Q10(b)'s (scale, rule) pairs, the pooled coverage of PD_g by the HAC interval for PD is between 0.90 and 0.985 and nearer 0.95 than the profile interval's: it sees part of the negative autocorrelation the iid likelihood misses (the Bartlett weights at lag 4 see only part of it); (b) under TTC in each arm, where periods are independent, its pooled PD coverage of the generating PD is at least the profile interval's less 0.03 and at most the profile interval's plus 0.01 (the cost of estimating S with T = 100 periods). The HAC interval for ρ is reported, not predicted: S-48 found it covers worse than the iid Wald interval (its run, #66) |
| Q12 | **Arm I, idiosyncratic momentum lowers the fitted correlation:** at T = 100, at θ_S = 0.25 and each λ ∈ {0.5, 1}, in every scale and every informative grade, median ρ̂ falls as θ_I falls from 1 to 0.5 to 0.25 (ρ_g falls by at least 5% at each step) |
| Q13 | **Arm I, above the asset correlation:** at T = 100, at θ_S = 0.25, θ_I = 1 and each λ ∈ {0.5, 1}, in every scale and every informative grade, median ρ̂ is above the generating ρ (ρ_g is 1.20–1.28 times it); and pooled, the ρ profile interval covers the generating ρ in under 0.90 of the fits |

**Reported without a prediction:** everything at T = 20 and 40 (where ρ̂'s small-T bias and the floor
dominate); the edge and thin grades and the excluded grade–rule pairs; the share of empty grade-periods by
rule; the generating-ρ distance per grade; the realised autocorrelation of each grade's default rates
against the closed form; arm I at θ_S = 1 with θ_I < 1 (where a step of θ_I moves ρ_g by −3% to +16%, not
always in one direction); the HAC interval at T = 20 and 40, its width against the profile interval's,
its ρ coverage under every rule, and its PD coverage on scale H; and, for capital, each grade's 99.9%
conditional PD from its fitted (PD̂, ρ̂) beside the one from the generating ρ and the master-scale PD
(S-23's interval), the practical consequence of reading a grade's ρ̂ as an asset correlation.

## Verdicts it could change

None of the pinned verdicts: it is a new DGP with its own reference values.

## Monitoring implication

To be filled in when the study finishes (D-152). The candidates: a check that the correlation a calibration
uses matches both the rating philosophy and the rating speed of the grades it is applied to, and a check of a
grade's realised default rate against its label, which a lagging rating moves.

## Registration

Approved by the owner on 2026-10-12 as it stands (D-401): the two arms (the smoothed-factor rule for
systematic momentum; an observable AR(1) driver with κ = 0.5 and ψ = 0.8 for idiosyncratic momentum;
obligors redrawn each period), the informative-grade rule, Q1–Q13 with their bands, and the HAC interval
with S-48's settings (D-209). Before the run: the generator's D-149 artifacts (prose specification, Python
mirror, reference-panel hash), the tool and scorer merged, and the runtime projection (D-187, D-193).
