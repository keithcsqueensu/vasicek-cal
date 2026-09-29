# S-4 Bartlett-corrected profile threshold (S-4a oracle, S-4b usable): predictions made before the run

Pre-registration for study S-4 (`studies/README.md`, D-148, D-158, D-168), both parts. This file is
committed **before** any likelihood-ratio statistic at the true value or any Bartlett factor has been
computed on a recovery panel. The results will be compared with it prediction by prediction, and misses
will be reported, not explained away. S-4 shares its run with S-10
([`../parametric-bootstrap/PREDICTION.md`](../parametric-bootstrap/PREDICTION.md)): S-4b's factor is
computed from S-10's bootstrap panels.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a new,
dated section appended at the end, committed before the run it concerns.

**Evidence used in writing it:** pinned results only, through
[`../parametric-bootstrap/reference.py`](../parametric-bootstrap/reference.py): the profile coverage of
the recovery summary and the Bartlett factor it implies, S-13's R = 10,000 coverage (D-158), and S-9's
Jeffreys coverage on the same replicates (D-175).

## The question

The profile interval is {ψ : W(ψ) ≤ c}, W the profile likelihood-ratio statistic and c = 3.8415, the
χ²₁ 0.95 quantile. At T = 20 it covers about 94% (S-13: a systematic shortfall of about one point in
group B, for PD, ρ and q; not six isolated findings). A Bartlett correction rescales the threshold to
c·k, k = E[W] at the truth, so that W/k is closer to χ²₁.

- **S-4a (oracle):** k from the truth, which no user has. It bounds what any Bartlett correction can do.
  **It is in-sample by design:** k_a is the mean of W₀ over the very replicates whose coverage it then
  corrects, so its coverage is an upper bound for S-4b, not a test of anything. H1 and H2 describe
  that bound; the tests are H3–H6.
- **S-4b (usable):** k estimated at θ̂ by parametric bootstrap.
- **Against Jeffreys:** the Jeffreys equal-tailed interval already closes part of the small-T gap at no
  extra width (D-175). If it closes as much as S-4b does, a Bartlett-corrected profile interval is not
  worth its machinery (a bootstrap of profile maximisations per fit).

## What is run

In S-10's run, on the recovery panels (81 × 1,000), for each panel and each parameter ψ ∈ {PD, ρ}:

- **W at the truth,** W₀ = 2(l_max − P_ψ(ψ₀)), P_ψ the profile log-likelihood (the inner maximum over the
  other parameter on the box, as in the profile code). The pinned profile interval covers exactly when
  W₀ ≤ c; the run checks that this reproduces the pinned coverage replicate by replicate and reports
  any replicate where it does not (an interval not computed, D-131).
- **S-4a:** k_a(scenario, ψ) = the mean of W₀ over the scenario's 1,000 replicates. Its corrected
  interval covers when W₀ ≤ c·k_a. This uses the same replicates that define the factor, which is part
  of why it is an oracle.
- **S-4b:** on the first 199 of the panel's S-10 bootstrap panels (generated at θ̂), W*_b =
  2(l*_max − P*_ψ(ψ̂)), the statistic at the value the bootstrap panel was generated from; k_b = mean of
  W*_b. Its corrected interval covers when W₀ ≤ c·k_b. Where ψ̂ is on the box's bound, or more than 10%
  of the W*_b are not finite, k_b = 1 (no correction), counted and reported.
- **Widths:** the corrected intervals are solved (the profile code with threshold c·k/2 in log-likelihood
  units) on replicates 0–199 of every scenario, for median logit widths and their ratio to the profile
  interval's; coverage uses all 1,000 replicates through W₀.
- **The verdicts it is judged against:** the pinned profile verdicts, with the group B shortfall as the
  target (S-13), in particular the six small-T findings (29 PD, 29 ρ, 55 ρ, 68 PD, 72 ρ, 74 ρ).

## Mechanism behind the predictions

- **The oracle factor.** If W at the truth were exactly k·χ²₁, the profile interval would cover with
  probability F(c/k). Inverting the pinned pooled coverage gives the factor an oracle would apply
  (`reference.py`): k ≈ 1.11 for ρ and 1.09 for PD at T = 20, 1.03–1.05 at T = 40, 1.02 at T = 100. W is
  not exactly a scaled χ²₁ (skew at T = 20 makes one tail heavier), so matching the mean fixes most of
  the shortfall, not all: the oracle should lift group B to about 0.945–0.952, not exactly 0.95.
- **The usable factor loses some of that.** k_b is evaluated at θ̂, not the truth. ρ̂ is biased low at
  T = 20, and the factor grows as ρ falls towards the floor (the likelihood is more skewed there), so
  k_b should on average be a little larger than k_a where ρ is small and smaller where ρ is large; with
  199 bootstrap panels its Monte Carlo spread is about ±0.1. Averaged over replicates, noise in k
  costs little coverage (the coverage is nearly linear in k near 1.1), so S-4b should keep most of the
  oracle's gain.
- **Against Jeffreys:** Jeffreys lifted group B ρ from 0.9371 to 0.9430 without widening the interval.
  A Bartlett correction lifts it by widening: raising the threshold by 11% widens a symmetric interval
  by about √1.11 − 1 ≈ 5%.

## Predictions

"Gain" is the pooled coverage of the corrected interval minus the profile interval's, over group B
(12 scenarios × 1,000), per parameter.

| # | Prediction |
|---|---|
| H1 | **The oracle bound (in-sample; not a test):** S-4a's pooled group B coverage is 0.944–0.955 for both PD and ρ, and all six small-T findings are in the band |
| H2 | **The oracle factor:** the median k_a over group B is 1.05–1.16 for ρ and 1.03–1.14 for PD; over group D it is 0.98–1.05 for both |
| H3 | **The usable factor keeps most of it:** S-4b's gain is at least half of S-4a's, for ρ and for PD; at least 4 of the six small-T findings are in the band under S-4b |
| H4 | **No harm where the profile works:** in groups C and D, S-4b's pooled coverage is within 0.005 of the profile's for each parameter, and at most 2 of the 76 C–D verdicts that are in the band under the profile leave it |
| H5 | **Against Jeffreys, coverage:** in group B, S-4b's pooled ρ coverage is within ±0.006 of Jeffreys' (0.9430) |
| H6 | **Against Jeffreys, width:** in group B, S-4b's median ρ width ratio to the profile interval is 1.02–1.08 in at least 10 of 12 scenarios. (Jeffreys' is already known: 0.994–1.015, median 1.008, D-175.) So where H5 holds, Jeffreys reaches the same coverage at smaller width |

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not; a prediction that holds in direction but
misses in size is not held. The rows are computed by a script committed before the full run's results
exist. The results go into the interval comparison in `studies/README.md` beside S-10's, and answer the
question stated above: whether a Bartlett-corrected profile interval is worth adding given the Jeffreys
interval. With a D-entry and a `MANIFEST.json`.
