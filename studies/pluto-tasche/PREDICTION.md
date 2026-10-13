# S-6 Pluto–Tasche most-prudent upper bounds: predictions made before the run

Pre-registration for study S-6 (`studies/README.md`, D-148, D-214). It is committed **before** any panel
is simulated for the study (D-148). The results will be compared with it prediction by prediction, and
misses will be reported, not explained away.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a new,
dated section at the end, committed before the run it concerns. The serially correlated version (an
AR(1) factor) is out of scope and stays with M6.

## Evidence used, and what was seen before this was written

- **For the numbers:** the model only, through [`reference.py`](reference.py) (scipy, no simulation),
  which writes [`reference.csv`](reference.csv). Every coverage and median ratio in it is exact, up to
  quadrature (checked against `scipy.integrate.quad` at 1.6e-13 relative) and the root solve (logit
  tolerance 1e-13).
- **Seen, and published:** the pinned recovery summary (`tests/golden/recovery/summary.csv`), for the
  scenario definitions and as context for the comparators: in these 27 scenarios the profile interval's
  two-sided 95% PD coverage is 0.930–0.986. It is not used for any prediction below.
- **No panel of this study's seed** has been drawn.

## The question

How conservative are Pluto and Tasche's most-prudent upper bounds in low-default portfolios, and do they
cover the true PD at their stated confidence? Does the one-factor version fix the independent version's
undercoverage, and what does a wrong assumed correlation cost?

## What is run

### The bounds

Per panel, at confidence γ ∈ {0.5, 0.75, 0.9, 0.95, 0.99}, with N = Σ n_t and D = Σ d_t:

- **(a) Independent:** the most prudent PD p with P(Bin(N, p) ≤ D) = 1 − γ, the one-sided
  Clopper–Pearson bound. At D = 0 it is 1 − (1 − γ)^(1/N).
- **(b) One factor, independent period factors:** the p with P(Σ_t D_t ≤ D | p, ρ_a) = 1 − γ.
  - Each period's count is the Vasicek mixed binomial. Its PMF on 0..D comes from the existing
    quadrature (`log_contrib`, with `lbinom`), entirely in log space. The T per-period PMFs are convolved,
    on counts 0..D only.
  - ρ_a is assumed, not estimated. There are two arms:
    - **oracle:** the generating ρ;
    - **fixed:** ρ_a = 0.12. It is misspecified on the 0.02 and 0.24 scales and equals the oracle on the
      0.12 scale.
- **The solve:** bisection on logit p for both bounds.
  - P(S ≤ D | p) falls as p rises, because the conditional PD p(z) rises with p at every z. So the root
    is unique and the bound rises with D.
  - The tolerance, planned as 1e-10 in logit p, is registered in `docs/methodology/tolerances.md` with
    the code and its tests (PR 2). The register's sync check needs the id in `tests/tolerances.toml`,
    and that comes with the code.
  - The reference shows the closest any bound comes to a true PD, at the D where coverage switches, is
    8.4e-5 in logit. That is far above the tolerance, so solver error cannot flip a coverage flag.
- **Both bounds depend on the panel only through D**, because n is constant over the periods. So their
  exact coverage is a sum over D (below), and on each scenario's panels the bound takes only as many
  distinct values as the panels have distinct totals.

### The comparators, on the same panels

- **The profile upper limit:** the binomial-mixture MLE's PD profile upper limit, one-sided at γ.
  - It is the PD where the signed root of the profile likelihood ratio equals Φ⁻¹(γ), i.e. the
    threshold (Φ⁻¹(γ))²/2 on the upper side. At γ = 0.5 it is the polished MLE.
  - Conventions are D-131's: a truncated end is still an end, and a not-computed interval counts as not
    covering.
- **The Jeffreys upper limit:** the γ quantile of the Jeffreys posterior's PD marginal, with the
  resolution rule on (D-176).

### The scenarios

- **The cells:** the recovery grid's low-default cells (`tests/recovery/recovery.hpp`): expected defaults
  per period n·PD ≤ 1. That gives 27 of the 81 scenarios, at every ρ ∈ {0.02, 0.12, 0.24} and every
  T ∈ {20, 40, 100}:
  - **n·PD = 0.1** (PD 0.1%, n = 100), 9 scenarios: 0, 3, 6, 9, 12, 15, 18, 21, 24;
  - **n·PD = 1**, 18 scenarios:
    - PD 0.1%, n = 1,000: 1, 4, 7, 10, 13, 16, 19, 22, 25;
    - PD 1%, n = 100: 27, 30, 33, 36, 39, 42, 45, 48, 51.
- **The subset's PD 0.1% cells,** scenarios 4 and 7 (D-150), are among the n·PD = 1 cells, so the list
  already holds them.
- **Replicates:** R = 1,000 per scenario. Replicate r of scenario id is the DGP panel (seed, id, r, PD,
  ρ, n, T) of `dgp/` (D-109), with n constant over the T periods. The panels are drawn from the base DGP
  (`dgp::simulate_panel`).
- **The seed:** `0x533036504C55544F` ("S06PLUTO"), new, reserved by D-214. It collides with no seed in
  the code or the decisions. The recovery seed is not reused, so the panels are independent of the
  pinned run's.

### What is recorded, per panel and γ

- D, and the bounds (a), (b) oracle and (b) fixed;
- the profile and Jeffreys upper limits, with their flags;
- whether each is ≥ the true PD.

### Monte Carlo bands

- **Coverage:** around the exact coverage c, the band c ± 3.29·√(c(1 − c)/1,000) (D-131's two-sided 0.1%
  band, centred on c rather than 0.95).
  - Where c = 1 (D_min = 0, so every panel's bound is at or above the truth), the band is the point 1.
  - This is the case in 36 of the 405 cells.
- **Median bound/true ratio:** the bound is increasing in D, so its median over panels is bound(D_med),
  where D_med is the median of the true total. The simulated median of 1,000 lies between bound(D_lo)
  and bound(D_hi) with probability about 0.999. Here D_lo and D_hi are the true total's quantiles at
  1/2 ∓ 3.29·√(0.25/1,000). Their ratios are `ratio_lo` and `ratio_hi` in `reference.csv`.

## Exact coverage, without simulation

bound(D) ≥ PD exactly when P_a(S ≤ D | PD) ≥ 1 − γ, where P_a is the bound's assumed law: binomial for
(a), the mixed binomial at ρ_a for (b). So:

- the bound covers exactly when D ≥ D_min, the smallest D with P_a(S ≤ D | PD) ≥ 1 − γ;
- coverage = P_true(S ≥ D_min) = 1 − F_true(D_min − 1), where F_true is the CDF of the true total (the
  T-fold convolution at the generating PD and ρ).

`reference.csv` holds, for each scenario × method × γ:
- D_min and the exact coverage, with the band;
- the bound's distance from the truth in logit at D_min and D_min − 1;
- D_med and the exact median ratio, with `ratio_lo` and `ratio_hi`.

**A note on the brief's wording.** The task described coverage as summing up to "the largest D whose bound
is still ≥ the true PD". The bound rises with D, so the covering totals are the large ones, D ≥ D_min, and
the sum runs over the non-covering totals 0..D_min − 1. The identity above is what `reference.py` uses.

**For the oracle, coverage ≥ γ holds by construction:** P_a = P_true, so F_true(D_min − 1) < 1 − γ. The
excess over γ is discreteness.

## Headline exact numbers (`reference.csv`)

Exact coverage of (a), (b) oracle and (b) fixed (ρ_a = 0.12). Each entry is the mean over the scenarios
of the tier (6 for n·PD = 1, 3 for n·PD = 0.1, per ρ), with the range in brackets. The last column is the
median over those scenarios of the exact median bound/true ratio, for (a) / oracle / fixed.

| n·PD | ρ | γ | (a) independent | (b) oracle | (b) ρ_a = 0.12 | median ratio |
|---|---:|---:|---|---|---|---|
| 1 | 0.02 | 0.5 | 0.516 (0.509–0.524) | 0.516 (0.509–0.524) | 0.574 (0.544–0.608) | 1.02 / 1.02 / 1.04 |
| 1 | 0.02 | 0.75 | 0.748 (0.731–0.762) | 0.768 (0.750–0.790) | 0.879 (0.837–0.911) | 1.13 / 1.14 / 1.23 |
| 1 | 0.02 | 0.9 | 0.902 (0.889–0.921) | 0.912 (0.900–0.921) | 0.975 (0.962–0.986) | 1.23 / 1.26 / 1.42 |
| 1 | 0.02 | 0.95 | 0.943 (0.934–0.951) | 0.957 (0.951–0.966) | 0.993 (0.989–0.997) | 1.30 / 1.34 / 1.55 |
| 1 | 0.02 | 0.99 | 0.988 (0.982–0.993) | 0.992 (0.990–0.993) | 1.000 (0.999–1.000) | 1.43 / 1.49 / 1.81 |
| 1 | 0.12 | 0.5 | 0.484 (0.465–0.497) | 0.525 (0.504–0.556) | = oracle | 0.99 / 1.01 / 1.01 |
| 1 | 0.12 | 0.75 | 0.652 (0.618–0.685) | 0.773 (0.757–0.795) | = oracle | 1.10 / 1.20 / 1.20 |
| 1 | 0.12 | 0.9 | 0.797 (0.757–0.844) | 0.909 (0.902–0.920) | = oracle | 1.21 / 1.39 / 1.39 |
| 1 | 0.12 | 0.95 | 0.850 (0.814–0.886) | 0.958 (0.954–0.967) | = oracle | 1.27 / 1.52 / 1.52 |
| 1 | 0.12 | 0.99 | 0.936 (0.899–0.967) | 0.992 (0.990–0.995) | = oracle | 1.40 / 1.78 / 1.78 |
| 1 | 0.24 | 0.5 | 0.447 (0.403–0.481) | 0.506 (0.500–0.517) | 0.475 (0.436–0.505) | 0.97 / 1.02 / 0.99 |
| 1 | 0.24 | 0.75 | 0.566 (0.510–0.621) | 0.764 (0.752–0.795) | 0.661 (0.620–0.701) | 1.05 / 1.30 / 1.15 |
| 1 | 0.24 | 0.9 | 0.683 (0.620–0.749) | 0.909 (0.906–0.917) | 0.798 (0.765–0.837) | 1.15 / 1.59 / 1.34 |
| 1 | 0.24 | 0.95 | 0.733 (0.673–0.795) | 0.956 (0.950–0.966) | 0.869 (0.837–0.908) | 1.23 / 1.79 / 1.46 |
| 1 | 0.24 | 0.99 | 0.836 (0.761–0.908) | 0.992 (0.991–0.993) | 0.950 (0.923–0.974) | 1.36 / 2.22 / 1.71 |
| 0.1 | 0.02 | 0.5 | 0.565 (0.541–0.591) | 0.565 (0.541–0.591) | 0.565 (0.541–0.591) | 1.17 / 1.17 / 1.19 |
| 0.1 | 0.02 | 0.75 | 0.799 (0.758–0.861) | 0.799 (0.758–0.861) | 0.848 (0.777–0.906) | 1.57 / 1.58 / 1.66 |
| 0.1 | 0.02 | 0.9 | 0.946 (0.906–1.000) | 0.946 (0.906–1.000) | 0.970 (0.931–1.000) | 2.00 / 2.02 / 2.18 |
| 0.1 | 0.02 | 0.95 | 0.983 (0.969–1.000) | 0.983 (0.969–1.000) | 0.983 (0.969–1.000) | 2.29 / 2.32 / 2.55 |
| 0.1 | 0.02 | 0.99 | 0.999 (0.997–1.000) | 0.999 (0.997–1.000) | 0.999 (0.997–1.000) | 2.90 / 2.95 / 3.35 |
| 0.1 | 0.12 | 0.5 | 0.546 (0.528–0.565) | 0.546 (0.528–0.565) | = oracle | 1.17 / 1.19 / 1.19 |
| 0.1 | 0.12 | 0.75 | 0.773 (0.731–0.836) | 0.824 (0.751–0.884) | = oracle | 1.57 / 1.66 / 1.66 |
| 0.1 | 0.12 | 0.9 | 0.932 (0.884–1.000) | 0.962 (0.912–1.000) | = oracle | 2.00 / 2.18 / 2.18 |
| 0.1 | 0.12 | 0.95 | 0.977 (0.957–1.000) | 0.977 (0.957–1.000) | = oracle | 2.29 / 2.55 / 2.55 |
| 0.1 | 0.12 | 0.99 | 0.998 (0.995–1.000) | 0.998 (0.995–1.000) | = oracle | 2.90 / 3.35 / 3.35 |
| 0.1 | 0.24 | 0.5 | 0.507 (0.499–0.516) | 0.540 (0.505–0.599) | 0.507 (0.499–0.516) | 1.17 / 1.26 / 1.19 |
| 0.1 | 0.24 | 0.75 | 0.720 (0.675–0.785) | 0.804 (0.785–0.838) | 0.774 (0.698–0.838) | 1.57 / 1.84 / 1.66 |
| 0.1 | 0.24 | 0.9 | 0.902 (0.838–1.000) | 0.960 (0.927–1.000) | 0.940 (0.868–1.000) | 2.00 / 2.53 / 2.18 |
| 0.1 | 0.24 | 0.95 | 0.960 (0.927–1.000) | 0.973 (0.954–1.000) | 0.960 (0.927–1.000) | 2.29 / 3.02 / 2.55 |
| 0.1 | 0.24 | 0.99 | 0.996 (0.987–1.000) | 0.999 (0.997–1.000) | 0.996 (0.987–1.000) | 2.90 / 4.13 / 3.35 |

**What the reference shows:**
- **At n·PD = 0.1, discreteness dominates.** D is 0, 1 or 2 in most panels, so each bound jumps between
  a few values. Coverage sits well above γ for every method, and the independent bound's undercoverage
  is mostly masked. So coverage there is reported, not predicted, and the bound/true ratio is the
  headline: at γ = 0.95 the median bound is 2.3–3.0 times the true PD.
- **At n·PD = 1, the independent bound's undercoverage is clear** once ρ ≥ 0.12.

## Mechanism behind the predictions

- **Mixing fattens the lower tail of the total.** Under a common period factor the total S is
  over-dispersed relative to Bin(N, p). Its lower tail is heavier, so at a fixed p a small D is more
  probable than the binomial says.
- **So (b) rises with ρ_a.** Solving P(S ≤ D | p) = 1 − γ under (b) needs a larger p than under (a), and
  a larger one still as ρ_a rises.
- **The oracle is exact, so its coverage is ≥ γ.** It inverts the true law of the statistic it uses.
- **(a), and (b) with ρ_a below the truth, invert a law that is too narrow.** They undercover by an
  amount that grows with the true ρ and with the information per period, n·PD.

## Predictions

Here "tier 1" means the 18 n·PD = 1 scenarios. Coverage is predicted in tier 1 only; at n·PD = 0.1 the ratio is predicted and coverage is reported. "Pooled" coverage is covered panels over all panels of the
cells named. "Paired" means on the same panels, where a dominance in D_min makes the comparison exact, with
no Monte Carlo noise.

| # | Prediction |
|---|---|
| P1 | **The oracle covers at ≥ γ.** Exact: coverage ≥ γ in all 135 oracle cells, by construction. Simulated, in tier 1: in each of the 90 tier-1 oracle cells (18 scenarios × 5 γ), coverage ≥ γ − 3.29·√(γ(1 − γ)/1,000) |
| P2 | **The independent bound undercovers once ρ ≥ 0.12, by the amounts below.** In the tier-1 scenarios at ρ ∈ {0.12, 0.24} (12 scenarios × 5 γ = 60 cells), the exact coverage is below γ in all 60. **(i)** The pooled simulated shortfall γ − coverage, per (ρ, γ), is within ±0.02 of the exact: at ρ = 0.12, 0.016 / 0.098 / 0.103 / 0.100 / 0.054; at ρ = 0.24, 0.053 / 0.184 / 0.217 / 0.217 / 0.154 (γ = 0.5 / 0.75 / 0.9 / 0.95 / 0.99). **(ii)** The simulated coverage is below γ in at least 55 of the 60 cells. The exact law expects 58.5, with 55 its 0.1% quantile, taking cells as independent |
| P3 | **How the ρ_a = 0.12 arm errs on each scale (tier 1).** **(i) ρ = 0.02 (ρ_a too high):** in the 30 tier-1 cells its D_min is ≤ the oracle's. So, paired, it covers every panel the oracle covers, and its simulated coverage is ≥ the oracle's in all 30. Its pooled excess over the oracle, per γ, is within ±0.02 of the exact 0.059 / 0.111 / 0.063 / 0.036 / 0.008. **(ii) ρ = 0.24 (ρ_a too low):** in the 30 tier-1 cells its D_min is ≥ the oracle's, so its simulated coverage is ≤ the oracle's in all 30. Its exact coverage is below γ in 27 of the 30. The pooled shortfall γ − coverage, per γ, is within ±0.02 of the exact 0.025 / 0.089 / 0.102 / 0.081 / 0.040, and the simulated coverage is below γ in at least 26 of the 30 (expected 28.3). **(iii) ρ = 0.12:** identical to the oracle, bound for bound, on every panel of all 9 ρ = 0.12 scenarios |
| P4 | **Each method's simulated coverage matches its exact value (tier 1).** In at least 268 of the 270 tier-1 cells (18 scenarios × 3 methods × 5 γ), it lies inside the band around the exact coverage. Expected out of band by chance: about 0.3 |
| P5 | **The bound/true ratio by γ: the headline at n·PD = 0.1.** **(i)** In at least 403 of all 405 cells, both tiers, the simulated median ratio lies in [`ratio_lo`, `ratio_hi`]. **(ii)** The exact medians rise with γ in every scenario and method, as in the table's last column. **(iii) Tier 1:** at γ = 0.95 the oracle's median ratio is 1.34, 1.52 and 1.79 (the median over scenarios) at ρ = 0.02, 0.12 and 0.24, against the independent bound's 1.30, 1.27 and 1.23. **(iv) n·PD = 0.1:** at γ = 0.95 the oracle's median ratio is 2.32, 2.55 and 3.02, against the independent bound's 2.29 at every ρ. Coverage in these cells is driven by small totals (P(D = 0) up to 0.22 at T = 20) and is reported, not predicted. At γ = 0.99 the oracle's median ratio reaches 6.81 there (scenario 18) |
| P6 | **The ordering of the Pluto–Tasche bounds.** On every panel, at every γ: independent ≤ (b, ρ_a = 0.02) ≤ (b, 0.12) ≤ (b, 0.24). The reference checks every pair of these laws at every total D from 0 to each (n, T)'s largest 0.9999 quantile of the true total, at all five γ, γ = 0.5 included: **0 violations and 0 ties in 32,700 pairs** (`reference_ordering.csv`, which also gives each pair's closest approach). The closest is 4.3e-4 in logit, independent against ρ_a = 0.02 at γ = 0.5 (n = 100, T = 100, D = 197): 4 million times the solver's planned tolerance. At γ = 0.5 the gaps narrow as D grows and are smallest at the top of the checked range, so a total beyond the 0.9999 quantile (probability below 1e-4 per panel) is not covered by the check; the scorer checks the order on every simulated panel directly. So no exception narrows P6. Paired on the same panels, with no Monte Carlo noise: **(i)** in the 90 tier-1 cells per method pair, the simulated coverage follows the same order (the oracle's ≥ the independent bound's; the fixed arm's ≥ the oracle's at ρ = 0.02 and ≤ it at ρ = 0.24); **(ii)** in all 135 cells per method pair, both tiers, the simulated median ratio follows it. Held if no cell breaks the order |

**Held** means held as stated. A prediction that holds in direction but misses in size is not held.

**Reported, not predicted:**

- **Coverage in the 9 n·PD = 0.1 scenarios,** for every method. It is driven by small totals, mostly
  D = 0, so it is reported beside `reference.csv`'s exact values, without a verdict. The bound/true ratio
  is the headline there (P5).
- **The profile and Jeffreys upper limits:** their coverage, median ratio and ordering against the
  Pluto–Tasche bounds, per scenario and γ.
  - They depend on the whole count vector, not only on D. Their exact coverage would need a sum over count
    vectors, which the reference cannot do in these cells, and no pinned result gives them at one-sided γ.
  - Their ordering against the bounds is therefore reported, paired on the same panels, with a paired z
    as in D-176.
- **The not-computed and truncated profile ends,** by scenario.
- **The bounds' spread:** the 5th and 95th percentiles of bound/true, per cell.
- **Group-level summaries** by T, which the predictions do not use.

## How the comparison will be reported

- **One row per prediction:** the prediction, the result, held or NOT HELD. The rows come from
  `score.py`, which reads only `rows.csv` and `reference.csv` and is merged before the run (PR 2).
- **Where the results go:** `studies/README.md`, a D-entry and `MANIFEST.json`. The per-panel rows are
  attached to a release, not committed (D-163).
- **Runtime:** projected before the run (D-187, D-193); expected minutes.
- **Verdicts it could change:** none pinned. The bounds are a `native` benchmark, not parity and not in
  the C ABI.
- **Monitoring implication (D-152):** to be confirmed by the result. If P2 and P3 hold, an independent
  Pluto–Tasche bound in a low-default portfolio understates the PD a stated confidence requires once
  correlation is material. Its coverage at γ = 0.95 in tier 1 is 0.85 at ρ = 0.12 and 0.73 at ρ = 0.24.
  The one-factor bound restores coverage only if ρ_a is not below the truth, so a conservative ρ_a is
  the safe choice.
