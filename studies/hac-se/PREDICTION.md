# S-48 HAC standard errors under serial dependence: predictions made before the run

Pre-registration for study S-48 (`studies/README.md`, D-206, D-209), the positive-dependence and
independence parts. It is committed **before** any HAC SE is computed on a study panel. The results will
be compared with it prediction by prediction, and misses will be reported, not explained away.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a new,
dated section at the end, committed before the run it concerns. The negative-dependence part (S-36's
momentum rules) is registered with S-36, and is not part of this file.

## Evidence used, and what was seen before this was written

- **For the numbers:** the model and pinned results only, through [`reference.py`](reference.py), which
  writes `reference.csv`. It uses S-35's per-period score structure (`studies/lrpd-serial-correlation/
  reference.py`) and the pinned φ = 0 Wald coverage (`tests/golden/recovery/summary.csv`). No panel was
  simulated or fitted for it.
- **Seen, and published:**
  - S-35's results at every φ, `studies/lrpd-serial-correlation/summary.csv`: the iid Wald and profile
    coverage on the very panels this study refits;
  - the M6 arm's, `summary_m6.csv`: the AR(1) MLE's coverage on the same panels.
- **A calibration check against them, disclosed.** The reference's coverage mapping, applied to the iid
  Wald interval, reproduces S-35's pooled iid Wald coverage, pooled the same way:
  - at φ = 0.3, within 0.002 in groups B, C and D, for PD and for ρ;
  - at φ = 0.6, within 0.007 for ρ in every group and for PD in C and D, but **0.027 too high for PD
    in B** (0.680 against 0.653).

  The same mapping gives this file's HAC numbers, so group B's PD number at φ = 0.6 may carry the same
  bias. The predictions were not adjusted for it.
- **The HAC SE itself** was computed only on its test panels: the M1.4 panel and the six M1.7 validation
  panels (`tests/unit/hac_test.cpp`, `validation/scipy/hac_se.py`). None is a study panel.

## The question

When periods are dependent, does a Wald interval with the HAC sandwich SE recover the coverage that the
iid intervals lose? How does it compare with the model-based fix, the AR(1)-aware MLE? And what does it
cost when the periods are independent?

## What is run

- **The panels:** S-35's, bit for bit. These are the recovery matrix (81 scenarios × 1,000 replicates),
  the recovery seed, scenario ids and replicates, under the AR(1) factor at φ ∈ {0, 0.3, 0.6}
  (`dgp::simulate_panel_ar1`). φ = 0 is the pinned recovery run's panels.
- **The fit:** the parity binomial-mixture MLE, fitted exactly as S-35 fitted it (`recovery::fit_counts`
  with the profile arm): estimate, iid SEs, iid Wald interval, profile intervals.
- **The HAC SE** (`engine/hac.hpp`, D-209, `docs/methodology/hac_se.md`), at the estimate:
  - V = H⁻¹ S H⁻¹ in logit units, with S the Bartlett-weighted long-run variance of the centred
    per-period scores;
  - lag L = ⌊4 (T/100)^(2/9)⌋: 2, 3 and 4 at T = 20, 40 and 100. It is also computed at L = 1 and at
    twice the rule, which are reported and not predicted;
  - the 95% Wald interval û ± 1.96 SE in logit units;
  - none on fits flagged at or near a bound (`kFlagGridEdge`, `kFlagFlatSurface`, `kFlagNearBound`),
    which have no iid Wald interval either.
- **Coverage** is the HAC Wald interval's, conditional on the fits that have an iid Wald interval, as
  the recovery harness reports the Wald interval's. A HAC SE missing on such a fit (H not negative
  definite, or V not positive) counts as a miss and is reported. Pooled coverage is covered fits over
  unflagged fits, summed over a group's scenarios.
- **The comparators, on the same panels:**
  - the iid Wald and profile intervals, from this run's own fits;
  - the AR(1)-aware MLE's profile interval, as published by S-35's M6 arm (`summary_m6.csv`, D-208).
    It is not refitted.
- **The tool and scorer** (`study_hac_se`, `score.py` with a self-test on a synthetic fixture) are
  committed and merged before the run, and the runtime projection (D-187) is recorded before it starts.
- **Checks before scoring:**
  1. **The fits are S-35's.** Per scenario and φ, the flagged count, and the iid Wald and profile
     coverage of PD and ρ, equal S-35's `summary.csv` to its printed digits. If they do not, the
     scoring stops.
  2. **The HAC SE is missing only where it should be.** Its count on unflagged fits is reported, and
     every such fit is named.

## Mechanism behind the predictions

- **The estimate's error is a mean of dependent scores.** To first order, the error in logit units is
  (T I)⁻¹ Σ_t s_t, so for each axis it is a sum of a scalar series x_t with a Toeplitz covariance Σ.
  S-35 derived Σ from the Hermite expansion of E[s | Z] under the AR(1) factor. The iid SE assumes Σ is
  diagonal and misses the factor f = 1′Σ1 / (T Σ₀₀), S-35's inflation.
- **The HAC SE catches only part of it at these T.** Its expectation is E[Q] = tr(M K M Σ), with M the
  centring and K the Bartlett weights. Two losses remain, g = 1′Σ1 / E[Q]:
  - **Truncation:** lags beyond L are dropped and the kept ones are down-weighted. At φ = 0.6 and L = 2,
    the Bartlett sum of the AR(1) autocorrelations is 2.04, against the full 4.
  - **Centring:** estimating the mean removes about f/T of the variance, the most at T = 20.
- **Its own variability matters at T = 20.** Its effective degrees of freedom are ν = 2 E[Q]² / Var[Q]:
  about 9, 14 and 29 at T = 20, 40 and 100 under independence. Dependence lowers them: at φ = 0.6, the
  median over a group's scenarios in `reference.csv` (`pd_nu`) is 5.4, 9.2 and 20.1 for PD, and
  (`rho_nu`) 6.9, 11.6 and 24.5 for ρ. The Wald interval with a 1.96 critical value then acts like a t
  interval with too few degrees of freedom.
- **Under independence it is pure cost.** g = 1.17, 1.11 and 1.05 at T = 20, 40 and 100 (the centring
  and the Bartlett weights on lags whose true autocovariance is 0), on top of the lost degrees of
  freedom.
- **The mapping**, a heuristic held loosely as S-35's was. z₀ is the iid Wald's effective critical value
  from its pinned φ = 0 coverage. The HAC interval's coverage is F_ν(z₀/√g − δ) − F_ν(−z₀/√g − δ), with
  F_ν Student's t, and δ S-35's extra downward bias of ρ̂ (ρ only). Pooled per group with the pinned
  unflagged counts as weights.
- **The expected answer.** The HAC SE helps where the history is long enough to estimate S: at T = 100
  and φ = 0.6, PD coverage about 0.86 against the iid Wald's 0.72. At T = 20 it barely helps (0.71 against
  0.68) and costs about 0.06 under independence. It stays short of nominal everywhere, while the AR(1)
  MLE's profile interval was 0.91–0.94 at φ = 0.6 (D-208).

## Predictions

Groups B, C and D are the informative scenarios at T = 20, 40 and 100. Group A is defined as in S-35 and
is reported only. "Paired" compares two intervals on the same unflagged fits. Every reference number is
from `reference.csv`, pooled as above.

| # | Prediction |
|---|---|
| H1 | **The HAC interval's PD coverage follows its bias and variability:** pooled per group (B, C, D) and φ (0, 0.3, 0.6), within ±0.04 of the reference in at least 8 of the 9 cells. Reference: B 0.864 / 0.815 / 0.707; C 0.901 / 0.871 / 0.799; D 0.927 / 0.909 / 0.862 |
| H2 | **The same for ρ:** within ±0.04 in at least 8 of the 9 cells. Reference: B 0.881 / 0.867 / 0.804; C 0.918 / 0.910 / 0.876; D 0.935 / 0.931 / 0.911 |
| H3 | **The gain under dependence grows with T:** at φ = 0.6, the paired gain in pooled PD coverage (HAC minus iid Wald) is ordered D > C > B, and D's is at least +0.08 (reference: D +0.14, C +0.10, B +0.03) |
| H4 | **It costs coverage under independence, most at T = 20:** at φ = 0, the HAC interval's pooled PD coverage is below the iid Wald's in each of B, C and D. The shortfall is ordered B > C > D, and B's is 0.03–0.10 (reference: B 0.059, C 0.035, D 0.017) |
| H5 | **It is not a fix at these T:** at φ = 0.6, the HAC interval's pooled PD coverage is below the band's lower edge (0.9273) in each of B, C and D, and below the AR(1) MLE's pooled profile PD coverage on the same panels in each (D-208: B 0.910, C 0.929, D 0.942) |
| H6 | **Its width follows E[Q]:** the root mean square of the per-fit ratio HAC SE / iid SE in logit PD, pooled per group, is within ±10% of the reference at φ = 0 (B 0.926, C 0.950, D 0.975: narrower than the iid SE) and at φ = 0.6 (B 1.170, C 1.318, D 1.425). The ratio is taken over the fits that have both SEs |

**Held** means held as stated. A prediction that holds in direction but misses in size is not held.

**Reported without a prediction:**

- the lag sensitivity, coverage and the SE ratio at L = 1 and at twice the rule, per group and φ (the
  reference has them, in `reference.csv`'s `_l1` and `_l2x` columns);
- group A;
- the HAC interval's verdicts against the recovery band, per scenario;
- ρ's SE ratio;
- the side of the misses (true value above or below the interval);
- the count of unflagged fits without a HAC SE;
- **the HAC interval with a fixed-b critical value**, per group and φ, from the same SE: û ± c(b) SE
  in logit units, with c(b) Kiefer and Vogelsang's (2005) two-sided 95% critical value for the Bartlett
  kernel, c(b) = 1.9600 + 2.9694 b + 0.4160 b² − 0.5324 b³, at b = L/T (the rule's L: b = 0.1, 0.075
  and 0.04 at T = 20, 40 and 100). It tests whether a critical-value change alone rescues the HAC
  interval at short T. Being unpredicted, it changes no verdict.

## How the comparison will be reported

- **One row per prediction:** prediction, result, held or not. The rows come from `score.py`, committed
  and merged before the run.
- **Where the results go:** `studies/README.md`, with a D-entry and a `MANIFEST.json`. The per-replicate
  rows are attached to a release, not committed (D-163).
- **Verdicts it could change:** none pinned. The HAC SE is a `native` option and changes no parity
  result.
- **The recommendation it informs:** if H3 to H5 hold, the HAC SE is a cheap partial correction at long
  histories, not a substitute for modelling the dependence. It would then be reported beside the iid SE
  as a diagnostic, a ratio far from 1 flagging dependence, rather than recommended as the interval.

**Monitoring implication (to be confirmed by the result, D-152):** a HAC/iid SE ratio well above 1 on a
long history is evidence of positive serial dependence. At T = 20 the ratio is too noisy to read on its
own: its effective degrees of freedom are about 9 under independence and 5.4 for PD at φ = 0.6
(`reference.csv`, `pd_nu`, the median over group B).

## Addendum, 2026-10-12: the runtime projection, recorded before the run (D-187)

- **The tool and scorer** (`study_hac_se`, `score.py`; D-210) were merged in cbdd234 (PR #65) before this
  projection was sampled. The run uses that commit's source, built MSVC Release.
- **The projection:** `perf/project_runtime.py --run "study_hac_se"` from 2 to 256 replicates per
  scenario, after the 120 s warm-up. At 256 replicates: 7,512 thread-seconds, 0.087 h (314 s) on 24
  threads at efficiency 0.997; 90% range ±1%; with the machine-state allowance (±5%, D-193), 298–330 s.
  The subset-based figure, for comparison, is 0.07 h.
- **The stopping rule:** the range (±1%) and the agreement with the previous sample size (1% apart) are
  within ±10%. The sample length is not: the 256-replicate sample ran 84 s against the 120 s minimum,
  and the tool's cap is 256. A 512-replicate sample would be about half the full run, so the projection
  is recorded as it stands.
- The sampler's rows were deleted unread; only its timing was used.
