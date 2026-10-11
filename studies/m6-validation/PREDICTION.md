# M6 validation: the AR(1)-aware MLE on its own model, predictions

Status: **registered before the run** (owner's design, D-197, with three changes; D-200). Committed with its scoring
script, [`score.py`](score.py), before any row of the run exists. Not edited after the run; results go in an
appended section.

## What this is, and what it is not

A correctness check of the estimator (`benchmark_models/ar1/`, D-195 to D-199) where theory says it should work,
before S-35's M6 arm uses it. It is not a preview of the arm's research question, and it cannot change it: the
S-35 M6 registration is frozen, and no amendment to it may follow this run (D-197). If a prediction below fails,
the cause is diagnosed before the arm runs; a fix is validated again on a fresh seed.

## The run

- **Panels:** the AR(1)-factor DGP (`dgp::simulate_panel_ar1`) on seed **0x7A11DA7E** (reserved in D-196, disjoint
  from S-35's, the information's, the comparator's, the recovery tests' and the scratch seeds), every recovery
  scenario (81: PD, ρ, T, n as `tests/recovery/recovery.hpp`) at φ ∈ {0, 0.3, 0.6}, replicates 0–99.
- **Tool:** `study_lrpd_serial_m6 --seed 0x7A11DA7E --replicates 100 --phi-intervals all` (D-198, D-200): per
  panel the binomial MLE, the M6 estimate with its checks, and 3-D profile intervals for PD, ρ and φ. No
  comparator (it is S-35's only).
- **Runtime:** projected from a timing sample before the run (D-187, D-193), on another seed (0x5AFE), so no
  validation panel is touched: 303,708 thread-seconds, **3.5 h** on 24 threads (90% range 3.3–3.7 h, 3.2–3.9 h with
  the machine-state allowance), accepted at 4 replicates per scenario, 1% from 2. The run is compared with it after.
- **Credible pairs:** a (scenario, φ) pair is *credible* where `studies/lrpd-serial-correlation/information_m6.csv`
  (the expected information at the true parameters, model-based, registered with S-35's M6 arm) marks all three
  parameters credible: 51 pairs at T = 100 (15 at φ = 0, 16 at 0.3, 20 at 0.6) and 13 at T = 40 (1 at 0.3, 12 at
  0.6); none at T = 20. Asymptotic arguments are applied only there.
- **SEs for standardising:** the information's expected SEs at the true parameters (`se_pd`, `se_rho`, `se_phi` in
  that file), not each fit's own, so a prediction does not lean on the estimator being checked. For the binomial MLE,
  the inverse of the information's (PD, ρ) block.

## Predictions

Strict pass/fail: W1–W5 and W7–W9. W6 is an expected range, reported against it (owner: small-T bias in φ̂ and ρ̂
can cause a genuine shortfall at T = 40 with a correct estimator).

| # | Prediction |
|---|---|
| W1 | **Checks before scoring** (as S-35's M6 arm; any failure stops scoring): every fit returns finite; its φ = 0 scan point's log-likelihood is at least the binomial MLE's on the same panel minus 1e-8; its maximum is at least the φ = 0 scan point's minus 1e-10 |
| W2 | **φ = 0 against the binomial MLE, paired** (owner: the strongest correctness check). In each of the 15 credible T = 100 cells at φ = 0, the 95th percentile over replicates of \|logit PD̂ − logit PD̂_bin\| is at most 0.25 of the information's SE of logit PD̂; the same for ρ. PD̂_bin is the binomial MLE polished (the scan's φ = 0 point). *Theory:* at φ = 0 the expected information makes φ nearly orthogonal to (PD, ρ); the asymptotic SD of the paired difference, sqrt(v₃ − v₂), is 0.3–10% of the SE in these cells, so its 95th percentile is at most about 0.21 SE |
| W3 | **The paired difference is the size theory says** (a Hausman-type comparison): pooled over the same 15 cells, the SD of (PD̂ − PD̂_bin)/sqrt(v₃ − v₂) is between 0.7 and 1.5, with v₃ and v₂ PD's variance from the full information and from its (PD, ρ) block; the same for ρ |
| W4 | **Kendall's bias at φ = 0, credible cells only:** pooled over the 15 credible T = 100 cells at φ = 0 (no T = 40 cell is credible at φ = 0), the median φ̂ is between −0.05 and +0.01: Kendall's −(1 + 3φ)/T = −0.01 for the latent factor, plus the noise of estimating it through counts. The other cells are excluded on purpose (owner): where φ is weakly identified the median φ̂ was near −0.4 on the exploratory panels (D-196), a property of the data, not a fault of the estimator, so pooling them would fail a correct estimator |
| W5 | **Coverage where asymptotics hold (strict):** pooled over the 51 credible T = 100 pairs (5,100 intervals per parameter: SE 0.3 points), the 95% profile intervals cover **PD and ρ** each with frequency in [0.93, 0.97], with no pair below 0.86 (four Monte Carlo SDs at 100 replicates); and **φ** with frequency in [0.92, 0.97], with no pair below 0.82. φ's band is computed, not assumed (owner), below. Per-pair coverages are reported beside the pooled ones (owner) |
| W6 | **Coverage at T = 40 (expected range, not pass/fail):** pooled over the 13 credible T = 40 pairs, coverage of PD and ρ in [0.90, 0.97], and of φ in [0.80, 0.97] (expected 0.920 at Kendall's bias, 0.821 at twice it; below) |
| W7 | **PD's median bias:** in every credible T = 100 pair, \|median over replicates of (PD̂ − PD)/se_pd\| ≤ 0.5 |
| W8 | **ρ's bias, judged against the iid MLE's (owner):** ρ̂ is biased low at small T in the iid binomial MLE already (M1). In at least 90% of the credible T = 100 pairs, the median standardised error of ρ̂, (ρ̂ − ρ)/se_rho, is within ±0.35 of the median standardised error of the binomial MLE's ρ̂ at φ = 0 in the same scenario (the correctly specified iid case, on the same factor draws) |
| W9 | **Flags:** pooled over every fit, at most 1% carry the placement-check flag and at most 1% the doubled-rule flag; and no (scenario, φ) pair has more than 5% of its fits with either |

**φ's coverage band, computed (owner).** φ̂ is biased low at finite T, which pulls the interval off the truth. For an
interval of φ's information SE s centred at φ + b, coverage is P(|Z + b/s| ≤ 1.96). With Kendall's b = −(1 + 3φ)/T
and s from the information file, the expected pooled coverage over the 51 credible T = 100 pairs is **0.945** (per
φ: 0.949 at φ = 0, 0.947 at 0.3, 0.940 at 0.6; worst pair 0.935). Kendall's formula is for an observed series; here
the factor is seen through counts, which can only add noise, so the band allows a bias of up to twice Kendall's:
pooled expectation **0.929**, worst pair 0.889. The lower bound 0.92 is three Monte Carlo SEs below that; the floor
0.82 is three per-pair SDs below the worst pair's. At T = 40 the same computation gives 0.920 (Kendall) and 0.821
(twice). PD and ρ keep the strict band, as the owner's design sets it; ρ̂'s known small-T bias is judged separately,
against the iid MLE's, in W8.

**Reported without a prediction:** coverage outside the credible pairs; how often φ̂ lies at ±0.95; the number of
local maxima along the scan and the gap to the second; profile intervals truncated at a bound; φ̂'s mean and median
per scenario and φ; and every per-pair number behind the pooled ones.

## Scoring

`score.py rows.csv` reads the tool's rows and the information file, runs W1 first and stops on a failure, then prints
each prediction's number and verdict and the per-pair table. Coverage counts a not-computed interval as missing
(D-131); a truncated interval covers when the truth lies inside it.

## Disclosure: the scorer's dry run (owner), after the predictions above were final

The predictions above were committed in c22c6a5. Then, as the owner asked, the scorer was run once on real tool
output to check its plumbing: `study_lrpd_serial_m6 --seed 0x5AFE --replicates 3 --phi-intervals all` (the smoke
seed, never the validation seed; 729 fits), scored with `score.py --dry-run`. It parsed every column, read the flags
from `tests/engine_flags.json`, joined the credible pairs and printed the per-pair table; without `--dry-run` it
refused the file (no validation seed). At three replicates per pair the numbers are not results, but they were seen,
so all of them are listed: W1 729 of 729; W2 worst 95th percentile 0.038 (PD), 0.071 (ρ); **W3 SD 1.597 (PD),
2.045 (ρ), against the band 0.7–1.5, over 45 paired differences**; W4 median φ̂ −0.0144; W5 pooled 0.922 (PD),
0.954 (ρ), 0.967 (φ), 22 pairs below a floor (three intervals per pair); W6 0.923, 0.872, 0.974; W7 worst +1.05
(medians of three); W8 66.7%; W9 0%. No prediction was changed after seeing them.

**After the dry run (owner): W3 kept as registered, and a descriptive diagnostic added (not a prediction).** If W3
fails, the likely cause is its denominator: sqrt(v₃ − v₂) is 0.3–10% of an SE in these cells, so finite-sample
effects that the asymptotic Hausman variance ignores can dominate the ratio, while W2 still checks the same paired
difference in absolute terms. `score.py` therefore also prints, per credible φ = 0 cell, the paired difference's
empirical SD and the asymptotic sqrt(v₃ − v₂) (both in SE units), their ratio, and the correlation across cells of
the log ratio with the log denominator, so the write-up can say how much of any excess the denominator explains.
It changes no verdict. Its plumbing was checked on the same smoke rows, where it showed correlations of −0.84 (PD)
and −0.54 (ρ) at three replicates per pair; those are disclosed here and are not results.

## Results (appended after the run, 2026-10-10; nothing above was edited)

The run: main at 1a97b4d, MSVC Release, 24 threads, 24,300 fits, one uninterrupted process; `rows.csv` attached to
the release `m6-validation-data-2026-10-10` (SHA-256 in [`MANIFEST.json`](MANIFEST.json)). Scored by `score.py`
unchanged; its full output, with the per-pair table, is [`results.txt`](results.txt).

**7 of 8 strict predictions held; W3 NOT HELD. W6 is inside its expected range.** The diagnosis below explains
W3's miss; it does not change the verdict.

| # | Result | Verdict |
|---|---|---|
| W1 | 24,300 of 24,300 fits pass the checks | held |
| W2 | worst 95th percentile of \|logit difference\| / SE: PD 0.056, ρ 0.099 (≤ 0.25) | held |
| W3 | SD of the paired difference / sqrt(v₃ − v₂): PD 1.689, ρ 1.592 (0.7–1.5) | **not held** |
| W4 | median φ̂ at φ = 0 over the 15 credible cells: −0.0188 (−0.05 to +0.01) | held |
| W5 | pooled coverage over 51 pairs: PD 0.9465, ρ 0.9418 (0.93–0.97), φ 0.9435 (0.92–0.97); no pair below its floor | held |
| W6 | pooled coverage over 13 T = 40 pairs: PD 0.9323, ρ 0.9315 (0.90–0.97), φ 0.9292 (0.80–0.97) | inside |
| W7 | worst PD median bias: −0.378 SE (\|·\| ≤ 0.5) | held |
| W8 | 98.0% of pairs within 0.35 SE of the iid MLE's ρ bias (≥ 90%) | held |
| W9 | placement flag 0.00%, doubled-rule flag 0.00%; worst pair 0% | held |

**W3, with the descriptive diagnostic registered before the run.** Per cell, the ratio of the paired difference's
empirical SD to sqrt(v₃ − v₂) ranges from 0.14 to 5.05 for PD and from 0.03 to 5.70 for ρ, and is below 1 in 12 of
15 cells for each. The pooled SD exceeds the band because of the cells with the smallest
denominators (scenario 52: 0.0045 SE for PD, ratio 5.05; scenario 61: 0.0034 SE, ratio 3.58). The log ratio's
correlation with the log denominator across cells is −0.96 for PD and −0.60 for ρ. Meanwhile the empirical SD
itself is nearly constant, 0.010–0.032 SE for PD, whatever the denominator.

**Diagnosis, after the run (not a prediction; changes no verdict; [`posthoc_w3.py`](posthoc_w3.py), output
[`posthoc_w3.txt`](posthoc_w3.txt)).** W3's premise does not hold at φ = 0, and the estimator behaves as theory says:

- *Theory:* at φ = 0 the periods are independent, so the φ score is a sum of products of adjacent periods'
  posterior factor means, each with mean zero. Its covariance with the PD and ρ scores and its expected derivative
  vanish, so the cross-information is exactly zero and v₃ − v₂ = 0. The Hausman first-order term is absent.
- *The denominator was noise:* the information is a Monte Carlo average over 500 panels per cell. At φ = 0 its
  cross-correlations average +0.005 (PD) and +0.010 (ρ), with a spread of 1.25 and 0.96 sampling SEs: zero plus
  noise. At φ = 0.6, by contrast, ρ's is −0.49, with every cell beyond two SEs.
- *The difference is second order:* within cells it does not move with φ̂ (R² 0.00–0.23; slopes near zero against
  the noise-based β), nor with φ̂², but its size does: corr(|difference|, |φ̂|) averages +0.63 for PD and +0.62 for
  ρ. That is φ̂ times a mean-zero coupling, the second-order term, and not a numerical floor, which would not
  scale with |φ̂|. φ̂'s own spread matches the information's SE (0.83–1.02 of it).

W2, which checks the same paired difference in absolute terms, held with room: its worst cell's 95th percentile is 0.1 SE (logit scale), against 0.25.

**Runtime: the projection held.** 12,596 s = 3.499 h against 3.5 h (3.3–3.7 h), projection / actual 1.005. The
machine was otherwise idle apart from about 20 minutes of light single-threaded work (`MANIFEST.json`).
