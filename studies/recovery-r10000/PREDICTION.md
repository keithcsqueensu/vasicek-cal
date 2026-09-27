# S-13 targeted: the borderline profile verdicts at R = 10,000 — predictions made before the run

Pre-registration for study S-13, targeted part (`studies/README.md`, D-148, D-151). This file is
committed **before** any replicate beyond the pinned 1,000 per scenario has been fitted. The
results will be compared with it prediction by prediction, and misses will be reported, not
explained away.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a
new, dated section appended at the end, committed before the run it concerns.

**Evidence used in writing it:** only the pinned recovery summary
(`tests/golden/recovery/summary.csv`: M1.8, M2a, M2b; D-121, D-131, D-137).

## What the study asks

D-151 puts this run first so that S-4 is judged against findings known to be real. The six
small-T profile findings are 29/PD, 29/ρ, 55/ρ, 68/PD, 72/ρ and 74/ρ, with coverage 0.920–0.927
against a band from 0.927. They and every other verdict near a band edge rest on 1,000
replicates, where the Monte Carlo SE of a coverage near 0.94 is 0.0075. At R = 10,000 it is
0.0024, and the band narrows from 0.927–0.973 to **0.9428–0.9572**. The question is which
verdicts survive.

## What is run

- **Scope (targeted):** every profile-likelihood verdict for PD or ρ whose R = 1,000 coverage lies
  within 0.010 of a band edge (0.9273 or 0.9727), on either side. That is 38 verdicts in 31
  scenarios: 1, 2, 6, 7, 9, 12, 15, 18, 21, 22, 24, 27, 28, 29, 30, 31, 37, 41, 45, 47, 49, 54, 55,
  59, 60, 64, 65, 68, 72, 73, 74. The six findings are all in it.
  - **Why the profile interval only:** it is the recommended interval, S-4 is judged against its
    findings, and the Wald and bootstrap verdicts near an edge would take in 66 of the 81
    scenarios, which is the whole-matrix run that stays after M4.
  - **The other families come for free:** each fit also gives the Wald, bootstrap and, once S-23
    is pinned, the conditional-PD verdicts of the same scenarios. They are reported at
    R = 10,000, without predictions and without changing any pinned verdict.
- **Replicates:** 0–9,999 of each scenario with the recovery seed, fitted exactly as in M1.8/M2a.
  Replicates 0–999 are the pinned ones, reproduced bitwise (same platform as the goldens).
- **Verdicts:** the D-131 policy with the band at m = 10,000. Coverage over all replicates.
- **Cost:** about 280,000 new fits; an estimated 3–4 h on 24 threads.
- **Output:** a summary with the same columns as `summary.csv`, for these scenarios at R = 10,000,
  committed in `studies/recovery-r10000/` with its provenance. The pinned R = 1,000 goldens and
  verdicts are not changed by this study; any change to them is a separate decision.

## Method behind the predictions: shrink each coverage towards its group

The 38 verdicts were selected *because* their coverage is near a band edge, so the observed
values are biased away from the centre: the six findings are the lowest of 162. The prediction
corrects for that selection with the simplest empirical-Bayes estimate, computed within each
group (D-136: A near a bound, B T = 20, C T = 40, D T = 100) and parameter:

- the group's mean coverage m, its observed variance v, and the mean binomial variance
  s² = mean(c(1 − c)/1,000), which is the Monte Carlo part of v;
- the between-scenario variance τ² = max(v − s², 0) and the shrinkage factor k = τ²/(τ² + s²);
- the estimate of a scenario's true coverage, m + k·(c − m); its sampling SD √(k·s²);
- the prediction at R = 10,000: 0.1·c + 0.9·(the shrunk estimate), since replicates 0–999 are
  included, with an SD combining the estimate's uncertainty and the binomial SD of 9,000 new
  replicates.

| Group | Parameter | Scenarios | Mean coverage | Monte Carlo SD | Between-scenario SD τ | k |
|---|---|---|---|---|---|---|
| A | PD | 31 | 0.9548 | 0.0065 | 0.0157 | 0.85 |
| A | ρ | 31 | 0.9645 | 0.0058 | 0.0221 | 0.94 |
| B | PD | 12 | 0.9401 | 0.0075 | 0.0055 | 0.35 |
| B | ρ | 12 | 0.9371 | 0.0077 | 0.0026 | 0.10 |
| C | PD | 17 | 0.9468 | 0.0071 | 0.0067 | 0.47 |
| C | ρ | 17 | 0.9448 | 0.0072 | 0 | 0 |
| D | PD | 21 | 0.9478 | 0.0070 | 0 | 0 |
| D | ρ | 21 | 0.9479 | 0.0070 | 0 | 0 |

**What this says before any run.** In group B (T = 20) almost all of the spread between scenarios
is Monte Carlo noise, and the group's mean coverage, 0.940 for PD and 0.937 for ρ, is itself
below the R = 10,000 band. The six findings are therefore expected to be real, but not special:
they are the scenarios where noise pushed a general T = 20 undercoverage of about one point below
the wider R = 1,000 band. Many of today's group B PASS verdicts are expected to fall below the
narrower band too. In group C the undercoverage is about half a point for PD and absent for ρ;
in group A the scenarios differ genuinely (k near 1), so little shrinkage applies.

That is not a contradiction of the pinned verdicts: a narrower band detects smaller shortfalls.
If the predictions hold, the conclusion is sharper than today's list of six findings: profile
intervals at T = 20 undercover by about one point, systematically, and not only in the flagged
scenarios. That is the evidence S-4a (the Bartlett-corrected threshold) should be judged against.
The pinned R = 1,000 verdicts stay as they are; S-13's results are its own.

## Predictions

| # | Prediction |
|---|---|
| Q1 | **The six findings:** all six below the band at R = 10,000, with coverage 0.915–0.942 |
| Q2 | **Group B:** every targeted group B verdict (PD and ρ, 11) below the band, with coverage 0.930–0.942 |
| Q3 | **Group C:** 68/PD and 41/PD below the band; 31/ρ and 49/ρ inside it; 59/PD inside it or above |
| Q4 | **Group A, near the lower edge:** 72/ρ below the band (0.915–0.938); of 2/ρ, 45/PD, 45/ρ, 22/ρ, 60/PD, 72/PD, 2/PD and 28/PD, at least 6 of 8 below the band |
| Q5 | **Group A, near the upper edge:** all seven CONSERVATIVE verdicts stay above the band; of the six PASS verdicts at 0.963–0.972 (1/PD, 7/ρ, 27/PD, 15/PD, 24/PD, 54/ρ), at least 4 above the band |
| Q6 | **Totals over the 38:** below the band 19–25; above the band 10–16; PASS 1–6 |
| Q7 | **Direction for ρ:** in every targeted group B ρ verdict below the band, more than two-thirds of the non-covering intervals lie entirely below the true ρ (ρ̂'s downward small-T bias, D-131) |

The table below gives each verdict's predicted coverage (±2 SD); the verdict predictions follow
from it, and Q1–Q7 are what is scored.

| Scenario | Group | Parameter | R = 1,000 coverage (verdict) | Shrunk estimate of the true coverage | Predicted coverage at R = 10,000 | Predicted verdict |
|---|---|---|---|---|---|---|
| 72 | A | ρ | 0.924 (KNOWN FINDING) | 0.9266 | 0.926 ± 0.011 | below |
| 29 | B | PD | 0.920 (KNOWN FINDING) | 0.9331 | 0.932 ± 0.009 | below |
| 2 | A | ρ | 0.931 (PASS) | 0.9332 | 0.933 ± 0.011 | below |
| 45 | A | PD | 0.930 (PASS) | 0.9337 | 0.933 ± 0.012 | below |
| 74 | B | ρ | 0.921 (KNOWN FINDING) | 0.9354 | 0.934 ± 0.006 | below |
| 22 | A | ρ | 0.933 (PASS) | 0.9350 | 0.935 ± 0.011 | below |
| 45 | A | ρ | 0.933 (PASS) | 0.9350 | 0.935 ± 0.011 | below |
| 68 | C | PD | 0.924 (KNOWN FINDING) | 0.9360 | 0.935 ± 0.010 | below |
| 55 | B | ρ | 0.926 (KNOWN FINDING) | 0.9359 | 0.935 ± 0.006 | below |
| 29 | B | ρ | 0.927 (KNOWN FINDING) | 0.9360 | 0.935 ± 0.006 | below |
| 73 | B | PD | 0.931 (PASS) | 0.9369 | 0.936 ± 0.009 | below |
| 64 | B | ρ | 0.935 (PASS) | 0.9369 | 0.937 ± 0.006 | below |
| 37 | B | ρ | 0.937 (PASS) | 0.9371 | 0.937 ± 0.006 | below |
| 74 | B | PD | 0.933 (PASS) | 0.9376 | 0.937 ± 0.009 | below |
| 47 | B | PD | 0.936 (PASS) | 0.9387 | 0.938 ± 0.009 | below |
| 55 | B | PD | 0.936 (PASS) | 0.9387 | 0.938 ± 0.009 | below |
| 65 | B | PD | 0.936 (PASS) | 0.9387 | 0.938 ± 0.009 | below |
| 60 | A | PD | 0.936 (PASS) | 0.9388 | 0.939 ± 0.012 | below |
| 72 | A | PD | 0.936 (PASS) | 0.9388 | 0.939 ± 0.012 | below |
| 2 | A | PD | 0.937 (PASS) | 0.9396 | 0.939 ± 0.012 | below |
| 28 | A | PD | 0.937 (PASS) | 0.9396 | 0.939 ± 0.012 | below |
| 41 | C | PD | 0.935 (PASS) | 0.9412 | 0.941 ± 0.010 | below |
| 49 | C | ρ | 0.934 (PASS) | 0.9448 | 0.944 ± 0.004 | PASS |
| 31 | C | ρ | 0.937 (PASS) | 0.9448 | 0.944 ± 0.004 | PASS |
| 59 | C | PD | 0.964 (PASS) | 0.9549 | 0.956 ± 0.010 | PASS |
| 1 | A | PD | 0.963 (PASS) | 0.9618 | 0.962 ± 0.011 | above |
| 7 | A | ρ | 0.964 (PASS) | 0.9640 | 0.964 ± 0.011 | above |
| 27 | A | PD | 0.966 (PASS) | 0.9643 | 0.964 ± 0.011 | above |
| 15 | A | PD | 0.970 (PASS) | 0.9677 | 0.968 ± 0.011 | above |
| 24 | A | PD | 0.971 (PASS) | 0.9686 | 0.969 ± 0.011 | above |
| 54 | A | ρ | 0.972 (PASS) | 0.9715 | 0.972 ± 0.011 | above |
| 6 | A | PD | 0.977 (CONSERVATIVE) | 0.9737 | 0.974 ± 0.011 | above |
| 21 | A | PD | 0.978 (CONSERVATIVE) | 0.9745 | 0.975 ± 0.011 | above |
| 9 | A | PD | 0.980 (CONSERVATIVE) | 0.9762 | 0.977 ± 0.011 | above |
| 12 | A | PD | 0.980 (CONSERVATIVE) | 0.9762 | 0.977 ± 0.011 | above |
| 18 | A | PD | 0.981 (CONSERVATIVE) | 0.9771 | 0.977 ± 0.011 | above |
| 27 | A | ρ | 0.979 (CONSERVATIVE) | 0.9781 | 0.978 ± 0.010 | above |
| 30 | A | ρ | 0.982 (CONSERVATIVE) | 0.9809 | 0.981 ± 0.010 | above |

## How the comparison will be reported

Each prediction Q1–Q7 gets a row: prediction, result, held or not. A prediction that holds in
direction but misses in size is recorded as not held, with the size of the miss. Results and the
comparison go in `studies/README.md` (S-13's entry) and in a D-entry. Whether any pinned verdict
or reviewed list changes as a result is a separate owner decision, recorded there.
