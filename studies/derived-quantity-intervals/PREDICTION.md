# S-23 Intervals for the 99.9% conditional PD: predictions made before the run

Pre-registration for study S-23 (`studies/README.md`, D-148). This file is committed
**before** any S-23 quantity has been computed: no interval, estimate or coverage of the
conditional PD has been run, on the subset or anywhere else. The results will be compared with it
prediction by prediction, and misses will be reported, not explained away.

**This file is not edited after the run.** If exploration on the study subset leads to a change
of method or to a sharper prediction for the full matrix, that goes in a new, dated section
appended at the end, committed before the full-matrix run. Nothing above it changes.

**Evidence used in writing it:** only results already pinned before this file:

- the recovery summary: profile, Wald and bootstrap coverage of PD and ρ per scenario, and ρ̂'s
  bias (M1.8, M2a, M2b; D-121, D-131, D-137);
- the true value of the quantity in each scenario, which is arithmetic on the truth (below).

## The quantity

The conditional PD at the 99.9% adverse factor level:

    q(PD, ρ) = Φ( (Φ⁻¹(PD) + √ρ · Φ⁻¹(0.999)) / √(1 − ρ) )

Φ⁻¹(0.999) = 3.0902. With the engine's convention (higher Z = better conditions), this is the PD
conditional on the factor at its 0.1% quantile, the adverse side.

True values (%), which depend only on PD and ρ:

| PD / ρ | 0.02 | 0.12 | 0.24 |
|---|---|---|---|
| 0.1% | 0.368 | 1.566 | 3.529 |
| 1% | 2.816 | 9.033 | 17.568 |
| 5% | 11.122 | 27.018 | 44.030 |

Over the estimation box (PD ∈ [1e-4, 0.2], ρ ∈ [1e-3, 0.5]) q ranges from 1.455e-4 to 0.9713.

## What is run

- **Panels and fits:** the recovery panels and fits, unchanged: same seeds, same replicates, same
  box and grid, the parity integrator and its check ([`docs/methodology/recovery.md`](../../docs/methodology/recovery.md)). The model is correctly
  specified, so the target is the true q above; there is no pseudo-true value.
- **Order:** the study subset first (scenarios 29, 37, 72, 4, 68, 49, 7, 43, 51, R = 1,000), then the
  full matrix (81 × 1,000) to pin verdicts.
- **Three 95% intervals for q per replicate:**
  1. **Profile likelihood.**
     - The profile is P_q(c) = max ℓ(PD, ρ) over the box subject to q(PD, ρ) = c. On that curve,
       PD = Φ(√(1 − ρ)·Φ⁻¹(c) − √ρ·Φ⁻¹(0.999)), so the inner maximisation is over ρ alone, restricted
       to the values for which that PD lies in the box.
     - The interval is {c : P_q(c) ≥ ℓ_max − 1.9207}, the threshold of D-129.
     - Each end point is solved against the objective, with the same endpoint residual,
       `TOL_PROFILE_ENDPOINT_RESIDUAL_LL` = 1e-7, asserted on every fit.
     - An end point is **box-limited** if the maximiser there lies on a bound of the box, or the end
       point is the limit of q attainable in the box. Box-limited end points are flagged and never
       extrapolated.
     - Coverage is over all replicates, and an interval that could not be computed counts as not
       covering (D-131).
  2. **Delta-method Wald** in logit(q), from the covariance of (logit PD, logit ρ) given by the
     Hessian at the estimate (D-119). Coverage is conditional on unflagged replicates, and the
     verdict is DEFERRED at 5% or more flagged (D-121, D-124). The flag is per replicate, so the
     deferred scenarios are, by construction, the same 31 as for PD and ρ today; that is a fact,
     not a prediction.
  3. **Bootstrap percentile:** q evaluated at each of the existing B = 999 iid bootstrap replicates'
     estimates (D-136), type-7 quantiles, over all replicates.
- **Also reported:** the point estimate q̂ = q(PD̂, ρ̂), with its median relative error
  median(q̂/q − 1) and its RMSE, and the width of the profile interval as the ratio of its upper
  end to its lower end.
- **Verdicts:** the band 0.95 ± 3.29·√(0.95·0.05/m), 0.927–0.973 at m = 1,000 (D-131).
  - **PASS:** inside the band.
  - **CONSERVATIVE:** above the band, reviewed, with its box-limited fraction reported.
  - **KNOWN FINDING:** below the band, reviewed, with a diagnosis.
  - **UNREVIEWED:** anything else out of the band. It fails CI until it is reviewed.
  - The pinned verdicts are recorded as the existing families are, in the recovery goldens.
- **Checks that are acceptance conditions, not predictions:**
  - the profile interval contains q̂ in every replicate;
  - the endpoint residual holds on every fit;
  - a scipy script reproduces the profile end points on the replay panels, within a tolerance
    entered in the register (`tests/tolerances.toml`) before the full-matrix run.

## Mechanism behind the predictions

q moves much more with ρ than with PD at the values studied: at PD = 1%, moving ρ from 0.07 to 0.20
takes q from about 5.9% to 14.6%. So the interval for q should inherit mostly ρ's behaviour:

- **Direction:** ρ̂ is biased downwards at small T (D-131), so q̂ should be too, and q's misses
  should mostly lie below the truth.
- **Size:** a profile interval for a smooth scalar function has the same first-order coverage as
  one for a parameter. q's coverage should therefore track that of PD and ρ, closer to the lower
  of the two.
- **Near a bound:** where ρ's profile is truncated at its lower bound (group A of D-136), q's lower
  end is held up by the ρ floor through the nuisance. Those intervals should be wide on the safe
  side and box-limited.
- **Bootstrap:** percentile intervals for ρ undercover almost everywhere (D-137), and q inherits
  that.

## Predictions, full matrix (81 scenarios, profile interval unless stated)

Groups as in D-136: A near-uninformative or near a bound (31), B informative T = 20 (12),
C informative T = 40 (17), D informative T = 100 (21).

| # | Prediction |
|---|---|
| P1 | **Group D:** all 21 PASS; at most 1 out of the band |
| P2 | **Group C:** at least 14 of 17 PASS; at most 3 below the band, and any below are among 31, 32, 41, 49, 58, 66, 68 and 75 (the C scenarios whose PD or ρ profile coverage is 0.940 or less today), with 68 the likeliest |
| P3 | **Group B:** coverage 0.910–0.950 throughout; 2–6 of 12 below the band; 29, 55 and 74 (today's ρ profile findings in B) are among them |
| P4 | **Group A:** at most 3 below the band, with 72 below or within 0.005 of the lower edge; 8–16 above the band, all of them in scenarios where PD or ρ is CONSERVATIVE today (0, 1, 3, 4, 6, 9, 12, 15, 18, 21, 27, 30, 33), each with box-limited end points in at least 25% of replicates, mostly at the lower end |
| P5 | **Totals:** PASS 55–72; above the band 8–18; below the band 3–10 |
| P6 | **Direction of misses:** in every below-band scenario, more than two-thirds of the non-covering intervals lie entirely below the true q |
| P7 | **Tracking:** in groups B–D (50 scenarios), q's coverage lies within [min(PD, ρ) − 0.015, max(PD, ρ) + 0.015] of today's PD and ρ profile coverages in at least 40 of 50 |
| P8 | **Point estimate:** median(q̂/q − 1) is negative in at least 9 of the 12 group B scenarios, between −1% and −10%; its magnitude at T = 100 is below its magnitude at T = 20 in at least 80% of the (PD, ρ, n) cells where neither is in group A |
| P9 | **Width:** in group B, the median ratio of upper to lower end is between 1.5 and 5; in cells where T = 20 and T = 100 are both informative, the median log-ratio at T = 100 is 0.35–0.60 of that at T = 20 (√(20/100) = 0.45) |
| P10 | **Wald (delta method):** among the 50 assessed scenarios, at least as many below the band as for the profile, concentrated in B, with misses mostly below the truth |
| P11 | **Bootstrap percentile:** below the band in at least 55 of 81; none above; at most 20 PASS |

**What would be a surprise worth its own diagnosis:** q's profile coverage 0.02 or more below the
lower of PD's and ρ's in an informative scenario. That would mean the curvature of q along the
likelihood surface makes the profile for a derived quantity worse than for either parameter.

## Predictions, study subset (explored first)

| Scenario | Group | PD, ρ, T, n | Predicted profile verdict for q | Predicted coverage |
|---|---|---|---|---|
| 29 | B | 1%, 0.02, 20, 10,000 | below the band | 0.910–0.935 |
| 37 | B | 1%, 0.12, 20, 1,000 | PASS, or just below | 0.925–0.950 |
| 72 | A | 5%, 0.24, 20, 100 | below the band, or within 0.005 above its lower edge | 0.915–0.940 |
| 4 | A | 0.1%, 0.02, 40, 1,000 | above the band, lower end box-limited in at least half the replicates | ≥ 0.974 |
| 68 | C | 5%, 0.12, 40, 10,000 | PASS, or just below | 0.920–0.950 |
| 49 | C | 1%, 0.24, 40, 1,000 | PASS | 0.930–0.960 |
| 7 | A | 0.1%, 0.02, 100, 1,000 | PASS or above the band | 0.945–0.985 |
| 43 | D | 1%, 0.12, 100, 1,000 | PASS | 0.935–0.965 |
| 51 | D | 1%, 0.24, 100, 100 | PASS | 0.935–0.965 |

**S1 (subset):** the predicted verdict holds in at least 7 of the 9 scenarios.

## How the comparison will be reported

Each prediction P1–P11 and S1 gets a row: prediction, result, held or not. A prediction that holds
in direction but misses in size is recorded as not held, with the size of the miss. Results and the
comparison go in `studies/README.md` (S-23's entry) and in a D-entry. The pinned verdicts get
diagnoses under D-131's policy.
