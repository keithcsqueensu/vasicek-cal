# Recovery harness: does the estimator recover the truth? (M1.8)

The unit and cross-reference tests show that the engine computes the binomial-mixture
likelihood and its maximum correctly (M1.2–M1.7). They do not show that the maximum-likelihood
estimator is any good: that its estimates centre on the truth, that its error shrinks as data
accumulate, and that its standard errors produce intervals with the coverage they claim. The
recovery harness answers those questions on synthetic panels whose true parameters are known
(D-121). The generated tables are in [recovery_results.md](recovery_results.md).

## Design

**Scenarios.** Every combination of

| Quantity | Values |
|---|---|
| PD | 0.1%, 1%, 5% |
| ρ | 0.02, 0.12, 0.24 |
| periods T | 20, 40, 100 |
| obligors per period n | 100, 1,000, 10,000 (the same in every period) |

81 scenarios, numbered id = ((i_PD·3 + i_ρ)·3 + i_T)·3 + i_n with each index 0, 1, 2 in the
order listed. The low corner (PD 0.1%, n = 100: about one default every ten periods) is there
on purpose. It is the low-default-portfolio case, where the likelihood carries little
information and estimates pile up at the edges of the parameter space.

**Data.** Replicate r of scenario id is the panel that the DGP produces for seed
`0x4D31385245434F56`, scenario id, replicate r, the scenario's PD and ρ, and n obligors in
each of T periods ([dgp.md](dgp.md)). Each replicate is therefore a pure function of those
numbers, bitwise identical on every platform. There are R = 1,000 replicates per scenario.

**Fitting.** Each replicate is fitted exactly as a user would fit it:

- `engine::calibrate` with the parity integrator and its doubled check (D-118);
- the estimation box PD ∈ [1e-4, 0.2] and ρ ∈ [1e-3, 0.5], on logit axes with 61 × 41 points
  (the M1.7 box, D-115);
- sub-grid refinement (D-095);
- standard errors from the Hessian of the log-likelihood at the estimate (D-119).

## What is reported, per scenario and parameter

Three numbers, kept apart so that no exclusion quietly improves the result (D-121).

1. **Bias and RMSE over all replicates.** Bias = mean(θ̂) − θ, with its Monte Carlo standard
   error sd(θ̂)/√R, and RMSE = √(mean((θ̂ − θ)²)). These do not depend on intervals, so nothing is
   excluded. An estimate on the edge of the grid counts at its grid value.
2. **The fraction of replicates flagged.** A replicate is flagged when it has no reliable Wald
   interval:
   - its estimate is on the edge of the grid (no interior curvature, so no SE);
   - its estimate is within 2 SEs of a bound of the box, on either axis (`kFlagNearBound`);
   - or its surface is flat.

   The flag is joint: an estimate of ρ near its upper bound also flags the replicate's PD.
3. **Coverage among the unflagged replicates, which is conditional.** The 95% Wald interval is
   built in the axis's logit coordinate u, where the Hessian is taken: û ± 1.959964·se_u, with
   se_u = se / (dθ/du) at the estimate. It covers when it contains the true u. Coverage is the
   fraction of unflagged replicates whose interval covers.

   This is **conditional** coverage. The flagged replicates are the ones whose estimates fell
   near the edges, so leaving them out selects on the outcome. Conditional coverage therefore
   overstates how well the intervals work across all replicates, and it is never reported
   without the flagged fraction next to it.

## Two kinds of criterion (D-124)

The harness tests two different things, and they are judged differently.

**Engine correctness must hold.** Its failure would mean the engine is wrong:

- RMSE falls strictly as T grows (20 → 40 → 100), for each (PD, ρ, n) and both parameters;
- no replicate is flagged for an unconverged quadrature or a non-finite likelihood;
- the replay replicates agree across platforms (below);
- in the cross-reference suite (M1.7), estimates and SEs agree with the independent reference.

Bias is reported with its Monte Carlo SE and is not tested against zero: ρ̂ has a known
small-T downward bias.

**Interval coverage is reported as a finding.** It is a property of the interval method (here
Wald), not of the engine. For each scenario and parameter the verdict is:

- **DEFERRED** if at least 5% of the replicates are flagged. The Wald interval is the wrong tool
  there: near a bound the likelihood is far from quadratic, and the right tool is a
  profile-likelihood interval, due in M2. The verdict is neither passed nor failed now. It is
  listed, with its M2 dependency, as part of parity sign-off (D-123).
- **PASS** if the conditional coverage lies in the Monte Carlo band
  0.95 ± 3.29·√(0.95 × 0.05 / m), where m is the number of unflagged replicates. 3.29 is the
  two-sided 0.1% normal point. With 162 verdicts, about 0.16 would fall outside by chance alone.
- **KNOWN FINDING** if it lies outside the band and has been reviewed. Each is listed with its
  diagnosis in `kKnownFindings` (`tests/recovery/recovery.hpp`) and in the results.
- **UNREVIEWED** if it lies outside the band and is not listed. The tests allow none: a new
  out-of-band verdict, or a listed one returning inside the band, fails CI until it is
  reviewed.

**Evidence for the diagnoses.** For each scenario and parameter, among the unflagged
replicates, in the logit coordinate: the *SE ratio*, meaning the standard deviation of û
divided by the root-mean-square Hessian SE. The root mean square is the average on the
variance scale. A ratio above 1 means the SEs understate the actual spread of the estimates, and
the interval is too short.

## Profile-likelihood intervals (M2, D-128–D-131)

Each replicate also gets its 95% profile-likelihood interval (`engine/profile.hpp`), reported side
by side with Wald:

- **Profile coverage is over all replicates.** A truncated interval is still an interval, and a
  not-computed one counts as not covering. Nothing is deferred.
- **The verdict** is PASS inside the band. Outside it, the verdict is one of two reviewed labels,
  each valid only on its own side of the band:
  - **CONSERVATIVE**, above the band: intervals contain the truth more often than advertised,
    typically because truncation at a bound of the box widens them. A safe-side property.
  - **KNOWN FINDING**, below the band: a genuine limitation.
- **UNREVIEWED**, meaning anything else outside the band, fails CI.
- **The t(T − 1) Wald interval** is a `native` comparison on the same replicates as Wald,
  reported without a verdict.

**Results.** 136 PASS, 20 CONSERVATIVE, 6 KNOWN FINDING, 0 UNREVIEWED.

| Wald \\ profile | PASS | CONSERVATIVE | KNOWN FINDING | total |
|---|---|---|---|---|
| PASS | 80 | 0 | 1 | 81 |
| DEFERRED | 41 | 20 | 1 | 62 |
| KNOWN FINDING | 15 | 0 | 4 | 19 |

- **What the profile interval fixes:** 15 of the 19 Wald findings, and 41 of the 62 verdicts Wald
  had to defer.
- **The 20 conservative verdicts** (0.977–0.999) are all near-uninformative settings, with 25–100%
  of intervals truncated at a bound.
- **The 6 known findings** (0.920–0.927 against a band from 0.927, at T = 20–40) are the
  small-sample undercoverage of the likelihood-ratio interval.
  - **For ρ:** most misses lie below the truth (61–65 against 9–16), which is ρ̂'s downward bias.
  - **For PD:** the misses are symmetric.
- **Engine correctness:** the worst endpoint residual over all 81,000 fits is 1.7e-8 in
  log-likelihood, and every fit is asserted against 1e-7.
- **Backlog, as `native` options:** a Bartlett-type threshold correction and a bias-corrected ρ̂,
  each to be evaluated against these pinned verdicts.

## Bootstrap percentile intervals (M2b, D-135–D-137)

- **Setup:** each panel also gets an iid bootstrap of periods, with B = 999 and a 95% percentile
  interval.
- **Stream:** seed = recovery seed XOR (scenario << 32 | replicate), in the resampling key domain.
- **Coverage:** over all replicates, with the same verdict rules as the profile interval.
- **Predictions:** the expected verdicts were committed before the run, in
  [bootstrap_predictions.md](bootstrap_predictions.md). The ordering was verified prior to public
  release, and D-146 records the file's commit time and SHA-256.
- **Outcome:** 37 PASS and 125 KNOWN FINDING, all below the band and none above. Every verdict is
  pinned, and D-137 records the full prediction-against-result comparison, including the
  predictions that failed.

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

Two diagnoses:

- **Boundary breakdown in group A.** The bootstrap is inconsistent near a bound; 3,500 intervals
  had zero width.
- **No bias or skewness correction elsewhere.** This is why ρ undercovers at every T.

The conclusion for users: profile-likelihood intervals are the recommended method for inference;
bootstrap percentile intervals are not recommended for ρ. BCa is a backlog `native` option, and
would not fix group A.

## Results (Wald, M1.8)

The full tables are in [recovery_results.md](recovery_results.md).

**Engine correctness holds.** RMSE falls with T in every (PD, ρ, n) and for both parameters.
None of the 81,000 fits is flagged for quadrature or numerics.

**Coverage verdicts (162 = 81 scenarios × 2 parameters):** 81 PASS, 62 DEFERRED, 19 KNOWN
FINDING, 0 UNREVIEWED. The known findings are of two kinds.

- **Small-T Wald undercoverage (14: 12 of PD, 2 of ρ), at T = 20 and 40.** Conditional coverage
  is 0.910–0.927, against a band starting at 0.927. The SE ratio is 1.01–1.10: the observed
  information understates the spread of the estimates by a few percent at small T.
  - A diagnostic re-fit of nine of these scenarios (not committed) showed coverage of
    0.924–0.943 with the t quantile on T − 1 degrees of freedom in place of z, and 0.93–0.955
    with the actual spread in place of the SE. This is the textbook behaviour of Wald intervals
    from few periods, not an engine error. The engine's SEs agree with the independent
    reference's (M1.7, D-119).
  - A t(T − 1) interval will be evaluated in M2 as a `native` option alongside profile
    likelihood. It is not a parity change, because changing the interval after seeing it fail
    is not a defensible way to pass.
- **Overcoverage of ρ on skewed estimates (5), coverage 0.973–0.998.** The SE ratio is at most 1,
  so the SEs do not understate uncertainty. But û is skewed in the logit coordinate, towards the
  lower bound: a diagnostic re-fit of all five measured skewness from −0.43 to −0.67. A symmetric
  interval is therefore miscalibrated. Where a few replicates are flagged (up
  to 3.4%), leaving out that tail adds to it.

**DEFERRED (62).** Most are low-information scenarios: n = 100 with PD 0.1% (92–100% flagged),
and ρ = 0.02 at small T or n, where ρ̂ sits within 2 SEs of the 0.001 lower bound in most
replicates. The upper bound matters too: at ρ = 0.24 and T = 20, ρ̂ often lands within 2 SEs
of 0.5. Because the flag is joint, the PD verdict is then deferred as well. This is the case
profile-likelihood intervals exist for, and it is not a defect of the estimator.

## Reproducing

    cmake --preset cpu-release && cmake --build --preset cpu-release --target recovery_harness
    build/cpu-release/tests/recovery_harness --write

This refits all 81,000 replicates and rewrites:

- `tests/golden/recovery/summary.csv`, one row per scenario with exact hexadecimal values;
- `tests/golden/recovery/replay.csv`, replicates 0 and 1 of every scenario;
- `tests/golden/recovery/replay_panels.csv`, those replicates' default counts, which the scipy
  cross-check of S-23 reads (`unit_recovery` checks them against `dgp/`);
- `tests/golden/recovery/MANIFEST.json`, the provenance;
- `docs/methodology/recovery_results.md`.

It then checks that every file exists and is non-empty. The run takes about 48 minutes on 24
threads, and results do not depend on the thread count. `recovery_harness --replicates 50`
prints a quick summary without writing anything. `--scenarios 29,37,72,4,68,49,7,43,51` restricts
a run to the study subset (D-150), and `--replay-dir DIR` writes the two replay files to DIR from
any run over every scenario.

`recovery_harness --check` refits everything and compares against the committed summary. On
the platform that wrote it (the same compiler string as `MANIFEST.json`), every count, verdict
and value must match exactly. Elsewhere:

- flagged and covered counts may differ by up to 2 replicates (`TOL_RECOVERY_CROSS_PLATFORM_COUNT`).
  A replicate within about 1e-12 of a threshold, such as the 2-SE near-bound test or an interval
  end, can fall either side under last-digit libm differences;
- values are compared at the replay tolerances (`TOL_RECOVERY_REPLAY_REL`, `TOL_RECOVERY_REPLAY_SE_REL`);
- verdicts are compared wherever the counts they rest on agree exactly.

In CI, `unit_recovery` re-fits the replay replicates on every platform, comparing estimates at
1e-12 and SEs at 5e-10, with flags compared exactly. This is a comparison of its own under D-107; the SE
tolerance is wider than M1.7's because these panels reach n = 10⁴. It also checks that the committed summary
meets the criteria above, recomputing each verdict from its recorded counts under the current
tolerances and the reviewed list. The full `--check` takes hours on a CI runner, so it is not
part of CI.

Anyone can reproduce any single replicate from its DGP inputs and a fit by any
maximum-likelihood method. Summary statistics then follow from the definitions above.
