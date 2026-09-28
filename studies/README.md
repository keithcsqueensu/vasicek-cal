# Studies

The index of this project's studies of its own estimators, and how each came out. Every research
study and every `native` variant gets one entry here, which states the study's own research
question. A study's files live beside this index, in `studies/<slug>/` (D-153). The rules for
running studies are in D-148:

- **A variant never changes parity.** It is a `native` option or a standalone study.
- **Every study is pre-registered.** Before any run, `studies/<slug>/PREDICTION.md` is committed,
  stating the expected results, their direction and their rough size. Results are compared with
  it, and misses are reported, not explained away. Studies that share one run share one
  registration (D-155): it lives in the first study's directory, and each other study's
  `PREDICTION.md` points to it.
- **Explore on the study subset; pin on the full matrix.** Exploration runs use the fixed study
  subset below. The full recovery matrix (81 scenarios × 1,000 replicates) is run only to pin a
  study's verdicts.
- **Verdicts follow the M1.8/M2a policy:** PASS, CONSERVATIVE, KNOWN FINDING or DEFERRED, each with
  a diagnosis.
- **A variant is judged on every verdict it could change** (D-149), not only on the findings it
  targets. A fix that moves the targeted findings into the band but pushes PASS verdicts out of it
  is reported as both.
- **New DGP variants are specified like the base DGP** (D-149): a prose description, a line-for-line
  Python mirror and a reference-panel hash (D-110–D-112).
- **Misspecification studies report two targets** (D-150). Where the fitted model is not the
  generating one:
  - **the primary target is the pseudo-true value**, the value the fit converges to on a very
    long panel, since that is what the estimator estimates under misspecification. Coverage is
    judged against it;
  - **the distance from the generating parameter is reported beside it**, because that gap is the
    misspecification bias a practitioner cares about.

  This applies to S-11, S-24, S-25 and S-29 to S-33.
  - **The pseudo-true value, defined (D-151).** A "very long panel" is well defined when
    every period has the same distribution (S-11, S-24). It is not when the design changes over
    the sample (a trend, a break, an inflow, a recording change: S-25, S-29, S-31, S-32, S-33),
    where a longer panel changes the design itself. The definition, which covers both: the
    pseudo-true value is the maximiser of the expected log-likelihood of the study's exact
    design (the same T, n_t and PD path), computed by summing over d with the quadrature where the
    observed counts have a computable distribution, and by a long simulation otherwise. For a
    stationary design it equals the very-long-panel value.

- **Every entry has a monitoring implication** (D-152), filled in when the study finishes: the
  monitoring metric, threshold or data check the result supports, if any, or "none". M5's
  monitoring design draws its metrics from finished entries and is not frozen until the studies it
  relies on are done.

Each entry records the study's research question, the experiment, the prediction, the result, the status of any
mitigation and the monitoring implication. An entry's prediction is "not registered" until its `PREDICTION.md` is committed;
nothing is run before then.

## The study subset (approved, D-150)

Nine recovery scenarios, one for each combination of period count T ∈ {20, 40, 100} and
correlation ρ ∈ {0.02, 0.12, 0.24}. Within each cell, PD and n were chosen so that the subset holds
some of the pinned findings and some of the pinned PASS verdicts, and touches all four bootstrap
groups (D-136).

| Scenario | PD | ρ | T | n | Group (D-136) | Pinned verdicts it exercises |
|---|---|---|---|---|---|---|
| 29 | 1% | 0.02 | 20 | 10,000 | B | Wald small-T (PD, ρ); profile small-T (PD, ρ); bootstrap bias/skew (PD, ρ) |
| 37 | 1% | 0.12 | 20 | 1,000 | B | Wald small-T (PD); profile PASS; bootstrap bias/skew (PD, ρ) |
| 72 | 5% | 0.24 | 20 | 100 | A | Profile small-T (ρ); Wald deferred; bootstrap (PD, ρ) |
| 4 | 0.1% | 0.02 | 40 | 1,000 | A | Profile truncation-conservative (ρ); bootstrap boundary breakdown (ρ) |
| 68 | 5% | 0.12 | 40 | 10,000 | C | Wald small-T (PD); profile small-T (PD); bootstrap (PD, ρ) |
| 49 | 1% | 0.24 | 40 | 1,000 | C | Wald small-T (PD); profile PASS; bootstrap bias/skew (PD, ρ) |
| 7 | 0.1% | 0.02 | 100 | 1,000 | A | Wald deferred (34% flagged); profile PASS; bootstrap boundary breakdown (ρ) |
| 43 | 1% | 0.12 | 100 | 1,000 | D | Wald PASS; profile PASS; bootstrap bias/skew (PD, ρ) |
| 51 | 1% | 0.24 | 100 | 100 | D | Wald skewed overcoverage (ρ); profile PASS; bootstrap (ρ) |

- **Replicates:** R = 1,000 per scenario, with the recovery seed and replicate indices. Every subset
  panel is therefore the same panel as in the pinned full-matrix run, so a variant is compared with
  the baseline pairwise, on identical data, and the Monte Carlo band (0.927–0.973) is unchanged.
- **Cost (measured):** 16.7 s of single-thread CPU per replicate across the nine scenarios, for fit,
  profile intervals and the B = 999 bootstrap (10.5 s fit, 4.2 s profile, 2.0 s bootstrap). For
  R = 1,000 that is 4.6 thread-hours: about 12 minutes on 24 threads at the speed of the
  measuring machine, and an estimated 7 minutes on the development machine, whose threads run the
  M1.8 fits about 1.7 times faster. On a 4-core machine it is about 70 minutes.
- **The full matrix, for comparison:** 50.6 thread-hours for the same three steps, about 2 hours on
  24 threads at the measuring machine's speed (an estimated 1.3 hours on the development machine).
- **The 72 other scenarios are held out.** Nothing tuned during exploration has seen them, which is
  what makes the full-matrix run a test rather than a confirmation.

## Priority order (D-151)

Placement says *when* a study can run; this order says which run first, so that the studies do
not delay M3.

1. **First batch, the highest value per hour, before M3:**
   1. S-23: intervals for the 99.9% conditional PD, the quantity capital and stress calculations
      use;
   2. S-13, targeted: R = 10,000 on the borderline scenarios, to settle whether the 6 small-T
      profile findings are real;
   3. the shared jackknife run: S-3, S-5 and S-21;
   4. S-1 and S-2: fast, and they give this index its first finished entries;
   5. S-34: sensitivity to severe new periods. It reuses S-23's interval for the 99.9% conditional
      PD and costs little (D-152).
2. **Then M3 itself,** with S-8, S-9, S-15, S-27 and S-28 folded into it (S-15 moves from after M4).
3. **Everything else** keeps its placement and is ordered when M3 ends.

## Index

| # | Slug | Question (short) | Milestone | Status |
|---|---|---|---|---|
| S-1 | `z-sign-invariance` | Is calibration invariant under z → −z? | now / M3 | finished: finding, reviewed (D-301) |
| S-2 | `sample-size-table` | How many years are needed for a given accuracy? | now / M3 | finished (D-304) |
| S-3 | `jackknife-bias-rho` | Does jackknife bias correction fix ρ̂'s small-T bias? | now / M3 | finished (D-159) |
| S-4 | `bartlett-profile` | Does a Bartlett-corrected threshold fix the 6 small-T profile findings? | now / M3 | not registered |
| S-5 | `bca-intervals` | Do BCa intervals fix the 125 percentile findings? | now / M3 | finished (D-159) |
| S-6 | `pluto-tasche` | How prudent are Pluto–Tasche upper bounds? | now / M3 | not registered |
| S-7 | `grid-resolution` | How do accuracy and runtime depend on grid resolution? | now / M3 | not registered |
| S-8 | `mle-vs-mom` | How efficient is MoM relative to MLE? | M3 | not registered |
| S-9 | `bayes-coverage` | Do grid-Bayesian credible intervals have frequentist coverage? | M3 | not registered |
| S-10 | `parametric-bootstrap` | Do parametric bootstrap intervals cover? | subset now, full after M4 | not registered |
| S-11 | `misspecification` | How wrong is standard Vasicek under a misspecified DGP? | subset now, full after M4 | not registered |
| S-12 | `double-bootstrap` | Does an iterated bootstrap calibrate interval coverage? | after M4 | not registered |
| S-13 | `recovery-r10000` | Do borderline verdicts survive R = 10,000? | now, targeted (first batch); full matrix after M4 | targeted part finished (D-158) |
| S-14 | `backtest-power` | How many years detect a misstated PD? | M5 | not registered |
| S-15 | `bayes-sbc` | Is the Bayesian estimator calibrated (SBC)? | M3 (D-151) | not registered |
| S-16 | `fp32-search` | Does FP32 search with FP64 finalisation match pure FP64? | after M4 | not registered |
| S-17 | `gpu-scaling` | How does performance scale across GPU generations? | after M4 | not registered |
| S-18 | `z-sign-macro` | Are macro sign filters mapped to the Z convention correctly? | M7 (deferred) | deferred; not started |
| S-19 | `z-extraction` | Z_t extraction, E[Z_t given d_t], as a standard output | M7 (deferred) | deferred; not started |
| S-20 | `bsf-apply` | A reference Belkin–Suchower–Forest apply function | M7 (deferred) | deferred; not started |
| S-21 | `period-influence` | How much do one or two extreme periods drive ρ̂? | now | finished (D-159) |
| S-22 | `box-sensitivity` | How much of the CONSERVATIVE group does the box create? | now | not registered |
| S-23 | `derived-quantity-intervals` | Are intervals for the 99.9% conditional PD reliable? | first batch | finished (D-156) |
| S-24 | `pd-heterogeneity` | How much does pooled PD heterogeneity inflate ρ̂? | now (after its DGP variant) | not registered |
| S-25 | `pd-trend` | How much does a PD trend inflate ρ̂, and does detrending fix it? | now (after its DGP variant) | not registered |
| S-26 | `varying-n` | Does anything assume a stable n? | now | not registered |
| S-27 | `large-portfolio` | When is Vasicek-rate MLE indistinguishable from binomial MLE? | M3 | not registered |
| S-28 | `zero-default-rates` | Refuse, drop or censor zero-default periods in rate-based estimators? | M3 | not registered |
| S-29 | `scale-version-change` | How much crosswalk error before ρ̂ inflation is material? | now (two-grade part); M6 (full) | not registered |
| S-30 | `grade-granularity` | Which number of grades K minimises error, given T and scale stability? | now (per-bucket part); M6 (shared ρ) | not registered |
| S-31 | `composition-shock` | How biased are PD̂ and ρ̂ when a riskier segment joins in a stress year? | now (bias, exclusion, indicator); M6 (separate segment) | not registered |
| S-32 | `default-misrecording` | How do misrecorded defaults bias the fit, and which treatment helps? | now | not registered |
| S-33 | `survivorship-backfill` | How biased is a backfilled, survivor-only history, and does truncation fix it? | now (portfolio part); M6 (late rating assignment) | not registered |
| S-34 | `severe-period-sensitivity` | How much do the estimates and the 99.9% conditional PD move after one or two severe periods? | first batch (after S-1 and S-2) | finished (D-162) |

Costs below are estimates from the subset measurement above unless marked measured, and are
refined in each study's `PREDICTION.md`.

---

## Now or alongside M3 (cheap: reuse the per-period surfaces)

### S-1 Z-sign invariance, calibration half (`z-sign-invariance`)

- **Question:** the engine writes p(z) = Φ((Φ⁻¹(PD) − √ρ·z)/√(1 − ρ)), so a higher Z means better
  conditions (R-2). Z is integrated out of the likelihood, and its distribution is symmetric, so
  the calibration should not depend on that choice. Are PD̂, ρ̂, the log-likelihood and every
  interval (Wald, profile, bootstrap percentile) identical under z → −z, within rounding?
- **Experiment (redesigned, D-300):** an exact invariance, so a pass/fail property check, not a
  statistical study. A test-only objective with the opposite sign (+√ρ·z) and its own hint is
  compared with parity on 34 fixed panels chosen for the hard cases (zero- and all-default
  periods, high ρ, estimates at or near a bound, n up to 10⁶): first the per-period surfaces cell
  by cell, then estimates, SEs, profile and bootstrap intervals and flags. The integration rules
  must mirror as well (adaptive GH centred on the mirrored mode; the Gauss–Legendre sinh map for
  d ∈ {0, n} centred on the mirrored half-point), so this also tests that the rules have no hidden
  asymmetry. A wrong-hint control proves the comparison can see an asymmetry. No coverage verdicts;
  the recovery-matrix version runs only if the check finds a difference that needs statistical
  characterisation.
- **Scope:** the calibration half only. The half where the sign matters, extracting Z_t and
  mapping macro effects to signs, is S-18 and S-19.
- **Cost:** a few minutes on 4 cores.
- **Prediction:** registered in [`studies/z-sign-invariance/PREDICTION.md`](z-sign-invariance/PREDICTION.md) before any run: agreement at the ε level (bounds derived from the rules' exact node symmetry and reversed summation order), identical flags, and any larger difference a finding.
- **Result ([`RESULTS.md`](z-sign-invariance/RESULTS.md), D-301): FINDING, reviewed; the invariance
  holds.**
  - **The surfaces agree at the ε level on all 34 panels.** The hints are exact mirrors at all
    912,865 cells. The worst surface cell is at 0.24 of its bound (5.7e-14), and every flag is
    identical. The wrong-hint control is caught (17 fits flagged by the quadrature check).
  - **Downstream, on 33 of 34 panels,** estimates agree to 8.5e-14, SEs to 2.7e-10, profile end
    points to 2.7e-12 logit and bootstrap ends to 2.4e-12. All are inside their bounds, though
    three sizes were underestimated (C3, C8, C10).
  - **Panel 26 (n = 1 in every period) is the finding.** ρ is not identified there, since a single
    obligor's likelihood does not involve ρ. The surface is flat along ρ to 2.8e-14, so rounding
    picks ρ̂: 0.262 under parity, 0.001 mirrored, and 767 of its 999 bootstrap argmaxes differ. The
    flags and the profile interval (the whole box, truncated at both ends) agree. This is not a
    sign asymmetry: any rounding perturbation would move ρ̂ the same way.
  - **The recovery-matrix version is not run (owner, D-301):** the mechanism is structural (every
    n_t = 1), every recovery scenario has n ≥ 100, and subset panels 1–9 pass every check.
- **Platform note (D-157):** the pinned result is from glibc. With the Windows UCRT libm (MSVC and
  MinGW GCC alike), panel 20 (n = 10⁵) exceeded two bounds at rounding level: C8 at 6.95·10⁻¹⁰
  against 5·10⁻¹⁰, and C11 at 3.6·10⁻¹² against an absolute 10⁻¹². C6 also differed on panel 26,
  the reviewed finding's own panel. By the prediction's rule these are findings. They are reviewed
  as rounding, not a convention error. C8 is now `TOL_ZSIGN_SE_REL` = 2·10⁻⁹; C11 is scaled by the
  engine's own rounding threshold; C6 joins C5, C7 and C10 in the panel-26 review. The prediction
  file is unchanged.
- **Mitigation: adopted (D-302).** A fit with no period of n_t ≥ 2 is flagged ρ not identified
  (`kFlagRhoNotIdentified`, `VCAL_FLAG_RHO_NOT_IDENTIFIED`), with its numbers still reported. A
  "flat to rounding" test is not added; that condition depends on the platform and stays with the
  existing flat-surface flag. `study_z_sign_invariance`
  (CTest, slow) guards the invariance and pins the finding.
- **Monitoring implication (D-152):**
  - **No metric for the sign itself.** PD̂, ρ̂ and every interval are the same under either
    convention, so a calibration run need not record it. Any output expressed on the factor scale
    (Z_t extraction, conditional PD by factor level, reverse stress) must state it, in adverse and
    benign terms (D-152, S-18, S-19).
  - **A data check:** ρ is identified only by periods with n_t ≥ 2. A monitoring report shows ρ
    as not identified, never as a number, when no period has n_t ≥ 2, or when ρ's profile interval
    is truncated at both ends of the box.

### S-2 Sample-size planning table (`sample-size-table`)

- **Question:** how many years of data are needed to estimate ρ within ±0.05, and PD within a
  stated relative error, as a function of PD and n?
- **Experiment:** no new fits. From the committed recovery summary, fit RMSE ∝ T^(−1/2) across
  T ∈ {20, 40, 100} for each (PD, ρ, n), check that the scaling holds on the three points, and
  solve for the T at which 1.96 × RMSE meets the target. Report the relative-error target for PD
  (proposed: ±25% relative) in the prediction. Mark every entry outside T ∈ [20, 100] as an
  extrapolation, and every scenario whose estimates are mostly at a bound as "not estimable at
  any T studied".
- **Cost:** seconds. An optional check refits the subset at one implied T (minutes); not planned
  (D-303).
- **Prediction:** registered in [`studies/sample-size-table/PREDICTION.md`](sample-size-table/PREDICTION.md) before any run (Q1–Q10), with the method fixed (D-303). The PD target is ±25% relative. A scenario is "mostly at a bound" when more than half its replicates are on the grid edge. The scaling holds when the log–log slope is in [−0.75, −0.30] with residuals ≤ 0.15 in log RMSE. The predictions come from large-n theory (SD(ρ̂) ≈ √2·ρ(1 − ρ)/√T), and no RMSE value was read in writing them.
- **Result ([`RESULTS.md`](sample-size-table/RESULTS.md), D-304, D-305):** the planning grids by
  PD, n and ρ. Figures marked † (beyond 100 years) or ‡ (below 20) extrapolate the T^(−½) law,
  which was verified only at T = 20, 40 and 100.
  - **The T^(−½) law holds** in all 24 cells where it could be checked (slope −0.40 to −0.59). One
    cell (PD 0.1%, ρ 0.02, n = 100) is not estimable at any T studied.
  - **ρ within ±0.05** takes about 35–40 years at ρ = 0.12 with n·PD ≥ 50, about 56 with
    n·PD = 10, and 180–225† with n·PD = 1. At ρ = 0.24 it is not met by T = 100 anywhere
    (110–120† for large portfolios).
  - **ρ within ±25% relative** (a post-hoc view, D-305) reverses the order: for large portfolios
    it takes 76–83 years at ρ = 0.24, 95–111 at 0.12 and 124–183† at 0.02. Much of the absolute
    target's high-ρ difficulty is its shape, since ±0.05 is ±20% of ρ at 0.24 but ±40% at 0.12.
    Relative to its size, ρ is estimated best when it is large (relative SD √2·(1 − ρ)/√T).
  - **PD within ±25%** is limited by the factor cycle: at 5% PD it takes 6–7, 38–39 and 82–86 years
    at ρ = 0.02, 0.12 and 0.24; at ρ = 0.12, 1% PD takes 71–81 years.
  - **Finite n** raises the years 1.6–2× at n·PD = 10 and 5–7× at n·PD = 1, against the large-n
    formula, which is close (T\* within about 16%) at n·PD ≥ 50.
  - **Six of the ten predictions held.** Q4, Q6, Q8 and Q9 missed: estimates with few defaults do
    better than predicted when ρ is low (PD is then close to the pooled binomial,
    T\* ≈ 61/(n·PD)), and (1%, 0.24, 10,000) needs 188† years against a predicted ≤ 170.
- **Mitigation:** n/a (a planning aid). The optional refit at an implied T was not run (D-303).
- **Monitoring implication (D-152):**
  - **A data-sufficiency check,** to be printed beside every calibration. It gives the expected 95%
    half-width 1.96·C/√T for the portfolio's PD, ρ and n, from the fitted C in `RESULTS.md`, or the
    large-n formula with the finite-n factor above, and the years T\* each target needs, marking any
    figure outside the verified 20–100 years as extrapolated, with ρ given against both an
    absolute and a relative target. When the history is shorter than T\*, the report says the
    data do not determine ρ (or PD) to that accuracy, so the value used rests partly on judgement
    or a floor and is documented as such.
  - **No threshold for year-on-year change.** The half-width is the scale of sampling error for
    one calibration. Successive calibrations share most of their periods, so how far one new year
    should move the estimates is S-34's question, and how far one period drives them is S-21's.

### S-3 Jackknife bias correction for ρ̂ (`jackknife-bias-rho`)

- **Question:** ρ̂ is biased downwards at small T (0.0002–0.014 where identified). Does the
  delete-one-period jackknife correction ρ̃ = T·ρ̂ − (T − 1)·mean(ρ̂₍₋ₜ₎) remove it, at what cost in
  RMSE, and do profile intervals shifted to the corrected estimate cover better?
- **Experiment:** a `native` option. The jackknife is already a W matrix (M2b), so it costs T grid
  reductions per replicate. Report bias, RMSE and the coverage of the shifted interval next to the
  parity profile interval, pairwise on the same panels.
- **Judged against:** the ρ small-T profile findings (29, 55, 72, 74) and the ρ PASS verdicts it
  could break.
- **Cost:** about the subset baseline (≈ 7–12 min); pinning ≈ 1.3–2 h.
- **Prediction:** registered, with S-5 and S-21, in [`jackknife-bias-rho/PREDICTION.md`](jackknife-bias-rho/PREDICTION.md) before any run: J1–J6 for S-3; J15–J18, a bias-corrected Wald interval for q, in its first addendum; J19–J23, a polished arm on the subset, in its second, merged before that arm ran.
- **Status:** finished 2026-09-28 (D-159). The summaries, per-replicate Parquet and provenance are in [`jackknife-bias-rho/`](jackknife-bias-rho/) (`MANIFEST.json`); the reviewed labels are kept apart in `reviewed.csv` (D-155), checked by `unit_study_jackknife`.
- **Result:** the jackknife removes ρ̂'s small-T bias, and **computed from exact delete-one refits it costs almost nothing** in RMSE or coverage. Computed from the grid refinement, as registered, it costs a great deal of coverage, and that cost is the grid's, not the jackknife's.
  - **Registered arm** (delete-one estimates from the 3 × 3 grid refinement; full matrix): the bias shrinks in 43 of the 50 informative scenarios. But the RMSE rises in all 12 of group B and by more than 3% in all 21 of group D, and the shifted interval's coverage moves by up to 0.143 from the parity profile's: 24 ρ profile PASS verdicts leave the band, all below. There are two patterns, and both come from the same source:
    - **ρ = 0.02** (13 scenarios in B–D): coverage 0.80–0.91, 0.04–0.14 below parity, and RMSE 1.3–6.2 × ρ̂'s. The refined delete-one estimates are off by up to 0.97 SE (0.27 SE where the refinement is accepted; S-21's J14), and ρ̃ = T·ρ̂ − (T − 1)·mean ρ̂₍₋ₜ₎ multiplies that error by T − 1. ρ̃ also falls to the box's floor, in up to 5% of replicates.
    - **ρ ≥ 0.12:** within 0.018 of parity at T = 20 and 40, but 0.009–0.117 below it at T = 100, where the factor T − 1 is largest, mostly from gaps of only 0.02–0.06 SE (0.38 SE at scenario 26, where some delete-one refinements are rejected).
    - **Group A:** ρ̃ is set to the floor in at least 5% of replicates in 23 of the 31; 16 shifted intervals fall below the band and 5 stay above it.
  - **Polished arm** (exact delete-one maxima; the 9 subset scenarios): in all 6 informative ones the shifted interval is within 0.008 of parity. At 29 it is 0.935 (registered 0.873, parity 0.927); at 51, 0.952 (registered 0.925, parity 0.953); at 43, 0.946 (0.930, 0.944). ρ̃'s bias is within ±0.0011 of zero, where ρ̂'s is −0.0011 to −0.0061, at an RMSE 0–5% above ρ̂'s. **So the registered arm's losses at T = 100 with ρ ≥ 0.12 were resolution as well, not the jackknife's own variance.** Group A keeps its losses: 4 and 7 stay 0.066 and 0.069 below parity, because setting ρ̃ to the floor is an effect of the box, not of resolution.
  - **What that means for J2, J3 and J5:** their mechanism (a small bias removed at a small cost in variance, with coverage barely moving) is what the polished arm shows. They missed because the registered estimator was not the one the mechanism described. They stand as not held; the polished arm explains the misses and does not rescore them.
- **q, a bias-corrected Wald interval (J15–J18, J23; S-23's mitigation 3):** neither sub-arm rescues the Wald interval for q.
  - **(a), q̃ with the delta-method SE:** above S-23's Wald coverage in 10 of its 17 below-band scenarios, in the band in 2. Polishing lifts it (29: 0.906, against 0.849 registered and 0.887 for S-23's Wald) but not into the band at T = 20.
  - **(b), q̃ with the jackknife SE:** covers at least as well as (a) in only 6 of the 17, and is in the band in none. It is never above the band.
  - **Diagnosis** ([`jackknife-bias-rho/diagnose.py`](jackknife-bias-rho/diagnose.py), written after the results). Outside ρ = 0.02, a symmetric interval q̃ ± 1.96 × (the error's standard deviation) covers 0.944–0.957: the corrected centre is nearly unbiased and nearly symmetric (skewness −0.22 to +0.02). Adding each SE's defects one at a time:
    - **SE_Δ:** slightly small on average (coverage 0.929–0.951 at its root-mean-square size), noisy (0.917–0.940 with its own values shuffled across replicates) and tied to the error (0.906–0.936 as it is; correlation +0.63 to +0.88);
    - **SE_J:** the right size on average (0.940–0.968), but twice as noisy (coefficient of variation 0.20–0.30 against SE_Δ's 0.10–0.15). The noise alone costs 2–4.5 points (0.907–0.928), because coverage is concave in the SE near 1.96: the replicates where the SE happens to be small lose more than those where it is large gain. Its link to the error is weaker (correlation +0.27 to +0.50; J17 held) and costs little.
    - **At ρ = 0.02** the corrected centre's error is 1.3–2.9 × wider than either SE, and skewed (up to +2.8), as the amplified refinement gaps would make it.

    In one line: **bias correction plus Wald fails for q because a symmetric interval needs a steady SE of the right size, and neither the delta-method SE (slightly small, and tied to the error) nor the jackknife SE (right on average but noisy) is one. Profile intervals need no SE.**
- **Comparison with the predictions** (scored by [`jackknife-bias-rho/compare.py`](jackknife-bias-rho/compare.py), committed before any result of the full run existed; D-159). Of S-3's J1–J6, J1 and J4 held; of the polished arm's J19–J23, 4 of 5 held; of the q sub-arms' J15–J18, J17 held. J20 missed favourably: polishing changed 43 and 51 by more than predicted, because the registered losses there were resolution too.

  | # | Prediction | Result | Held |
  |---|---|---|---|
  | J1 | \|bias(ρ̃)\| < \|bias(ρ̂)\| in at least 40 of 50 (B-D) | 43 of 50 | held |
  | J2 | RMSE(ρ̃) > RMSE(ρ̂) in at least 9 of 12 (B); within ±3% in at least 15 of 21 (D) | B: 12 of 12; D: 0 of 21 | not held |
  | J3 | B-D: shifted coverage within 0.015 of parity in every scenario, higher in at least 30 of 50 | largest \|difference\| 0.143 (26 beyond 0.015); higher in 9 | not held |
  | J4 | At most 2 of the ρ small-T findings (29, 55, 74, 72) move into the band | moved: 74 | held |
  | J5 | At most 2 B-D ρ profile PASS verdicts leave the band with the shift | 24: 56, 5, 13, 31, 32, 49, 58, 59, 8, 16, 17, 25, 26, 34, 35, 42, 44, 51, 52, 53, 61, 62, 71, 80 | not held |
  | J6 | ρ̃ clamped in at least 5% of replicates in at least 10 of 31 (A); under 1% in every C and D | A: 23 of 31; C/D at 1% or more: 5, 13, 31, 32, 58, 59, 8, 26, 34, 35, 61, 62 | not held |
  | J19 | 29: polished shifted coverage at least 0.91; ρ̃_p's RMSE at least 20% below the registered ρ̃'s | 0.935; RMSE 0.0070 vs 0.0099 | held |
  | J20 | 37, 68, 49, 43, 51: polished and registered shifted coverages within 0.01 | 37: +0.004; 68: +0.007; 49: +0.008; 43: +0.016; 51: +0.027 | not held |
  | J21 | 4, 7: polished shifted coverage still at least 0.03 below parity | 4: -0.066; 7: -0.069 | held |
  | J22 | At least 5 of the 6 informative: polished shifted coverage within 0.015 of parity | 6 of 6 (29: +0.008; 37: -0.004; 68: -0.003; 49: +0.000; 43: +0.002; 51: -0.001) | held |
  | J15 | (a) above S-23's Wald in at least 14 of 17; (a) in the band in 5-12 of them | above in 10; in the band in 2 | not held |
  | J16 | (b) covers at least as well as (a) in at least 14 of the 17; in the band at least as often | (b) ≥ (a) in 6; in the band: (b) 0, (a) 2 | not held |
  | J17 | error-SE correlation lower for SE_J than SE_Δ in at least 12 of the 17 | 17 of 17 | held |
  | J18 | of S-23's 32 Wald PASS at most 3 leave under (a) and at most 3 under (b); (b) above the band in at most 5 of 50 | (a) 11; (b) 16; (b) above 0 | not held |
  | J23 | q (a): polished at least registered - 0.005 in each of the 6, higher at 29 | 29: 0.906 vs 0.849; 37: 0.924 vs 0.916; 68: 0.933 vs 0.924; 49: 0.934 vs 0.929; 43: 0.932 vs 0.928; 51: 0.964 vs 0.944 | held |

- **Mitigation:**
  1. **Jackknife bias correction requires exact delete-one refits.** A `native` option must maximise each delete-one panel off the grid (T profile maximisations: seconds for one production fit), since ρ̃ multiplies the grid refinement's error by T − 1. *Supported by S-3's polished arm* on 6 informative scenarios; the full matrix was not polished.
  2. **Bias correction is a point-estimate adjustment, not an interval method.** Even exact, shifting the profile interval to ρ̃ gains nothing in coverage (within 0.008 of parity). Report ρ̃ beside ρ̂ where the small-T bias matters, and keep the profile interval around the MLE.
  3. **Not near a bound:** in group A, ρ̃ is set to the floor often, and the shifted interval undercovers by about 0.07 even when exact. Do not apply it where the profile interval is box-limited.
  4. **q:** bias correction plus Wald is not a remedy (above). S-23's mitigation 3 is updated accordingly.
- **Monitoring implication:** ρ̃ is a useful check on ρ̂ only when computed from exact delete-one fits: a large gap between them flags small-T bias or an influential period (S-21). From the grid refinement it mostly measures the grid. Intervals and thresholds stay on the profile likelihood.

### S-4 Bartlett-corrected profile threshold (`bartlett-profile`)

- **Question:** does scaling the χ²₁ threshold by a Bartlett factor fix the 6 small-T profile
  findings (scenario/parameter 29 PD, 29 ρ, 55 ρ, 68 PD, 72 ρ, 74 ρ; coverage 0.920–0.927) without
  breaking the 136 PASS verdicts?
- **Experiment, in two parts:**
  - **S-4a, oracle factor (cheap):** the mean likelihood-ratio statistic at the true value, from
    the recovery replicates themselves. It is not usable in practice, since it uses the truth, but
    it bounds what any Bartlett correction can achieve. It costs one profile maximisation per
    replicate.
  - **S-4b, feasible factor:** E[LR] estimated at θ̂ by parametric bootstrap. That needs new data
    per replicate, so it shares S-10's cost and belongs with it.
- **Cost:** S-4a ≈ subset baseline; S-4b as S-10.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** backlog item since D-131.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-5 BCa intervals (`bca-intervals`)

- **Question:** do bias-corrected and accelerated intervals fix the 125 pinned percentile findings?
  The two diagnoses in D-137 give the expected split: they should help with "no bias or skewness
  correction" (68 verdicts) and not with "boundary breakdown" (57).
- **Experiment:** a `native` option. z₀ comes from the existing B = 999 replicates; the
  acceleration from the jackknife (as S-3). Report BCa next to percentile, pairwise.
- **Cost:** about the subset baseline; pinning ≈ 1.3–2 h.
- **Prediction:** registered in the shared run's [`jackknife-bias-rho/PREDICTION.md`](jackknife-bias-rho/PREDICTION.md) (J7–J10 for S-5) before any run.
- **Status:** finished 2026-09-28, with S-3 and S-21 (D-159); data in [`jackknife-bias-rho/`](jackknife-bias-rho/).
- **Result:** **for ρ, BCa is a clear improvement on the percentile interval; for PD it changes nothing.**
  - **ρ:** BCa covers better than the percentile interval in all 50 informative scenarios, by 0.043–0.069 at T = 20, 0.012–0.053 at T = 40 and 0.001–0.031 at T = 100. In the band: 20 of 21 at T = 100 (percentile: 6), 8 of 17 at T = 40 (none), 1 of 12 at T = 20 (none). Group B stays below the band (0.861–0.930).
  - **PD:** BCa minus percentile has a median of +0.006, +0.003 and 0.000 at T = 20, 40 and 100. In the band: 26 of 50, the same number as the percentile interval.
  - **Group A (boundary breakdown, 57 verdicts):** 48 stay below the band, as predicted, although BCa is not computed in at least 2% of replicates in only 7 scenarios, not 10.
  - **Totals over 162:** 97 below the band (percentile: 125), none above.
  - **Why PD gains nothing:** z₀ and a correct the interval's location and skew; PD's percentile interval is short on width. A percentile interval for a mean of T periods is too narrow twice over: it uses the normal quantile where Student's t applies, and resampling T periods gives a variance (T − 1)/T of the truth. That alone predicts coverage of 0.929, 0.940 and 0.946 at T = 20, 40 and 100. PD's percentile coverage averages 0.901, 0.924 and 0.935 (BCa: 0.905, 0.925, 0.934), so width accounts for about half the shortfall at T = 20, and BCa for none of it. The rest is not diagnosed here. For ρ the shortfall is largely median bias and skew (D-137), which is what BCa corrects.
- **Comparison with the predictions:** J7 and J10 held; J8 and J9 did not. J8's miss is an effect that did not appear (BCa does not improve PD at all); J9's is in size.

  | # | Prediction | Result | Held |
  |---|---|---|---|
  | J7 | ρ B-D: BCa > percentile in at least 45 of 50; BCa below in at least 8 of 12 (B); at least 14 of 21 PASS (D) | 50 of 50; B below 11; D PASS 20 | held |
  | J8 | PD B-D: BCa in the band in at least 35 of 50 | 26 of 50 | not held |
  | J9 | Of the 57 boundary-breakdown verdicts at least 45 stay below; BCa not computed in at least 2% of replicates in at least 10 scenarios | 48 stay below; not computed at 2% or more in 7 scenarios | not held |
  | J10 | Totals over 162: BCa below 60-100; above at most 5 | below 97; above 0 | held |

- **Mitigation:**
  1. **Profile-likelihood intervals first.** They cover everywhere except about one point short at T = 20 (S-13).
  2. **If a bootstrap interval for ρ is needed, use BCa, never the percentile interval.** *Supported by S-5.* At T = 20 even BCa undercovers (0.86–0.93).
  3. **For PD, BCa is no remedy.** A correction of width (a studentised or parametric bootstrap, S-10; an iterated bootstrap, S-12) is the candidate, untested.
  4. **Near a bound, no bootstrap interval:** BCa does not repair the boundary breakdown.
- **Monitoring implication:** a bootstrap band for ρ should be BCa. For PD, do not build monitoring bands on the iid bootstrap at T ≤ 40: they are 2–5 points short of nominal, in either form.

### S-6 Pluto–Tasche most-prudent upper bounds (`pluto-tasche`)

- **Question:** how conservative are Pluto–Tasche upper bounds, and do they cover the true PD at
  their stated confidence, in the low-default scenarios?
- **Experiment:** the closed-form independent-period bound, and the bound under the one-factor
  model with independent period factors (the per-period mixture distributions convolved over T,
  which the existing quadrature supports). Report the frequency with which the bound is at or above
  the true PD, and its ratio to the truth.
- **Scope:** the serially correlated version needs the AR(1) factor and stays in M6, which already
  lists Pluto–Tasche as a benchmark.
- **Cost:** minutes for the full matrix.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-7 Grid resolution against accuracy and runtime (`grid-resolution`)

- **Question:** how do estimate accuracy (after sub-grid refinement), SE accuracy, interval end
  points and runtime depend on the grid? Is 61 × 41 (D-115) near the knee?
- **Experiment:** the subset at 31 × 21, 61 × 41, 121 × 81 and 241 × 161, against the finest grid as
  reference, pairwise per replicate.
- **Cost:** the surface cost grows with the number of points, so about 1 h at R = 1,000 or 15 min
  at R = 200 (enough, since this compares per-replicate differences, not coverage).
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

## With M3

### S-8 MLE against method of moments (`mle-vs-mom`)

- **Question:** what is MoM's efficiency relative to the MLE (relative RMSE, and bias), across the
  full matrix?
- **Experiment:** MoM (joint-default-probability form) on every recovery panel, compared pairwise
  with the MLE.
- **Cost:** minutes for MoM itself; ≈ 1 h if the MLE is refitted for pairing rather than read from
  the goldens.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-9 Frequentist coverage of grid-Bayesian credible intervals (`bayes-coverage`)

- **Question:** do 95% credible intervals under flat and Jeffreys priors cover at 95% in repeated
  sampling, and where do they differ from the profile intervals?
- **Experiment:** the grid posterior on the recovery panels; equal-tailed and HPD intervals. Grid
  resolution matters for a discretised posterior, so S-7 informs the grid.
- **Cost:** about the fit cost (≈ 5 min subset, ≈ 1 h full matrix), plus the Jeffreys prior once
  per n.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

## Subset now, full matrix after M4 (new data per replicate)

### S-10 Parametric bootstrap intervals (`parametric-bootstrap`)

- **Question:** do parametric bootstrap intervals, simulating from the DGP at θ̂, cover better than
  the iid percentile intervals (D-137)?
- **Experiment:** for each replicate, B panels from the DGP at θ̂ in their own key domain, each
  refitted; percentile and studentised intervals.
- **Cost:** with a full refit per bootstrap panel, B = 999 on the subset is ≈ 70 h on 24 threads;
  R = 200 and B = 199 is ≈ 5 h. If surfaces are cached by (n, d) (see the note in D-148), this
  becomes a W × L job at about the cost of the iid bootstrap.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-11 Misspecification (`misspecification`)

- **Question:** what bias and coverage does the standard Vasicek fit have when the data come from a
  t-copula (a common scale mixture on the factor), an AR(1) factor, or a beta mixture?
- **Experiment:** one parity fit per replicate from each alternative DGP; the same summaries as the
  recovery harness. Each DGP needs a prose description, a Python mirror and a hash, as `dgp/` does
  (D-110–D-112).
- **Targets:** the pseudo-true value (primary) and the distance from the generating parameters
  (D-150).
- **Cost:** one recovery run per DGP: ≈ 7–12 min on the subset, ≈ 1.3–2 h for the full matrix. The
  blocker is the DGP extensions, not compute.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

## After M4 (GPU scale)

### S-12 Double (iterated) bootstrap (`double-bootstrap`)

- **Question:** does calibrating the nominal level with an inner bootstrap bring iid bootstrap
  intervals into the band, away from the bounds?
- **Experiment:** nested iid bootstraps. An inner resample of an outer resample is again a weight
  vector over the original periods, so the whole study is W × L with no refits.
- **Cost:** B × C = 999 × 199 grid reductions per replicate: ≈ 3 h on the subset on CPU, ≈ 1 day
  for the full matrix on CPU, which is the GPU's fused-reduce kernel's job.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-13 R = 10,000 recovery re-run (`recovery-r10000`)

- **Question:** which borderline verdicts survive a Monte Carlo band about three times narrower?
- **Experiment:** the recovery harness at R = 10,000 for the scenarios whose verdicts lie within
  about 0.01 of a band edge, with the replicates 0–999 unchanged.
- **Placement (D-151):** the targeted run is in the first batch, before S-4 is judged against the
  6 small-T profile findings; the whole matrix at R = 10,000 stays after M4.
- **Cost:** the whole matrix ≈ 12–21 h on 24 threads on CPU; the borderline scenarios alone
  ≈ 3–4 h.
- **Targeted scope (D-154):** the 38 profile-likelihood verdicts for PD or ρ within 0.010 of a band
  edge at R = 1,000, in 31 scenarios; the other families are reported for those scenarios without
  predictions. Results are S-13's own: the pinned R = 1,000 verdicts are not overwritten.
- **Prediction:** the targeted part is registered in
  [`recovery-r10000/PREDICTION.md`](recovery-r10000/PREDICTION.md) before any run: predictions
  Q1–Q7, from an empirical-Bayes shrinkage of each coverage towards its group. The whole-matrix
  run is not registered.
- **Result (targeted part, D-158):** at R = 10,000 (band 0.9428–0.9572), **profile-likelihood
  intervals at T = 20 undercover by about one point, systematically**, not only in the six flagged
  scenarios.
  - **Group B:** all 11 targeted verdicts are below the band, with coverage 0.932–0.942. Seven of
    them were PASS at R = 1,000.
  - **The six findings:** all six stay below the band (0.932–0.942). They are real, but not special.
  - **q (S-23):** its profile interval shows the same at T = 20. Every informative T = 20 scenario
    here is at 0.936–0.941, below the band, where R = 1,000 had flagged only 72 and 74.
  - **Group C (T = 40):** a smaller shortfall for PD. 68/PD is just below the band (0.9422); 41/PD
    is inside it (0.9471), as are ρ's two borderline verdicts.
  - **Group A:** 6 of the 8 lower-edge verdicts are below the band. All 13 upper-edge verdicts are
    above it, including the six that were PASS at R = 1,000: the conservative, box-limited ones.
  - **Direction:** for ρ, 75–80% of misses lie entirely below the truth (ρ̂'s downward small-T
    bias).
  - **The other families of these scenarios, at R = 10,000, reported without predictions:**
    - Wald: 18 below the band, 1 above, 3 PASS, 40 DEFERRED;
    - bootstrap: 61 below the band, 1 PASS;
    - q's delta-method Wald: 10 below, 1 above, 20 DEFERRED;
    - q's bootstrap: all 31 below.
  - **Provenance:** replicates 0–999 reproduce the pinned fits field for field. The per-replicate
    estimates and profile ends are committed as
    [`recovery-r10000/fits_r10000.parquet`](recovery-r10000/fits_r10000.parquet), and the summary
    as `summary_r10000.csv`, with a `MANIFEST.json`. The pinned R = 1,000 verdicts are unchanged
    (D-154).
- **Comparison with the predictions** (scored by
  [`recovery-r10000/compare.py`](recovery-r10000/compare.py), committed before the run's results):
  4 held and 3 did not. The misses are in size, not direction.
  - **Q1** missed by 0.0002: 68/PD came in at 0.9422 against a predicted ceiling of 0.942, though
    all six are below the band as predicted.
  - **Q3:** 41/PD is inside the band, not below it.
  - **Q4:** 72/ρ came in at 0.9399, above its predicted range of 0.915–0.938, though below the band
    as predicted.

  | # | Prediction | Result | Held |
  |---|---|---|---|
  | Q1 | The six findings below the band, coverage 0.915-0.942 | 29/PD 0.9416, 29/ρ 0.9408, 55/ρ 0.9323, 68/PD 0.9422, 72/ρ 0.9399, 74/ρ 0.9369 | not held |
  | Q2 | Every targeted group B verdict (11) below the band, coverage 0.930-0.942 | 29/PD 0.9416, 74/ρ 0.9369, 55/ρ 0.9323, 29/ρ 0.9408, 73/PD 0.9420, 64/ρ 0.9376, 37/ρ 0.9414, 74/PD 0.9405, 47/PD 0.9417, 55/PD 0.9388, 65/PD 0.9418 | held |
  | Q3 | 68/PD and 41/PD below; 31/ρ and 49/ρ inside; 59/PD inside or above | 68/PD 0.9422, 41/PD 0.9471, 31/ρ 0.9443, 49/ρ 0.9461, 59/PD 0.9468 | not held |
  | Q4 | 72/ρ below (0.915-0.938); at least 6 of the other 8 group A lower-edge verdicts below | 72/ρ 0.9399; 6 of 8 below (2/ρ 0.9399, 45/PD 0.9352, 45/ρ 0.9387, 22/ρ 0.9399, 60/PD 0.9446, 72/PD 0.9422, 2/PD 0.9419, 28/PD 0.9450) | not held |
  | Q5 | All seven CONSERVATIVE above; at least 4 of the six upper-edge PASS verdicts above | CONSERVATIVE above: 7 of 7; PASS above: 6 of 6 (1/PD 0.9593, 7/ρ 0.9610, 27/PD 0.9615, 15/PD 0.9762, 24/PD 0.9748, 54/ρ 0.9765) | held |
  | Q6 | Totals over the 38: below 19-25; above 10-16; PASS 1-6 | below 19; above 13; PASS 6 | held |
  | Q7 | Group B ρ verdicts below the band: more than 2/3 of non-covering intervals entirely below the truth | 74/ρ 0.77; 55/ρ 0.80; 29/ρ 0.80; 64/ρ 0.78; 37/ρ 0.75 | held |

- **Mitigation:** the evidence S-4a (a Bartlett-corrected threshold) is to be judged against is
  now a systematic shortfall of about one point at T = 20, in PD, ρ and q, and about half a point
  for PD at T = 40. It is not six isolated findings. A correction should move group B as a whole,
  not only the flagged scenarios.
- **Monitoring implication:** at T = 20 a nominal 95% profile interval is about a 94% interval.
  A monitoring threshold built on it at T ≈ 20 should either widen it (S-4) or state its
  coverage as about 94%. Long histories (T ≥ 100) show no such shortfall.

### S-14 Backtest power (`backtest-power`)

- **Question:** how many years of data does each M5 backtest need to detect a PD misstated by 20%
  or 50%, at a given ρ and n?
- **Experiment:** power under the one-factor model. For tests on per-period or summed defaults it
  is exact by quadrature and convolution; simulation is the cross-check.
- **Extension (D-152): the operating characteristics of M5's threshold table.** For each tier
  boundary of the table that maps next-period default counts to statuses:
  - **false-alarm rate:** the probability, per period and over a horizon, that a correctly
    calibrated model reaches "warning" or "threshold exceeded". Under the one-factor model with
    independent periods it is exact by quadrature over the factor;
  - **detection delay:** the distribution of the number of periods until a miscalibrated model (PD
    misstated by 20% or 50%) first reaches each tier. With independent periods it is geometric in
    the per-period exceedance probability; with an AR(1) factor (M6) it needs simulation;
  - both with the table built from the true parameters and from estimated ones, so that the effect
    of estimation error on the false-alarm rate is visible.
- **Placement:** M5, where the backtests and the threshold table it needs are built (moved from
  after M4).
- **Cost:** minutes.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-15 Simulation-based calibration of the Bayesian estimator (`bayes-sbc`)

- **Question:** are the rank statistics of the true parameter within the grid posterior uniform
  (Talts et al. 2018)?
- **Experiment:** draw (PD, ρ) from the prior, simulate a panel, compute the posterior rank; a few
  (n, T) settings. Ties on a discrete grid need a stated tie-breaking rule.
- **Cost:** thousands of calibrate-only fits: minutes on CPU.
- **Placement (D-151):** folded into M3 with the grid-Bayesian estimator and S-9; it needs no GPU.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-16 FP32 search, FP64 finalisation (`fp32-search`)

- **Question:** does evaluating the grid in FP32 (with FP64 accumulation, D-039) and refining in
  FP64 give the same estimates, SEs and intervals as pure FP64, within the documented tolerances?
- **Experiment:** the subset under both precision policies, on the CPU first (the policies are
  host and device code), then on the GPU, so that precision and device-math effects are separated.
- **Cost:** about the subset baseline per policy per backend.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-17 Performance across GPU generations (`gpu-scaling`)

- **Question:** how do surface and fused-reduce throughput scale across GPU generations (at least
  sm_89 and sm_120)?
- **Experiment:** the `perf/` suite on each architecture. A performance benchmark rather than a
  methodology study: its results live in `perf/`, with this entry as the pointer.
- **Cost:** hours per GPU; needs access to hardware beyond the development machine's sm_120.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

## Addendum: now (reuse existing machinery; D-149)

### S-21 Period influence (`period-influence`)

- **Question:** how much do one or two extreme periods drive ρ̂ (and PD̂)? Report the influence of
  every period and the most influential periods in each scenario.
- **Experiment:** leave-one-period-out estimates from the existing jackknife weights (M2b), and
  leave-two-out from the T(T − 1)/2 pairs, which is also a W matrix (190, 780 and 4,950 rows at
  T = 20, 40, 100). Influence is reported in SE units, against the period's default rate and the
  realised factor Z_t (known from the DGP).
- **Resolution:** jackknife replicates use the 3 × 3 quadratic refinement on the grid surface, not
  the polished off-grid maximum. For a few replicates per scenario the leave-one-out estimates are
  also refitted exactly, so the refinement error is measured rather than assumed small.
- **Shared run (approved, D-150):** the same jackknife pass serves S-3 and S-5, and the three share
  one run, with the exact refits. It is also the first slice of M6's window and influence analysis, which then extends it to windows.
- **Verdicts it could change:** none; it is descriptive.
- **Cost:** ≈ 10–15 min on the subset including the fits; ≈ 3–5 min on top of an S-3/S-5 run
  (leave-two-out at T = 100 costs about five iid bootstraps). Full matrix ≈ 1.5–2.5 h. The exact
  refits add ≈ 10 min.
- **Prediction:** registered in the shared run's [`jackknife-bias-rho/PREDICTION.md`](jackknife-bias-rho/PREDICTION.md) (J11–J14 for S-21) before any run.
- **Status:** finished 2026-09-28, with S-3 and S-5 (D-159); data in [`jackknife-bias-rho/`](jackknife-bias-rho/).
- **Result:**
  - **Size:** in groups B–D the median of the largest leave-one-out change in ρ̂ is 0.53–0.69 SE at T = 20, 0.49–0.61 at T = 40 and 0.40–0.55 at T = 100, falling with T in every cell; its 90th percentile is 0.99–1.11, 0.79–1.21 and 0.62–0.98. For PD̂ the medians are 0.56–0.61, 0.46–0.58 and 0.33–0.55.
  - **Pairs:** the largest leave-two-out change is 1.65–1.87 times the largest leave-one-out change (median).
  - **Which period:** **the most influential period is often not the one with the most extreme factor.** It is in 25–80% of replicates, and under 60% in 33 of the 50 scenarios. The share rises with the expected defaults per period, n·PD: its median is 0.34 at n·PD ≤ 5, 0.45 at 10, 0.60 at 50–100 and 0.75 at 500. With few defaults a period's count is a noisy reading of its factor, so the realised count and the period's place relative to the fit decide which period moves ρ̂, not Z_t.
  - **Resolution:** the refined delete-one estimates differ from the exact maxima by up to 0.97 SE (scenario 35) and by 0.27 SE where the refinement is accepted; 18 scenarios reach 0.05 SE, mostly at ρ = 0.02. Against influences of 0.4–0.7 SE that is small outside ρ = 0.02 (gaps of at most 0.06 SE), but S-3 multiplies it by T − 1.
- **Comparison with the predictions:** J11 and J13 held; J12 and J14 did not, both in size.

  | # | Prediction | Result | Held |
  |---|---|---|---|
  | J11 | B-D median largest \|Δρ\| in range by T (0.45-0.85, 0.4-0.75, 0.3-0.65) in every scenario, falling with T in every cell outside A | outside range: none; cells not falling (by T = 20 scenario): none | held |
  | J12 | B-D: the most influential period has the most extreme Z in at least 60% of replicates, every scenario | range 0.25-0.80; below 60%: 33 | not held |
  | J13 | B-D: median leave-two/leave-one ratio 1.3-1.9 in every scenario | range 1.65-1.87 | held |
  | J14 | Refined vs exact delete-one estimates under 0.05 SE in every checked replicate outside A | largest 0.968 SE (scenario 35); over unflagged delete-one fits 0.267; scenarios at 0.05 or more: 18 | not held |

- **Mitigation:** compute influence from the data, with exact delete-one refits where ρ is small; do not infer it from the factor or from the default rate alone.
- **Monitoring implication:** "the crisis year drives ρ̂" is a hypothesis to check, not an assumption: with a few defaults a period, the period with the most extreme factor is the most influential in only about a third of histories. Report each period's leave-one-out influence beside ρ̂. At T = 20, a largest influence above about 1.1 SE lies beyond the 90th percentile of what a correctly specified model produces, and is worth a review.

### S-22 Box sensitivity (`box-sensitivity`)

- **Question:** how much of the CONSERVATIVE group (20 profile verdicts, coverage 0.977–0.999) does
  the estimation box itself create?
- **Prior evidence, from the pinned run and not from this study:** in a diagnostic refit of
  replicates 0–39 of those 20 verdicts, done while planning, truncated intervals were truncated
  almost only at the *lower* bound (ρ ≥ 1e-3, PD ≥ 1e-4). Upper-bound truncation occurred mostly
  together with lower (intervals spanning the whole axis, in near-uninformative panels); upper-only
  truncation was at most 6 of 40 in any verdict. Raising the ρ cap alone therefore tests the
  smaller part of the question. This is recorded here so that the prediction is written knowing
  it.
- **Experiment, two arms, each against parity on the same panels (both approved, D-150):**
  - **Upper cap:** ρ ≤ 0.9 instead of 0.5 (as asked).
  - **Lower bounds:** ρ ≥ 1e-5 and PD ≥ 1e-6 instead of 1e-3 and 1e-4. As ρ → 0 the model reaches
    the binomial, a true boundary of the parameter space, so this arm separates conservativeness
    caused by the box from conservativeness the boundary causes (D-129's boundary mixture).
  - In both arms the grid keeps its logit spacing and gains points, so the box effect is not
    confounded with grid resolution (S-7).
  - Above ρ = 0.5 the quadrature's precision rests on D-118 and the per-run check (D-092); any
    flagged fit is reported, not dropped.
- **Verdicts it could change:** every one. The box moves the grid-edge and near-bound flags, so
  Wald deferrals, profile truncation and bootstrap boundary breakdown can all move; hence the full
  matrix.
- **Cost:** each arm ≈ 1.1–1.4 × the baseline: ≈ 20–30 min on the subset for both arms,
  ≈ 3–5 h for the full matrix.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-23 Intervals for derived quantities (`derived-quantity-intervals`)

- **Question:** are intervals for the quantity a stress or capital calculation actually uses as
  reliable as those for PD and ρ? The primary quantity is the 99.9% conditional PD,
  q = Φ((Φ⁻¹(PD) + √ρ·Φ⁻¹(0.999))/√(1 − ρ)); other quantiles are optional.
- **Experiment:** a profile-likelihood interval for q, profiling along the curve on which q is
  fixed (reparametrise to (q, ρ) with PD = Φ(√(1 − ρ)·Φ⁻¹(q) − √ρ·Φ⁻¹(0.999))). q is not a grid
  axis, so the crossing is bracketed by search rather than from the grid, and the endpoint residual
  is asserted as in D-128. The delta-method Wald interval and the bootstrap percentile interval
  of q come for free and are reported beside it. A scipy script cross-checks the profile for q.
- **Verdicts it could change:** none of the existing ones; it adds a verdict family (81 per
  quantile) under the same band and policy.
- **Cost:** ≈ 1.15 × the baseline (one more profile per replicate): ≈ 8–14 min on the subset,
  ≈ 1.5–2.3 h for the full matrix.
- **Prediction:** registered in [`derived-quantity-intervals/PREDICTION.md`](derived-quantity-intervals/PREDICTION.md) before any run. It also defines a *box-limited* end point (the inner maximiser on a bound of the box), because q's interval can be held by the ρ floor through the nuisance without q itself reaching its range.
- **Status:** finished, pinned 2026-09-27 (D-156). Method: `engine/conditional_pd.hpp`; methodology note §6c; the tables are in `docs/methodology/recovery_results.md`.
- **Result:** profile-likelihood intervals for q are as reliable as those for PD and ρ, and in fact track them in every informative scenario. The other two intervals are not.
  - **Profile likelihood** (over all replicates): 64 PASS, 15 CONSERVATIVE, 2 KNOWN FINDING (72 at 0.918, 74 at 0.925, both T = 20; misses mostly below the truth, 80% and 73%). In groups B–D (50 informative scenarios) q's coverage lies within 0.015 of PD's and ρ's in all 50, and 49 are PASS. No scenario came 0.02 or more below the lower of PD and ρ: the pre-registered surprise did not occur.
  - **Delta-method Wald in logit(q)** (unflagged replicates): 32 PASS, 18 KNOWN FINDING, 31 DEFERRED. 17 are below the band: 11 of 12 at T = 20, 5 of 17 at T = 40, and 1 of 21 at T = 100 (35, n = 10⁴). Across them the misses are 1,373 below the truth and 76 above. The diagnosis is low estimates with narrow intervals, bias and width together:
    - q̂ is biased low (median q̂/q − 1 from −1% to −11%; mean −0.13 to −0.36 SE), carrying ρ̂'s downward bias;
    - its delta-method SE moves with the estimate: across replicates the error and the SE correlate at +0.64 to +0.88, so the lowest estimates get the narrowest intervals. The SE is right on average (sd/rms ratio 0.94–1.07) but not where it matters;
    - so the bias alone does not explain the shortfall: removing the mean bias exactly (an oracle shift) would lift coverage only to 0.916–0.947.

    The shortfall shrinks with T, but a longer history does not reliably remove it. Where n is large the SE is small enough that a 1.4% bias still matters at T = 100. **Scenario 13** is the one exception, the only scenario where the error goes the other way: above the band (0.983) because the delta-method SE overstates the spread there (ratio 0.94).
  - **Bootstrap percentile of q** (over all replicates): 8 PASS, 73 below the band, none above. Each finding inherits ρ's diagnosis in the same scenario: boundary breakdown in group A (31), no bias or skewness correction elsewhere (42).
  - **Point estimate:** median q̂/q − 1 is negative in every one of the 81 scenarios; −3% to −11% in group B, −0.7% to −4.5% at T = 100 with n ≥ 1,000.
  - **Acceptance checks:** every profile interval contains q̂ and every end-point residual is within 10⁻⁷ (worst 9.4·10⁻⁹). The first full-matrix run failed the first check in 11 of 81,000 fits, by 1 ulp at the box's corner; that was fixed before the pinning run (D-156). The scipy cross-check agrees within 1.2·10⁻⁹ in logit(q) on the 162 replay panels, and on every truncation and box-limited flag.
- **Comparison with the predictions** (scored by [`derived-quantity-intervals/compare.py`](derived-quantity-intervals/compare.py), committed before any result was seen; D-156). 7 held, 4 did not. Every miss is recorded as a miss, including the favourable ones: **the predictions underestimated how well profile intervals perform for q.** P3 and P5 missed because only 2 scenarios fell below the band, not 3–10; P4 missed on two conservative scenarios outside the predicted set; P9 missed on group B's widths (up to 9.7, in scenario 11, where PD is 0.1% and defaults are few), although its second half (the shrinkage from T = 20 to 100, 0.445–0.464 against √(20/100) = 0.45) came out as predicted.

  | # | Prediction | Result | Held |
  |---|---|---|---|
  | P1 | Group D: all 21 PASS; at most 1 out of the band | 21 of 21 PASS; out of the band: none | held |
  | P2 | Group C: at least 14 of 17 PASS; at most 3 below, all among 31, 32, 41, 49, 58, 66, 68, 75 | 17 of 17 PASS; below: none | held |
  | P3 | Group B: coverage 0.910-0.950 throughout; 2-6 of 12 below; 29, 55, 74 among them | coverage 0.925-0.952; below: 74 | not held |
  | P4 | Group A: at most 3 below, 72 below or within 0.005 of the lower edge; 8-16 above, all among 0, 1, 3, 4, 6, 9, 12, 15, 18, 21, 27, 30, 33, each box-limited in at least 25% of replicates, mostly at the lower end | below: 72; 72 at 0.918; above: 0, 1, 3, 4, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 54; box-limited condition met | not held |
  | P5 | Totals: PASS 55-72; above 8-18; below 3-10 | PASS 64; above 15; below 2 | not held |
  | P6 | Every below-band scenario: more than 2/3 of non-covering intervals entirely below the true q | 72: 0.80; 74: 0.73 | held |
  | P7 | Groups B-D: q's coverage within [min(PD, rho) - 0.015, max(PD, rho) + 0.015] in at least 40 of 50 | 50 of 50 | held |
  | P8 | Group B: median(q-hat/q - 1) between -10% and -1% in at least 9 of 12; its magnitude smaller at T = 100 than at T = 20 in at least 80% of cells with neither in group A | 10 of 12 in [-10%, -1%] (range -0.113 to -0.034); smaller at T = 100 in 12 of 12 cells | held |
  | P9 | Group B: median width ratio (hi/lo) 1.5-5 in every scenario; median log-ratio at T = 100 0.35-0.60 of that at T = 20 in every cell where both are informative | group B ratios 1.69-9.72; T = 100 / T = 20 log-ratios 0.445-0.464 (median 0.452) over 12 cells | not held |
  | P10 | Wald: among the assessed scenarios, at least as many below the band as the profile, concentrated in B (more than half), misses mostly below the truth | 50 assessed; Wald below 17 (11 in B) vs profile below 1; misses below/above 1373/76 | held |
  | P11 | Bootstrap: below the band in at least 55 of 81; none above; at most 20 PASS | below 73; above 0; PASS 8 | held |

  - **P4's two unexpected conservative verdicts** fit the same box-limited mechanism as the other 13:
    - **24** (PD 0.1%, ρ 0.24, T = 100, n = 100; coverage 0.973): the data are nearly uninformative, about 0.1 defaults a period. 95.7% of intervals are box-limited, at both ends: the ρ floor holds the lower end and the ρ cap of 0.5 the upper. PD's own profile coverage, 0.971, is just inside the upper edge, which is why the prediction left it out.
    - **54** (PD 5%, ρ 0.02, T = 20, n = 100; 0.977): 82.5% of intervals are box-limited, all at the lower end, where the ρ floor holds q's lower end through the nuisance. ρ's own coverage, 0.972, is 0.0007 inside the upper edge.
  - **Subset (S1, explored first):** held, 8 of 9 predicted verdicts; scenario 29 came in at 0.929, PASS, where below the band was predicted.
- **Mitigation for intervals on the 99.9% conditional PD:**
  1. **Use profile-likelihood intervals for q.** *Supported by S-23.* They cover correctly except in 2 scenarios at T = 20, and are conservative where the data are nearly uninformative. This is the recommended method.
  2. **Do not use delta-method Wald intervals for q.** *Supported by S-23.* The failure is low estimates with narrow intervals (q̂ inherits ρ̂'s downward bias, and its SE shrinks with it): 11 of 12 below the band at T = 20, 5 of 17 at T = 40, and still one at T = 100 with n = 10⁴, so a longer history reduces it without reliably removing it. Moving to the logit(q) scale is no remedy: the interval tested here is already symmetric in logit(q).
  3. **If a Wald-type interval is unavoidable** (a downstream system that takes only an estimate ± SE), candidate fixes, each still to be tested:
     - **bias-correct first:** *tested and failed* (S-3's q sub-arms, D-159). Neither the delta-method SE nor the jackknife SE rescues it, even with exact delete-one fits: a symmetric interval needs a steady SE of the right size, and the delta-method SE is slightly small and tied to the error, the jackknife SE right on average but noisy.
     - **parametric bootstrap:** *pending S-10*, the one alternative still untested.

     Profile-likelihood intervals are therefore the only supported method for q.
  4. **Treat short histories explicitly.** At T = 20 even profile intervals can undercover slightly, and S-13 is testing whether that is systematic. Report the history length alongside q's interval, and prefer the upper end of the profile interval when the estimate feeds a stress or capital figure.
- **Monitoring implication:** compare the production q against the *profile* interval of each fresh estimate. Wald intervals for q sit too low, so a check built on them fails in the costly direction: a production q that understates risk looks consistent with the data (a missed alarm), while a correctly set one is flagged as too high more often than the nominal 5%. Keep the history length beside every interval reported (see mitigation 4).

## Addendum: synthetic data variants (D-149)

S-24 and S-25 are DGP variants: each is built from the existing binomial DGP and gets its own key
domain, a prose description, a Python mirror and a hash. S-26 is a scenario definition on the
existing DGP, which already takes an n per period, so it needs no new generator code, mirror or
hash. S-24 and S-25 are misspecification studies of the same
kind as S-11 and report both targets (D-150).

### S-24 PD heterogeneity (`pd-heterogeneity`)

- **Question:** each panel pools two sub-segments with different PDs, loading on the same factor
  with the same ρ (no extra factor correlation). How much does the unmodelled heterogeneity
  inflate ρ̂?
- **DGP change:** the second segment's Bernoulli draws need their own stream but the same Z_t.
  Two calls of the existing `simulate_panel` with the same key would share the uniforms as well
  (common random numbers, an artificial dependence), so a segment index joins the key (approved,
  D-150). S-29, S-30 and S-31 reuse it for K segments.
- **Design:** the scenario's n split between the segments and their PDs set so that the pooled
  PD equals the scenario's; PD ratios of about 2, 5 and 10.
- **Verdicts it could change:** none pinned; it reports bias and coverage against the pseudo-true
  value, and the distance from the generating ρ.
- **Cost:** ≈ the baseline per ratio: ≈ 20–35 min on the subset for three ratios, ≈ 4–6 h for the
  full matrix. Plus the DGP extension and its mirror.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-25 PD trend (`pd-trend`)

- **Question:** PD drifts smoothly over the sample. How much does the trend inflate ρ̂, and does
  detrending fix it?
- **DGP change:** a PD per period (today's `PanelSpec` takes one PD).
- **Design:** a trend linear in probit(PD) with the scenario's PD at mid-sample, two slopes. Three
  fits per panel: naive; **oracle-detrended** (the true trend as a per-period probit offset, an
  upper bound on what detrending can do); and **estimated-detrended** (the trend estimated first,
  then ρ fitted with the offsets fixed, which ignores the trend's own estimation error).
- **Needs:** the detrended fits need an objective with per-period probit offsets, a `native`
  objective. S-29 and S-31 reuse it for a break and for a joining indicator.
- **Placement (approved, D-150):** the oracle and two-stage arms now; the joint fit of trend, PD
  and ρ, a three-parameter problem with the structure of M7's specification fits, in M7.
- **Verdicts it could change:** none pinned; bias and coverage against the pseudo-true value, and
  the distance from the generating ρ.
- **Cost:** ≈ the baseline per slope and fit arm: ≈ 45–70 min on the subset, ≈ 8–12 h for the full
  matrix (CPU, overnight). Plus the DGP extension and the offset objective.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-26 Varying portfolio size (`varying-n`)

- **Question:** n_t changes over time, for example a steady decline as in bank-count data. Does
  anything in the estimator or intervals implicitly assume a stable n?
- **No DGP change needed:** the DGP already takes an n per period, and its stream addressing
  already guarantees that changing one period's n leaves the others unchanged (M1.6). The variant
  is a scenario definition, not new generator code, so no new mirror or hash is required.
- **Design (approved, D-150):** for each scenario, n_t declining linearly from 1.5·n to 0.5·n, so the
  mean n is the scenario's and results pair with the constant-n scenario. Repeated in M3 for
  method of moments and the Vasicek-rate MLE.
- **Where a stable n is assumed:** D-122's deduplication only saves time when n repeats, so fits
  cost more but must be bitwise unaffected; the iid bootstrap resamples periods of different sizes;
  walk-forward windows differ in information. Among M3's estimators, method of moments needs a
  weighting across n and the Vasicek-rate MLE treats every period's rate alike, so they are the
  likelier places for the assumption; S-26 is repeated for them in M3.
- **Verdicts it could change:** every verdict of the paired constant-n scenario.
- **Cost:** without deduplication, fits cost ≈ 1.5–3 × the baseline: ≈ 15–35 min on the subset,
  ≈ 3–6 h for the full matrix.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

## Addendum: with M3 (D-149)

### S-27 Large-portfolio approximation (`large-portfolio`)

- **Question:** at what n does the Vasicek-rate MLE become indistinguishable from the binomial MLE,
  in estimates and in coverage? Deliver a rule of thumb by PD and ρ.
- **Experiment:** both estimators on the same panels; the matrix's n ∈ {100, 10³, 10⁴} is extended
  with n ∈ {10⁵, 10⁶}. "Indistinguishable" is defined in the prediction, for example a mean
  difference below 0.1 SE and coverage that differs by less than the Monte Carlo band allows.
- **Depends on S-28:** at low PD and moderate n many periods have no defaults, where the rate model
  is undefined, so S-27's answer depends on the zero-default treatment. Run S-28 first, or together.
- **Shared run (approved, D-150):** S-8, S-27 and S-28 run as one estimator-comparison pass.
- **Cost:** the Vasicek-rate MLE is closed form (seconds). The binomial fits at the new n, without
  bootstrap, ≈ 1.5–2 h for the full matrix and ≈ 10–20 min for the subset's 18 new scenarios; for
  n ≤ 10⁴ the goldens are reused.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-28 Zero-default treatments for rate-based estimators (`zero-default-rates`)

- **Question:** for the Vasicek-rate MLE and method of moments on rates, compare refusing (parity,
  D-044), dropping zero-default periods, and the censored likelihood (native, D-044) by bias, RMSE
  and coverage.
- **Design notes:** refusal produces no estimate, so its row is the refusal share (approved,
  D-150), not a bias. Runs in the S-8 estimator-comparison pass with S-27.
  Dropping is a data edit, so it can only be an explicit `native` option, never silent (D-044).
  The censored likelihood uses the Vasicek rate CDF, which is closed form, at a censoring point stated
  in the prediction. The binomial MLE, which handles zero defaults exactly, is the benchmark row.
- **Verdicts it could change:** none of the binomial verdicts; it adds the rate-based estimators'
  verdicts for each treatment.
- **Cost:** closed-form likelihoods: minutes for the full matrix.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

## Theme: data and population instability (S-29 to S-33; D-150)

These studies ask what happens to the calibration when the data describe a population, or a
scale, that changed over the sample. They extend three earlier studies:

- **S-24** (PD heterogeneity): pooling obligors with different PDs on the same factor. S-29's
  crosswalk error, S-30's bucketing and S-31's inflow all create this heterogeneity; S-24 is the
  stationary special case of each.
- **S-26** (varying n): S-31's inflow and S-33's survivor-only history both change n_t over the
  sample.
- **S-21** (period influence): S-31's joining year and S-32's misrecorded periods are single
  periods with outsized influence. Whether S-21's influence measure flags them is a check reported
  in both studies.

All five report the two targets of the misspecification rule (D-150).

### What runs on today's engine, and what needs multi-grade calibration

| Part | Single-segment engine and synthetic data variants: now | Needs multi-grade calibration (M6) |
|---|---|---|
| S-29 | Two-grade version with a crosswalk; per-grade separate fits on the remapped history; portfolio-level fit with allocation to grades; the structural break as a two-stage fit with S-25's offset objective | Probabilistic mapping in the likelihood; the structural break fitted jointly with a shared ρ; the full multi-grade scale |
| S-30 | Per-bucket separate fits (each bucket its own ρ); per-obligor PD error; the 99.9% loss quantile via S-23's conditional PD | Calibration with a shared ρ across buckets; the cross with S-29's crosswalk error |
| S-31 | Bias of the pooled fit; exclusion of the inflow; a joining-period indicator as a two-stage fit with S-25's offset objective | The inflow as a separate segment |
| S-32 | All of it: exclusion, the sensitivity band and the misclassification-aware likelihood (a `native` objective) | nothing |
| S-33 | Portfolio-level survivorship and backfill; truncation at the recording change, at the true and at a misjudged date | Late rating assignment, which acts per grade |

**Prerequisites (decided, D-151):**

- **Multi-grade calibration is an M6 item:** a single-factor model with K grades, each with its
  own PD, and a shared ρ (the model under the ASRF). All grades share the same Z_t, so each
  period's likelihood is still a one-dimensional integral,
  ∫ φ(z) ∏_k Binom(d_kt; n_kt, p_k(z)) dz. The quadrature stays cheap and only the optimiser
  grows. The approach: a one-dimensional grid in ρ, with the K PDs maximised for each fixed ρ
  (Newton or coordinate ascent), which keeps the profile-likelihood machinery for ρ intact. The
  adaptive rule's mode hint is revisited for the product integrand.
- **One DGP observation layer.** Every study here takes the true panel from the base DGP (with
  S-24's segments) and alters what is *recorded*: grades relabelled through a crosswalk (S-29),
  defaults missed or added (S-32), defaulters dropped from a backfilled history (S-33), a segment
  added from a date on (S-31). These are binomial and multinomial thinnings of counts, built once
  as a single layer with one prose description, one Python mirror and one hash, shared by S-29 and
  S-31 to S-33.
- **The grade scenario matrix**, below. The study subset and the 81-scenario matrix are single-PD,
  and S-29 and S-30 need grade scales.

### The grade scenario matrix (D-151)

Nine scenarios: three grade scales × T ∈ {20, 40, 100}, numbered id = i_scale·3 + i_T with the
scales in the order L, M, H.

- **The master scale:** 17 grades with PD_k = 0.03% × 1.5^(k−1), k = 1, …, 17, from 0.03% to 19.7%.
- **Obligors:** 10,000 per period, the same in every period. Grade k's share is proportional to
  exp(−(k − c)²/32), a discretised normal of standard deviation 4 grades centred on grade c. The
  counts n_k are 10,000 × share rounded down, with the remainder given one each to the largest
  fractional parts (ties to the lower k), so they sum to 10,000 exactly.
- **The three scales** (ρ shared by all grades of a scale; the values are the recovery matrix's,
  so results can be set beside the recovery scenarios of similar PD):

  | Scale | Centre c | ρ | Portfolio PD | Smallest grade |
  |---|---|---|---|---|
  | L (low default) | 5 | 0.24 | 0.61% | 13 obligors (grade 17) |
  | M (mid) | 8 | 0.12 | 1.54% | 83 obligors (grade 17) |
  | H (high default) | 11 | 0.02 | 3.42% | 46 obligors (grade 1) |

  Obligors per grade, k = 1 … 17:
  - L: 696, 866, 1012, 1112, 1147, 1112, 1012, 866, 696, 525, 373, 248, 155, 91, 50, 26, 13
  - M: 224, 337, 475, 629, 783, 916, 1006, 1037, 1006, 916, 783, 629, 475, 337, 224, 140, 83
  - H: 46, 84, 143, 229, 343, 484, 641, 797, 932, 1024, 1056, 1024, 932, 797, 641, 484, 343
- **Data:** replicate r of grade scenario id is the panel the DGP produces for seed
  `0x475241444553434E` ("GRADESCN"), scenario id, replicate r, with grade k drawn under S-24's
  segment index k on one factor Z_t per period. The observation layer then alters what is recorded.
  R = 1,000 replicates, so the Monte Carlo band is the recovery band.
- **Study-specific settings**, such as S-29's two-grade version (derived from these scales) and its
  crosswalk error levels, or S-30's bucketings, are stated in each study's `PREDICTION.md`.
- **Cost:** the costs in S-29 and S-30 below assume these 9 scenarios at R = 1,000.

### S-29 Rating-scale version changes (`scale-version-change`)

- **Question:** the master scale changes version partway through the history, and the crosswalk
  between versions is imperfect, including many-to-many mappings. How much crosswalk error can
  there be before ρ̂ inflation becomes material? "Material" is defined in the prediction.
- **DGP:** K grades on one factor (S-24's segment key); before the change date, each true grade's
  obligors and defaults are recorded under the old scale's grades through a crosswalk matrix with
  controlled error. Obligors and defaults are split separately by multinomial draws, in the
  observation layer.
- **Treatments compared:**
  1. remapping the history to the new scale, then fitting per grade;
  2. a probabilistic mapping in the likelihood (M6);
  3. a structural break: separate PDs per version, shared ρ (two-stage now, joint in M6);
  4. portfolio-level calibration with allocation to grades.
- **Now:** a simplified two-grade version with treatments 1, 3 (two-stage) and 4, at crosswalk
  error levels from 0 to about 20% plus a many-to-many case. **M6:** the full study.
- **Cost (now):** about 5 fits per replicate per error level, on smaller per-grade panels: about
  30–50 min at R = 200 on 9 grade scenarios, 2.5–4 h at R = 1,000. Full grade matrix about 8–12 h on
  CPU. The M6 part is estimated once the multi-grade optimiser exists.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-30 Grade granularity (`grade-granularity`)

- **Question:** given a fine true scale (15–20 grades with known PDs and a shared factor),
  calibrated at K = 4, 7, 10 and 20 buckets, which K minimises error, as a function of history
  length T and of scale stability?
- **Measures:** per-obligor PD error against the true grade PD; ρ̂ bias and coverage per bucket;
  the error in the 99.9% loss quantile (under one factor, the portfolio's conditional loss at unit LGD is
  the exposure-weighted sum of the buckets' conditional PDs, S-23); all as a function of T.
- **Design:** obligors bucketed by true PD (perfect ranking), so any error is due to granularity
  and data, not to the rank ordering. With K = 20 the buckets are finer than the truth in places;
  that arm measures over-splitting.
- **Now:** per-bucket separate fits, each with its own ρ. Many low-default buckets will sit near a
  bound, which is part of the answer. **M6:** a shared ρ across buckets, and the cross of K with
  S-29's crosswalk error.
- **Cost (now):** 41 bucket fits per replicate (4 + 7 + 10 + 20), mostly cheap low-default fits:
  about 35–60 min at R = 200, 3–5 h at R = 1,000. The cross with S-29 multiplies by the number of
  error levels.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-31 Composition shock (`composition-shock`)

- **Question:** a riskier sub-portfolio joins in a stress year, so the inflow coincides with an
  adverse factor. How biased are PD̂ and ρ̂?
- **DGP:** a second segment (S-24's key) with a higher PD on the same factor, present from the
  first period whose realised Z_t falls below a threshold. The joining date is selected on the
  factor, which is the point.
- **Treatments:**
  1. the pooled fit, which measures the bias (now);
  2. excluding the inflow's obligors (now; a single-segment fit on the original segment);
  3. an indicator for the joining periods as a probit offset, two-stage (now, with S-25's offset
     objective);
  4. the inflow as a separate segment (M6).
- **Cross-checks:** S-26's varying n (n jumps at joining) and S-21's influence of the joining year.
- **Cost (now):** 2 inflow sizes × 3 arms, with fit and profile: about 40–70 min on the subset,
  8–12 h for the full matrix.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-32 Default measurement error (`default-misrecording`)

- **Question:** in some periods defaults are misrecorded, missed (a detection rate below 1) or
  spurious (a false-default rate above 0), at controlled rates. How much does this bias PD̂ and ρ̂,
  and which treatment helps?
- **DGP:** in the observation layer, affected periods chosen at a set rate; recorded defaults are
  the true defaults thinned by the detection rate, plus false defaults drawn from the survivors.
- **Treatments:**
  1. exclusion of the affected periods, assuming they are known;
  2. a sensitivity band: the misclassification-aware fit repeated over a range of plausible rates;
  3. a misclassification-aware likelihood with an assumed detection rate δ and false-default rate
     ε: the conditional PD inside the binomial term becomes δ·p(z) + ε·(1 − p(z)).
- **Identification:** δ and ε are not identified from the counts, so they are assumed, not
  estimated. The question is how wrong the assumption can be before the fit is worse than the
  naive one.
- **Numerics:** the aware likelihood changes the integrand, so the adaptive rule's mode hint (built
  for the binomial-mixture integrand) may no longer fit it. The per-run check catches a poor rule.
  Before any study result is trusted, the objective is cross-checked against two oracles (D-151):
  a scipy script and `ref/`, whose adaptive Gauss–Kronrod integration shares no quadrature code
  with core.
- **All of it runs now** on the single-segment engine.
- **Cost:** 3 error settings × about 5 fits (naive, exclusion, the aware fit at 3 assumed rates,
  which also form the band): about 20–35 min at R = 200 on the subset, 1.5–3 h at R = 1,000; the
  full matrix about 20–30 h on CPU, so pin on the subset scenarios plus the settings that exploration
  shows matter.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-33 Survivorship and backfilled history (`survivorship-backfill`)

- **Question:** part of the early history is recorded only for obligors that survived to a later
  date, or ratings were assigned late rather than at origination. How large is the bias, and does
  truncating the history at the recording change remove it?
- **A design point (approved, D-151):** in its pure form the mechanism is trivial. An obligor recorded only if it
  survived to date s cannot have defaulted before s, so every backfilled period shows zero
  defaults, and truncating at s removes the problem by construction. The study is informative only
  where that is not true:
  - **partial survivorship:** a fraction of the defaulters in the backfilled periods is still
    recorded;
  - **a recording-change date that is misjudged**, so the truncation cuts too early or too late.
- **Overlap with S-32:** at portfolio level, partial survivorship is S-32's missed-default thinning
  concentrated in the early block of periods. They share the observation layer, and S-33's distinct
  questions are the truncation date and the late rating assignment.
- **Now:** portfolio-level partial survivorship at 3 retention levels, fitted naive, truncated at
  the true date and truncated at a misjudged date. **M6:** late rating assignment, which selects
  per grade.
- **Cost (now):** about 9 fits per replicate, the truncated ones on shorter panels: about 45–90 min
  on the subset, 10–15 h for the full matrix.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

## Addendum: sensitivity to severe new periods (D-152)

### S-34 Sensitivity to severe new periods (`severe-period-sensitivity`)

- **Question:** how much do PD̂, ρ̂ and the 99.9% conditional PD (S-23) move after one or two severe
  periods are added, as a function of T?
- **Experiment:** each recovery replicate is extended by one or two periods drawn at a fixed
  adverse factor level, a 1-in-100 and a 1-in-1,000 adverse year, with the period's default count
  set to the median of its binomial at that factor level, so the added periods are the same for
  every replicate of a scenario. The extended panel is refitted, and the shift in each estimate
  and in each interval (both end points) is reported in SE units and relative terms. Adverse
  means the direction of worse credit conditions, stated through the engine's convention in one
  place (S-1, S-18).
- **What it is not:** a bias study. The added periods are possible draws from the same model, so
  the truth does not change; the question is how strongly the estimates respond, which is what a
  monitoring process sees after a bad year.
- **Relations:**
  - S-21, which removes periods, is the mirror image. A severe added period should show up there
    as the most influential period.
  - S-23 provides the interval for the conditional PD.
  - M5's what-if recalibration is the same computation offered as a feature.
- **Cost:** a hypothetical period is one more surface row. Each variant costs a refinement and three
  profiles, not a new surface. Four variants per replicate (1 or 2 periods × 2 severities): about
  12–20 min on the subset and about 2–3.5 h for the full matrix. Measured: 35 min and 6 h 12 min,
  because the study tool refits each extended panel from a new surface rather than adding a row.
- **Verdicts it could change:** none; descriptive, with the coverage of the refitted intervals for
  the unchanged truth reported beside the originals.
- **Placement:** the first batch, after S-1 and S-2. It needs S-23's interval for the conditional PD
  and none of M3's estimators.
- **Prediction:** registered in [`severe-period-sensitivity/PREDICTION.md`](severe-period-sensitivity/PREDICTION.md) before any run: K1–K10, from the large-n limit of the model ([`large_n_reference.py`](severe-period-sensitivity/large_n_reference.py), which reads no recovery panel). The estimates are the exact off-grid maxima, not the grid refinement (D-159).
- **Status:** finished 2026-09-28 (D-162). The summary, per-replicate Parquet (28 MB, exact) and provenance are in [`severe-period-sensitivity/`](severe-period-sensitivity/) (`MANIFEST.json`).
- **Result:** **the size of the jump after a severe year is predictable from T alone, and matches the model's large-n limit almost exactly; what a refit does to the interval is not what a symmetric interval would suggest.**
  - **The jump in q̂** (medians, SE units of the original panel), groups B–D:

    | Added | T = 20 | T = 40 | T = 100 | Large-n reference, T = 20 / 40 / 100 |
    |---|---|---|---|---|
    | 1 × 1-in-100 | 0.64–0.84 | 0.41–0.61 | 0.24–0.40 | 0.81 / 0.59 / 0.38 |
    | 2 × 1-in-100 | 1.09–1.49 | 0.74–1.11 | 0.45–0.77 | 1.46 / 1.11 / 0.73 |
    | 1 × 1-in-1,000 | 1.21–1.45 | 0.79–1.12 | 0.53–0.89 | 1.40 / 1.02 / 0.67 |
    | 2 × 1-in-1,000 | 2.03–2.51 | 1.38–1.95 | 1.00–1.58 | 2.45 / 1.91 / 1.29 |

    Where n·PD ≥ 100 the shift is 0.93–1.02 × the reference for q and 0.93–1.06 × for ρ, in all 72 scenario-variants. In relative terms, one 1-in-100 year raises q̂ by 11–57% at T = 20 and 2–21% at T = 100; two 1-in-1,000 years by 35–248% at T = 20.
  - **Finite n does not dilute it.** Where n·PD ≤ 10 the shift is 0.61–1.33 × the reference (median 0.98): the original SE grows with the binomial noise as fast as the pull of the added period shrinks.
  - **The refitted interval moves up rather than widening.** Its lower end rises more than its upper end for ρ in all 50 B–D scenarios under every variant, and for q in most (33–40 of 50, by variant). The added periods make low dispersion implausible, so the lower end is pinned up. Hence:
    - **Coverage of the unchanged truth:** after one 1-in-100 year it barely moves (−0.020 to +0.034 from the pinned coverage). After one 1-in-1,000 year at T = 20 it falls to 0.81–0.89; after two, to 0.26–0.51 at T = 20, 0.54–0.79 at T = 40 and 0.69–0.88 at T = 100. The symmetric large-n model had predicted 0.72 at T = 20.
    - **Exceedance of the original interval:** the refitted q̂ lies above the original profile interval's upper end in at most 0.7% of replicates after one 1-in-100 year, at most 6.3% after one 1-in-1,000 year at T = 20, and 28–45% after two 1-in-1,000 years at T = 20; at T = 100 at most 8.3% under any variant.
  - **Group A:** in 10 of the 31 scenarios one "1-in-100" year *lowers* q̂ (by 1–12%): where n·PD ≤ 1 its median count (0–2 defaults) is milder than the dispersion binomial noise already implies. Three of them (0, 3, 6) were registered (K10).
  - **Flags:** the extended fits carry a Wald flag more often in B–D, from 328 of 50,000 originals to 482 (1 × 1-in-100) and 2,345 (2 × 1-in-1,000); every q profile interval was computed.
- **Comparison with the predictions** (scored by [`severe-period-sensitivity/compare.py`](severe-period-sensitivity/compare.py), committed before the full run started): **8 of 10 held.** K5 missed because the predicted effect did not appear (no dilution at small n·PD). K8 missed in its first part only, and in size: after one 1-in-100 year coverage stayed level rather than rising, and in 13 of the 50 it fell by 0.001–0.020. The subset, run first, had shown both tendencies; the predictions were not amended.

  | # | Prediction | Result | Held |
  |---|---|---|---|
  | K1 | B-D: median shifts of PD, rho, q positive in every scenario-variant (200 each; SE-unit and relative medians both) | pd: 200 of 200; rho: 200 of 200; q: 200 of 200 | held |
  | K2 | B-D: q shift (SE) falls with T in every cell outside A, each variant | 12 cells; not falling: none | held |
  | K3 | B-D: two/one q shift ratio 1.55-2.05 in at least 90 of 100 scenario-severity pairs | 100 of 100 (range 1.60-1.95) | held |
  | K4 | n*PD >= 100 (18 scenarios): q shift within ±25% of the reference in at least 58 of 72 | 72 of 72 (ratio to reference 0.93-1.02) | held |
  | K5 | n*PD <= 10 (23 scenarios): q shift below the reference in at least 69 of 92 | 49 of 92 | not held |
  | K6 | n*PD >= 100 (18 scenarios): rho shift within ±25% of the reference in at least 58 of 72 | 72 of 72 (ratio to reference 0.93-1.06) | held |
  | K7 | q-hat above the old upper end: 1x1-in-100 at most 5% in every B-D; 2x1-in-1000 at least 40% in at least 8 of 12 B; every D variant at most 10% | 1x1-in-100 over 5%: none; 2x1-in-1000 at 40% or more in B: 8 of 12 (range 0.279-0.446); D over 10%: none | held |
  | K8 | q coverage: 1x1-in-100 at least pinned in at least 40 of 50; 2x1-in-1000 at least 0.05 below pinned in at least 8 of 12 B, and below pinned in at least 15 of 21 D | 37 of 50; 12 of 12; 21 of 21 | not held |
  | K9 | T = 20: 1x1-in-1000 raises q-hat by a median of at least 15% in at least 10 of 12 B | 12 of 12 (range +0.189 to +1.167) | held |
  | K10 | Scenarios 0, 3, 6: q-hat's median relative shift at most 0 under all four variants | 0: -0.046, -0.089, -0.046, -0.089; 3: -0.024, -0.046, -0.024, -0.046; 6: -0.010, -0.019, -0.010, -0.019 | held |

- **Mitigation:** none needed for the estimator; the results are what a correctly specified model does. What needs care is how a refit is read after a severe year (below). The size of the jump is the M5 what-if recalibration's expected output, and it can be quoted from T alone: about 3.7/√T SE per 1-in-100 year and 6.5/√T SE per 1-in-1,000 year.
- **Monitoring implication:**
  - **A jump of the expected size is not evidence of a model change.** After one severe year, a rise in q̂ of up to about 3.7/√T SE (1-in-100) or 6.5/√T SE (1-in-1,000) is what the unchanged model produces. Compare the observed jump with this before reading it as a change in risk.
  - **The new q̂ above the previous interval is a strong signal.** One severe year, even a 1-in-1,000 one, almost never produces it (at most 6.3% at T = 20, 0.4% at T = 100). It takes two 1-in-1,000 years in a 20-year history. A breach therefore warrants a review.
  - **A refit right after a crisis overstates the long-run q.** Its interval is conditional on the bad years and covers the unchanged truth far less often than 95% (0.81–0.89 after one 1-in-1,000 year at T = 20, a quarter to a half after two). This is the conservative direction for capital, but a later run of ordinary years will then look like improvement. Report the number of severe years in the window beside every refit, and keep the pre-crisis fit alongside the post-crisis one.

## Deferred to M7 (recorded together, not started)

### S-18 Z sign convention and macro sign filters (`z-sign-macro`)

- **Question:** R-2. Expected macro effects are declared in economic terms ("worsens" or "improves
  credit conditions") and mapped to coefficient signs in one place, through the engine's
  convention (higher Z = better conditions). Is the mapping right?
- **Experiment:** a synthetic DGP with known macro effects of both signs; the sign filter must keep
  exactly the specifications with the right economic direction.
- **Status:** deferred to M7; not started.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-19 Z_t extraction (`z-extraction`)

- **Question:** E[Z_t | d_t] under the fitted parameters, as a standard output: is it unbiased for
  the simulated Z_t, and is its sign consistent with S-1 and S-18?
- **Experiment:** posterior mean of each period's factor on DGP panels whose factors are known
  (`simulate_panel` returns them).
- **Status:** deferred to M7; not started.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

### S-20 Belkin–Suchower–Forest apply function (`bsf-apply`)

- **Question:** a reference function taking a through-the-cycle migration matrix, ρ and z to the
  conditional migration matrix, with a test that fixes its conventions (the sign of z, the order of
  rating rows, the default column).
- **Status:** deferred to M7; not started. Migration matrices are outside the methodology scope
  listed in CLAUDE.md, so this needs a scope decision first.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

