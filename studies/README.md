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
| S-2 | `sample-size-table` | How many years are needed for a given accuracy? | now / M3 | not registered |
| S-3 | `jackknife-bias-rho` | Does jackknife bias correction fix ρ̂'s small-T bias? | now / M3 | registered (shared jackknife run) |
| S-4 | `bartlett-profile` | Does a Bartlett-corrected threshold fix the 6 small-T profile findings? | now / M3 | not registered |
| S-5 | `bca-intervals` | Do BCa intervals fix the 125 percentile findings? | now / M3 | registered (shared jackknife run) |
| S-6 | `pluto-tasche` | How prudent are Pluto–Tasche upper bounds? | now / M3 | not registered |
| S-7 | `grid-resolution` | How do accuracy and runtime depend on grid resolution? | now / M3 | not registered |
| S-8 | `mle-vs-mom` | How efficient is MoM relative to MLE? | M3 | not registered |
| S-9 | `bayes-coverage` | Do grid-Bayesian credible intervals have frequentist coverage? | M3 | not registered |
| S-10 | `parametric-bootstrap` | Do parametric bootstrap intervals cover? | subset now, full after M4 | not registered |
| S-11 | `misspecification` | How wrong is standard Vasicek under a misspecified DGP? | subset now, full after M4 | not registered |
| S-12 | `double-bootstrap` | Does an iterated bootstrap calibrate interval coverage? | after M4 | not registered |
| S-13 | `recovery-r10000` | Do borderline verdicts survive R = 10,000? | now, targeted (first batch); full matrix after M4 | targeted part registered |
| S-14 | `backtest-power` | How many years detect a misstated PD? | M5 | not registered |
| S-15 | `bayes-sbc` | Is the Bayesian estimator calibrated (SBC)? | M3 (D-151) | not registered |
| S-16 | `fp32-search` | Does FP32 search with FP64 finalisation match pure FP64? | after M4 | not registered |
| S-17 | `gpu-scaling` | How does performance scale across GPU generations? | after M4 | not registered |
| S-18 | `z-sign-macro` | Are macro sign filters mapped to the Z convention correctly? | M7 (deferred) | deferred; not started |
| S-19 | `z-extraction` | Z_t extraction, E[Z_t given d_t], as a standard output | M7 (deferred) | deferred; not started |
| S-20 | `bsf-apply` | A reference Belkin–Suchower–Forest apply function | M7 (deferred) | deferred; not started |
| S-21 | `period-influence` | How much do one or two extreme periods drive ρ̂? | now | registered (shared jackknife run) |
| S-22 | `box-sensitivity` | How much of the CONSERVATIVE group does the box create? | now | not registered |
| S-23 | `derived-quantity-intervals` | Are intervals for the 99.9% conditional PD reliable? | first batch | registered |
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
| S-34 | `severe-period-sensitivity` | How much do the estimates and the 99.9% conditional PD move after one or two severe periods? | first batch (after S-1 and S-2) | not registered |

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
  - **The recovery-matrix version is not run:** every recovery scenario has n ≥ 100, where ρ is
    identified, and subset panels 1–9 pass every check.
- **Mitigation:** proposed, not implemented (parity unchanged; the owner's decision, D-301). Either
  flag or refuse ρ when no period has n_t ≥ 2, or flag a surface that is flat to rounding along an
  axis at its argmax, and report no point estimate for that parameter. `study_z_sign_invariance`
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
- **Cost:** seconds. An optional check refits the subset at one implied T (minutes).
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a (a planning aid).
- **Monitoring implication:** to be filled in when the study finishes (D-152).

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
- **Prediction:** registered, with S-5 and S-21, in [`jackknife-bias-rho/PREDICTION.md`](jackknife-bias-rho/PREDICTION.md) (J1–J6 for S-3) before any run. **Result:** not run. **Mitigation:** backlog item since D-131.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

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
- **Prediction:** registered in the shared run's [`jackknife-bias-rho/PREDICTION.md`](jackknife-bias-rho/PREDICTION.md) (J7–J10 for S-5) before any run. **Result:** not run. **Mitigation:** backlog item since D-137.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

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
  run is not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

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
- **Prediction:** registered in the shared run's [`jackknife-bias-rho/PREDICTION.md`](jackknife-bias-rho/PREDICTION.md) (J11–J14 for S-21) before any run. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

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
- **Prediction:** registered in [`studies/derived-quantity-intervals/PREDICTION.md`](derived-quantity-intervals/PREDICTION.md) before any run. It also defines a *box-limited* end point (the inner maximiser on a bound of the box), because q's interval can be held by the ρ floor through the nuisance without q itself reaching its range. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

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
  12–20 min on the subset and about 2–3.5 h for the full matrix.
- **Verdicts it could change:** none; descriptive, with the coverage of the refitted intervals for
  the unchanged truth reported beside the originals.
- **Placement:** the first batch, after S-1 and S-2. It needs S-23's interval for the conditional PD
  and none of M3's estimators.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.
- **Monitoring implication:** to be filled in when the study finishes (D-152).

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

