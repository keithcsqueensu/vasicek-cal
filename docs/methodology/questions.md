# Methodology questions

An index of the methodology questions this project asks of its own estimators, and how each was
answered. Every research study and every `native` variant gets one entry here. The rules for
running them are in D-148:

- **A variant never changes parity.** It is a `native` option or a standalone study.
- **Every study is pre-registered.** Before any run, `studies/<slug>/PREDICTION.md` is committed,
  stating the expected results, their direction and their rough size. Results are compared with
  it, and misses are reported, not explained away.
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
- **Misspecification studies state their estimand.** Where the fitted model is not the generating
  one, the prediction names the value the estimate is judged against: the generating parameter,
  or the pseudo-true value the fit converges to on a very long panel. Bias against one is not bias
  against the other.

Each entry records the question, the experiment, the prediction, the result, and the status of
any mitigation. An entry's prediction is "not registered" until its `PREDICTION.md` is committed;
nothing is run before then.

## The study subset (proposed, awaiting owner approval)

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

## Index

| # | Slug | Question (short) | Milestone | Status |
|---|---|---|---|---|
| S-1 | `z-sign-invariance` | Is calibration invariant under z → −z? | now / M3 | not registered |
| S-2 | `sample-size-table` | How many years are needed for a given accuracy? | now / M3 | not registered |
| S-3 | `jackknife-bias-rho` | Does jackknife bias correction fix ρ̂'s small-T bias? | now / M3 | not registered |
| S-4 | `bartlett-profile` | Does a Bartlett-corrected threshold fix the 6 small-T profile findings? | now / M3 | not registered |
| S-5 | `bca-intervals` | Do BCa intervals fix the 125 percentile findings? | now / M3 | not registered |
| S-6 | `pluto-tasche` | How prudent are Pluto–Tasche upper bounds? | now / M3 | not registered |
| S-7 | `grid-resolution` | How do accuracy and runtime depend on grid resolution? | now / M3 | not registered |
| S-8 | `mle-vs-mom` | How efficient is MoM relative to MLE? | M3 | not registered |
| S-9 | `bayes-coverage` | Do grid-Bayesian credible intervals have frequentist coverage? | M3 | not registered |
| S-10 | `parametric-bootstrap` | Do parametric bootstrap intervals cover? | subset now, full after M4 | not registered |
| S-11 | `misspecification` | How wrong is standard Vasicek under a misspecified DGP? | subset now, full after M4 | not registered |
| S-12 | `double-bootstrap` | Does an iterated bootstrap calibrate interval coverage? | after M4 | not registered |
| S-13 | `recovery-r10000` | Do borderline verdicts survive R = 10,000? | after M4 | not registered |
| S-14 | `backtest-power` | How many years detect a misstated PD? | after M4 | not registered |
| S-15 | `bayes-sbc` | Is the Bayesian estimator calibrated (SBC)? | after M4 | not registered |
| S-16 | `fp32-search` | Does FP32 search with FP64 finalisation match pure FP64? | after M4 | not registered |
| S-17 | `gpu-scaling` | How does performance scale across GPU generations? | after M4 | not registered |
| S-18 | `z-sign-macro` | Are macro sign filters mapped to the Z convention correctly? | M7 (deferred) | deferred; not started |
| S-19 | `z-extraction` | Z_t extraction, E[Z_t given d_t], as a standard output | M7 (deferred) | deferred; not started |
| S-20 | `bsf-apply` | A reference Belkin–Suchower–Forest apply function | M7 (deferred) | deferred; not started |
| S-21 | `period-influence` | How much do one or two extreme periods drive ρ̂? | now | not registered |
| S-22 | `box-sensitivity` | How much of the CONSERVATIVE group does the box create? | now | not registered |
| S-23 | `derived-quantity-intervals` | Are intervals for the 99.9% conditional PD reliable? | now | not registered |
| S-24 | `pd-heterogeneity` | How much does pooled PD heterogeneity inflate ρ̂? | now (after its DGP variant) | not registered |
| S-25 | `pd-trend` | How much does a PD trend inflate ρ̂, and does detrending fix it? | now (after its DGP variant) | not registered |
| S-26 | `varying-n` | Does anything assume a stable n? | now | not registered |
| S-27 | `large-portfolio` | When is Vasicek-rate MLE indistinguishable from binomial MLE? | M3 | not registered |
| S-28 | `zero-default-rates` | Refuse, drop or censor zero-default periods in rate-based estimators? | M3 | not registered |

Costs below are estimates from the subset measurement above unless marked measured, and are
refined in each study's `PREDICTION.md`.

---

## Now or alongside M3 (cheap: reuse the per-period surfaces)

### S-1 Z-sign invariance, calibration half (`z-sign-invariance`)

- **Question:** the engine writes p(z) = Φ((Φ⁻¹(PD) − √ρ·z)/√(1 − ρ)), so a higher Z means better
  conditions (R-2). Z is integrated out of the likelihood, and its distribution is symmetric, so
  the calibration should not depend on that choice. Are PD̂, ρ̂, the log-likelihood and every
  interval (Wald, profile, bootstrap percentile) identical under z → −z, within rounding?
- **Experiment:** a test-only objective with the opposite sign (+√ρ·z), fitted to the 162 replay
  panels and to the subset. Compare every output with the parity fit. The integration rules must
  mirror as well (adaptive GH centred on the mirrored mode; the Gauss–Legendre sinh map for
  d ∈ {0, n} centred on the mirrored half-point), so this also tests that the rules have no
  hidden asymmetry.
- **Scope:** the calibration half only. The half where the sign matters, extracting Z_t and
  mapping macro effects to signs, is S-18 and S-19.
- **Cost:** under a minute.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** none needed unless it fails.

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
- **Prediction:** not registered. **Result:** not run. **Mitigation:** backlog item since D-131.

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

### S-5 BCa intervals (`bca-intervals`)

- **Question:** do bias-corrected and accelerated intervals fix the 125 pinned percentile findings?
  The two diagnoses in D-137 give the expected split: they should help with "no bias or skewness
  correction" (68 verdicts) and not with "boundary breakdown" (57).
- **Experiment:** a `native` option. z₀ comes from the existing B = 999 replicates; the
  acceleration from the jackknife (as S-3). Report BCa next to percentile, pairwise.
- **Cost:** about the subset baseline; pinning ≈ 1.3–2 h.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** backlog item since D-137.

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

### S-7 Grid resolution against accuracy and runtime (`grid-resolution`)

- **Question:** how do estimate accuracy (after sub-grid refinement), SE accuracy, interval end
  points and runtime depend on the grid? Is 61 × 41 (D-115) near the knee?
- **Experiment:** the subset at 31 × 21, 61 × 41, 121 × 81 and 241 × 161, against the finest grid as
  reference, pairwise per replicate.
- **Cost:** the surface cost grows with the number of points, so about 1 h at R = 1,000 or 15 min
  at R = 200 (enough, since this compares per-replicate differences, not coverage).
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

## With M3

### S-8 MLE against method of moments (`mle-vs-mom`)

- **Question:** what is MoM's efficiency relative to the MLE (relative RMSE, and bias), across the
  full matrix?
- **Experiment:** MoM (joint-default-probability form) on every recovery panel, compared pairwise
  with the MLE.
- **Cost:** minutes for MoM itself; ≈ 1 h if the MLE is refitted for pairing rather than read from
  the goldens.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-9 Frequentist coverage of grid-Bayesian credible intervals (`bayes-coverage`)

- **Question:** do 95% credible intervals under flat and Jeffreys priors cover at 95% in repeated
  sampling, and where do they differ from the profile intervals?
- **Experiment:** the grid posterior on the recovery panels; equal-tailed and HPD intervals. Grid
  resolution matters for a discretised posterior, so S-7 informs the grid.
- **Cost:** about the fit cost (≈ 5 min subset, ≈ 1 h full matrix), plus the Jeffreys prior once
  per n.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

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

### S-11 Misspecification (`misspecification`)

- **Question:** what bias and coverage does the standard Vasicek fit have when the data come from a
  t-copula (a common scale mixture on the factor), an AR(1) factor, or a beta mixture?
- **Experiment:** one parity fit per replicate from each alternative DGP; the same summaries as the
  recovery harness. Each DGP needs a prose description, a Python mirror and a hash, as `dgp/` does
  (D-110–D-112).
- **Cost:** one recovery run per DGP: ≈ 7–12 min on the subset, ≈ 1.3–2 h for the full matrix. The
  blocker is the DGP extensions, not compute.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

## After M4 (GPU scale)

### S-12 Double (iterated) bootstrap (`double-bootstrap`)

- **Question:** does calibrating the nominal level with an inner bootstrap bring iid bootstrap
  intervals into the band, away from the bounds?
- **Experiment:** nested iid bootstraps. An inner resample of an outer resample is again a weight
  vector over the original periods, so the whole study is W × L with no refits.
- **Cost:** B × C = 999 × 199 grid reductions per replicate: ≈ 3 h on the subset on CPU, ≈ 1 day
  for the full matrix on CPU, which is the GPU's fused-reduce kernel's job.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-13 R = 10,000 recovery re-run (`recovery-r10000`)

- **Question:** which borderline verdicts survive a Monte Carlo band about three times narrower?
- **Experiment:** the recovery harness at R = 10,000 for the scenarios whose verdicts lie within
  about 0.01 of a band edge, with the replicates 0–999 unchanged.
- **Cost:** the whole matrix ≈ 12–21 h on 24 threads on CPU; the borderline scenarios alone
  ≈ 3 h.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-14 Backtest power (`backtest-power`)

- **Question:** how many years of data does each M5 backtest need to detect a PD misstated by 20%
  or 50%, at a given ρ and n?
- **Experiment:** power under the one-factor model. For tests on per-period or summed defaults it
  is exact by quadrature and convolution; simulation is the cross-check.
- **Cost:** minutes. Needs the M5 backtests.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-15 Simulation-based calibration of the Bayesian estimator (`bayes-sbc`)

- **Question:** are the rank statistics of the true parameter within the grid posterior uniform
  (Talts et al. 2018)?
- **Experiment:** draw (PD, ρ) from the prior, simulate a panel, compute the posterior rank; a few
  (n, T) settings. Ties on a discrete grid need a stated tie-breaking rule.
- **Cost:** thousands of calibrate-only fits: minutes on CPU.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-16 FP32 search, FP64 finalisation (`fp32-search`)

- **Question:** does evaluating the grid in FP32 (with FP64 accumulation, D-039) and refining in
  FP64 give the same estimates, SEs and intervals as pure FP64, within the documented tolerances?
- **Experiment:** the subset under both precision policies, on the CPU first (the policies are
  host and device code), then on the GPU, so that precision and device-math effects are separated.
- **Cost:** about the subset baseline per policy per backend.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-17 Performance across GPU generations (`gpu-scaling`)

- **Question:** how do surface and fused-reduce throughput scale across GPU generations (at least
  sm_89 and sm_120)?
- **Experiment:** the `perf/` suite on each architecture. A performance benchmark rather than a
  methodology question: its results live in `perf/`, with this entry as the pointer.
- **Cost:** hours per GPU; needs access to hardware beyond the development machine's sm_120.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

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
- **Overlap:** the same jackknife pass serves S-3 and S-5, so the three should share one run. It is
  also the first slice of M6's window and influence analysis, which then extends it to windows.
- **Verdicts it could change:** none; it is descriptive.
- **Cost:** ≈ 10–15 min on the subset including the fits; ≈ 3–5 min on top of an S-3/S-5 run
  (leave-two-out at T = 100 costs about five iid bootstraps). Full matrix ≈ 1.5–2.5 h. The exact
  refits add ≈ 10 min.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

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
- **Experiment, two arms, each against parity on the same panels:**
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
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

## Addendum: small DGP variants (D-149)

Each variant is built from the existing binomial DGP and gets its own key domain, a prose
description, a Python mirror and a hash. S-24 and S-25 are misspecification studies of the same
kind as S-11 and share its estimand rule.

### S-24 PD heterogeneity (`pd-heterogeneity`)

- **Question:** each panel pools two sub-segments with different PDs, loading on the same factor
  with the same ρ (no extra factor correlation). How much does the unmodelled heterogeneity
  inflate ρ̂?
- **DGP change:** the second segment's Bernoulli draws need their own stream but the same Z_t.
  Two calls of the existing `simulate_panel` with the same key would share the uniforms as well
  (common random numbers, an artificial dependence), so a segment index joins the key.
- **Design:** the scenario's n split between the segments and their PDs set so that the pooled
  PD equals the scenario's; PD ratios of about 2, 5 and 10.
- **Verdicts it could change:** none pinned; it reports bias and coverage against the stated
  estimand.
- **Cost:** ≈ the baseline per ratio: ≈ 20–35 min on the subset for three ratios, ≈ 4–6 h for the
  full matrix. Plus the DGP extension and its mirror.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-25 PD trend (`pd-trend`)

- **Question:** PD drifts smoothly over the sample. How much does the trend inflate ρ̂, and does
  detrending fix it?
- **DGP change:** a PD per period (today's `PanelSpec` takes one PD).
- **Design:** a trend linear in probit(PD) with the scenario's PD at mid-sample, two slopes. Three
  fits per panel: naive; **oracle-detrended** (the true trend as a per-period probit offset, an
  upper bound on what detrending can do); and **estimated-detrended** (the trend estimated first,
  then ρ fitted with the offsets fixed, which ignores the trend's own estimation error).
- **Needs:** the detrended fits need an objective with per-period probit offsets, a `native`
  objective. A trend is a covariate, so this is the structure M7's specification fits need; a joint
  fit of trend, PD and ρ is a three-parameter problem and belongs there.
- **Verdicts it could change:** none pinned; bias and coverage against the stated estimand.
- **Cost:** ≈ the baseline per slope and fit arm: ≈ 45–70 min on the subset, ≈ 8–12 h for the full
  matrix (CPU, overnight). Plus the DGP extension and the offset objective.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-26 Varying portfolio size (`varying-n`)

- **Question:** n_t changes over time, for example a steady decline as in bank-count data. Does
  anything in the estimator or intervals implicitly assume a stable n?
- **No DGP change needed:** the DGP already takes an n per period, and its stream addressing
  already guarantees that changing one period's n leaves the others unchanged (M1.6). The variant
  is a scenario definition, not new generator code, so no new mirror or hash is required.
- **Design:** for each scenario, n_t declining linearly from 1.5·n to 0.5·n, so the mean n is the
  scenario's and results pair with the constant-n scenario.
- **Where a stable n is assumed:** D-122's deduplication only saves time when n repeats, so fits
  cost more but must be bitwise unaffected; the iid bootstrap resamples periods of different sizes;
  walk-forward windows differ in information. Among M3's estimators, method of moments needs a
  weighting across n and the Vasicek-rate MLE treats every period's rate alike, so they are the
  likelier places for the assumption; S-26 is repeated for them in M3.
- **Verdicts it could change:** every verdict of the paired constant-n scenario.
- **Cost:** without deduplication, fits cost ≈ 1.5–3 × the baseline: ≈ 15–35 min on the subset,
  ≈ 3–6 h for the full matrix.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

## Addendum: with M3 (D-149)

### S-27 Large-portfolio approximation (`large-portfolio`)

- **Question:** at what n does the Vasicek-rate MLE become indistinguishable from the binomial MLE,
  in estimates and in coverage? Deliver a rule of thumb by PD and ρ.
- **Experiment:** both estimators on the same panels; the matrix's n ∈ {100, 10³, 10⁴} is extended
  with n ∈ {10⁵, 10⁶}. "Indistinguishable" is defined in the prediction, for example a mean
  difference below 0.1 SE and coverage that differs by less than the Monte Carlo band allows.
- **Depends on S-28:** at low PD and moderate n many periods have no defaults, where the rate model
  is undefined, so S-27's answer depends on the zero-default treatment. Run S-28 first, or together.
- **Overlap:** the same estimator-comparison pass as S-8.
- **Cost:** the Vasicek-rate MLE is closed form (seconds). The binomial fits at the new n, without
  bootstrap, ≈ 1.5–2 h for the full matrix and ≈ 10–20 min for the subset's 18 new scenarios; for
  n ≤ 10⁴ the goldens are reused.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

### S-28 Zero-default treatments for rate-based estimators (`zero-default-rates`)

- **Question:** for the Vasicek-rate MLE and method of moments on rates, compare refusing (parity,
  D-044), dropping zero-default periods, and the censored likelihood (native, D-044) by bias, RMSE
  and coverage.
- **Design notes:** refusal produces no estimate, so its row is the refusal rate, not a bias.
  Dropping is a data edit, so it can only be an explicit `native` option, never silent (D-044).
  The censored likelihood uses the Vasicek rate CDF, which is closed form, at a censoring point stated
  in the prediction. The binomial MLE, which handles zero defaults exactly, is the benchmark row.
- **Verdicts it could change:** none of the binomial verdicts; it adds the rate-based estimators'
  verdicts for each treatment.
- **Cost:** closed-form likelihoods: minutes for the full matrix.
- **Prediction:** not registered. **Result:** not run. **Mitigation:** n/a.

## Deferred to M7 (recorded together, not started)

### S-18 Z sign convention and macro sign filters (`z-sign-macro`)

- **Question:** R-2. Expected macro effects are declared in economic terms ("worsens" or "improves
  credit conditions") and mapped to coefficient signs in one place, through the engine's
  convention (higher Z = better conditions). Is the mapping right?
- **Experiment:** a synthetic DGP with known macro effects of both signs; the sign filter must keep
  exactly the specifications with the right economic direction.
- **Status:** deferred to M7; not started.

### S-19 Z_t extraction (`z-extraction`)

- **Question:** E[Z_t | d_t] under the fitted parameters, as a standard output: is it unbiased for
  the simulated Z_t, and is its sign consistent with S-1 and S-18?
- **Experiment:** posterior mean of each period's factor on DGP panels whose factors are known
  (`simulate_panel` returns them).
- **Status:** deferred to M7; not started.

### S-20 Belkin–Suchower–Forest apply function (`bsf-apply`)

- **Question:** a reference function taking a through-the-cycle migration matrix, ρ and z to the
  conditional migration matrix, with a test that fixes its conventions (the sign of z, the order of
  rating rows, the default column).
- **Status:** deferred to M7; not started. Migration matrices are outside the methodology scope
  listed in CLAUDE.md, so this needs a scope decision first.
