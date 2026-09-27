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
