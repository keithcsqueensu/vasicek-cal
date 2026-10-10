# S-35 Long-run PD under serial correlation: predictions made before the run

Pre-registration for study S-35 (`studies/README.md`, D-148, D-188, D-190). This file is committed
**before** any panel of the AR(1)-factor DGP has been fitted. The results will be compared with it
prediction by prediction, and misses will be reported, not explained away.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a new,
dated section appended at the end, committed before the run it concerns. The AR(1)-aware MLE (M6) is not
part of this registration; its predictions are appended before it is run.

**Evidence used in writing it:** the model and pinned results only, through
[`reference.py`](reference.py) (outputs `reference.csv`, `reference_rate.csv`): the pinned φ = 0 coverage
of the recovery summary and of the Vasicek-rate recovery summary, and the variance inflation that an
AR(1) factor implies for each estimator. No AR(1) panel was simulated for it.

## The question

When the systematic factor is serially correlated, how far off are the long-run PD by simple averaging
of default rates and the MLEs that assume independent periods, and do their intervals still cover?
Miu and Ozdemir (2008) find that the simple average tends to understate the long-run PD and that an MLE
that models the serial correlation does better (`docs/references.md`).

## What is run

- **The DGP** (D-190, `docs/methodology/dgp.md`, "The AR(1) factor"): Z_t = φ Z_{t−1} + √(1 − φ²) e_t,
  Z_0 = e_0, with e_t and every Bernoulli draw the base DGP's. φ ∈ {0, 0.3, 0.6}.
- **The panels:** the recovery matrix, 81 scenarios × 1,000 replicates, with the recovery seed, scenario
  ids and replicates. Every φ uses the same random numbers, so comparisons across φ are paired. φ = 0
  is the pinned run's panels bit for bit.
- **Estimators, per panel:**
  - **The simple average** ā = (1/T) Σ d_t/n, with the t interval ā ± t_{T−1, 0.975} s/√T (s the sample
    standard deviation of the d_t/n, T − 1 denominator).
  - **The binomial-mixture MLE** (parity), fitted exactly as the recovery harness fits it: estimate, SEs,
    Wald and t intervals, profile intervals, the iid bootstrap and q.
  - **The Vasicek-rate MLE** on the n → ∞ rates, as its own recovery harness builds them: the rates
    r_t = p(Z_t) of the AR(1) factor of the n = 10,000 scenario with the same (PD, ρ, T), fitted with
    profile intervals. This isolates serial correlation from binomial noise. It is 27 (PD, ρ, T) cells ×
    1,000 replicates per φ.
- **Targets.** Under the AR(1) factor every period keeps its φ = 0 distribution. The expected
  log-likelihood of an iid model therefore depends only on that marginal, and its maximiser, the
  pseudo-true value (D-150, D-151), is the generating (PD, ρ). The simple average's mean is PD exactly. So
  coverage is judged against the generating values, which are the pseudo-true values here.
- **Checks before scoring:** at φ = 0, the binomial MLE's per-scenario counts (covered, profile covered,
  bootstrap covered, flagged, for PD and ρ) must equal the pinned recovery summary's exactly, and the
  rate MLE's covered counts its pinned summary's. Either failing stops the scoring.
- **Order:** the study subset (D-150) first, scored as a subset, then the full matrix. The runtime
  projection (D-187) is recorded before the full run starts.

## Mechanism behind the predictions

- **Nothing moves on average, but the variance of a mean over T dependent periods does.** Every
  estimator here is a function of period-by-period averages: of the default rates, of the per-period
  scores, of the probit rates. With lag-k correlation ρ_k between periods, the variance of such a mean is
  inflated over the iid value by f_T = 1 + 2 Σ_{k<T} (1 − k/T) ρ_k. Every interval here takes its width
  from the iid value, so it is too narrow by √f_T.
- **The inflation per estimator** (`reference.py`). For the binomial MLE the per-period score's dependence
  comes only through E[score | Z_t], whose Hermite expansion gives its lag-k covariance as
  Σ_j a_j a_j′ φ^{jk}. The binomial noise is independent across periods and dilutes the dependence where
  it dominates (small n·PD).
  - **PD-type quantities** are carried mostly by the first Hermite term, so f approaches (1 + φ)/(1 − φ)
    where the factor dominates. Over the informative groups B–D, f_PD is 1.35–1.81 at φ = 0.3 and
    2.31–3.82 at φ = 0.6.
  - **ρ** is carried mostly by the second term (correlation φ^{2k}): f_ρ is 1.04–1.19 and 1.26–2.09.
  - **Group A** (near-uninformative or near a bound) is diluted by binomial noise: f_PD is 1.02–1.56 and
    1.06–2.96, f_ρ 1.00–1.12 and 1.00–1.70.
- **Longer history does not repair it.** f_T tends to (1 + φ)/(1 − φ) as T grows. The interval narrows
  like 1/√T, and so does its true spread; their ratio does not change.
- **ρ is biased further down.** For the rate MLE, σ̂² = (1/T) Σ (x_t − x̄)² of the probit rates has mean
  σ²(1 − f_μ/T) exactly, against σ²(1 − 1/T) at φ = 0. The extra downward bias is (f_μ − 1)/T: 4% and 13%
  at T = 20 (φ = 0.3, 0.6), 0.8% and 2.9% at T = 100. In ρ that is about (1 − ρ) times as much.
- **The simple average** is unbiased in mean at every φ. Its distribution is right-skewed (the default rate
  is), and fewer effective periods leave more of that skew in the mean. So it falls below PD more often
  as φ grows, which is the sense in which it "understates" the long-run PD. L2's size is a rough
  heuristic, held loosely:
  - a period's default rate at PD 1% and ρ 0.12 has skewness of about 4;
  - the mean of T/f effective periods has skewness of about 4·√(f/T);
  - P(mean < its expectation) ≈ 0.5 + 0.066 × that skewness;
  - so about 0.56 at φ = 0 and 0.62 at φ = 0.6 for T = 20.
- **Coverage mapping.** A scenario's pinned φ = 0 coverage c₀ gives its effective critical value
  z₀ = Φ⁻¹((1 + c₀)/2). With inflation f and an extra bias of δ standard deviations, the predicted
  coverage is Φ(z₀/√f − δ) − Φ(−z₀/√f − δ). δ is applied to ρ only.

## Predictions

Coverage is pooled over all replicates of a group's scenarios (groups A–D as in
`docs/methodology/bootstrap_predictions.md`). "Fall" is a coverage at φ minus its coverage at φ = 0 on the
same panels. The band is the profile band of the recovery summary (lower edge 0.9273 at R = 1,000).

| # | Prediction |
|---|---|
| L1 | **The simple average is unbiased at every φ:** the mean of ā over the replicates is within 3 Monte Carlo SEs of PD in at least 160 of the 162 (scenario, φ > 0) cells |
| L2 | **It falls below PD more often as φ grows:** P(ā < PD) is higher at φ = 0.6 than at φ = 0 in at least 40 of the 50 informative scenarios (groups B–D); its mean over group B at φ = 0.6 is 0.56–0.70 |
| L3 | **The rate MLE's σ̂² bias is the closed form:** the mean of σ̂²/σ² (σ̂² = ρ̂/(1 − ρ̂)) is within 3 Monte Carlo SEs of 1 − f_μ/T in at least 50 of the 54 (cell, φ > 0) pairs, with box-limited fits counted and reported |
| L4 | **The rate MLE's profile coverage follows the inflation:** pooled per T and φ, PD and ρ coverage are each within ±0.03 of `reference_rate.csv`'s means: PD 0.86 (φ = 0.3) and 0.70 (φ = 0.6) at every T; ρ 0.91–0.92 and 0.79–0.82 |
| L5 | **The binomial MLE's profile coverage follows the inflation:** pooled per group (B, C, D) and φ, within ±0.03 of `reference.csv`'s group means: PD 0.86 / 0.87 / 0.88 at φ = 0.3 and 0.71 / 0.73 / 0.73 at φ = 0.6; ρ 0.92 / 0.93 / 0.93 and 0.81 / 0.84 / 0.86 |
| L6 | **Diluted where binomial noise dominates:** group A's pooled profile PD coverage falls less than group D's, at both φ |
| L7 | **Every iid interval falls alike:** for the Wald, t and iid-bootstrap intervals, each group's fall (B, C, D; PD and ρ; both φ) is within ±0.04 of the profile interval's fall. Resampling periods independently ignores the dependence as the Hessian does |
| L8 | **Longer history does not repair it:** at φ = 0.6, the binomial MLE's pooled profile PD coverage differs by less than 0.04 between groups B (T = 20) and D (T = 100), and the rate MLE's by less than 0.04 between T = 20 and T = 100 |
| L9 | **Verdicts:** at φ = 0.3, all 50 informative PD profile verdicts are below the band, and 15–40 of the 50 ρ profile verdicts. At φ = 0.6, all 50 PD and at least 45 of the 50 ρ verdicts are below the band |
| L10 | **ρ is biased further down:** the paired change in the binomial MLE's mean relative bias of ρ̂ from φ = 0 to φ = 0.6 is −0.06 to −0.17 pooled over group B and −0.01 to −0.05 over group D |
| L11 | **The simple average's t interval follows its inflation:** pooled per group (B, C, D) and φ, its coverage is within ±0.04 of the mapping applied to its own φ = 0 coverage in the same run, with `reference.csv`'s f_avg |

Reported without a prediction: the median relative error of the binomial MLE's PD̂ and of ā by φ, which
is the comparison Miu and Ozdemir draw; the q interval's coverage; and the flagged and box-limited
counts by φ.

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not; a prediction that holds in direction but
misses in size is not held. The rows are computed by a script committed before the full run's results
exist (`compare.py`). The results go into `studies/README.md` with a D-entry and a `MANIFEST.json`. The
per-replicate rows are attached to a release, not committed (D-163). The verdict tables per φ answer "the
verdicts it could change" (D-149) for every PD and ρ verdict of every interval.

**Monitoring implication (to be confirmed by the result, D-152):** if L5 to L9 hold, a PD interval
computed on default history is too narrow by up to about √((1 + φ)/(1 − φ)) when the history is serially
correlated, whatever its length. A residual autocorrelation check belongs before any interval is
reported.

## Addendum, 2026-10-09: a second rate-MLE arm on the observed rates (owner), registered before any run

Appended before any run of S-35; nothing above changes. The rate MLE on the n → ∞ rates stays the
primary rate arm: it isolates serial correlation. The owner asked for a second arm on the rates as a
practitioner observes them, as Miu and Ozdemir's (2009) first step applies the rate model. It shows how
serial correlation and granularity error combine.

- **The arm:** the Vasicek-rate MLE on the observed rates d_t/n of the same count panels, with S-28's
  censored treatment of zero-default periods (detection limit 1/(2n), D-044). It is fitted exactly as the S-8
  pass fitted it (`study_estimator_pass`): one refusing fit where no period has d = 0 or d = n, since the
  two likelihoods are the same there, and otherwise the censored likelihood. The estimate is the polished
  maximum, the interval is the profile interval, and coverage is of the generating PD and ρ.
- **Check before scoring:** at φ = 0 the arm's estimates and coverage indicators must equal the pinned S-8
  rows (`studies/mle-vs-mom/fits.parquet`, the `censor_*` columns) bit for bit, replicate by replicate.
- **Reference numbers** ([`reference_observed.py`](reference_observed.py), `reference_observed.csv`) come
  from those pinned rows only. Unlike the n → ∞ arm, this estimator is biased at φ = 0: S-28 found ρ̂
  biased up where n·PD is small, because binomial noise is read as factor variance. So the mapping
  carries each scenario's own φ = 0 bias b₀ and spread s₀ (logit) forward:
  - the interval's half-width w is the one that reproduces the pinned coverage c₀ for a normal estimate
    with bias b₀ and spread s₀;
  - at φ the spread becomes s₀√f, with f the binomial MLE's inflation for the parameter (`reference.csv`);
  - ρ's bias moves down by (f_PD − 1)/T in logit, the extra downward bias of a variance about the mean of T
    dependent periods;
  - PD's bias does not move.
- **What it implies:** where granularity error dominates (n = 100, ρ̂ biased up by about 2.7 spreads at
  φ = 0), serial correlation's downward pull on ρ̂ offsets its wider spread, and ρ coverage barely
  moves. Where granularity error is small (n = 10,000) ρ coverage falls as in the n → ∞ arm.

| # | Prediction |
|---|---|
| L12 | **PD coverage follows the inflation from its own starting point:** pooled per group and φ, within ±0.04 of `reference_observed.csv`'s means: A 0.84 / 0.77, B 0.87 / 0.73, C 0.89 / 0.74, D 0.88 / 0.74 (φ = 0.3 / 0.6) |
| L13 | **ρ coverage, pooled per n and φ,** is within ±0.04 of the reference means: n = 100: 0.51 / 0.51; n = 1,000: 0.69 / 0.65; n = 10,000: 0.84 / 0.76 |
| L14 | **The two errors partly cancel for ρ:** at n = 100 the pooled ρ coverage at φ = 0.6 is within 0.03 of its φ = 0 value, while at n = 10,000 it is at least 0.06 below it |
| L15 | **Serial correlation pulls ρ̂ down here too:** the paired mean of logit ρ̂ at φ = 0.6 is below its φ = 0 value in at least 45 of the 50 informative scenarios |

## Addendum, 2026-10-09: scoring fixed before the subset run

Appended before the subset run; no prediction, threshold or definition above changes.

- **The scoring script** is [`compare.py`](compare.py), committed with this note. It runs the checks first and stops on any failure:
  - the binomial and rate MLEs' φ = 0 counts against the pinned summaries;
  - the observed-rate arm's φ = 0 rows against S-8's, bit for bit.
- **Pooling:** per-scenario coverages are averaged over the scenarios pooled, R being equal. Wald and t coverage is among unflagged fits, as pinned; profile, bootstrap and every rate-MLE coverage is over all replicates, with a not-computed interval missing. "Below the band" is coverage below the scenario's pinned profile band (`profile_band_lo_hex`).
- **The subset run** (the study subset of D-150, 9 scenarios × 1,000) is scored by the same script and reported as a subset, not as the study's result.
- **Disclosure:** the tool was run on 3 replicates of every scenario at each φ, to check the code and the pipeline: the φ = 0 checks passed, and the script ran to completion. Only its exit code and row count were read; those rows are not results. Before that, replicates 0–2 at φ = 0 alone were compared with S-8's rows (243 of 243 identical), and 2 replicates at each φ ran for the timing sample, whose rows were not read.

## Addendum, 2026-10-10: the M6 arm, the AR(1)-aware binomial MLE, registered before any of its fits on S-35's panels

Appended before the M6 arm has fitted any of S-35's panels; nothing above changes. This arm tests one
mechanism behind Miu and Ozdemir's (2008) comparison, not the comparison as a whole. Their claim sets simple
averaging against their MLE as estimators of the long-run PD under correlated defaults, and part of
averaging's understatement comes from nonlinearity: the mean of default rates is not the PD implied by the
latent factor, whatever the dependence over time. The question here is narrower. Given a binomial MLE that
already handles that nonlinearity, does also modelling the serial correlation (AR(1)) stop the long-run PD
being understated, and do the intervals cover again?

- **The estimator** (D-195, D-196; `benchmark_models/ar1/`): the binomial-mixture likelihood under an AR(1)
  factor, by a forward filter over each period's parity nodes placed for the law of Z_t given its
  neighbours (forward, backward, forward). It is maximised over (PD, ρ, φ) in the recovery box, with
  |φ| ≤ 0.95: a φ scan {−0.6, −0.3, 0, 0.3, 0.6, 0.9} with 2-D quasi-Newton in (PD, ρ), then a 3-D
  polish from each local maximum. Profile intervals are 3-D, with D-129's threshold and truncation at
  bounds (D-130).
- **What is run:**
  - **Full matrix**, the same panels (81 scenarios × 1,000 replicates × φ ∈ {0, 0.3, 0.6}): the point
    estimate, its flags and checks, and profile intervals for PD and ρ.
  - **The study subset** (D-150: 29, 37, 72, 4, 68, 49, 7, 43, 51), the same replicates: also the profile
    interval for φ (owner: secondary to PD; if its χ² calibration fails, that becomes its own registered
    study).
  - **The parametric-bootstrap comparator** (owner; D-043), on the study subset, replicates 0–99, φ = 0.6
    only: B = 199 panels drawn from the AR(1) DGP at the panel's own estimate (PD̂, ρ̂, φ̂). Panel b of
    (scenario s, replicate r) uses seed 0xB0071A61, scenario s and replicate 256 r + b. Each is refitted for
    the point estimate only; the intervals are percentile intervals for PD and ρ.
- **Checks before scoring.** Any failure stops the scoring.
  - Every fit's log-likelihood at the φ = 0 scan point is at least the binomial MLE's on the same panel
    minus 1e-8, the same function being maximised continuously.
  - Every fit's maximum is at least its φ = 0 scan point's minus 1e-10.
  - The placement and doubled-rule checks are reported. The number of fits flagged by either is a
    result, not a filter.
- **Runtime:** projected from a timing sample before the full run (D-187, D-193). The projection is
  recorded here before the run starts; the estimate is about 18 h on 24 threads.

### Disclosure: what was seen before this was written

The AR(1) MLE was run during its development, on panels from seeds disjoint from S-35's panels, the
validation run and the information computation. Exactly:
- **The search trial:** seed 0xA11CE5EED, scenario 900, replicates 0–3. Fits at PD 1%, ρ 0.12 for n ∈
  {100, 1,000} × T ∈ {20, 100} × true φ ∈ {0, 0.6}. These are the parameters of S-35's scenarios 36 (A),
  37 (B), 42 (D) and 43 (D). Seen: the mean φ̂ per cell: at n = 100, T = 20, −0.40 (true 0) and −0.16
  (true 0.6), one fit in each at the −0.95 bound; at n = 100, T = 100, −0.19 and 0.39; at n = 1,000, T =
  20, −0.25 and 0.15; at n = 1,000, T = 100, −0.03 and 0.55. Also seen: run time and evaluation counts.
- **The profile trial:** seed 0xA11CE5EED, scenario 901, replicates 0–2, the same four cells at true
  φ = 0.6. Seen: profile intervals for PD, ρ and φ. The true value was inside 32 of 36 intervals; the φ
  intervals of several short panels ran to a bound.
- **The owner's check** of the n = 100, T = 20 cell: the 8 panels of the search trial against ref's
  independent profile. ref's profile maximum was at core's φ̂ in all 8, so the search is right. The profiles
  were nearly flat: 0.01 to 2.4 nats between φ = 0 and the maximum, and the 95% profile set spanned all of
  [−0.95, 0.95] in two panels. In 2 of the 4 panels at true φ = 0.6 the set excluded 0.6.
- **Likelihood only, no fits:** the filter's accuracy cases (adjacent periods with very different modes,
  D-195) and its timing.

The predictions below rest on the theory stated with them, not on these panels. Where a prediction concerns
a cell above (the φ̂ < 0 and bound-hit rates at short T, and the φ interval's calibration there), it was
written knowing what is listed here.

### Model-based reference numbers: where asymptotics can be trusted

[`information_m6.csv`](information_m6.csv) is the expected Fisher information at the true parameters of
every scenario and φ (`study_ar1_information`, committed before it was run). It is the mean of the score's
outer product over 500 panels per cell, from a fresh seed 0x1F020261010. It contains no fits and no
estimates. The credibility thresholds were stated in the tool before its first run:
- **φ:** credible where the information-based SE ≤ 0.15;
- **PD and ρ:** credible where 1.96 SE in logit units ≤ the logit distance to the nearer box bound.

Where a parameter is not credible, its predictions below are qualitative (owner).

- **PD** is credible in every scenario at every φ; **ρ** in 17 of group A's 31 and in all of B–D.
- **φ** is credible in no scenario at T = 20, whatever n. Even at n = 10,000 its SE is 0.18–0.23, about the
  observed-factor value √((1 − φ²)/T); at small n·PD it reaches 11. At T = 40 it is credible in 1 of 17
  (group C) at φ = 0.3 and 11 of 17 at φ = 0.6. At T = 100 it is credible in 15–20 of group D's 21.
- **The score's sample mean** is within 3.0 Monte Carlo SEs of 0 in all 729 entries, as at the truth it
  must be.

### Mechanism

- **PD coverage returns where φ is estimated well.** With the dependence modelled, the profile interval is
  asymptotically correct. Its width now carries the variance factor (1 + φ̂)/(1 − φ̂), which the binomial
  MLE's interval lacks; that factor was 4.0 at φ = 0.6.
- **At short T the return is partial.** φ̂ is biased down: Kendall's −(1 + 3φ)/T is −0.14 at T = 20 and
  φ = 0.6, before the extra attenuation of a factor seen through binomial noise. A smaller φ̂ gives a
  narrower PD interval. At φ̂ ≈ 0.46 the factor is 2.7 against 4.0, and the width ratio √(2.7/4) = 0.82
  maps a nominal 0.95 to about 0.89. φ̂'s spread at T = 20 (SE ≥ 0.18) lowers it further: φ̂ near 0 halves
  the width.
- **The point estimate barely moves (Grenander).** For a mean with AR(1) errors, the equal-weight mean is
  asymptotically as efficient as generalised least squares; the efficient weights differ only near the
  ends of the series. PD̂ is driven by an average of period-level information, so the AR(1) MLE's PD̂ stays
  close to the binomial MLE's. Its median understatement, which comes from the skewness of a mean over
  few effective periods, is a property of the data, not of the weights, and remains. The model fixes the
  interval, not the point.
- **φ̂ at small n·PD is weakly identified.** About one default per period carries almost no information
  about Z_t, so the profile in φ is flat (the information above), and φ̂ is often at a bound.
- **The χ² calibration of the φ interval** needs the information to be large; where it is not, the profile
  set is set by a flat profile and coverage is not nominal.

### Predictions

Pooling, groups and the band are as above. "Binomial" means the binomial MLE's result on the same panels
in S-35's full run.

| # | Prediction |
|---|---|
| L16 | **PD coverage returns at T = 100:** the AR(1) MLE's pooled PD profile coverage in group D is 0.92–0.96 at both φ = 0.3 and φ = 0.6 |
| L17 | **Partly at short T:** pooled PD profile coverage at φ = 0.6 is 0.82–0.92 in group B and 0.87–0.95 in group C; at φ = 0.3, 0.88–0.95 in B and 0.90–0.96 in C |
| L18 | **Above the binomial MLE's wherever it matters:** at φ = 0.6, pooled PD profile coverage exceeds the binomial MLE's by at least 0.10 in each of groups B, C and D |
| L19 | **Little cost when there is no correlation:** at φ = 0, pooled PD profile coverage is within −0.03 to +0.01 of the binomial MLE's in each of groups B, C and D |
| L20 | **ρ coverage returns too:** pooled ρ profile coverage at φ = 0.6 is 0.88–0.96 in each of groups B, C and D |
| L21 | **The point estimate does not rescue the median (Grenander):** at φ = 0.6, the AR(1) MLE's median relative error of PD̂ is closer to 0 than the binomial MLE's in between 19 and 31 of the 50 informative scenarios, and differs from it by less than 0.01 in at least 40 of the 50. This tests one mechanism only: whether modelling the AR(1) structure, on top of an MLE that already handles the nonlinearity, moves the median point estimate. If L21 holds, the conclusion is "modelling serial correlation fixes the interval but not the median bias", not "Miu and Ozdemir's claim fails". Their comparison of averaging against an MLE also rests on the nonlinearity, which S-35's main arms address separately. The AR(1) model would be shown to move the median by 32 or more of 50 (a one-sided sign test at about 3%) |
| L22 | **φ̂ < 0 at true φ = 0, where φ is credible:** over the credible scenarios of group D, the pooled share of φ̂ < 0 is 0.50–0.60 (Kendall's bias −1/T against an SE of about 0.11 gives 0.54) |
| L23 | **φ̂ < 0 at true φ = 0, short T (qualitative; φ not credible):** over group B the pooled share is above 0.5 (0.55–0.80), and over group A's T = 20 scenarios it is higher than over group B's |
| L24 | **φ̂ < 0 at true φ = 0.6:** in group D's credible scenarios under 0.5% (φ̂ ≈ 0.57 ± 0.10); in group B (qualitative) between 1% and 15% |
| L25 | **φ̂ at a bound (±0.95), reported apart from the sign:** in group D's credible scenarios under 0.5% at both φ = 0 and φ = 0.6. Qualitative where φ is not credible: in group B between 1% and 20% at each φ, and in group A's T = 20 scenarios at least 10% at each φ |
| L26 | **The φ interval (study subset), where φ is credible** (φ = 0.6: 68, 49, 43, 51; φ = 0 and 0.3: 43): coverage at least 0.92 in each (R = 1,000, Monte Carlo SE 0.007) |
| L27 | **The φ interval where φ is not credible (qualitative):** at φ = 0.6 its coverage is below 0.93 in at least 2 of the three T = 20 subset scenarios (29, 37, 72). The χ² calibration fails at short T |
| L28 | **The bootstrap comparator** (study subset, replicates 0–99, φ = 0.6): its percentile PD coverage is within ±0.05 of the profile interval's on the same replicates in at least 7 of the 9 scenarios. At 100 replicates a coverage near 0.95 has a Monte Carlo SE of 0.022 (about ±2.2 points), so differences under 0.044 are not read as results (owner) |

Reported without a prediction: the number of local maxima along the scan and the gap to the second (owner's
diagnostic); the placement- and doubled-check flag counts; the q interval; and φ̂'s mean and median per
scenario and φ.

## Addendum, 2026-10-10: the M6 arm's runtime projection (a fact, not a prediction)

Appended after the M6 registration above; no prediction, threshold, scope or definition changes. The
registration estimated the full M6 run at about 18 h on 24 threads. A timing sample (D-187, D-193) projects
it far higher:

- **Projected:** 3,968,311 thread-seconds, **46.1 h** on 24 threads at efficiency 0.997; 90% range
  42.4–50.0 h, 40.3–52.5 h with the ±5% machine-state allowance.
- **Method:** `perf/project_runtime.py --run` on `study_lrpd_serial_m6` (D-198) as of d68c364, with its
  warm-up and minimum sample wall time; accepted at 4 replicates per scenario (10% from the 2-replicate
  sample). The comparator was sampled at 8 refits per job, its seconds scaled to B = 199 (`--boot-sample`,
  which writes no rows). The sample's rows were never read.
- **Next:** P-18 (the filter's transition, with the shared in-house exp) is re-projected the same way. The
  arm's scale is decided after that number and before the validation run, since any reduction in replicates
  or scope would be an amendment, and amendments are allowed only before the validation run (D-197).
