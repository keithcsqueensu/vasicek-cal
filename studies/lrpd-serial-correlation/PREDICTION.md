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
