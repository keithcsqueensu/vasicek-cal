# The estimator-comparison pass (S-8, S-27, S-28): predictions made before the run

Pre-registration for studies S-8 (`mle-vs-mom`), S-27 (`large-portfolio`) and S-28
(`zero-default-rates`), which share one estimator-comparison pass (D-150). This file is committed
**before** M3's estimators exist and before any of them has run on a recovery panel. The results
will be compared with it prediction by prediction, and misses will be reported, not explained away.
The S-27 and S-28 entries point here.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a
new, dated section appended at the end, committed before the run it concerns. M3 implements the
estimators as defined below; if M3 has to define one differently, that change is appended here,
dated, before the pass runs.

**Evidence used in writing it:**

- the model alone: [`model_reference.py`](model_reference.py), which reads no recovery panel and
  runs no estimator of the pass on simulated panels, and gives every reference number below;
- results already pinned: the groups of D-136 and the binomial MLE's recovery summary.

## What is run

- **Panels:** the recovery panels, unchanged, for n ∈ {100, 10³, 10⁴} (81 scenarios). For S-27,
  54 new scenarios extend n to {10⁵, 10⁶} with the same PD, ρ and T: scenario id
  81 + 2·(9·i_PD + 3·i_ρ + i_T) + (i_n − 3), under the recovery seed and DGP (D-109). The study
  subset first (scenarios 29, 37, 72, 4, 68, 49, 7, 43, 51, and their n = 10⁵ and 10⁶
  counterparts), then the full set (135 × 1,000).
- **Estimators, on each panel:**
  - **Binomial MLE (the benchmark):** the parity fit and profile intervals, as pinned. At n ≥ 10⁵
    it is fitted the same way, without the bootstrap.
  - **Method of moments, joint-default-probability form (S-8):** PD̂ = Σd/Σn and
    PD̂₂ = Σd(d − 1)/Σn(n − 1), both exact in finite n; ρ̂ solves Φ₂(c, c; ρ) = PD̂₂ with
    c = Φ⁻¹(PD̂), within the box. When PD̂₂ ≤ PD̂², ρ̂ is set to the box's floor and flagged; when
    Σd = 0 it refuses. Zero-default periods need no treatment.
  - **Vasicek-rate MLE (S-27, S-28):** the likelihood of the observed rates d_t/n_t under the
    Vasicek rate distribution, maximised over the box off the grid. For data without zero or 100%
    rates its maximum is the closed form PD̂ = Φ(μ̂/√(1 + v̂)), ρ̂ = v̂/(1 + v̂), from the mean μ̂ and
    1/T variance v̂ of Φ⁻¹(d_t/n_t), when that lies in the box. Its 95% interval is the profile
    likelihood interval of the rate likelihood (threshold 1.9207, D-129). Three zero-default
    treatments (D-044):
    - **refuse** (parity): no estimate if any period has d_t = 0 or d_t = n_t;
    - **drop** (`native`, an explicit data edit): those periods are removed, and the fit refuses if
      fewer than 3 remain;
    - **censored** (`native`): a zero-default period contributes P(DR ≤ c_n) =
      Φ((√(1 − ρ)·Φ⁻¹(c_n) − Φ⁻¹(PD))/√ρ) with the censoring point **c_n = 1/(2n)** (half a
      default); a 100% period contributes P(DR ≥ 1 − c_n).
  - **Method of moments on rates** (for rate series, D-046): PD̂ the mean rate, and ρ̂ solving
    Φ₂(c, c; ρ) − PD̂² = var(d_t/n_t); zero rates enter as zeros. Reported, with the drop treatment
    beside it, without predictions.
- **Reported, per scenario and estimator:** bias and RMSE of PD̂ and ρ̂ over the replicates each
  estimator estimates, the refusal share, and, for the rate MLE, profile coverage under the D-131
  band and policy. Pairwise against the binomial MLE on the same panels: the mean difference of ρ̂
  in units of the binomial fit's root-mean-square Hessian SE, and the coverage difference.
- **"Indistinguishable" (S-27):** the rate MLE's mean ρ̂ is within 0.1 SE of the binomial MLE's (the
  SE as above), and its ρ profile coverage is within 0.0227 of the binomial's (the half-width of the
  D-131 band at R = 1,000).
- **One fitting run pins the results** (D-156's process change). Verdicts of the new estimators
  that fall outside the band are reviewed in a separate labels file, as in D-155.

## Mechanism behind the predictions

- **MoM against the MLE (large-n limit).** MoM matches the first two moments of the rate on the
  natural scale; the MLE uses the whole distribution on the probit scale. Rates are right-skewed,
  more so at low PD and high ρ, and a second moment on the natural scale is dominated by the few
  worst periods. So MoM loses efficiency for ρ, more as ρ rises and PD falls, and **the loss grows
  with T**, because MoM's variance does not approach the MLE's. It is also biased downwards. For PD
  both are the mean rate to first order. The large-n reference (2,000 draws per cell):

  | PD | ρ | T | RMSE ratio MoM/MLE, ρ | RMSE ratio, PD | bias of ρ̂: MoM | bias of ρ̂: MLE |
  |---|---|---|---|---|---|---|
  | 0.001 | 0.02 | 20 | 1.12 | 1.00 | -0.0017 | -0.0009 |
  | 0.001 | 0.02 | 40 | 1.18 | 1.00 | -0.0009 | -0.0004 |
  | 0.001 | 0.02 | 100 | 1.22 | 1.00 | -0.0004 | -0.0001 |
  | 0.001 | 0.12 | 20 | 1.44 | 1.02 | -0.0306 | -0.0060 |
  | 0.001 | 0.12 | 40 | 1.66 | 1.04 | -0.0199 | -0.0029 |
  | 0.001 | 0.12 | 100 | 1.96 | 1.03 | -0.0109 | -0.0010 |
  | 0.001 | 0.24 | 20 | 1.91 | 1.04 | -0.0935 | -0.0126 |
  | 0.001 | 0.24 | 40 | 2.20 | 1.11 | -0.0676 | -0.0061 |
  | 0.001 | 0.24 | 100 | 2.68 | 1.12 | -0.0419 | -0.0021 |
  | 0.01 | 0.02 | 20 | 1.08 | 1.00 | -0.0013 | -0.0009 |
  | 0.01 | 0.02 | 40 | 1.11 | 1.00 | -0.0007 | -0.0004 |
  | 0.01 | 0.02 | 100 | 1.13 | 1.00 | -0.0003 | -0.0001 |
  | 0.01 | 0.12 | 20 | 1.34 | 1.01 | -0.0191 | -0.0060 |
  | 0.01 | 0.12 | 40 | 1.49 | 1.02 | -0.0110 | -0.0029 |
  | 0.01 | 0.12 | 100 | 1.63 | 1.01 | -0.0053 | -0.0010 |
  | 0.01 | 0.24 | 20 | 1.63 | 1.03 | -0.0580 | -0.0126 |
  | 0.01 | 0.24 | 40 | 1.86 | 1.05 | -0.0360 | -0.0061 |
  | 0.01 | 0.24 | 100 | 2.12 | 1.04 | -0.0184 | -0.0021 |
  | 0.05 | 0.02 | 20 | 1.04 | 1.00 | -0.0011 | -0.0009 |
  | 0.05 | 0.02 | 40 | 1.06 | 1.00 | -0.0005 | -0.0004 |
  | 0.05 | 0.02 | 100 | 1.07 | 1.00 | -0.0002 | -0.0001 |
  | 0.05 | 0.12 | 20 | 1.23 | 1.00 | -0.0115 | -0.0060 |
  | 0.05 | 0.12 | 40 | 1.30 | 1.01 | -0.0060 | -0.0029 |
  | 0.05 | 0.12 | 100 | 1.33 | 1.00 | -0.0026 | -0.0010 |
  | 0.05 | 0.24 | 20 | 1.42 | 1.02 | -0.0314 | -0.0126 |
  | 0.05 | 0.24 | 40 | 1.55 | 1.02 | -0.0168 | -0.0061 |
  | 0.05 | 0.24 | 100 | 1.60 | 1.01 | -0.0073 | -0.0021 |

- **The rate MLE against the binomial MLE (S-27).** The rate model reads the binomial noise in d/n
  as factor variance, so ρ̂_rate is biased **upwards**, by about B/(s²√(2/T)) SE, where
  B = E[p(1 − p)/(n·φ(Φ⁻¹p)²)] is the binomial variance on the probit scale and s² = ρ/(1 − ρ).
  n·B does not depend on n, so the bias falls as 1/n; it grows with T, because the SE shrinks.
  Reference, bias in SE units at T = 20 / 40 / 100. The delta method is unreliable where n·p(z) is
  small over much of the factor's range, which is also where zero-default periods are common;
  PD 0.1% at ρ = 0.24 is the clearest case:

  | PD | ρ | n = 100 | n = 1,000 | n = 10,000 | n = 100,000 | n = 1,000,000 |
  |---|---|---|---|---|---|---|
  | 0.001 | 0.02 | 164.20 / 232.22 / 367.17 | 16.42 / 23.22 / 36.72 | 1.64 / 2.32 / 3.67 | 0.16 / 0.23 / 0.37 | 0.02 / 0.02 / 0.04 |
  | 0.001 | 0.12 | 83.36 / 117.89 / 186.40 | 8.34 / 11.79 / 18.64 | 0.83 / 1.18 / 1.86 | 0.08 / 0.12 / 0.19 | 0.01 / 0.01 / 0.02 |
  | 0.001 | 0.24 | 558.70 / 790.12 / 1249.29 | 55.87 / 79.01 / 124.93 | 5.59 / 7.90 / 12.49 | 0.56 / 0.79 / 1.25 | 0.06 / 0.08 / 0.12 |
  | 0.01 | 0.02 | 23.91 / 33.82 / 53.47 | 2.39 / 3.38 / 5.35 | 0.24 / 0.34 / 0.53 | 0.02 / 0.03 / 0.05 | 0.00 / 0.00 / 0.01 |
  | 0.01 | 0.12 | 7.01 / 9.92 / 15.69 | 0.70 / 0.99 / 1.57 | 0.07 / 0.10 / 0.16 | 0.01 / 0.01 / 0.02 | 0.00 / 0.00 / 0.00 |
  | 0.01 | 0.24 | 13.69 / 19.36 / 30.61 | 1.37 / 1.94 / 3.06 | 0.14 / 0.19 / 0.31 | 0.01 / 0.02 / 0.03 | 0.00 / 0.00 / 0.00 |
  | 0.05 | 0.02 | 7.28 / 10.30 / 16.29 | 0.73 / 1.03 / 1.63 | 0.07 / 0.10 / 0.16 | 0.01 / 0.01 / 0.02 | 0.00 / 0.00 / 0.00 |
  | 0.05 | 0.12 | 1.52 / 2.16 / 3.41 | 0.15 / 0.22 / 0.34 | 0.02 / 0.02 / 0.03 | 0.00 / 0.00 / 0.00 | 0.00 / 0.00 / 0.00 |
  | 0.05 | 0.24 | 1.37 / 1.94 / 3.07 | 0.14 / 0.19 / 0.31 | 0.01 / 0.02 / 0.03 | 0.00 / 0.00 / 0.00 | 0.00 / 0.00 / 0.00 |
  (each cell: rho-hat_rate's bias in SE units at T = 20 / 40 / 100)

- **Zero-default periods (S-28).** A period has no defaults with probability
  P₀ = E_z[(1 − p(z))ⁿ], so the parity rate estimator refuses a T-period panel with probability
  exactly 1 − (1 − P₀)^T. That is most of the matrix at n ≤ 1,000:

  | PD | ρ | n | P(d = 0) | refused, T = 20 | T = 40 | T = 100 |
  |---|---|---|---|---|---|---|
  | 0.001 | 0.02 | 100 | 0.9059 | 1.0000 | 1.0000 | 1.0000 |
  | 0.001 | 0.02 | 1,000 | 0.4072 | 1.0000 | 1.0000 | 1.0000 |
  | 0.001 | 0.02 | 10,000 | 0.0038 | 0.0729 | 0.1404 | 0.3149 |
  | 0.001 | 0.12 | 100 | 0.9135 | 1.0000 | 1.0000 | 1.0000 |
  | 0.001 | 0.12 | 1,000 | 0.5577 | 1.0000 | 1.0000 | 1.0000 |
  | 0.001 | 0.12 | 10,000 | 0.1171 | 0.9171 | 0.9931 | 1.0000 |
  | 0.001 | 0.24 | 100 | 0.9260 | 1.0000 | 1.0000 | 1.0000 |
  | 0.001 | 0.24 | 1,000 | 0.6828 | 1.0000 | 1.0000 | 1.0000 |
  | 0.001 | 0.24 | 10,000 | 0.3186 | 0.9995 | 1.0000 | 1.0000 |
  | 0.01 | 0.02 | 100 | 0.3914 | 1.0000 | 1.0000 | 1.0000 |
  | 0.01 | 0.02 | 1,000 | 0.0016 | 0.0309 | 0.0608 | 0.1452 |
  | 0.01 | 0.02 | 10,000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
  | 0.01 | 0.12 | 100 | 0.4984 | 1.0000 | 1.0000 | 1.0000 |
  | 0.01 | 0.12 | 1,000 | 0.0596 | 0.7075 | 0.9144 | 0.9979 |
  | 0.01 | 0.12 | 10,000 | 0.0010 | 0.0202 | 0.0401 | 0.0972 |
  | 0.01 | 0.24 | 100 | 0.5994 | 1.0000 | 1.0000 | 1.0000 |
  | 0.01 | 0.24 | 1,000 | 0.1936 | 0.9865 | 0.9998 | 1.0000 |
  | 0.01 | 0.24 | 10,000 | 0.0295 | 0.4510 | 0.6986 | 0.9501 |
  | 0.05 | 0.02 | 100 | 0.0143 | 0.2505 | 0.4382 | 0.7635 |
  | 0.05 | 0.02 | 1,000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
  | 0.05 | 0.02 | 10,000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
  | 0.05 | 0.12 | 100 | 0.0832 | 0.8241 | 0.9691 | 0.9998 |
  | 0.05 | 0.12 | 1,000 | 0.0007 | 0.0140 | 0.0279 | 0.0682 |
  | 0.05 | 0.12 | 10,000 | 0.0000 | 0.0000 | 0.0000 | 0.0001 |
  | 0.05 | 0.24 | 100 | 0.1898 | 0.9852 | 0.9998 | 1.0000 |
  | 0.05 | 0.24 | 1,000 | 0.0178 | 0.3017 | 0.5124 | 0.8339 |
  | 0.05 | 0.24 | 10,000 | 0.0007 | 0.0148 | 0.0295 | 0.0721 |

  **Dropping** zero periods keeps only the worse periods: PD̂ is biased upwards, and the kept
  periods' rates are dominated by binomial noise at small counts (d = 1 or 2), which the rate model
  reads as factor variance. **Censoring** at half a default keeps the zero periods' information
  about the mean, so its PD̂ should be much less biased; but it cannot remove the binomial noise in
  the non-zero periods, so where the S-27 bias is large neither treatment rescues the rate MLE.

## Predictions

Groups as in D-136 for n ≤ 10⁴: A (31), B (12), C (17), D (21). The 54 scenarios at n ≥ 10⁵ are
grouped by T and are not in A–D.

### S-8: MoM against the binomial MLE (n ≤ 10⁴)

| # | Prediction |
|---|---|
| E1 | **ρ, efficiency:** RMSE(MoM)/RMSE(MLE) for ρ exceeds 1 in at least 45 of the 50 B–D scenarios |
| E2 | **Pattern:** that ratio increases from ρ = 0.02 to 0.12 to 0.24 in at least 80% of the (PD, T, n) cells whose three ρ scenarios are all in B–D |
| E3 | **Size at n = 10⁴:** in the B–D scenarios with n = 10⁴, the ratio is within ±25% of the large-n reference (table above) in at least 75% |
| E4 | **ρ, bias:** MoM's ρ̂ has a more negative mean error than the MLE's in at least 40 of the 50 B–D scenarios |
| E5 | **PD:** RMSE(MoM)/RMSE(MLE) for PD is 0.95–1.15 in at least 45 of the 50 B–D scenarios |

### S-27: the rate MLE (censored) against the binomial MLE

| # | Prediction |
|---|---|
| E6 | **n = 10⁶:** the two are indistinguishable in at least 24 of the 27 scenarios |
| E7 | **Direction:** the rate MLE's mean ρ̂ exceeds the binomial MLE's in every scenario whose reference bias is at least 0.5 SE and where both estimate in at least half the replicates |
| E8 | **The threshold:** over the (PD, ρ, T, n) cells with n ≥ 10³ whose reference is unambiguous (below 0.05 or above 0.3 SE), "indistinguishable" agrees with "reference below 0.1" in at least 80% |

### S-28: zero-default treatments

| # | Prediction |
|---|---|
| E9 | **Refusal is exactly predictable:** the parity rate estimator's refusal share lies within ±3.29·√(q(1 − q)/1,000) of q = 1 − (1 − P₀)^T (table above) in at least 79 of the 81 scenarios with n ≤ 10⁴ |
| E10 | **Drop biases PD upwards:** where P₀ ≥ 0.2, the drop treatment's mean relative error of PD̂ exceeds +10% in every scenario where it estimates in at least half the replicates |
| E11 | **Censoring helps PD:** where P₀ ≥ 0.05, the censored treatment's absolute mean relative error of PD̂ is smaller than drop's in at least 90% of the scenarios where both estimate in at least half the replicates |
| E12 | **Neither rescues ρ where the rate model is wrong:** where the S-27 reference bias is at least 1 SE, the ρ profile coverage of both drop and censored is below the band in at least 80% of the scenarios where they estimate in at least half the replicates |

## How the comparison will be reported

Each prediction gets a row: prediction, result, held or not. A prediction that holds in direction
but misses in size is recorded as not held, with the size of the miss. The rows are computed by a
script committed before the full run's results exist. Results, the comparison, a mitigation and a
monitoring implication go in `studies/README.md` (the S-8, S-27 and S-28 entries) and in a D-entry,
with a `MANIFEST.json`. The pinned binomial verdicts are unchanged.
