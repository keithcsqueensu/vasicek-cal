# Parity sign-off: binomial-mixture MLE (M1 exit)

What D-123 requires, and where each item is. Tag: `v0.1.0-parity-binomial`, placed on a
commit whose CI is green.

| Requirement | Where |
|---|---|
| Methodology in plain prose | [binomial_mixture_mle.md](binomial_mixture_mle.md) |
| scipy replication script | `validation/scipy/binomial_mixture_mle.py`, run in CI by the `validation (scipy)` job |
| Tolerance register | [tolerances.md](tolerances.md), with values in `tests/tolerances.toml` |
| Recovery tables as goldens with provenance | `tests/golden/recovery/` (summary, replay, MANIFEST.json); [recovery_results.md](recovery_results.md) |
| Deferred and known-finding coverage verdicts, each with its M2 dependency | this page, below |

**Engine correctness holds** (D-124):

- RMSE falls with T in every setting;
- none of the 81,000 recovery fits has a quadrature or numeric flag;
- the replay agrees across all CI platforms;
- estimates and SEs agree with the independent reference (M1.7);
- the scipy replication passes all three checks.

**Interval coverage** is a property of the Wald interval method (D-124). Of 162 scenario ×
parameter verdicts, 81 PASS. The other 81 are listed below, and all of them are re-assessed when
profile-likelihood intervals exist (M2).

The list below is derived from `tests/golden/recovery/summary.csv`. The engine's near-bound flag
is joint, so a scenario is deferred for both parameters together.

## Known findings (19): re-assess with M2's profile-likelihood intervals and the t(T − 1) option

| Scenario | PD | ρ | T | n | Parameter | Flagged | Coverage | Diagnosis |
|---|---|---|---|---|---|---|---|---|
| 11 | 0.001 | 0.12 | 20 | 10000 | PD | 0.0% | 0.926 | small-T undercoverage (SE ratio 1.0069) |
| 13 | 0.001 | 0.12 | 40 | 1000 | ρ | 3.4% | 0.998 | overcoverage, skewed estimate (SE ratio 0.8919) |
| 29 | 0.01 | 0.02 | 20 | 10000 | PD | 0.0% | 0.910 | small-T undercoverage (SE ratio 1.0747) |
| 29 | 0.01 | 0.02 | 20 | 10000 | ρ | 0.0% | 0.911 | small-T undercoverage (SE ratio 1.0992) |
| 31 | 0.01 | 0.02 | 40 | 1000 | ρ | 2.4% | 0.984 | overcoverage, skewed estimate (SE ratio 0.9436) |
| 37 | 0.01 | 0.12 | 20 | 1000 | PD | 0.2% | 0.917 | small-T undercoverage (SE ratio 1.0289) |
| 38 | 0.01 | 0.12 | 20 | 10000 | PD | 0.0% | 0.927 | small-T undercoverage (SE ratio 1.0332) |
| 42 | 0.01 | 0.12 | 100 | 100 | ρ | 0.3% | 0.981 | overcoverage, skewed estimate (SE ratio 0.9988) |
| 47 | 0.01 | 0.24 | 20 | 10000 | PD | 4.7% | 0.920 | small-T undercoverage (SE ratio 1.0553) |
| 49 | 0.01 | 0.24 | 40 | 1000 | PD | 2.8% | 0.924 | small-T undercoverage (SE ratio 1.0266) |
| 51 | 0.01 | 0.24 | 100 | 100 | ρ | 2.6% | 0.973 | overcoverage, skewed estimate (SE ratio 0.9763) |
| 55 | 0.05 | 0.02 | 20 | 1000 | PD | 0.9% | 0.921 | small-T undercoverage (SE ratio 1.0370) |
| 58 | 0.05 | 0.02 | 40 | 1000 | PD | 0.0% | 0.910 | small-T undercoverage (SE ratio 1.0593) |
| 63 | 0.05 | 0.12 | 20 | 100 | ρ | 0.8% | 0.984 | overcoverage, skewed estimate (SE ratio 0.9855) |
| 65 | 0.05 | 0.12 | 20 | 10000 | PD | 0.0% | 0.920 | small-T undercoverage (SE ratio 1.0202) |
| 68 | 0.05 | 0.12 | 40 | 10000 | PD | 0.0% | 0.918 | small-T undercoverage (SE ratio 1.0461) |
| 73 | 0.05 | 0.24 | 20 | 1000 | PD | 4.9% | 0.911 | small-T undercoverage (SE ratio 1.0300) |
| 74 | 0.05 | 0.24 | 20 | 10000 | PD | 2.7% | 0.923 | small-T undercoverage (SE ratio 1.0172) |
| 74 | 0.05 | 0.24 | 20 | 10000 | ρ | 2.7% | 0.924 | small-T undercoverage (SE ratio 1.0497) |

## Deferred (31 scenarios × both parameters = 62): at least 5% of fits flagged; re-assess with M2's profile-likelihood intervals

| Scenario | PD | ρ | T | n | Flagged | of which on the grid edge |
|---|---|---|---|---|---|---|
| 0 | 0.001 | 0.02 | 20 | 100 | 100.0% | 91.8% |
| 1 | 0.001 | 0.02 | 20 | 1000 | 70.1% | 45.3% |
| 2 | 0.001 | 0.02 | 20 | 10000 | 6.0% | 0.5% |
| 3 | 0.001 | 0.02 | 40 | 100 | 99.9% | 81.7% |
| 4 | 0.001 | 0.02 | 40 | 1000 | 51.5% | 26.9% |
| 6 | 0.001 | 0.02 | 100 | 100 | 98.4% | 61.4% |
| 7 | 0.001 | 0.02 | 100 | 1000 | 33.6% | 11.5% |
| 9 | 0.001 | 0.12 | 20 | 100 | 100.0% | 84.4% |
| 10 | 0.001 | 0.12 | 20 | 1000 | 32.1% | 6.8% |
| 12 | 0.001 | 0.12 | 40 | 100 | 100.0% | 67.5% |
| 15 | 0.001 | 0.12 | 100 | 100 | 91.9% | 33.5% |
| 18 | 0.001 | 0.24 | 20 | 100 | 100.0% | 85.0% |
| 19 | 0.001 | 0.24 | 20 | 1000 | 68.3% | 8.1% |
| 20 | 0.001 | 0.24 | 20 | 10000 | 25.5% | 0.8% |
| 21 | 0.001 | 0.24 | 40 | 100 | 99.9% | 62.9% |
| 22 | 0.001 | 0.24 | 40 | 1000 | 33.2% | 1.9% |
| 23 | 0.001 | 0.24 | 40 | 10000 | 5.6% | 0.0% |
| 24 | 0.001 | 0.24 | 100 | 100 | 94.2% | 22.6% |
| 27 | 0.01 | 0.02 | 20 | 100 | 76.6% | 48.9% |
| 28 | 0.01 | 0.02 | 20 | 1000 | 13.4% | 2.3% |
| 30 | 0.01 | 0.02 | 40 | 100 | 62.5% | 37.7% |
| 33 | 0.01 | 0.02 | 100 | 100 | 50.4% | 24.9% |
| 36 | 0.01 | 0.12 | 20 | 100 | 46.9% | 10.8% |
| 39 | 0.01 | 0.12 | 40 | 100 | 8.1% | 2.0% |
| 45 | 0.01 | 0.24 | 20 | 100 | 71.7% | 8.4% |
| 46 | 0.01 | 0.24 | 20 | 1000 | 18.5% | 0.0% |
| 48 | 0.01 | 0.24 | 40 | 100 | 30.0% | 1.0% |
| 54 | 0.05 | 0.02 | 20 | 100 | 51.3% | 24.5% |
| 57 | 0.05 | 0.02 | 40 | 100 | 31.4% | 10.0% |
| 60 | 0.05 | 0.02 | 100 | 100 | 8.5% | 1.4% |
| 72 | 0.05 | 0.24 | 20 | 100 | 20.6% | 0.5% |

## M2a re-assessment with profile-likelihood intervals (D-128–D-131)

The profile interval re-assesses every verdict above. Its coverage counts all replicates, so
nothing is deferred. Band: 0.95 ± 3.29·√(0.95·0.05/1000) = 0.927–0.973, a two-sided 0.1% band,
so about 0.16 of 162 verdicts would fall outside by chance.

| Wald \ profile | PASS | CONSERVATIVE | KNOWN FINDING | total |
|---|---|---|---|---|
| PASS | 80 | 0 | 1 | 81 |
| DEFERRED | 41 | 20 | 1 | 62 |
| KNOWN FINDING | 15 | 0 | 4 | 19 |

**Headline:**

- The profile interval fixes 15 of the 19 Wald findings and 41 of the 62 deferred verdicts.
- **CONSERVATIVE (20):** near-uninformative settings, coverage 0.977–0.999, where truncation at a
  bound widens the intervals. This is a safe-side property.
- **KNOWN FINDING (6):** coverage 0.920–0.927 at T = 20–40, within a point of the band's lower
  edge. This is the genuine small-sample limitation of the likelihood-ratio interval. For ρ it is
  driven by ρ̂'s downward bias.
- Every verdict, with its diagnosis, is in [recovery_results.md](recovery_results.md).

**Backlog, as `native` options:**

- a Bartlett-type threshold correction;
- a bias-corrected ρ̂.

Each will be evaluated against these same pinned verdicts.

## M2b: bootstrap percentile intervals (D-135–D-137)

- **Undercoverage:** bootstrap percentile intervals undercover for ρ at every T and for PD at
  T ≤ 40.
- **Two reasons:** boundary breakdown near a bound of the box, and no correction for the bias or
  skewness of ρ̂.
- **Pinning:** the 125 per-scenario verdicts are pinned in the goldens, not listed here. The
  group-level coverage and the comparison with predictions made before the run are in
  [recovery.md](recovery.md) and D-137.
- **Recommendation:** profile-likelihood intervals are the recommended method for inference.
  Bootstrap percentile intervals are provided for comparison and are not recommended for ρ.
