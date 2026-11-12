# The rate MLE as a first-class estimator: the plan (D-177 item 3, D-416)

The Vasicek-rate MLE (`vasicek_rate_mle.md`, M3) is an exposed estimator (`vcal_calibrate_rate`, ABI 0.3,
D-169). It still runs on the engine's generic fallbacks: central-difference SEs (D-119, kept for it by
D-185), and nested Brent for the profile intervals and the polished maximum (which places a flat maximum
only to about 10⁻⁷ in logit units, D-164). D-177 ranked making it first-class third after M3. This note
plans that work in three PRs and records what each must show. It is rewritten as measured as each PR lands.

## PR A: analytic derivatives (done, D-417)

- **What:** `VasicekRate::log_contrib_derivs`, the score and Hessian of each period's contribution in
  (PD, ρ), for the density and for both native boundary treatments (censored and substituted). It returns
  `log_contrib`'s value, bit for bit, from the same call (D-170's contract).
- **What the engine now does, with no other change** (each path is chosen by `has_contrib_derivs`):
  - the SEs are the analytic observed information (D-185), not D-119's differences;
  - the profile intervals and the polished maximum are solved by safeguarded Newton (D-170), not Brent;
  - the HAC SE (D-209) becomes available to the rate MLE.
- **Checks, with tolerances stated before they ran:**

  | Check | Tolerance | Measured |
  |---|---|---|
  | Score vs Richardson differences, density | `TOL_SCORE_HESSIAN_FD_REL` (1e-8) | 3.2e-10 |
  | Hessian vs Richardson differences, density | 1e-8 | 9.4e-10 |
  | Score / Hessian, censored boundary terms | 1e-8 | 5.4e-10 / 3.0e-9 |
  | Score / Hessian, substituted boundary terms | 1e-8 | 2.9e-10 / 9.6e-10 |
  | Analytic SE vs the closed-form SE, 35 panels | 1e-10 relative | 5.6e-15 |
  | Polished maximum vs the closed form, unit panels | `TOL_VASICEK_RATE_POLISH_U` | 7.1e-10 (4.8e-8 under Brent) |
  | Polished maximum vs the closed form, 27,000 recovery fits | as above | 2.6e-9 (1.6e-7 under Brent) |

  `TOL_VASICEK_RATE_POLISH_U` was set from the recovery measurement to 1.2e-8 (worst × 4), from 4e-7.
- **Results change: within solver tolerance, and no verdict flips.**
  - The rate recovery summary was refreshed. Every count and class is unchanged, so 0 of the 81 verdicts
    flip. The doubles moved by at most 1.6e-7 (scaled), and the three closed-form gap maxima above
    tolerance shrank from about 1.6e-7 to 2.6e-9.
  - The replay panels for the scipy script were refreshed. scipy now agrees to 1.1e-7 in the estimates
    (scipy's own Brent is the looser side), 3.4e-13 in the log-likelihood and 2.5e-10 in the interval ends.
  - Pinned study results that used the rate MLE (S-8, S-27, S-28, S-35's rate arms) stay pinned to their
    runs (D-406), and are not rerun.

## PR B: two diagnostics from S-27

- **The adequacy flag.** S-27 found that for count data the rate MLE matches the binomial MLE only once
  n·PD reaches about 500–1,000, and that below it ρ̂ is biased up (D-174). So a fit of count data gets a
  flag when its n·PD is below a threshold.
  - **The threshold is the owner's decision.** It is proposed as 1,000, S-27's mitigation line. It is
    applied with n the panel's median n_t and PD the fit's PD̂.
  - It is a `native` diagnostic. Parity output is unchanged.
- **The rate-vs-count disagreement diagnostic.** For count data both fits are possible. The diagnostic
  reports the gap between the two ρ̂ in units of the binomial fit's SE, and S-27's predicted gap at that
  (n, PD, ρ, T) from `model_reference.py s27`. A gap well beyond the prediction points at something other
  than binomial noise.
- **Checks:**
  - unit tests at both sides of the threshold;
  - the predicted gap against S-27's pinned reference rows;
  - a scipy replication of the predicted gap.

## PR C: interval coverage and a data-sufficiency table

- **Coverage:** the rate MLE's Wald, profile and bootstrap intervals on its own model. The profile is
  measured already (section 5 of `vasicek_rate_mle.md`); after PR A it is re-measured with analytic SEs and
  bootstrap intervals (D-042's W × L, not yet validated for this objective).
- **The data-sufficiency table:** like the binomial MLE's sample-size table, the smallest T at which the
  rate MLE's intervals for PD and ρ are in the band and their widths are below stated targets, by
  (PD, ρ).
- **Pre-registered:** this PR runs recovery studies, so its predictions are committed and merged before
  the run, as for S-53.

## Order and what waits on whom

1. PR A, now.
2. PR B, after the owner sets the threshold.
3. PR C, registered after PR A merges, since its predictions use PR A's analytic SEs.

Nothing here touches the binomial MLE or the parity profile's other estimators.
