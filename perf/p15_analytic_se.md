# P-15: calibrate's SEs from the analytic Hessian (D-185)

Changes results within a stated tolerance: the SEs, their correlation, q's delta-method SE and the
near-bound flag. Every estimate, log-likelihood, profile interval and bootstrap interval is bit for bit
as before. Machine: Intel Core Ultra 9 275HX (8 P-cores, 16 E-cores, 24 threads), MSVC 19.51 Release.

## Speed

`study_parametric_bootstrap --replicates 40 --scenarios 51,29,80,63`, the binaries before (A) and after
(B) P-15, run in the order A, B, B, A to separate the change from drift. Thread-busy seconds summed over
fits:

| run | wall s | total | calibrate after the argmax | other phases |
|---|---|---|---|---|
| A | 135 | 2,123 | 641 | 1,482 |
| B | 134 | 1,877 | 248 | 1,629 |
| B | 128 | 1,843 | 240 | 1,603 |
| A | 149 | 2,293 | 686 | 1,607 |

The calibrate phase falls 2.7× (mean 664 → 244). In this sample Jeffreys, which P-15 does not touch, is 22%
of the work, so the total falls 16%. In the S-10/S-4 full run the calibrate phase was 43% of 555,004 CPU s
and Jeffreys 2%, so at 2.7× the same run would save about a quarter of its CPU.

On the recovery harness (81 × 1,000), whose fits reuse cached surfaces and spend most of their time in the
bootstrap and the profile, the wall time was 288 s before and 285 s after: no measurable change.

The rows of the S-10 tool are byte-identical before and after on every column except `se_calib_u_pd` and
`se_calib_u_rho`. That includes the Jeffreys ends, which are seeded by calibrate's SE.

## Accuracy

D-119's central differences (step 0.15 × the stencil SE) carried a truncation error that was not small:

- **What was checked:** at chosen recovery fits, the observed information was recomputed independently,
  by Richardson-extrapolated central differences of `validation/scipy/binomial_mixture_mle.py`'s
  log-likelihood in logit coordinates.
- **The analytic SEs** agree with it to 1e-10 relative, and to 1.6e-7 in the worst near-bound case.
- **D-119's SEs** were off by up to 9e-3 in interior fits (ρ SE at ρ ≈ 0.017) and by up to 96% near a
  bound: fit 39/888's ρ SE was 23 times too small.

Over the 81,000 recovery fits, the relative change in the SE was:

| | above 2e-4 | above 1e-3 | above 1e-2 | above 0.1 | above 1 |
|---|---|---|---|---|---|
| PD SE, interior fits | 2,531 | 0 | 0 | 0 | 0 |
| PD SE, near-bound fits | 3,638 | 469 | 40 | 0 | 0 |
| ρ SE, interior fits | 8,326 | 1,776 | 0 | 0 | 0 |
| ρ SE, near-bound fits | 6,148 | 3,881 | 2,211 | 573 | 118 |

## Flip check (every recovery fit)

Fields that differ: `se_pd`, `se_rho`, `q_se_s` in 72,341 fits, and `flags` in 16 (near-bound: 15 on,
1 off). Everything else is bit for bit the same.

Coverage flips per replicate, gained / lost:

| Interval | Gained | Lost |
|---|---|---|
| Wald PD | 0 | 0 |
| Wald ρ | 1 | 0 |
| t PD | 0 | 1 |
| t ρ | 1 | 1 |
| q Wald | 32 | 0 |
| Profile (PD, ρ, q) | 0 | 0 |
| Bootstrap (PD, ρ, q) | 0 | 0 |

The flagged set gains 15 fits and loses 1. No verdict changes in any of the nine verdict columns.

Reproduce: build `recovery_harness` before and after P-15, run each with `--save-fits FILE`, and compare
the two files field by field (the `Fit` struct of `tests/recovery/recovery.hpp`).
