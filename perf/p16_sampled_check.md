# P-16: a sampled quadrature check on bootstrap fits (D-180)

No effect on results other than the check's own fields (bit for bit otherwise). Measured with
`study_parametric_bootstrap --replicates 100 --scenarios 51,29` before and after P-16 (both with P-14),
24 threads, MSVC 19.51 Release; CPU seconds summed over fits:

| scenario | T | n | calibrate after the argmax: before | after | change |
|---|---|---|---|---|---|
| 51 | 100 | 100 | 176.3 | 156.0 | -11% |
| 29 | 20 | 10,000 | 251.4 | 238.7 | -5% |

Rows byte-identical on every shared column. The saving is small because the check integrates each distinct
observation once, against the D-119 Hessian's eight: P-15 (analytic Hessian SEs) is the lever for this phase.
P-16's purpose is the evidence: the flagged counts on 5% of bootstrap fits and every original fit.
