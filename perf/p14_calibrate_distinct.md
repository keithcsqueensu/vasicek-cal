# P-14: calibrate's quadrature check and Hessian per distinct observation (D-178)

Bitwise identical. Measured with `study_parametric_bootstrap` (the S-10/S-4 tool) on two subset
scenarios, 1,000 replicates each, B = 999, 24 threads, MSVC 19.51 Release, before and after P-14.
CPU seconds summed over fits, by phase:

| scenario | T | n | calibrate after the argmax: before | after | ratio | rows |
|---|---|---|---|---|---|---|
| 51 | 100 | 100 | 30,062 | 2,780 | 10.8× | 1,000 byte-identical |
| 29 | 20 | 10,000 | 3,074 | 3,100 | 1.0× | 1,000 byte-identical |

Before P-14, the subset (9 scenarios) spent 81% of 124,109 CPU s in this phase. The gain follows the
number of distinct observations per panel: largest where T is long and n·PD small (few distinct
default counts), none where every period's count differs (n = 10⁴).

Reproduce: build `study_parametric_bootstrap` at the commits before and after P-14 and run
`study_parametric_bootstrap --scenarios 51,29 --out rows.csv`; the per-scenario phase table is printed.
