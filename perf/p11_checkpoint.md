# P-11: checkpoint and resume (D-179)

No effect on results by construction: a unit on disk is the rows the run holds in memory, and the
output is built from the units. Checked by `study_parametric_bootstrap_resume` (stop after one scenario,
resume, byte-identical CSV; a different configuration refused) and once with a real kill:

| run | scenarios | units on disk after `taskkill /F` | resumed from the checkpoint | output |
|---|---|---|---|---|
| `study_parametric_bootstrap --replicates 3 --scenarios 72,51,4,29 --checkpoint DIR` | 4 | 2 (72, 51) | 2 of 4 | SHA-256 equal to the uninterrupted run's |

Overhead: one file write per scenario (1,000 rows, a few hundred kB); negligible against a scenario's
minutes of fitting. Cost saved: a killed or crashed long run loses at most the scenario in progress.
