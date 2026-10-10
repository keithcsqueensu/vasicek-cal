# P-9: the repeat rate of identical panels, measured (D-189)

P-9 would memoise complete fits of identical panels (STATE.md, the P-items table). It was to be built only
if the repeat rate projected a real saving. It does not: memoising exactly repeated panels would save about
1% of the S-10/S-4 run, and nothing measurable in the recovery harness. **P-9 is not built.**

Measured with `perf/p9_repeat_rate.cpp` (target `p9_repeat_rate`), which regenerates the panels exactly as
the runs draw them and counts repeats; it fits nothing beyond what it needs to regenerate them. Machine:
Intel Core Ultra 9 275HX, MSVC 19.51 Release; about 22 minutes for the whole matrix.

## What can be shared, and under which key

- **The ordered key** (the default counts in period order) is the only key under which a memoised fit is
  bit for bit the fit it replaces. The log-likelihood sums the periods in order, so a fit of the same counts
  in another order differs in the last bits.
- **The sorted key,** the one STATE.md's P-9 entry named, would share fits between permutations of the same
  counts. That changes results within rounding. It is also wrong for any estimator that uses the order of
  the periods, such as an AR(1)-aware MLE (M6, S-35), so it could not carry over.
- **What a memo could share,** by scope:
  - an outer (recovery) panel's estimate, SEs, profile and q intervals, across replicates and scenarios of
    a (T, n) cell. Its iid-bootstrap arm is drawn from a stream seeded by (scenario, replicate), so it
    cannot be shared;
  - an S-10 bootstrap panel's estimate and SE, across the cell; S-4b's W* only within a replicate, since it
    is evaluated at that replicate's estimate.

## Repeat rates

Both runs' panels, 81 scenarios × 1,000 replicates. S-10's bootstrap panels are 999 per replicate,
simulated at that replicate's estimate (recomputed by the tool, and equal bit for bit to S-10's rows: checked
on 200 replicates). For the n = 10,000 cells the bootstrap panels were counted on replicates 0–99: there,
not one panel repeated under the ordered key, and 306 of 8.7 million under the sorted key.

Fraction of S-10 bootstrap panels that repeat an earlier panel of the same (T, n) cell:

| T | n | ordered | sorted | S-4b's first 199, within a replicate (ordered) |
|---|---|---|---|---|
| 20 | 100 | 0.319 | 0.686 | 0.155 |
| 20 | 1,000 | 0.019 | 0.305 | 0.001 |
| 20 | 10,000 | 0 | 0.000 | 0 |
| 40 | 100 | 0.209 | 0.606 | 0.055 |
| 40 | 1,000 | 0.000 | 0.256 | 0.000 |
| 40 | 10,000 | 0 | 0 | 0 |
| 100 | 100 | 0.017 | 0.475 | 0.002 |
| 100 | 1,000 | 0 | 0.110 | 0 |
| 100 | 10,000 | 0 | 0 | 0 |

The repeats are panels with almost no defaults: at PD 0.001 with n = 100 the expected count per period
is 0.1. In the (20, 100) cell, 11% of all bootstrap panels have no defaults at all.

The outer panels: 1,974 of 81,000 repeat within their scenario under the ordered key (10,561 sorted), all
in the n = 100 cells at PD 0.001.

## The projected saving

Each cell's repeat fraction times its per-replicate cost of the phases a memo would skip (the estimate and
SE, and W* for the within-replicate repeats), from a timing sample of the current code (2 replicates per
scenario; `perf/README.md`, "Projecting a full run"). The simulation itself is not saved, since a panel
must be drawn to be looked up.

| | thread-s | share of the run |
|---|---|---|
| The S-10/S-4 run, projected | 436,676 | |
| Saved, ordered key | 4,313 | 1.0% |
| Saved, sorted key (not bit for bit) | 27,570 | 6.3% |

The repeats fall in the cheapest cells: n = 100 is 10% of the run's cost, n = 10,000 is 64%, and the latter
never repeats. In the recovery harness (288 s for the whole matrix) the shareable part of 2.4% of the
outer fits is a few seconds.

## Where the time goes instead

In the same timing sample, the simulation is 31% of the run's projected thread-seconds, and 38–46% in the
n = 10,000 cells, against 17% of the S-10 run's total before P-14 and P-15. That bears on SIMD Philox (STATE.md,
next steps).

## The AR(1) DGP

Measured with S-35's own tool on its exact panels (`study_lrpd_serial --count-repeats`, D-190), the recovery
matrix's outer panels at each φ. Under an AR(1) factor the periods are no longer exchangeable, so only the
ordered key is valid. Persistence makes runs of quiet periods more alike, but the effect on the rate is
small.

Fraction of the 81,000 panels that repeat an earlier panel of the same (T, n) cell, by φ:

| φ | repeats | of which at n = 100, T ≤ 40 |
|---|---|---|
| 0 | 2,443 (3.0%) | 2,443 |
| 0.3 | 2,478 (3.1%) | 2,471 |
| 0.6 | 2,624 (3.2%) | 2,594 |

The φ = 0 count equals this tool's outer-panel count exactly (2,443 per cell; 1,974 per scenario), a check
that the two tools draw the same panels.

S-35 fits no parametric bootstrap. Its binomial fits are 97% of its projected run, and a memo could share
their estimate, SEs, profile and q intervals but not their iid-bootstrap arm. So the saving is at most about
0.6% of S-35's run, from that run's own timing sample. **This decision stands: P-9 is not built.**

Reproduce: `p9_repeat_rate --scenarios <the n <= 1,000 ids> --out small.csv`, then `p9_repeat_rate
--scenarios <the n = 10,000 ids> --boot-replicates 100 --out large.csv` (ids with id % 3 == 2 have
n = 10,000).
