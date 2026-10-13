# P-19: the analytic SE shared with calibrate (D-406)

## Expected, written before the build (2026-10-13)

- **What is computed twice (D-199):** calibrate takes the analytic observed information (P-15, D-185) at its
  refined estimate's scaled coordinates `r.scaled`. S-10's tool then calls `analytic_se_scaled` on the same panel
  at `to_scaled(value)`, the estimate after a round trip through the natural scale. It does this for the panel
  and for each of its B = 999 bootstrap refits, to studentise them. The round trip is inexact in about 13% of PD
  and 25% of ρ coordinates, so the two Hessians agree bit for bit in only about 75% of panels. The second one is
  the "analytic SE" phase: 13.5% of the S-10 run in D-199's profile, and 8-11% of per-replicate time in today's
  timing sample (scenarios 37 and 4).
- **The change:** calibrate keeps its Hessian in scaled units, and its scaled coordinates, in `Estimate2`. The
  studentised SE is taken from that Hessian with `analytic_se_scaled`'s own rule (−H positive definite, det
  finite). Where calibrate did not take its Hessian (a grid-edge estimate; P-15 skips those), the SE is computed
  as before, but at `r.scaled` instead of `to_scaled(value)`. The SE is then one computation, at one point.
- **Expected:** the "analytic SE" phase falls to its grid-edge share, so **S-10's tool 8-13% faster.**
  `study_parametric_bootstrap` is the only caller.
- **Changes results: within rounding.** About 25% of the studentised SEs move by an ulp or so, and the studentised
  interval's ends with them. By the P-item rules (D-177) that needs a golden refresh with a flip check of every
  fit, its counts reported even when zero, and its own commit. The full S-10 run (3.6 h on 24 threads) is the
  owner's laptop's. Here: a flip check on a timing-sized sample, and the identity tests.
- **Stop rule:** under 5% on S-10's tool: recorded as measured, not merged, since it would change results for less.

## What was built

- **`Estimate2::scaled` and `Estimate2::hess_u`** (`engine/calibrate.hpp`): the refined estimate's scaled
  coordinates, and the analytic Hessian calibrate takes there (P-15), in scaled units. Both are NaN where
  calibrate takes none, which is the grid-edge and flat estimates and the finite-difference path.
- **`engine/profile.hpp`:** `analytic_se_from_hessian`, the one validity rule (−H positive definite, det finite);
  `analytic_se_at_scaled`, the SE computed at scaled coordinates; `analytic_se_of(est)`, which uses calibrate's
  Hessian where it exists and otherwise computes the SE at `est.scaled`. `analytic_se_scaled` (natural values,
  through `to_scaled`) is kept for callers without an estimate.
- **`study_parametric_bootstrap`:** the panel's SE and every bootstrap refit's SE come from `analytic_se_of`.

## Identity and the change in results

- **`analytic_se_of_is_calibrates_hessian_bit_for_bit`** (`unit_profile`): on the profile test's panels and a
  zero-default panel at the grid edge, the SE from calibrate's Hessian equals `analytic_se_at_scaled` at the
  scaled estimate, bit for bit, in both branches. The P-15 test against the old, round-tripped computation still
  holds at `TOL_CALIBRATE_SE_VS_ANALYTIC_REL`. The non-slow suite passes (79 of 79, GCC 13, `-Werror`).
- **A sample flip check** (subset scenarios 37, 4 and 29, 8 replicates, 24 rows; `main` against this branch, row
  by row): the only columns that change are the analytic SEs (8 of 24 for PD, 5 for ρ) and the studentised ends
  (8-9 of 24). The largest end-point shift is 1.2e-14 in logit units. **Coverage flips: 0.** Every other column is
  byte-identical.
- **Still owed before merge:** the full flip check, every fit of the S-10/S-4 run against the run that wrote
  `studies/parametric-bootstrap/summary.csv`, with the summary refreshed from it (D-177). That run is the owner's
  laptop's.

## Measured (2026-10-13)

Cloud Xeon 2.8 GHz, GCC 13, one thread, Release; `main` with P-6.

- **The "analytic SE" phase:** 2.9 → 0.2 thread-s over the sample (−93%). The rest is grid-edge refits.
- **Per-replicate work** (every timed phase but the once-per-scenario cache fill): 40.4 → 37.6 thread-s (−6.9%).
  That is the run-level saving once fixed costs are spread over R = 1,000.
- **Wall time of the sample:** 500, 500 s against 498, 497 s. Each run has about 450 s of fixed set-up before the
  first replicate (outside the timed phases; probably the Jeffreys tables for n = 10,000, not checked), so the
  sample's wall time cannot show a per-replicate saving.
- **Against the expectation:** missed. 8-13% was expected; the measured saving is 6.9%. The phase was 7.2% of
  per-replicate work in this sample, against D-199's 13.5% of the whole run. D-199's profile spans every cell; this
  sample has three.
- **Against the stop rule (under 5%):** not triggered.
