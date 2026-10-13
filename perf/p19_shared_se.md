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
