# P-24, P-25 and P-26: three bit-for-bit savings in the S-10 fit (D-407, D-408, D-409)

From the October scouting note's candidates C1, C3 (tier A) and C2's first step
([scouting_2026-10.md](scouting_2026-10.md); D-500). None changes a result: every row of every tool is to be
byte-identical to `main`'s, and no golden is refreshed.

## Expected, written before the build (2026-11-06)

Shares are of S-10's projected thread-seconds after P-19 (the scouting note's re-profile, with P-19's analytic SE
phase taken out: W* about 36%, argmax + refine + SE about 27%, simulate about 28%, surface about 5%).

### P-24: calibrate's log-likelihood from its Hessian pass (D-407; C1)

- **The repeat:** `calibrate_from_surface` integrates every distinct observation at the refined estimate twice. The
  value pass (`log_contrib`) gives `l_primary`, which feeds the log-likelihood and the quadrature check; then the
  Hessian pass (`panel_derivs`, P-15) integrates the same observations at the same θ again. The two θ are the same
  doubles: `out.value` is `from_scaled(r.scaled)`, which is what `panel_derivs` evaluates at. And
  `log_contrib_derivs` returns `log_contrib`'s value bit for bit (the objective's contract, D-170).
- **The change:** where calibrate takes the analytic Hessian (no grid-edge or flat flag, an objective with
  derivatives), the Hessian pass runs first and hands back each distinct observation's value; `l_primary` is filled
  from those in period order, as now. The doubled rule's `l_check` pass, flagged fits and the finite-difference path
  are unchanged.
- **Expected:** the value pass is about 40% of the argmax + refine + SE phase at T = 100, n = 10,000 and about 26% at
  T = 100, n = 100 (scouting note, callgrind). **The fit phase −25% to −40%; S-10 about 8-11%.**
- **Stop rule:** under 4% of per-replicate work on the S-10 tool: recorded, not merged.

### P-25: W* without repeated evaluation points (D-408; C3 tier A)

- **The repeat:** in `NewtonProfile`, about one evaluation per W* panel (7-8% of them) is an exact repeat of a point
  already evaluated: `start_axis1` re-evaluates the polished maximum, whose derivatives an inner solve computed.
- **The change:** `NewtonProfile` keeps the full derivatives of its last few evaluated points and returns the memo
  when asked for exactly the same pair (compared bit for bit). `evaluations` still counts every point requested, so
  every output, the C ABI's counts included, is unchanged. `profile_intervals_newton` uses the same class and gains
  the same.
- **Tier A′ (seeding the memo with calibrate's own evaluation) is not built.** The memo would then hold a value from
  outside the `NewtonProfile`, which is right only if the caller profiles with the integrator it calibrated with.
  The C ABI calls calibrate and the profile separately and cannot enforce that, and a wrong memo would be a silent
  wrong result. It would also hit only where `to_scaled(value)` round-trips exactly (about 75% of panels, D-199), for
  about 1.5% of S-10.
- **Expected:** the W* phase −6% to −8%; **S-10 about 2-3%.**
- **Stop rule:** under 3% of the W* phase: recorded, not merged.

### P-26: no zero fill before the surface copy (D-409; C2, first step)

- **The waste:** before every cached fit, `L.assign(T × K, 0.0)` zero-fills the surface, and
  `evaluate_surface_cached` then overwrites every element (each period's row is a full `memcpy`). At T = 100 that is
  2 MB written twice per bootstrap fit, 999 times per replicate.
- **The change:** `resize` instead of `assign`, in the cached paths (`calibrate_cached`, the S-10 tool's
  `fit_panel`). A reused vector of the same size is not touched; a new one is still value-initialised.
- **Expected:** the zero fill is about half of the surface phase's instructions (3.1% of 6.3% at T = 100,
  n = 10,000; 17% of 35% at T = 100, n = 100), memory-bound, so time is less certain. **The surface phase −30% to
  −50%; S-10 about 1.5-2.5%.**
- **Stop rule:** under 15% of the surface phase: recorded, not merged. The rest of C2 (a row view in place of the
  copy, and calibrate's argmax bounded by cached tile maxima) waits on what is left of the surface phase after this.

### The checks, for all three

- **Identity:** each item's own unit test against the old path, bit for bit; the non-slow suite with the committed
  goldens; the S-10 tool's rows on a sample against `main`'s, byte for byte.
- **Timing:** the S-10 tool on scenarios 42 (T = 100, n = 100) and 44 (T = 100, n = 10,000), one replicate, `main`
  and this branch run A B B A on one machine state, per phase (`--timing-out`). Each item acts on its own phase (fit,
  W*, surface), so one pair of builds attributes all three.

## What was built

- **P-24:** `panel_derivs` and `panel_derivs_natural` (`engine/derivs.hpp`) take an optional `values` array and write
  each distinct observation's `log_contrib_derivs` value to it. `calibrate_from_surface` (`engine/calibrate.hpp`)
  takes the Hessian pass first where it takes the analytic Hessian, and fills `l_primary` from those values in period
  order; the value pass runs only for flagged fits and objectives without derivatives.
- **P-25:** `NewtonProfile::derivs` (`engine/profile.hpp`) keeps the last 8 evaluated points and their `PanelDerivs`,
  and answers a pair equal bit for bit from it. `memo_hits` counts the answers; `evaluations` is unchanged in
  meaning. `memo_enabled()` turns it off, for the identity test only.
- **P-26:** `calibrate_cached` and the S-10 tool's `fit_panel` `resize` the surface instead of `assign`ing zeros.
  `evaluate_surface_cached`'s comment now states that it writes every element.

## Identity

- **Unit tests** (all bit for bit):
  - `calibrate_takes_its_value_pass_from_the_hessian_pass_bit_for_bit` (`unit_profile`): on the profile test's
    panels, an n = 100 panel with one-sided periods and a grid-edge panel, with the quadrature check on and off,
    calibrate's log-likelihood, the check's max, total and flagged count, and its flag equal the value pass computed
    as before P-24. 8 fits took the Hessian pass, 2 the value pass; 0 differ.
  - `newton_profile_memo_changes_no_number` (`unit_profile`): with the memo on and off, every field of
    `profile_intervals_newton` and of `profile_lr_statistic` (at a truth-like point and at two profile ends) is the
    same, the evaluation counts included: 16 calls, 0 differ. W*'s sequence asked for 61 points, of which 4 (6.6%)
    came from the memo.
  - `cached_surface_overwrites_every_element` (`unit_surface_cache`): over a NaN-filled buffer and over the previous
    panel's surface, `evaluate_surface_cached` gives the zeroed fill's result.
- **The non-slow suite** passes with the committed goldens (81 of 81, GCC 13, warnings as errors).
- **The recovery harness:** `recovery_harness --replicates 50 --save-fits` on `main` and this branch (81 scenarios
  × 50, so calibrate, `profile_intervals_newton` with the memo, the bootstrap and q's interval on every scenario):
  every field of all 4,050 fits identical. The two files are not byte-identical, because `--save-fits` writes
  `rc::Fit` raw and its two 4-byte padding slots are uninitialised; the comparison skips them.
- **The S-10 tool's rows**, all four timing runs below (scenarios 42 and 44, 10 replicates each): byte-identical to
  `main`'s.

## Measured (2026-11-06)

Cloud Xeon, 2 vCPU (AVX2), GCC 13, Release; `study_parametric_bootstrap --replicates 10 --scenarios 42,44
--threads 2`, `main` (f54db25) and this branch run A B B A; thread-seconds per phase summed over the 10 replicates.

| Phase | (100, 100): `main` | branch | change | (100, 10,000): `main` | branch | change |
|---|---|---|---|---|---|---|
| surface | 1.66, 1.60 | 0.89, 0.91 | −44.6% | 3.42, 3.79 | 2.31, 2.47 | −33.7% |
| argmax + refine + SE | 4.45, 4.32 | 3.40, 3.60 | −20.3% | 15.47, 15.77 | 10.54, 10.43 | −32.8% |
| W* (S-4b) | 3.79, 3.66 | 3.46, 3.39 | −8.0% | 16.35, 16.26 | 15.16, 15.14 | −7.1% |
| simulate | 1.16, 1.15 | 1.15, 1.15 | −0.6% | 20.85, 21.16 | 21.25, 21.14 | +0.9% |
| Jeffreys (unchanged code) | 5.95, 5.98 | 5.98, 5.85 | −0.8% | 39.56, 40.13 | 40.28, 39.81 | +0.5% |

- **Per-replicate work, with the full run's schedule** (Jeffreys on 5% of replicates, the S-4 widths on 20%): 11.4 →
  9.4 thread-s at (100, 100), **−17.0%**; 59.2 → 51.9 at (100, 10,000), **−12.4%**. The (100, 10,000) cell is 36%
  of the run.
- **Against the expectations, by phase:**
  - **P-24, the fit phase:** expected −25% to −40%; measured −20% at n = 100 and −33% at n = 10,000. Held at
    n = 10,000, a little under at n = 100, where the value pass is a smaller share of the phase (the scouting note's
    26%) and the argmax's K × T reduction a larger one.
  - **P-25, W*:** expected −6% to −8%; measured −7% to −8%. Held.
  - **P-26, the surface phase:** expected −30% to −50%; measured −34% to −45%. Held.
- **Against the stop rules:** none triggered (P-24's 4% of per-replicate work, P-25's 3% of W*, P-26's 15% of the
  surface phase).
- **Run level:** on the scouting note's shares after P-19, about 12-14% of the S-10/S-4 run. The owner's laptop
  measures it on the next full run. **Measured (D-414):** P-24 to P-27 together save about 16% of the laptop's
  full run, normalised for machine state ([p24_p27_full_run.md](p24_p27_full_run.md)).
- **Still left in the surface phase:** after P-26 it is about 10% of per-replicate work at (100, 100) and 4% at
  (100, 10,000). The rest of the scouting note's C2 (a row view in place of the copy, and calibrate's argmax bounded
  by cached tile maxima) would act on that and on the K × T reduction in the fit phase. It is not proposed until
  the laptop's full-run phase split shows how much is left.
