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
