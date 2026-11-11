# P-24 to P-27 on the full S-10/S-4 run, the owner's laptop (D-414)

P-24 to P-27 (D-407 to D-410) were measured on a cloud machine with GCC, on three to nine scenarios at 10
replicates: about 19-21% of a run together. This is the full run on the owner's laptop, as `STATE.md` asked:
its total, its phase split, and what the split says comes next. No change in results was expected, and none
was found.

## Method

- **Machine:** Intel Core Ultra 9 275HX, 24 threads, Windows 11, on AC power, the Balanced power plan, sleep
  off. MSVC 19.51.36260, the `cpu-release` preset. Every kernel was AVX2 (`dgp_kernel=avx2`).
- **Binaries:**
  - **A, the baseline:** `f54db25`, `main` just before P-24 to P-26 merged (#82). The tool's SHA-256 is
    `f437fc00…6461`.
  - **B:** `34e55bb`, `main` with P-24 to P-27, the audit fixes (D-218, D-219) and the multi-grade work, which
    S-10 does not call. The tool's SHA-256 is `ad7fb7ca…cf84`.
- **The runs:** `study_parametric_bootstrap --out ROWS --timing-out TIMING`, all 81 scenarios × 1,000
  replicates, A and then B, on 2026-11-11. That is the order D-406 used for P-19, and the same machine. A ran
  08:17-10:25 on an otherwise idle machine. B ran 10:25-12:07; from about 11:15 the machine was also in use
  for other work.
- **The projection,** recorded before either run (D-187, D-193): `python perf/project_runtime.py --run
  "build/cpu-release/tests/study_parametric_bootstrap.exe"` on B. It stopped at 16 replicates per scenario
  (a 193 s sample) with **1.62 h** (90% range 1.59-1.65 h; 1.51-1.73 h with the machine-state allowance).
- **The split:** `python perf/full_run_split.py A_timing.csv B_timing.csv`, and the same with `--log
  B_stderr.txt --cut 2981` for the scenarios B finished by 11:15. The timing files and rows are kept outside
  the repository; they are not needed to replicate a result, since P-items change no result.

## Results

**Rows:** B's 81,000 rows are byte-identical to A's (47,866,886 bytes). A's are also byte-identical to the
P-19 full run's of 2026-10-14 (`13072d8`): S-10's code did not change between the two.

**Raw totals:**

| | A (baseline) | B (P-24 to P-27) | Change |
|---|---|---|---|
| Wall | 7,701 s | 6,129 s | −20.4% |
| CPU (timed phases) | 183,930 s | 146,173 s | −20.5% |
| Efficiency, CPU / (24 × wall) | 0.995 | 0.994 | |

**The raw change overstates the gain.** The phases whose code P-24 to P-27 leave alone (simulate, profile,
W0, Jeffreys, the S-4 widths) also ran faster in B: 5% over the full run, and 8-18% by cell over the
scenarios B finished before 11:15. The machine itself was faster in B than in A. A was also slower than the identical
code's run on 2026-10-14, by 15-21% in every phase (17.6% in total). So the change is normalised: within each (T, n)
cell, B's cost is divided by B / A on those control phases (`perf/full_run_split.py`).

**Normalised, per phase** (A's thread-seconds and share; B's change):

| Phase | A | Share in A | B, full run | B, to 11:15 (48 scenarios) | Share in B |
|---|---|---|---|---|---|
| simulate (control) | 52,383 | 28.5% | −0.4% | +0.5% | 33.9% |
| surface | 8,939 | 4.9% | −66.0% | −73.8% | 2.0% |
| argmax + refine + SE | 48,628 | 26.4% | −44.1% | −46.3% | 17.7% |
| analytic SE | 896 | 0.5% | −2.2% | +5.9% | 0.6% |
| W* (S-4b) | 64,996 | 35.3% | −4.1% | −6.6% | 40.5% |
| Jeffreys (control) | 5,335 | 2.9% | +3.3% | −4.7% | 3.6% |
| profile, W0, S-4 widths, cache fill | 2,754 | 1.5% | −4% to +3% | −7% to +4% | 1.7% |
| **Total** | **183,930** | | **−16.3%** | **−16.6%** | |

- **The run-level saving is about 16%,** against 19-21% expected from the cloud machine. The two cuts agree
  to 0.3 points, although the machine factor moved by up to 18% between them. The factor's assumption (that
  machine state scales every phase of a cell alike) is not tested here. The control phases scatter by about
  ±4-7% around zero after normalisation, so read the per-phase changes to about that.
- **By item, against the cloud measurement:**
  - **The surface phase, −66% to −74%** (P-26 and P-27 together). The cloud measured −34% to −45% for P-26
    and −60% on what was left for P-27, about −75% together.
  - **Argmax + refine + SE, −44% to −46%** (P-24's log-likelihood from the Hessian pass, P-27's bounded argmax).
    The cloud measured P-24 alone at −20% to −33% of this phase, and P-27 removed about 60% of calibrate's argmax
    inside it; no combined expectation for the whole matrix was written down.
  - **W\*, −4% to −7%** (P-25's memo). Expected −7% to −8%; at the low end or a little under it.
- **The shortfall** against 19-21% is these phases' own readings, weighted by A's shares: the surface phase gives
  3.2-3.6 points, the fit phase 11.7-12.2 and W* 1.4-2.3, 16.3-18.2 in all. If the fit phase had fallen by half
  and W* by 7.5%, the same shares would give 19.5. That is within the ±4-7% the per-phase readings carry, so the
  laptop does not show a different effect from the cloud's, only a slightly smaller one. The n = 10,000 cells,
  where simulate is 36-50% of the cost, are 59% of the run.

**The projection: not a test.** B took 1.70 h against 1.62 h projected: inside the allowance (1.51-1.73 h),
outside the 90% range (1.59-1.65 h). About the last third of B shared the machine with other work, so this
run neither validates nor counts against the method (D-187).

## The phase split after P-27, and what it says next

B's normalised split: **W\* 40.5%, simulate 33.9%, argmax + refine + SE 17.7%,** Jeffreys 3.6%, surface 2.0%,
everything else 2.4%.

- **W\* is the largest phase.** P-20 (warm-started W* refits, the scouting note's C3 tier B, within solver
  tolerance, with its own flip check) targets it, so P-20's ceiling has risen to about 40% of the run.
- **The surface phase is down to 2%.** The P-27 follow-up (each panel's cached-row lookups resolved once, and
  the bootstrap's panels looked up in a flat table) can save at most 2% of the run. That is under P-27's stop
  rule of 3%.
- **Simulate is a third of the run** and on no candidate list. P-17 (SIMD Philox) is its last item, and it is
  half of the T = 100, n = 10,000 cell. Any further item there would be a new candidate, not on the owner's list.
- **Argmax + refine + SE** is 18%. What is left after P-24 and P-27 is the bounded search's own cost and the
  refinement; no candidate is open there.
