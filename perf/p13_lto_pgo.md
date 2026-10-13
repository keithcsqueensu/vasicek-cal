# P-13: link-time and profile-guided optimisation for local study builds (D-403): rejected

**Rejected under the stop rule (2026-10-13; D-403).** Measured on GCC (here) and MSVC (the owner's laptop): LTO
gains 0-2.4%, PGO 0-4% on GCC and nothing or worse on MSVC, all bit for bit. The options were built and are kept
on the branch `p13-lto-pgo`'s history (`cmake/VcalLtoPgo.cmake`, `perf/p13_pgo.py`), not merged.

## Expected, written before the build (2026-10-13)

- **Scope (D-192 addendum, D-199):** local study builds only. Options, off by default; CI's presets, the shared
  library and every golden are unchanged.
- **Where LTO can act:** almost nowhere. `core/`, `engine/`, `resample/` and `backend_cpu/` are header-only, so
  each study tool's hot code (the quadrature, the surface fill, the searches) is already one translation unit.
  The only calls across units go to the compiled kernels (`vcal_core_kernels`: erfcx; `vcal_dgp`: the DGP's
  draws; `vcal_benchmark_models`: the AR(1) transition) and the test harness. They are called once per batch or
  per panel, not per node. **Expected: 0-3%.**
- **Where PGO can act:** branch layout and inlining choices in code with data-dependent branches, such as the
  integrand's one-sided and interior paths, the hint's bracket and Newton, the composite rule's extent bisection,
  the argmax's tile visits and Newton's polish. **Expected: 3-10%** on the study tools, from a training run on
  the recovery subset.
- **Stop rule:** if LTO + PGO together give under 3% on both measured tools, P-13 is recorded as measured and
  the options are not added.
- **Bit for bit:** neither changes floating-point semantics. Contraction stays off (no `-mfma`; GCC and Clang
  in ISO mode; MSVC `/fp:precise`, which does not contract), and no reassociation is allowed. So every output must
  be byte-identical to the plain Release build's. The check compares each measured tool's rows byte for byte, and
  runs the non-slow test suite under the LTO + PGO build.
- **Measured on:** this cloud machine with GCC and Clang. The owner's laptop builds with MSVC, so its LTO
  (`/GL`, `/LTCG`) and PGO (`/GENPROFILE`, `/USEPROFILE`) are built here but measured and checked there.

## What was built

- **`cmake/VcalLtoPgo.cmake`**, included from the root `CMakeLists.txt`. It adds `VCAL_LTO` (CMake's IPO, checked
  with `check_ipo_supported`) and `VCAL_PGO` = `OFF`, `GENERATE` or `USE` with `VCAL_PGO_DIR`, for GCC
  (`-fprofile-generate`/`-fprofile-use`, atomic counter updates), Clang (the same, with `llvm-profdata merge`)
  and MSVC (`/GL`, `/LTCG`, `/GENPROFILE` and `/USEPROFILE` with one `.pgd` per binary). They are off by default,
  and refused with `VCAL_ENABLE_CUDA`.
- **`perf/p13_pgo.py`** runs the cycle in one build directory: GENERATE, build, train, (Clang) merge, USE,
  rebuild. With `--check` it also compares each tool's output against a plain Release build, byte for byte.
- **Training and check workloads,** split so the check is never trained on. Training: the bootstrap tool on subset
  scenarios 29 and 72, the M6 tool on scenario 49 at φ = 0.3. Check: the bootstrap tool on 37 and 4, the M6 tool on
  37 at φ = 0.6. One replicate each, one thread.

## Measured (2026-10-13, GCC 13 only)

Cloud Xeon 2.8 GHz, 2 cores, Release, one thread. The three builds were run one at a time, in the order
A B C C B A; times are in seconds.

| Tool | Release | LTO | LTO + PGO |
|---|---|---|---|
| `study_parametric_bootstrap` (S-10's tool) | 51.3, 52.7 | 52.8, 52.6 (+1.3%) | 51.5, 53.6 (+1.1%) |
| `study_lrpd_serial_m6` (S-35's M6 tool) | 162.1, 167.4 | 164.2, 163.9 (−0.5%) | 159.4, 158.1 (−3.6%) |

- **Identity:** every output of the LTO and LTO + PGO builds is byte-identical to the Release build's.
- **Tests:** the full non-slow suite passes under the LTO build (71 of 71, the committed goldens included).
- **LTO:** nothing, as expected. The hot code is header-only and already lives in one translation unit.
- **PGO:** nothing on the bootstrap tool, and about 3-4% on the M6 tool. The two Release runs of the M6 tool are
  3% apart, so that gain is at the edge of the noise.
- **Cost:** training with the instrumented build took 28.5 min for the bootstrap tool, mostly the atomic counters
  in single-threaded hot loops, and 71 s for the M6 tool. The profiles must be retrained whenever the code
  changes.
- **Against the expectation:** LTO at 0-3% held. PGO at 3-10% held only at its floor, and only on one tool.
- **Against the stop rule (under 3% on both tools):** borderline, not cleanly met, because the M6 tool shows
  −3.6% with a noise floor of about 3%.
- **Short of the plan written above:** only GCC was measured here, not Clang. The non-slow suite ran under the LTO
  build only (below), not under LTO + PGO, whose build covered just the two measured tools.
- **Not measured: MSVC,** which is what the owner's runs use. MSVC's inliner and code layout without profiles
  differ from GCC's, and its PGO is often worth more. Its LTO and PGO paths are built into the CMake module, but
  untested.

## Measured on MSVC (2026-10-13, the owner's laptop)

MSVC 14.51, Ninja, 24 threads; the check workloads above, single-threaded, A B C C B A. Repeat runs of one build
agree to about 1%.

| Tool | Release | LTO | LTO + PGO |
|---|---|---|---|
| `study_parametric_bootstrap` | 24.04, 24.07 s | 23.48, 23.47 s (−2.4%) | 23.91, 23.95 s (−0.5%) |
| `study_lrpd_serial_m6` | 78.89, 77.95 s | 78.63, 77.69 s (−0.3%) | 83.32, 82.41 s (+5.7%, slower) |

- **Identity:** all 12 runs byte-identical to Release.
- **LTO:** a real 2.4% on the bootstrap tool, nothing on M6.
- **PGO:** slows the M6 tool. The training is small (one replicate of scenario 49 at φ = 0.3) and the check runs
  another case (37 at φ = 0.6). A likely reason, not checked: code the training never ran is laid out as cold.
- **A defect found in the MSVC path:** the instrumented binaries wrote their `.pgc` files next to the executables,
  not next to the `.pgd` files in `VCAL_PGO_DIR` as `VcalLtoPgo.cmake` assumed. So the USE link merged no counts
  (`pgomgr /summary`: 0 dynamic instructions), and `p13_pgo.py --check` timed an untrained build without noticing.
  The PGO column above is a corrected build, with the `.pgc` files moved by hand and the binaries relinked (3.6e12
  and 9.3e10 dynamic instructions). `--check` also timed each build once, which cannot resolve 3%. Neither is
  fixed, since the options are dropped.

## Decision

Under the stop rule written before the build (under 3% on both tools: drop the options), met on both compilers:
P-13 is rejected and the options are not merged. If it is ever revisited, the `.pgc` placement (`pgomgr /merge`
into `VCAL_PGO_DIR`), a training workload covering the checked cases, and repeated, interleaved check timings
come first.
