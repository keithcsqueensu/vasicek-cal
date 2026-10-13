# P-13: link-time and profile-guided optimisation for local study builds (D-403)

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

## Status: the owner's decision, after the MSVC measurement

On GCC, LTO is worth nothing and PGO about 0-4%, inside the noise: not worth adding there. The runs that matter,
the M6-class runs, are built with MSVC, where PGO is often worth more and where these options are untested. So
the next step is the laptop: `python perf/p13_pgo.py --build-dir build/pgo --generator Ninja --check
build/cpu-release/tests` in a developer prompt. Then apply the same stop rule (under 3% on both tools: record and
drop the options). Until then this is a draft and is not merged.
