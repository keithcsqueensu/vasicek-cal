# State

Last updated: 2026-09-27 (pre-publication pass, D-146)

## Current status

- **M0 (design): done.** Q1–Q13 answered and recorded as D-025–D-046. ARCHITECTURE.md updated to match.
- **M1.1 (build skeleton): done.** Verified locally on Windows and in CI on Linux + Windows.
  - Git repository on `main`. Contains LICENSE (the Apache-2.0 text, verbatim from apache.org), NOTICE, `.gitattributes` (LF) and `.gitignore`.
  - CMake ≥ 3.25 with presets `cpu-debug` and `cpu-release` (Ninja), warnings-as-errors in presets, FP contraction off (D-048), and the CUDA option/arch logic (D-026, D-030, D-040).
  - In-house test runner under CTest, plus self-tests: passing semantics, and a deliberately failing suite whose exit code and failure count are both asserted.
  - `layering_check` enforces the ARCHITECTURE §2 rules; `layering_check_fires` proves it catches all 5 fixture violations.
  - `core/precision.hpp` holds the precision policies and the `Precision` concept, which enforces FP64 accumulation (D-039).
  - **Verified locally:** 7/7 tests pass on MSVC 19.51 (both presets) and the MSVC 14.44 toolset (Release) with `/WX`. Verified by assembly that MSVC `/fp:precise` emits no FMA even with `/arch:AVX2`, and `/fp:contract` does.
  - **CI:** `.github/workflows/ci.yml` (GCC 11 floor, GCC 14, Clang 18, MSVC 2025 runner; both presets). The first run was green on all 8 jobs. This machine has no Linux toolchain, so CI is the only GCC/Clang check.
  - **Not exercised:** `cmake/VcalCuda.cmake`. The development machine has an NVIDIA driver and a Blackwell (sm_120) GPU, but **no CUDA toolkit installed**. That's fine for M1–M3; install CUDA ≥ 12.8 before M4.
- **Q14 resolved** (D-052–D-057): `dgp/` gets its own minimal math and a cross-OS panel-hash check. **C++17 zone** for device-visible code decided (D-058).
- **M1.1b (C++17 zone): done** (D-061–D-065).
  - `core/precision.hpp` now uses traits plus per-clause `static_assert` messages instead of concepts.
  - `cxx17_header_check` compiles every header in core/, engine/ and backends/cuda/ as C++17. The `vcal_core` target now requires only C++17.
  - Two compile-fail tests prove that the check and the messages fire.
  - `VcalCuda.cmake`: CUDA is C++17, with an nvcc floor of 11.8 and the host-compiler range documented.
  - **Verified locally:** 9/9 tests pass on MSVC 19.51 (both presets) and 14.44 (Release).
- **M1.2 (core/special): done** (D-066–D-079).
  - Functions: `log_phi`, `erfcx`, `probit` / `probit_upper`, `log_add_exp`, `lbinom` and `stirling_error`, all header-only and C++17.
  - `tools/gen_special_tables.py` (mpmath 1.4.1, 50 digits) writes every coefficient, constant and table into `core/special/generated/`, plus the golden CSVs and `MANIFEST.json`. `--check` proves it regenerates them exactly.
  - Worst observed error against mpmath (MSVC 19.51):

    | Function | Worst error |
    |---|---|
    | `erfcx` | 2 ulp |
    | `log_phi` | 4 ulp relative, both tails |
    | `probit` | 6 ulp |
    | `lbinom` | 2 ulp (lgamma difference: 14,823 ulp at n = 10⁶) |
    | `stirling_error` | 2 ulp |
    | `log_add_exp` | 1 ulp |

  - Tolerance register started (`docs/methodology/tolerances.md` ↔ `tests/tolerances.hpp`, sync-checked).
  - **Verified locally:** 11/11 tests pass on MSVC 19.51 (Debug + Release, identical errors) and 14.44.
- **M1.3 (quadrature): done** (D-080–D-087).
  - GH tables for N ∈ {8, …, 128} generated from mpmath. They are checked at 50 digits for symmetry, all moments to degree 2N−1, and the closed-form weights, and stored as log-weights.
  - `GaussHermiteFixed` and `GaussHermiteAdaptive` are C++17, device-visible and deterministic.
  - The binomial-mixture log-integrand and its deterministic mode/scale hint are in core/model; `inverse_mills` joins core/special.
  - Shared generator machinery lives in `tools/tablegen.py`.
  - **Results against mpmath:**
    - Moments are relative-exact to 2.04·(m+1)·ε.
    - The adaptive rule is exact for Gaussian log-integrands (1.6 ε).
    - Adaptive GH is at full precision (≤ 4 ε) for N ≥ 64 and ρ ≤ 0.5, from n = 1 to 10⁶. ρ = 0.9 with n ≤ 2 converges slowly (5.6e-14 at N = 128).
    - Fixed GH is wrong by ~12 in log-likelihood at n = 10⁶.
  - Review fixes: the stale ARCHITECTURE header is gone, `VCAL_CHECK_ULP_STATS` reports the caller's line, and device math is recorded as known M4 risk K-1.
- **M1.3b: done** (D-088–D-091).
  - Parity N = 128 with a 256-node table for the twice-the-nodes check; the moment checks now run in log space.
  - Golden cases extended to ρ = 0.95/0.99 and to low-default periods.
  - That exposed a correction to D-086: zero-default periods need N = 128 at ρ = 0.5 (4.4e-13). Above ρ = 0.5 they are not at full precision at any tabulated N (7e-5 at ρ = 0.9, n = 50, N = 128).
  - The claims are stated as measured regimes A and B, and the rest is open item Q15.
  - A split-at-mode half-range rule was prototyped and rejected (D-091).
- **Q15 and O-1 resolved** (D-092–D-094). Parity claims precision for ρ ≤ 0.5, with a separate materiality statement above it, and the 128-vs-256 check is proven never to miss a golden period. Composite Gauss–Legendre for d ∈ {0, n} is a native enhancement for later. CUDA floor 11.8 + MSVC 14.39 is confirmed, so the C++17 device zone is required.
- **M1.4 (engine, CPU): done** (D-095–D-098).
  - `core/grid`: four axis scales, with the ρ bound configurable (default 0.5).
  - The `BinomialMixture` objective, and the `ArgMax` reducer, which merges identically in any order.
  - `engine/`: surface, compensated weighted sums, fixed-tile reduction, 3×3 quadratic refinement with delta-method SEs, and calibrate with the per-run quadrature check.
  - `backends/cpu`: OpenMP `parallel for` only.
  - **Measured:**
    - Surface and estimate are bitwise identical for 1, 2, 3 and 8 threads.
    - The quadratic stencil recovers the vertex and covariance to 9e-16.
    - On the test panel the coarse-grid refinement is within 0.010 SE of a 10×-finer grid, and the SEs agree to within 1%.
    - Edge estimates (no defaults; extreme clustering above ρ = 0.5) are flagged, never clamped.
  - 15/15 tests pass on MSVC 19.51 (Debug, Release) and 14.44.
- **K-1 retired** (D-100). nvcc 11.8 and 12.8 compile-only CI jobs (Linux containers, GCC 11) compile every device-visible function as real device code: cubins for sm_89 (+ sm_120 on 12.8) and compute_80 PTX, with host-only calls from device code as errors, proven by a compile-fail fixture. The Windows nvcc/14.39 job is deferred (D-101).
- Owner follow-ups recorded: 2-D refinement kept (D-103); golden-estimate tolerance 1e-12 relative (D-102); the ×5 quadrature-check factor is labelled as an empirical bound.
- **Publication hygiene** (D-104). Committed docs use neutral wording; machine specifics live in gitignored `*.local.md` notes. The public repository is published from a single squashed snapshot.
- **M1.5 (ref/): done.** Independent FP64 reference, standard library only, no core code or algorithm.
  - Measured against the mpmath goldens: `log_ncdf` 11 ulp, `ncdf_inv` 4 ulp, `lchoose` 1 ulp.
  - `log_mixture` is at 7.5 ε on every golden case, including ρ = 0.99 and zero-default periods, so it is a valid oracle exactly where core's adaptive GH is weakest.
  - `fit` uses nested golden-section search in logit coordinates, with SEs from a finite-difference Hessian. About 8 s on a 20-period, 1000-obligor panel.
  - Plain-prose description for replication: `docs/methodology/reference_implementation.md`.
- **Build compatibility** (D-105, D-106).
  - CMake minimum ≤ 3.26 and presets schema ≤ 6, with CI pinning CMake 3.25.x and 3.26.x in the nvcc jobs.
  - CMake versions in CI (verified prior to public release): Linux runners 3.31.6, Windows runner 4.4.3, nvcc jobs 3.25.2 and 3.26.4.
  - CUDA architectures must be explicit numbers (all/all-major/native rejected at configure).
  - CUDA builds require Ninja (a Visual Studio generator is rejected).
  - Machine-specific notes moved to gitignored `*.local.md` notes.
- **M1.6 (dgp/): done** (D-110–D-112).
  - Philox4x32-10 matches Random123's known answers.
  - Own `log`/`exp` (1 ulp), Φ (≤ 1.1e-13) and AS241 Φ⁻¹ (6 ulp), from exact or correctly rounded operations only.
  - Polar factor on blocks 0 … 2¹⁶−1 and Bernoulli sums from block 2¹⁶; uniforms (k + ½)·2⁻⁵².
  - A line-for-line Python mirror produces six reference panels and their SHA-256. C++ matches the mirror bit for bit, and CI checks the hash on every platform.
  - Stream addressing is proven: more periods or a changed n elsewhere leaves every other period unchanged.
  - The moments agree with the model over 40,000 panels (batch-means SEs).
  - Prose description: `docs/methodology/dgp.md`.
  - Follow-ups done: Φ uses fixed iteration counts (128 Lentz, 32 series; D-113); the reference-panel hash is unchanged (recorded in D-113); every CI job prints the hash (D-114).
- **M1.7 (cross-reference): done** (D-116, D-117).
  - Core vs ref, per period, 0 < d < n: within 35 ε (term-scaled) over 1,704 points up to ρ = 0.9.
  - Estimates within 0.024 SE; SEs within 10.4% (Q17). Cross-platform regression values at 1e-12.
  - **Finding:** zero-default periods with few expected defaults are imprecise in core at any ρ (up to 1.9e-5 absolute at ρ ≤ 0.5, 6e-3 at ρ = 0.9). The precision claim is restated to 0 < d < n; the 128-vs-256 check catches every period above threshold. Q16 proposes the dedicated rule in parity.
- **M1.7b (Q16 and Q17): done** (D-118–D-120; D-107 amended).
  - **Parity rule (D-118):** adaptive GH (N = 128) for 0 < d < n; composite Gauss–Legendre, 16 panels × 16 points under a sinh map centred on the survival factor's half-point, for d ∈ {0, n}. Fixed panel count and placement. The per-run check is the same rule doubled.
    - Against mpmath: every golden case is within 11.7 ε, including ρ = 0.99.
    - Against ref: 4,296 points (all periods) within 35.2 ε, term-scaled.
    - Zero-default surfaces are smooth: second differences of core vs ref at 3.7e-14.
    - The D-116 exclusion and its materiality statement are withdrawn.
  - **SEs (D-119):** the objective's own Hessian at the refined estimate, step 0.15 × the stencil SE (sensitivity 0.10–0.20: 4.7e-5). New flag `kFlagNearBound` within 2 SEs of a bound; profile-likelihood intervals are due in M2. Core vs ref SEs: within 1.56% (was 10.4%), near-bound panels reported only.
  - **Check threshold (D-120):** max(1e-10, 64·ε·term size), so rounding at n = 10⁶ is not flagged. Detection is tested on a deliberately coarse rule of the same family (ratio 0.99998).
  - **Fourth cross-platform tolerance:** SEs at 1e-10 (D-107). Regression values were regenerated.
- **M1.8 (recovery harness): done** (D-121–D-124).
  - **Setup:** 81 scenarios × R = 1000 DGP replicates, each fitted exactly as a user would. Goldens in `tests/golden/recovery/` with provenance; generated tables in `docs/methodology/recovery_results.md`; design and interpretation in `recovery.md`.
  - **Engine correctness holds:** RMSE falls with T everywhere, and no fit of the 81,000 is flagged for quadrature or numerics.
  - **Coverage (Wald, conditional on unflagged):** 81 PASS, 62 DEFERRED (at least 5% flagged; waits for M2's profile intervals), 19 KNOWN FINDING, 0 UNREVIEWED. The findings are 14 small-T undercoverage and 5 ρ overcoverage on skewed estimates, each listed with its diagnosis and the SE-ratio evidence.
  - **D-122:** the surface evaluates repeated observations once, with bit-identical results. It speeds fits 3–4× on low-default panels.
  - **CI:** replays 162 fits on every platform, and enforces the exit criteria and the pinned findings on the committed summary. `recovery_harness --check` does the full comparison (exact on the writing platform, counts ±2 elsewhere).
- **M1.9 (methodology and scipy replication): done** (D-125–D-127).
  - `docs/methodology/binomial_mixture_mle.md`: everything needed to replicate the fit, in one place.
  - `validation/scipy/binomial_mixture_mle.py`: 125 lines of plain scipy. It agrees with the engine to 4.1e-13 on 702 surface points, its own optimum is within 3.9e-4 of log-likelihood of the engine's estimate, and the estimates agree within 0.026 SE. It runs in CI.
  - Tolerances now live in `tests/tolerances.toml`, from which the C++ header is generated; the TOL_SCIPY_* entries are a fifth, named kind of comparison.
  - `docs/methodology/parity_signoff.md` lists the D-123 items, the 19 known findings and the 31 deferred scenarios.
- **M1 exit: parity sign-off**, verified in CI prior to public release: all 12 jobs green, and the scipy job on Linux reproduced the Windows figures exactly. In the public repository the tag `v0.1.0-parity-binomial` marks the initial release snapshot, which contains M1.
- **M2a (profile-likelihood intervals): done** (D-128–D-131).
  - `engine/profile.hpp`:
    - polished off-grid maximum;
    - crossing bracketed from the grid, then solved against the objective (residual ≤ 1e-7, asserted on every recovery fit; worst 1.7e-8);
    - truncated at the box with flags, never extrapolated;
    - χ²₁ threshold with truncation (D-129), with the boundary mixture for testing ρ = 0 documented.
  - An independent `ref/` implementation agrees to 1.5e-10 in logit units, and on every truncation.
  - **Recovery re-run:** profile, Wald and t(T−1) reported side by side. Profile: 136 PASS, 20 CONSERVATIVE (safe side, truncation in near-uninformative settings), 6 KNOWN FINDING (small-T undercoverage 0.920–0.927). It fixes 15 of the 19 Wald findings and 41 of the 62 deferred verdicts.
- **M2b (resampling engine): done** (D-132–D-137).
  - `resample/`:
    - Philox in its own key domain, checked bit for bit against a Python mirror;
    - iid and moving-block bootstrap, jackknife, walk-forward, and external index matrices;
    - a fused reduction, never B × K, bitwise equal to a materialised one;
    - percentile intervals.
  - **Coverage study** (full matrix, iid, B = 999; predictions committed before the run): 37 PASS and 125 below the band, all pinned.
    - Boundary breakdown near bounds; no bias or skewness correction elsewhere; ρ undercovers at every T.
    - **Conclusion:** profile-likelihood intervals are the recommended method for inference; bootstrap percentile intervals are not recommended for ρ.
    - Moving-block coverage moves to M6; BCa stays on the backlog.
- **M2c (C ABI v0.1): done** (D-138–D-144). M2 is complete.
  - `include/vcal/vcal.h`, plain C99, and `abi/`, the shared library `vcal`.
  - **Functions:** calibrate (with profile intervals), surface, resample (six W sources, percentile intervals, the implied W), DGP simulate, grid helpers, build info, errors and the context.
  - **Conventions:**
    - a caller-set `struct_size` on every struct; unknown nonzero fields are refused;
    - no implicit padding;
    - caller-owned buffers with size queries;
    - a status from every call and a thread-local message naming the field;
    - a catch-all at every entry point;
    - `__cdecl`, and `vcal_*` exports only (version script on Linux, checked against the header on Linux and MSVC).
  - **Tests:**
    - `abi_c` (C11): 866 checks, including the DGP reference panels bit for bit;
    - a C99 `-pedantic-errors` header check;
    - `abi_engine_equivalence`: every number bitwise equal to direct engine calls, all six schemes;
    - `abi_exports` and `abi_exports_fires`.
  - **Verified locally:** MSVC 19.51 Release 32/32 and Debug 29/29. The README's C example compiles with `/W4 /WX` and reproduces the C++ example.
- **CI economy (D-145):**
  - the full 12 jobs run on pushes to `main` and on pull requests;
  - other branches get a reduced set of 3 jobs;
  - the slow tests run nightly (only if `main` moved) and on `workflow_dispatch`;
  - prose-only commits start no run. The Markdown files the code reads are re-included, and `ci_path_filter_sync` derives and checks that list: today `tolerances.md`, `recovery_results.md` and README.md, whose C example `readme_c_example` compiles, runs and compares with its documented output;
  - newer pushes cancel older runs, except on `main`.
  - **Why:** faster feedback and fewer redundant runs. Each Release job had spent 7–8 minutes on the slow tests, and the accuracy step ran `recovery_replay` a second time.
- **Pre-publication pass (D-146):**
  - `ci-economy` (D-145) merged;
  - the project is described as an independent credit-risk research library for practitioners, researchers and students, designed so every result can be independently replicated, and a personal open-source project; NOTICE names Keith Yip, and LICENSE is the unmodified Apache-2.0 text;
  - toolchain floors restated as portability goals, with every technical detail kept;
  - ARCHITECTURE.md §3.0 gives each layer a language-agnostic contract;
  - `ci_path_filter_sync` also scans CUDA, TOML and JSON test inputs and `CMakePresets.json`.
  - CI run numbers and commit hashes from before the public release are replaced by neutral wording ("verified in CI prior to public release");
  - the public repository is published from a single squashed snapshot, tagged `v0.1.0-parity-binomial`.
  - **Verified locally (Linux, GCC 13.3, cpu-release):** 35/35 tests pass, slow tests included.
- **Open:** R-1 is deferred to M6 and R-2 to M7. The toolchain floor is a portability goal: CUDA ≥ 11.8 and older host compilers (GCC 11, MSVC 14.39), with device-visible code in C++17 (D-094).

## Next steps

1. M3: the remaining single-factor estimators. Awaiting the owner's go.
2. Backlog (`native` options, each evaluated against the pinned recovery verdicts; D-131, D-137):
   - a Bartlett-type correction of the profile threshold;
   - a bias-corrected ρ̂, the small-sample bias benchmark;
   - BCa bootstrap intervals (they will not fix the boundary breakdown near a bound);
   - the stationary bootstrap and a data-driven block length (Politis–White); moving-block coverage is studied in M6 with AR(1) data.

## Milestones

Ordering principle: the CPU path is complete and validated before any GPU development, since the
GPU is checked against it (D-010, D-011). Each milestone ends with its methodology doc, its
tolerance-register entries and its golden values.

### M1 — Core math + ref + synthetic DGP + recovery harness (single-factor binomial MLE, CPU only)

| # | Item | Done when |
|---|---|---|
| M1.1 ✅ | CMake skeleton + presets (cpu-debug, cpu-release), warnings-as-errors, per-target FP-contraction flags, in-house test runner under CTest, CI on Linux + Windows (D-025–D-030, D-047–D-051) | empty targets build and one test runs on both OSes |
| M1.2 ✅ | `core/special`: `log_phi` (upper tail via `log1p(−Φ(−x))`, D-068), `erfcx`, `probit` + `probit_upper(q)` (D-069), `log_add_exp`, `lbinom` via Stirling differences (D-070); mpmath table generator + manifest with version and `mp.dps` (D-067) | matches mpmath tables (D-033) over documented ranges (e.g. `log_phi` on [−40, 10]); **relative** error asserted on both tails of `log_phi`; `probit` tested near 0 and near 1; `lbinom` tested at n up to 10⁶ and beyond; `log_add_exp` edge cases for −∞ and equal arguments (D-071); every tolerance entered in the register |
| M1.3 ✅ | `core/quadrature`: mpmath-generated GH tables (D-080, D-081), `GaussHermiteFixed`, `GaussHermiteAdaptive`, binomial-mixture integrand + hint (D-084) | exact to degree 2N−1 on polynomials; `E[p(Z)] = PD` and `E[p(Z)²] = Φ₂(·;ρ)` closed forms reproduced; spike integrands (n = 10⁶) converge under the adaptive rule |
| M1.4 ✅ | `core/objectives/BinomialMixture`, `core/grid`, `engine` + `backends/cpu`, `ArgMax` reducer, sub-grid refinement (D-038), curvature SEs, edge/flat/convergence flags | calibrate runs end-to-end through the engine with W = one row of ones; results bitwise identical for 1, 2 and N threads |
| M1.5 ✅ | `ref/`: independent FP64 Φ/logΦ (std::erfc + continued-fraction tail), adaptive Gauss–Kronrod on a truncated/transformed z range, a different optimiser (e.g. nested golden-section) | shares no code with core (include check passes) |
| M1.6 ✅ | `dgp/`: Philox4x32-10, deterministic normal and binomial samplers, scenario spec — CPU only (D-034). Own `log`/`exp`/`sqrt`/`erfc` → Φ, Φ⁻¹ via AS241, Marsaglia polar normals, Bernoulli-sum binomials, `/fp:strict` (D-052–D-056) | Philox matches the published known-answer vectors; sampler moments match theory; outputs are reproducible from (seed, stream) alone; DGP functions match mpmath within 1e-12 relative; the reference-panel hash is identical on every CI OS (D-057) |
| M1.7 ✅ | Crossref suite: core vs ref | per-period ℓ_t agrees within tolerance on ~10⁴ fixed (PD, ρ, n, d) points including tails (PD 1e-6…0.5, ρ 1e-4…0.9, n up to 10⁶, d = 0 and d = n); estimates agree within the grid-refinement tolerance |
| M1.7b ✅ | Parity rule for d ∈ {0, n}: composite Gauss–Legendre (Q16, D-118); SEs from the objective's Hessian at the estimate (Q17, D-119); term-scaled check threshold (D-120) | every golden case within tolerance at every ρ; every cross-reference period within tolerance; zero-default surfaces smooth; SE step sensitivity bounded; SE cross-platform tolerance in CI |
| M1.8 ✅ | Recovery harness: scenario matrix PD ∈ {0.1%, 1%, 5%} × ρ ∈ {0.02, 0.12, 0.24} × T ∈ {20, 40, 100} × n ∈ {100, 10³, 10⁴}, R = 1000 seeded replicates each; bias, RMSE, flagged fraction, conditional coverage (D-121, D-124) | engine correctness holds (RMSE falls with T, no quadrature/numeric flags, cross-platform replay); coverage PASS inside the Monte Carlo band, DEFERRED at ≥ 5% flagged, or a reviewed KNOWN FINDING, with none unreviewed; tables committed as goldens with provenance |
| M1.9 ✅ | Docs: `docs/methodology/binomial_mixture_mle.md` (plain prose), tolerance register, `validation/scipy/binomial_mixture_mle.py` (D-032); the note explains why fixed-node GH fails at large n (D-037) | the scipy script reproduces the M1.7 surface cells within the documented tolerance |

**M1 exit:** all of the above green on Linux and Windows. The
ρ-MLE's known small-T downward bias is documented, not hidden.

### M2 — Profile intervals, resampling engine + C ABI v0 (CPU)

| # | Item | Done when |
|---|---|---|
| M2a ✅ | Profile-likelihood intervals (D-128–D-130): polished maximum, grid-bracketed crossing solved against the objective, truncation flagged; an independent `ref/` implementation; recovery re-run with profile, Wald and t(T−1) side by side (D-131) | endpoint residual and core-vs-ref tolerances met; every DEFERRED verdict and M1.8 KNOWN FINDING re-assessed |
| M2b ✅ | Resampling engine: the fused, tiled `reduce_weighted` builds W for iid bootstrap, block bootstrap, jackknife, walk-forward and custom schemes (D-042, D-045) | bootstrap CIs join the recovery harness; W × L reuse bitwise identical for any thread count |
| M2c ✅ | C ABI v0: `vcal_calibrate` (with the profile intervals), `vcal_resample`, `vcal_surface`, `vcal_dgp_simulate`, build info and errors, with a pure-C ABI test | a `.c` test calibrates, resamples and simulates through the ABI |

### M3 — Remaining single-factor estimators (CPU)

Vasicek-rate MLE, method of moments (joint-default-probability form, exact in finite n) and
grid Bayesian (`LogSumExpPosterior`, documented priors). Parity and native profiles become
real registry whitelists (D-035, D-036, D-044). Each estimator gets recovery tests, goldens, a
methodology note and a scipy script.

### M4 — CUDA backend

Surface and fused-reduce kernels, a fat binary (sm_89 + sm_120 when nvcc ≥ 12.8 + compute_80 PTX, D-026/D-040), precision
policies (D-039), and a CPU↔GPU parity suite with tolerance definitions (D-041). Install the CUDA toolkit first. Performance goes
in `perf/` on both architectures. The nvcc compile-only CI jobs already exist (D-100). Install the development machine's toolkit (11.8 + MSVC 14.39, D-094) before starting.

### M5 — Monitoring

Backtests (exact binomial, Jeffreys, correlation-adjusted binomial, traffic light). Also a
run ledger (canonical serialisation + SHA-256, links calibration and backtest runs) and
report data plus rendering outside the core (D-045).

### M6 — Benchmark models

Multi-segment correlated factors (two-stage + pairwise composite likelihood, PSD repair
documented), AR(1) factor (forward filter, parametric bootstrap (D-043), revisit GPU-side DGP (R-1)), alternative mixing
distributions, Pluto–Tasche and window/influence analysis.

### M7 — Macro pipeline

Transform library with lineage, pre-screen filters, batched specification fits (fits are
the B dimension), pluggable scorers and a config-driven ranker (config parsed outside the core, D-045).

Deferred item R-2 (DECISIONS.md): the Z sign convention and macro sign filters. Expected effects are declared in economic terms (worsens or improves credit conditions) and mapped to coefficient signs in one place through the engine's Z convention (higher Z = better conditions). The mapping is tested on a synthetic DGP. Revisit when the macro pipeline starts interacting with the estimator.

### M8 — Public-data demonstrations + replication pack + ABI 1.0

FDIC failed-bank counts as the binomial showcase; FRED/StatCan rate series via Vasicek-rate/MoM (D-046), end-to-end worked examples, a complete scipy replication
set, and an ABI freeze at 1.0.
