# Decisions

Append-only. A decision that is later reversed gets a new entry that names the one it
supersedes; the old line stays. Format: **ID — decision.** rationale.

## Adopted — from the project charter (CLAUDE.md), recorded 2026-09-26

- **D-001 — Core is C++ + CUDA; no Rust in the core.** One toolchain compiles both host and device code; Rust stays a consumer via the C ABI.
- **D-002 — Core dependencies: CUDA toolkit + CMake only.** Smallest supply-chain surface for anyone reviewing, auditing or packaging the code.
- **D-003 — Stable, flat C ABI (`extern "C"`, plain structs, pointers + scalars).** Any language can bind it, and engine releases don't break consumers.
- **D-004 — Optional Rust CLI in a separate repo, linking the C ABI.** Keeps the core single-language and proves the ABI is sufficient.
- **D-005 — Apache-2.0.** Permissive with an explicit patent grant, and widely accepted for commercial, academic and personal use.
- **D-006 — Public or synthetic data only (DGP fixtures, FDIC, FRED, StatCan).** Every result can be reproduced by anyone, with no confidentiality barrier to independent replication.
- **D-007 — Policies composed as templates: Objective × Integrator × Reducer × Precision.** Each axis is tested on its own and swapped without duplicating code.
- **D-008 — Surface engine: evaluate a scalar surface over a parameter grid, then reduce.** Global view of the objective, no local-optimiser convergence failures, and it parallelises trivially.
- **D-009 — Per-period surfaces L (T × grid) computed once; bootstrap, walk-forward, jackknife are weight matrices W, surfaces = W × L.** The expensive quadrature runs once, and every resampling scheme becomes a cheap weighted reduction.
- **D-010 — CPU backend (OpenMP) is first-class for correctness and CI.** Correctness and CI must not depend on having a GPU.
- **D-011 — GPU must match CPU within tolerance.** The CPU backend is the engine's reference; the GPU is an accelerator, not a second source of truth.
- **D-012 — `ref/` is an independent FP64 reference that shares no quadrature code.** A bug in shared code would be reproduced by both sides and pass unnoticed.
- **D-013 — Counter-based RNG (Philox), deterministic across CPU/GPU.** Random streams are addressed by (seed, counter), so results don't depend on thread scheduling or device.
- **D-014 — Two profiles: `parity` (textbook, replicable) and `native` (enhancements).** Anyone replicating a result gets a path they can reproduce exactly; enhancements are opt-in and labelled.
- **D-015 — Numerics in log space.** Tail PDs and large-n binomial likelihoods underflow otherwise.
- **D-016 — Every tolerance is documented.** Anyone replicating a result has to know what "matches" means.
- **D-017 — Every feature ships with synthetic recovery tests and golden values.** Recovery shows the estimator recovers the truth; goldens catch silent regressions.
- **D-018 — Anything needed to replicate a result is describable in plain prose and a short scipy script.** Independent replication is the project's central requirement.
- **D-019 — Fat binaries for sm_89 and sm_120 plus PTX; don't tune for one GPU.** Runs natively on current Ada and Blackwell parts, and PTX gives forward compatibility.
- **D-020 — Methodology scope for 1.0 is fixed:** binomial-mixture MLE, Vasicek-rate MLE, MoM, grid Bayesian; multi-segment, AR(1), alternative mixing, Pluto–Tasche, window/influence benchmarks; macro pipeline; backtests, ledger, report. A bounded scope keeps 1.0 reviewable.
- **D-021 — STATE.md and DECISIONS.md are updated at the end of every session.** Records how and why the engine reached its current state.

## Derived — consequences of the charter, not new choices

- ~~**D-022 — CUDA toolkit ≥ 12.8.**~~ *Superseded by D-026.* It is the first toolkit that can target sm_120 (follows from D-019).
- **D-023 — The public C ABI is the only installed header.** Follows from D-003: templates can't cross a C ABI, so supported policy combinations are pre-instantiated and chosen by enum at runtime.
- **D-024 — `ref/` also shares no special-function code with `core/`.** Follows from the intent of D-012; if Φ were shared, a Φ bug would pass crossref. *(Confirmed by D-031.)*

## Resolved — owner's answers to Q1–Q13, 2026-09-26

- **D-025 (Q1) — C++20 throughout.** *Narrowed by D-058.* Concepts state the policy contracts in code and give readable errors.
- **D-026 (Q1; supersedes D-022) — No hard CUDA floor. sm_120 SASS is emitted only when the toolkit is ≥ 12.8; an older toolkit builds sm_89 + PTX.** Ada (sm_89) hardware is widely installed, so reaching it matters more than Blackwell support. The real floor is set once the toolchain-floor question is answered (O-1; settled by D-094).
- **D-027 (Q1) — Linux and Windows are both first-class. Host compilers: GCC ≥ 11, MSVC 2022+, Clang ≥ 14.** Provisional until O-1.
- **D-028 (Q2a) — OpenMP allowed, only as `parallel for` over independent outputs, with no reduction clauses.** Compatible with MSVC's default `/openmp` (2.0); `/openmp:llvm` (3.x) exists but isn't needed. Deterministic by construction.
- **D-029 (Q2b) — In-house test runner registered with CTest; no GoogleTest or Catch2.** Keeps D-002 intact.
- **D-030 (Q2c) — `VCAL_ENABLE_CUDA` defaults ON only if a CUDA compiler is found; the CPU-only build needs no toolkit.** CI runners without a GPU can build and test everything else.
- **D-031 (Q3) — `ref/` is C++, depends only on the standard library, and shares no code with core: no quadrature, special functions or optimiser.** Independence has to cover every place a shared bug could hide.
- **D-032 (Q3) — scipy replication scripts live in `validation/scipy/`, outside the build, and run in an optional CI job.** Anyone replicating a result gets runnable prose-equivalents, and the core build stays Python-free.
- **D-033 (Q3) — Goldens are generated by `ref/`, cross-checked by scipy, and checked in with a provenance manifest (generator, version, command, date, content hash). Special-function tables come from mpmath, generated once and checked in.** The engine never certifies itself, and anyone can trace every expected value.
- **D-034 (Q4) — The DGP runs on the CPU only; GPU code consumes copied panels.** It's simpler and reproducible without a device port. Revisit at M6 (R-1).
- **D-035 (Q5) — `parity` = FP64 only, mode-centred adaptive GH with a fixed documented N, fixed-order reductions, argmax + parabolic refinement. `native` = a superset, with each enhancement labelled.** Every parity step fits in a short scipy script.
- **D-036 (Q5) — The public `parity` profile means "textbook and replicable in a short scipy script", not conformance to any particular organisation's method. A user who needs one defines their own conformance profile as a whitelist over the registry.** Keeps the public profile neutral and makes the registry the extension point.
- **D-037 (Q5) — The large-n failure of fixed-node GH is documented in the methodology note (M1.9).** It's the reason parity uses the adaptive rule, and anyone replicating the method will ask.
- **D-038 (Q6) — Per-axis parabolic sub-grid refinement; SEs from finite-difference curvature of the surface; grid-edge estimates are flagged, never clamped.** Stays inside `W × L` and fits in one sentence of prose.
- **D-039 (Q7) — On the GPU, `native` may evaluate the integrand in FP32. Accumulation stays FP64, `parity` stays FP64 on every device, and achieved agreement is reported per run.** sm_89/sm_120 run FP64 at ~1/64 rate.
- **D-040 (Q8; refines D-019) — SASS for sm_89 (+ sm_120 per D-026) plus PTX for compute_80. A100/H100 "should run" via JIT; tested hardware is Ada and Blackwell only.** PTX JITs forward only, so compute_80 gives the widest reach; we test only what we have.
- **D-041 (Q9) — CPU/GPU agreement is defined at two levels: surface (relative, F64 vs F64) and estimate (absolute, curvature-derived, for mixed precision).** Surface checks catch kernel bugs; estimate checks bound what users see.
- **D-042 (Q10a) — L generalises to "per-period contributions": grid values or sufficient statistics. MoM reduces by moment inversion.** Keeps one resampling mechanism for all separable estimators.
- **D-043 (Q10b; exception to D-009) — AR(1)-factor models are resampled by parametric bootstrap through `dgp/`.** Their likelihood isn't separable across periods, and reweighting destroys the dependence.
- **D-044 (Q11) — Zero or 100% default rates: `parity` refuses and names the periods; `native` offers a censored likelihood.** No silent data edits in the replicable path.
- **D-045 (Q12) — Opaque context handle allowed in the ABI. Core emits ledger records as canonical bytes + in-house SHA-256; config parsing and report rendering live outside the core, which accepts structs only.** Hashes are identical whichever front-end writes them, and the core stays dependency-free.
- **D-046 (Q13) — FDIC failed-bank data is count data (d_t = failures per year, n_t = insured institutions) and is the primary real-data showcase for the binomial estimators. FRED/StatCan rate series use Vasicek-rate/MoM. Dated snapshots + checksums go under `data/`, with fetch scripts outside the core.** Correction from the owner: the original Q13 treated FDIC as rate data.

## Build mechanics — M1.1, 2026-09-26

- **D-047 — CMake ≥ 3.25 (presets schema v6).** Oldest CMake with the preset features used; provisional until O-1.
- **D-048 — FP contraction off on every host target (`-ffp-contract=off`; MSVC `/fp:precise`); fast-math is never used.** Results must not depend on whether a compiler chose to fuse a multiply-add.
- **D-049 — `VCAL_WARNINGS_AS_ERRORS` is an option: ON in presets and CI, OFF by default.** Development stays strict, and a user building from source with a newer compiler isn't blocked by a new warning.
- **D-050 — The layering rules in ARCHITECTURE.md §2 are enforced by a CMake-script check that runs under CTest, with a negative fixture proving it fires.** Runs identically on every OS with no Python or bash.
- **D-051 — Every source file carries `SPDX-License-Identifier: Apache-2.0`.** Automated licence scanners read it.

## Resolved — owner's answers on Q14 and the nvcc/C++20 interaction, 2026-09-26

- **D-052 (Q14) — `dgp/` has its own elementary functions. The requirement is bitwise identity across OSes; accuracy only needs to be ≤ 1e-12 relative, documented.** Determinism matters more than accuracy here, which makes these functions much simpler than a real libm.
- **D-053 (Q14) — The DGP function set is minimal: `log`, `exp`, `sqrt` (IEEE correctly rounded, so the compiler's is fine), `erfc` by rational approximation → Φ, and Φ⁻¹ by AS241.** Nothing else, so there is little to verify. *Confirmed by D-061.*
- **D-054 (Q14) — Normals by the Marsaglia polar method.** It needs only `log` and `sqrt`, avoiding Box–Muller's `sin`/`cos`. Rejection uses a variable but deterministic number of counter values per stream.
- **D-055 (Q14) — Binomial draws sum n Bernoulli comparisons `u < p(z_t)` against one p per period.** Slow at n = 10⁶ but trivially deterministic, and fine for fixtures. Move to inversion (which needs an in-house `lgamma`) only if generation time actually hurts.
- **D-056 (Q14; refines D-048 for `dgp/`) — The `dgp/` target builds with `/fp:strict` on MSVC and `-ffp-contract=off` (never `-ffast-math`) on GCC/Clang.** The strictest setting goes where bitwise identity is the contract, and only there.
- **D-057 (Q14) — DGP functions are tested against mpmath tables (D-033), and a cross-OS CI check hashes a reference panel and compares it with a checked-in hash.** Bitwise identity is asserted, not assumed.
- **D-058 (narrows D-025) — Device-visible code is C++17-compatible: no concepts, and policy contracts checked by `static_assert`/SFINAE traits. Everything else is C++20. A test TU compiles every device-visible header with `-std=c++17`.** Cheap now; retrofitting after O-1 showed nvcc 11.x would be painful. If O-1 shows nvcc ≥ 12, a one-line entry can loosen this. *Scope fixed by D-062.*
- **D-059 (version from O-1; fulfilled by D-100) — A compile-only CI job builds with the exact floor nvcc version, targeting sm_89.** Runners can install a specific toolkit without a GPU, so this catches version drift without a GPU or a local install of that toolkit.
- **D-060 — The development machine's CUDA toolkit is installed only after O-1 settles the floor version.** Avoids developing against a toolkit newer than the supported floor.

## Resolved — N-1, N-2 and M1.1b notes, 2026-09-26

- **D-061 (N-1; confirms D-053) — The DGP function set is `log`, `exp`, `sqrt`, an in-house `erfc` for Φ, and AS241 for Φ⁻¹.** `exp` is needed by the `exp(−x²)` factor in rational `erfc` approximations. AS241 needs only `log` and `sqrt` and is accurate to about 1e-16, well within the 1e-12 target.
- **D-062 (N-2; scopes D-058) — The C++17 zone is `core/`, `engine/` and `backends/cuda/`, and the C++17 header check covers all three.** `engine/` templates are instantiated inside `.cu` files, so they are device-visible too.
- **D-063 — Each policy contract is a boolean trait (`is_precision_v<P>`) plus a `check_<policy><P>()` function with one `static_assert` and message per clause. A compile-fail test asserts the message text actually reaches the user.** Recovers most of the readable errors concepts would give.
- **D-064 — The C++17 header check is a generated TU including every header in the zone (globbed, so new headers are covered automatically), compiled with the host compiler for now. Once the D-059 job exists, nvcc compiles the same TU; `.cuh` files join the check then.** The nvcc run is the one that matters, and the host run keeps the rule honest until then.
- **D-065 (replaces the nvcc ≥ 12.0 check implied by D-025) — CUDA floor is nvcc 11.8, the first release with sm_89. This floor comes with its host-compiler range, per NVIDIA's CUDA 11.8 installation guides: GCC ≤ 11, Clang ≤ 14.0, MSVC 19.1x–19.3x (VS 2017 15.x – VS 2022 up to 17.9; 17.10+ reports `_MSC_VER` ≥ 1940). Combined with D-027, a CUDA build on nvcc 11.8 means exactly GCC 11, Clang 14, or MSVC 19.30–19.39.** An nvcc version caps the host compiler, so a bare number would be misleading. Theoretical until O-1. Note: this machine's MSVC (19.44, 19.51) is outside that range.

## M1.2 design, from owner's notes, 2026-09-26

- **D-066 (fallback for D-065) — If O-1 confirms nvcc 11.8, reproduce the floor CUDA build (nvcc 11.8 + MSVC 14.39) locally by installing the MSVC v143 14.39 toolset (VS 2022 17.9) as a side-by-side individual component, selected with `vcvars_ver=14.39`.** No downgrade needed, and no reliance on CI alone. The `vcvars_ver` mechanism is already used here to test the 14.44 toolset.
- **D-067 (extends D-033) — The golden provenance manifest records the generator command, script hash, Python version, mpmath version and working precision (`mp.dps`).** Tables must be regenerable exactly, not just approximately.
- **D-068 — `log_phi` computes logΦ(x) as `log1p(−Φ(−x))` for large positive x, where logΦ → 0. Its tests assert *relative* error on both tails.** Absolute error hides total loss of relative accuracy near 0.
- **D-069 — A complement entry point takes q = 1 − p directly: `probit_upper(q) = Φ⁻¹(1 − q)`, computed as `−probit(q)` by symmetry, so it's exact whenever `probit` is accurate for small q. Tests cover p near both 0 and 1.** Forming 1 − q at the call site loses q's digits.
- **D-070 — `lbinom` uses a cancellation-free Stirling-difference form, not `lgamma(n+1) − lgamma(d+1) − lgamma(n−d+1)`:**
  `log C(n,k) = [δ(n) − δ(k) − δ(n−k)] + k·log(n/k) − (n−k)·log1p(−k/n) + ½·log(n / (2π·k·(n−k)))`, where δ is the Stirling error term (Loader, 2000). δ comes from mpmath-generated constants for small arguments and its asymptotic series beyond; the edge cases k = 0 and k = n are exactly 0. It's tested against mpmath at n up to 10⁶ and beyond. The lgamma difference cancels at n ~ 10⁶, and while the argmax is unaffected, absolute ℓ, AIC/BIC and goldens are not. The in-house form also avoids two portability traps: glibc's `lgamma` writes the global `signgam` (not thread-safe under OpenMP), and host vs CUDA `lgamma` differ, which D-011 cares about.
- **D-071 — `log_add_exp(a, b)` is tested with −∞ inputs (either, and both → −∞, never NaN), and with equal arguments (→ a + log 2).** Zero-probability quadrature terms produce −∞ routinely, and `−∞ − (−∞)` is NaN.

## M1.2 implementation, 2026-09-26

- **D-072 — Where mpmath can't compute a reference directly, the generator uses a documented identity and cross-checks it where both are usable:** `erfcx` by its asymptotic series for x ≥ 100, and logΦ via `erfcx` for x ≤ −200. The self-checks agree to 1e-45. mpmath's `erfc` overflows near x = 1e300 and its `ncdf` near x = −1e155.
- **D-073 — Golden tables store exact hex floats, with decimal columns only for readers. Generated C++ uses hex literals (C++17).** No decimal round-trip ambiguity across `strtod` implementations or compilers.
- **D-074 — Tolerances are set at twice the worst error observed on the golden table, rounded up, and stated in ulp (relative). Each has a register entry with value, observed figure, rationale and enforcing test; a CTest check keeps `tests/tolerances.hpp` and `docs/methodology/tolerances.md` in sync. Raising one means editing the register with a reason.** Headroom covers libm differences across platforms; ulp can't hide relative-accuracy loss the way an absolute bound can.
- **D-075 (amends D-033) — The provenance manifest has no generation date; git history dates the files. `--check` compares everything except `python_version` and `mpmath_backend`, which are recorded but don't affect results.** A date would make exact regeneration impossible.
- **D-076 — `core/special` is binary64-only for now. Float variants for `PrecisionMixed` (D-039) come with M4, when their tolerances can be measured on a GPU.** Nothing uses them before then.
- **D-077 — `probit` is plain AS241 (6 ulp worst observed, about 7e-16 relative), with no Newton refinement.** One step through `log_phi` would reach about 2 ulp, but the textbook algorithm is what someone replicating it would implement. Revisit if a consumer needs it.
- **D-078 — `erfcx` on x ≥ 0 uses four mpmath-fitted Chebyshev pieces ([0,1], [1,3], [3,8], and a tail in 8/x). x² is split so its large part is exact inside `exp`/`log`.** 2 ulp observed; the split removes an error that would otherwise grow as ~x² ulp.
- **D-079 — A `generated-files` CI job reruns the generator with `--check` at the pinned mpmath, in a venv on ubuntu-24.04.** It proves every generated header and golden table can be regenerated exactly.

## M1.3 quadrature, 2026-09-26

- **D-080 (replaces the runtime Golub–Welsch plan in ARCHITECTURE §3.4) — Gauss–Hermite nodes and log-weights are mpmath-generated tables for N ∈ {8, 16, 24, 32, 48, 64, 96, 128}, through `tools/gen_quadrature_tables.py` and its manifest.** Golub–Welsch in double gives outer weights with good absolute but poor relative accuracy, the wrong error model for a log-space engine. Tables also keep an eigensolver out of the code a validator must trust. The N set covers the 2N convergence check up to 64. Parity uses one fixed N, chosen in M1.8.
- **D-081 — The generator checks the *converted* standard-normal rule (z = √2·x, w = W/√π, converted in one function) at 50 digits before writing anything:**
  - node and weight symmetry;
  - every moment up to degree 2N−1 (E[Z²ᵏ] = (2k−1)!!, odd moments 0);
  - each weight against the closed form 2^(N−1) N! / (N² H_(N−1)(xᵢ)²).

  In C++, the same moments are checked in double, relative, with a degree-proportional tolerance and terms summed smallest first. No Golub–Welsch cross-check is kept: in double it could only confirm the outer weights to an absolute tolerance, which is the weakness being avoided. mpmath 1.4.1's own weights were measured at ≤ 4e-48 relative even for the 1.8e-102 outer weight at N = 128.
- **D-082 — `tools/tablegen.py` holds the shared generator machinery and the mpmath reference functions. Each manifest records the sha256 of every generator source file (`generator_files`). The refactor left every existing special-function output byte-identical.** One implementation of provenance, and a manifest that can't miss a shared dependency.
- **D-083 — Mixture references: log E[pᵈ(1−p)ⁿ⁻ᵈ] at n ≤ 2 by closed forms (E[p] = PD; E[p²] = Φ₂(c, c; ρ) by Plackett's identity, an integral over ρ unrelated to the z-integral). For n up to 10⁶, by tanh-sinh quadrature at 50 digits over breakpoints around the mode, with estimated error ≤ 1e-35. The generator cross-checks Plackett against the generic integral.** The references share no method with the code under test.
- **D-084 — The binomial-mixture log-integrand and its adaptive hint are in `core/model/binomial_mixture.hpp`, landed with M1.3 because the integrator can't be tested on the real integrand without them. The objective (M1.4) reuses them.** h is strictly concave, so bracketed Newton with a bisection fallback always finds its unique mode, deterministically.
- **D-085 — `inverse_mills(x)` = φ(x)/Φ(x) joins `core/special`. For x < 0 it uses the exact identity √(2/π)/erfcx(−x/√2), which never cancels. 4 ulp observed against mpmath.** The hint needs slopes deep in the tail, where φ/Φ computed naively cancels.
- **D-086 — Adaptive-GH accuracy is established by regime.**
  - N ≥ 64 and ρ ≤ 0.5: ≤ 4 ε observed across n = 1 to 10⁶ (tolerance 8).
  - ρ = 0.9 with n ≤ 2: slow convergence (252 ε ≈ 5.6e-14 at N = 128; 7.9e-9 at N = 64), a documented shortfall with its own tolerance.
  - Fixed-node GH is wrong by ~12 in log-likelihood at n = 10⁶, even at N = 128 (D-037).

  This is the input to the parity-N choice (M1.8). If the ρ grid must reach ~0.9 with small-n periods, parity needs N = 128 or accepts ~1e-8.
- **D-087 (amends D-070) — The worst intermediate in `lbinom`'s exact path (n ≤ 60) is 3.55×10¹⁸, at n = 60, k = 30, about 19% of 2⁶⁴, by exhaustive enumeration.** Recorded in the code comment.

## M1.3b: the 256-node rule and the ρ range, 2026-09-26

- **D-088 (owner; amends D-086) — Parity uses N = 128 everywhere, with a 256-node table (generated and checked like the others; smallest weight 3e-211) for the twice-the-nodes convergence check. `native` keeps N = 64 available.** N only affects building the per-period surfaces L once; every resampling scheme reuses L through W. **Measured correction:** D-086's "N ≥ 64 is full precision for ρ ≤ 0.5" was wrong. Zero-default periods at ρ = 0.5 (not in the golden set then) need N = 128 (4.4e-13). The claims now stand as:
  - regime A (ρ ≤ 0.24, N ≥ 64): full precision;
  - regime B (ρ ≤ 0.5, N ≥ 128): ≤ 4.4e-13;
  - beyond ρ = 0.5: zero-default and all-default periods are not at full precision at any tabulated N (Q15).
- **D-089 (owner) — The ρ upper bound is a configurable grid parameter, default 0.5. Any estimate or bootstrap replicate on the grid edge sets the edge flag.** Bounds are a methodology choice, and a user's own methodology may set its own. **With D-088, a bound above 0.5 carries a precision caveat until Q15 is resolved.** The twice-the-nodes check (128 vs 256) detects it per run: the two disagree on exactly the periods that aren't converged.
- **D-090 — Golden cases now include ρ ∈ {0.95, 0.99} closed forms and low-default spikes (d = 0 with n ∈ {10, 50, 200} at ρ ∈ {0.5, 0.9, 0.95}; d = 5 and d = n at ρ = 0.9).** A precision claim has to cover the range a user can configure. The low-default cases are also the M6 Pluto–Tasche regime.
- **D-091 — Tried and rejected: splitting at the mode with half-range Gauss–Hermite (generalised Gauss–Laguerre, α = −½) per side.** Prototyped in mpmath: errors fall only as a power of N, because the half-range integrand isn't smooth in the Laguerre variable. It is worse than plain adaptive GH on every golden case. The prototype also reproduced the C++ adaptive-GH errors independently.

## Q15 and O-1 resolved, 2026-09-26

- **D-092 (Q15, option c; owner) — (a) now: parity claims numerical precision for ρ ≤ 0.5, as measured. Above that, the tolerance register makes a *separate* statistical-materiality statement: the error in log-likelihood per period against the 1.92 LR threshold (worst golden period 4.9e-4, ~4,000× below). Each run reports it through the 128-vs-256 check, which flags any period whose difference exceeds 1e-10 and reports the per-period maximum and the total.** The cap is then a precision statement, not a correctness problem. A test shows the check never misses a golden period whose true error exceeds the threshold, and never reports less than 0.43× the true error (tolerance 0.2).
- **D-093 (Q15 b; owner) — The dedicated rule for d = 0 and d = n periods is composite Gauss–Legendre: breakpoints at the sigmoid's transition point (where the conditional PD crosses 0.5) plus or minus a few transition widths √(1−ρ)/√ρ, fixed panels over a truncated range. It is a `native` enhancement, promoted to parity only if it passes the same golden checks.** Deterministic, fixed cost, easy to describe in prose and to reproduce in scipy. Tanh-sinh would also be accurate but is harder to explain.
- **D-094 (O-1; supersedes the provisional parts of D-026, D-027, D-065) — The supported CUDA floor is nvcc 11.8 with the MSVC v143 14.39 toolset, the newest MSVC that nvcc 11.8 accepts: the project supports CUDA ≥ 11.8 and older host compilers for broad compatibility. With the floor of 11.8 confirmed, the C++17 scope for device-visible code (D-058, D-062) is required, not provisional. sm_120 is an optional target for development hardware only. The development machine builds with the floor, 11.8 + MSVC 14.39, as primary.** The CMake floor is still open, so D-047 stays provisional.

## M1.4 engine, 2026-09-26

- **D-095 (amends D-038) — Refinement fits the exact 2-D quadratic through the 3×3 stencil around the grid argmax, with the PD–ρ cross term, in the axes' scaled coordinates. The estimate is the vertex −H⁻¹g; SEs come from (−H)⁻¹, converted to natural scale by the delta method with dv/du at the estimate.**
  - H not negative definite: grid point, flat flag, SEs NaN.
  - Vertex outside the stencil: grid point, refinement-rejected flag, SEs kept, since the curvature is valid.
  - Argmax on an axis end: grid point, edge flag, SEs NaN.
  - Nothing is ever extrapolated.

  PD and ρ estimates are strongly correlated (0.30 on the test panel; far more along the ρ ridge at short T), so per-axis parabolas would bias the vertex. The SEs need the cross term anyway, to report `corr`.
- **D-096 — Engine structure.**
  - A Backend exposes `parallel_for(count, f)`; the engine includes no backend (layering, D-050).
  - `weighted_sum` (ascending t, Neumaier-compensated, zero weights skipped) is the single definition every reducer and refinement path uses, so all of them see identical values.
  - The reduction runs over fixed 1024-point k-tiles merged in tile order, whatever the thread count.
  - Panel validity is objective-specific (`Objective::panel_error`).
  - `calibrate` recomputes each l_t at the estimate with both the primary and the check integrator. It reports the maximum and the total of their differences, counts the periods above 1e-10, and sets the quadrature flag if any are.
- **D-097 — Estimate flags:** `grid_edge`, `flat_surface`, `quadrature_unconverged`, `refinement_rejected`, `numeric` (a NaN anywhere in the surface). Each is set, never used to repair a result silently. They map onto the ABI's `VCAL_FLAG_*` in M2.
- **D-098 — OpenMP is found by CMake and optional (`VCAL_ENABLE_OPENMP`, default ON). Without it the CPU backend runs serially with identical results. CI installs `libomp` for Clang so the threaded path is tested on every compiler. MSVC's OpenMP 2.0 accepts the 64-bit signed loop index.** Determinism never depends on threading, so a missing runtime costs speed, not correctness.
- **D-099 — Cross-platform scope of determinism, as measured in CI prior to public release.** GCC 11, GCC 14 and Clang 18 (Debug and Release) give bitwise-identical estimates. MSVC differs from them by ~3e-15 relative (PD, ρ) and ~4e-16 (log-likelihood), because glibc and the UCRT round `log`/`exp` differently in the last ulp. So bitwise identity is claimed per platform (§6.1), cross-platform only for the DGP (D-052), and golden estimates are compared with a documented tolerance, never bitwise.

## K-1 retired; owner follow-ups, 2026-09-26

- **D-100 (fulfils D-059; retires K-1) — nvcc compile-only CI on Linux containers: `nvidia/cuda:11.8.0-devel-ubuntu22.04` and `12.8.0-devel-ubuntu22.04`, both with GCC 11.4. Each job:**
  - compiles the C++17 header-check TU with nvcc (D-064);
  - compiles `tests/cuda/device_compile_check.cu`, a kernel that calls every device-visible function (special functions, both integrators on the binomial-mixture integrand, objective, grid, `weighted_sum`, `ArgMax`), so device code must actually be generated;
  - treats a host-only callee and any other nvcc warning as an error (`--Werror=cross-execution-space-call,all-warnings`), with a compile-fail fixture proving that check fires;
  - also runs the CPU suite.

  An evidence step prints the applied flags and `cuobjdump` listings. The CI run verified prior to public release shows cubins for sm_89 (11.8 and 12.8) and sm_120 (12.8), plus compute_80 PTX, all compiled clean. `std::log`, `std::log1p`, `std::exp`, `std::trunc`, `std::sqrt`, `std::isnan` and `std::isfinite` compile in device code on both toolkits. Headers alone would prove little, because `VCAL_HD` inline functions only get device code when device code calls them.
- **D-101 — A Windows nvcc + MSVC 14.39 CI job is deferred.** GitHub's Windows runners don't ship that toolset, and installing it per run is slow. Add it (possibly with a cached toolset) if Windows builds at the 11.8 + 14.39 floor need CI coverage; until then the development machine covers that combination (D-066, D-094).
- **D-102 (owner) — Golden estimates are compared with a relative tolerance of 1e-12: far above the ~3e-15 cross-platform libm spread (D-099), and far below statistical significance.** Enters the tolerance register with its first golden-estimate test (M1.7).
- **D-103 (owner; confirms D-095) — 2-D stencil refinement is kept.**

## Publication hygiene, 2026-09-26

- **D-104 (owner) — Committed files describe the build environment only in neutral terms: toolchain floors are portability goals, and the machine the owner develops on is "the development machine". Machine-specific facts (versions, drivers, packaging, hardware) live in gitignored `*.local.md` notes. To make that true of the current files, D-026, D-059, D-060, D-066, D-094 and D-101 were reworded in place, the one exception to this log's append-only rule (extended by D-146). The public repository is published from a single squashed snapshot.** The project is independent of any organisation.

## Build compatibility, 2026-09-26

- **D-105 (confirms D-047) — CMake compatibility.**
  - `cmake_minimum_required` stays at or below 3.26 (currently 3.25), and `CMakePresets.json` at or below schema 6, so the build configures with the older CMake releases that long-term-support Linux distributions still ship (3.26 among them).
  - CI proves both ends: the nvcc 11.8 job pins CMake 3.25.x (the declared floor), the nvcc 12.8 job pins 3.26.x, and the runner-image jobs report their CMake version, which covers 4.x.
  - CUDA architectures are explicit numbers, optionally `-real`/`-virtual` (currently `89-real;80-virtual`, plus `120-real` on nvcc ≥ 12.8). `all`, `all-major` and `native` are rejected at configure time, both in `VCAL_CUDA_ARCHITECTURES` and in a cached `CMAKE_CUDA_ARCHITECTURES`.

  A build must not depend on which GPU, if any, the build machine has, and must configure with the oldest CMake the project supports.
- **D-106 — CUDA builds use the Ninja generator on every platform, including Windows. Configuring CUDA with a Visual Studio generator is rejected at configure time with a pointer to the `cuda-release` preset, which names Ninja explicitly.** The CUDA toolkit's Visual Studio (MSBuild) integration is tied to particular Visual Studio releases and can be missing for the one installed. Ninja plus a developer shell with the required MSVC toolset (`vcvars_ver=14.39` for nvcc 11.8, D-066) does not depend on that integration.

## Comparison tolerances and M1.6 design, 2026-09-26

- **D-107 (owner; scopes D-102) — Three different comparisons, three separately named tolerances in the register, never interchangeable:**
  - **core vs ref, per-period l_t:** tight, in ε units, from the two measured error bounds (core's regime tolerances, D-086/D-088; ref's `TOL_REF_LOG_MIXTURE_EPS`). Zero-default and all-default periods with ρ > 0.5 are excluded, because there core is not at full precision (D-092).
  - **core vs ref, estimates:** in SE units (around 0.05 SE), justified by refinement accuracy (0.010 SE observed, D-095). The two use different optimisers (grid + 2-D stencil vs golden-section), so agreement is bounded by refinement and optimiser precision, not rounding.
  - **core across platforms, estimates:** 1e-12 relative (D-102), covering only last-digit libm differences between builds of the *same* engine (D-099).

  Naming them separately stops the 1e-12 cross-platform tolerance from being applied to core vs ref, which it was never meant for.
- **D-108 (owner) — ref is the oracle for the composite Gauss–Legendre enhancement (D-093).** Its `log_mixture` is at 7.5 ε on ρ = 0.99 zero-default periods, exactly the regime that enhancement targets.
- **D-109 (owner) — DGP design points for M1.6** (with D-052–D-057, D-061):
  - **Stream addressing:** Philox4x32-10 is keyed by the seed. The counter is (scenario, replicate, period, draw index), so adding replicates or scenarios never changes an existing panel, and goldens stay stable as the harness grows.
  - **Floating-point environment:** besides `/fp:strict` and `-ffp-contract=off`, nothing may enable flush-to-zero or denormals-are-zero. The DGP checks at start-up that subnormals survive, and refuses to generate otherwise, since some libraries flip those modes globally.
  - **Cross-OS check:** CI computes the hash of a reference panel on every platform and compares it with one checked-in value, the direct test of Q14's goal.
  - **Known answers:** Random123's published Philox4x32-10 vectors. The in-house `log`, `exp`, `erfc`, Φ and AS241 are tested against mpmath tables, like everything else.

## M1.6 DGP, 2026-09-26

- **D-110 (amends D-053/D-061) — The DGP computes Φ without a rational `erfc`. It uses the series ½ + φ(x)·Σ x^(2n+1)/(2n+1)!! for abs(x) ≤ 2, and φ(x)·R(abs(x)) with the Mills ratio R from its Lentz continued fraction beyond, with φ from the in-house `exp`. Its function set is therefore `log`, `exp`, `sqrt`, Φ and AS241 Φ⁻¹.** There are no coefficient tables to transcribe or generate, and the Python mirror is simpler. Accuracy is 496 ulp worst (≈ 1.1e-13), inside the 1e-12 target.
- **D-111 (owner's M1.6 notes, implemented) — The DGP contract is exact.**
  - **Blocks:** the factor draw uses Philox blocks 0 … 2¹⁶−1 (generation is refused if they are exhausted); Bernoulli draws start at block 2¹⁶, two obligors per block, so rejections can never shift them.
  - **Uniforms:** (k + ½)·2⁻⁵² from the top 26 bits of two words, exact and strictly inside (0, 1).
  - **The Python mirror** (`tools/gen_dgp_tables.py`) uses only exact or correctly rounded IEEE operations and line-for-line ports of the in-house maths, never `math.log`/`exp`/`erfc`. It reproduces Random123's known answers and produces the reference panels and their SHA-256.
  - **Bitwise checks:** C++ must match the mirror bit for bit, function by function and panel by panel, and match the checked-in hash on every CI platform.
  - **Known answers:** Random123's vectors are copied verbatim with provenance (commit, file hash) and a NOTICE attribution.
  - **Build:** the `dgp` library builds with `FP_STRICT` (MSVC `/fp:strict`; `-fno-fast-math` elsewhere), and the DGP refuses to run under flush-to-zero or denormals-are-zero.

  Full description: `docs/methodology/dgp.md`.
- **D-112 — Statistical tests take their SEs from batch means, not fourth-moment formulas.** For the right-skewed default rate, the fourth-moment SE of the variance came out 1.7× too small (measured over ten seeds), which would make a fixed-seed test fragile.

## M1.6 follow-ups and M1.7 design, 2026-09-26

- **D-113 (owner; done, amends D-110) — The DGP's Φ uses fixed iteration counts, chosen for the worst case: 128 Lentz iterations for the continued fraction (the worst measured is 106, at t = 2; 52 at t = 3; 6 at t = 38) and 32 series terms (the worst measured is 23, at abs(x) = 2).** A fixed count is simpler to state in `dgp.md` and makes the GPU port (R-1) trivially identical: no data-dependent loop exits. Results may move by ulps, so the reference panels and their hash were regenerated.
  **Reference-panel hash record:** before `367c26d845c07fe0b82cee31b3a9b509e30794e927ad9dad9cf7af52d1b5887c`; after `367c26d845c07fe0b82cee31b3a9b509e30794e927ad9dad9cf7af52d1b5887c`, **unchanged**.
  - Φ values did move: up to a few ulps in the tails (for example Φ(−8) and Φ(−37.5) changed in the last hex digits of `tests/golden/dgp/mirror_values.csv`), and the worst Φ error against mpmath went from 496 to 488 ulp.
  - The factor draws use only `log`, so they are unaffected.
  - A shift of a few ulps in p_t flips a draw U < p_t only if some uniform lands within those ulps. None did, in the reference panels or in the 40,000-panel moment test, whose mean and variance came out identical to every printed digit.
  - So a panel generated before D-113 still matches bit for bit, and any future hash change needs its own entry here.
- **D-114 (owner; done) — Every CI job prints the reference-panel hash, including the nvcc containers, even though the test already asserts it.** An auditor can then point at each platform's log without re-running anything.
- **D-115 (owner) — M1.7 design:**
  - **Identical domains:** core's grid bounds and ref's optimiser bounds are the same box, including the ρ cap, so estimates near a bound never differ for reasons unrelated to correctness.
  - **Flagged fits are compared separately:** estimates with edge, flat or refinement-rejected flags have no meaningful curvature SE, so they are excluded from the SE-unit comparison. Instead, both engines must put them on the same bound, or both must flag them.
  - **ref's cost is budgeted:** about 8 s per fit, so ref estimates run on a fixed, seeded subset of panels, with the selection rule recorded. The per-period l_t comparison (about 10⁴ points) is cheap and stays comprehensive.

## M1.7 cross-reference, 2026-09-26

- **D-116 (restates D-086/D-088/D-092; found by the M1.7 cross-reference) — Core's precision claim covers periods with 0 < d < n only. Zero-default and all-default periods carry no precision claim at any ρ.**
  - **What was found:** against ref, confirmed by mpmath at 50 digits, core (N = 128) misses zero-default periods where few defaults are expected (tiny PD, large n), even at low ρ: 9.2e-9 absolute at PD 1e-6, ρ 0.24, n 10⁶; 1.9e-5 at PD 1e-5, ρ 0.5, n 10⁶. ref is exact there (≤ 1e-15).
  - **Why the old claims held:** the golden set had no such cases. They are now added (five periods).
  - **What still holds:** the 128-vs-256 check flags every period whose error exceeds 1e-10, on the golden set and on 3,580 cross-reference points, and never reports less than 0.33× the true error.
  - **Also found:** the worst zero-default error at ρ = 0.9 is 6.0e-3 per period, not the 4.9e-4 previously stated. The materiality statement is updated; see Q16.
  - **A metric correction:** l_t = log C + log I can cancel heavily, so the cross-reference error scale is the size of the terms, max(1, abs(log C) + abs(log I)), not the result. The apparent regime-C failure (35,000 ε at n = 10⁶, d = 350,000) is that cancellation, identical in both engines and confirmed by mpmath, not an engine error.
- **D-117 — M1.7 results.**
  - **Per-period, 0 < d < n:** within 35 ε (term-scaled) over 1,704 points up to ρ = 0.9.
  - **Estimates:** core vs ref within 0.024 SE (PD) and 0.022 SE (ρ) on 6 unflagged DGP panels with an identical box. The no-defaults panel lands on the same bound in both engines.
  - **SEs:** agree within 10.4%; see Q17.
  - **Regression values:** core's estimates, written by MSVC 19.51, are compared at 1e-12 on every CI platform. Measured (CI): 0 on Windows. On Linux (GCC 11/14, Clang 18) the worst is **3.9e-13, in `se_pd` of scenario 5**. Point estimates and log-likelihoods agree far more tightly (about 3e-15, D-099). SEs come from second differences of the surface, which amplify last-ulp libm differences. 1e-12 holds, with a 2.5× margin.
  - **Cost:** ref estimates run in Release builds only (CTest label `slow`, excluded by the Debug preset), about 40 s.
  - **Panel set:** 6 fixed scenarios, replicate 0, T = 20, n = 1000, DGP seed `0x4D31375852454631`, written in `tests/xref/xref_test.cpp`.

## M1.7b: Q16 and Q17 resolved, 2026-09-26

- **D-118 (Q16; owner: yes, as M1.7b; amends D-093, restates D-116) — Parity integrates zero- and all-default periods (d ∈ {0, n}) with composite Gauss–Legendre, K = 16 panels × M = 16 points, and every other period with adaptive Gauss–Hermite, N = 128. The per-run check (D-092) is the same rule doubled: GH 256 and K = 32.** Together the two rules are one `SplitRule`, selected per period by the objective's hint. The rule:
  - **Domain:** where h(z) = g(z) − z²/2 has fallen by 50 from its value at the mode, on each side. The crossing is found by doubling a step, then exactly 64 bisections, so it is the root to rounding and a smooth function of the parameters.
  - **Placement (fixed relative to the transition point, owner's condition):** the map z = z_c + w·sinh(u), with K equal panels in u. z_c is where the survival factor (1−p)ⁿ (or pⁿ) equals one half, i.e. p = 1 − 2^(−1/n) (mirrored for d = n). w is that factor's width there, 1/(n·β·λ), with λ the inverse Mills ratio and β = √ρ/√(1−ρ). This differs from D-093's wording (conditional PD = 0.5, breakpoints a few widths √(1−ρ)/√ρ either side): for large n the cliff sits where p ≈ ln 2/n, not where p = 0.5, and one smooth map replaces hand-placed breakpoints.
  - **Fixed panel count (owner's condition):** K and M are the same at every parameter point, with no adaptive subdivision, so the surface has no panel count to jump.
  - **Why K = 16, M = 16:** at fixed K·M, more panels beat more points per panel. 16 × 16 is the smallest configuration under 12 ε in every golden cell (ρ ≤ 0.5 / ρ > 0.5: 11.7 / 7.6 ε), and every larger one stays at the rounding floor (table in `tolerances.md`).
  - **Measured:**
    - **Golden set:** every case within 11.7 ε of mpmath (interior 2.15 ε), including ρ = 0.99.
    - **Cross-reference against ref:** 4,296 points, term-scaled. 35.2 ε interior; 23.5 ε for d ∈ {0, n} at ρ ≤ 0.5; 8.5 ε at ρ > 0.5.
    - **Rule vs doubled rule:** worst difference 1.2e-10, on a period whose terms are about 6.5e5 in size. The engine would flag 0 periods (D-120).
    - **Smoothness (owner's condition):** zero-default surfaces on a 21 × 21 logit grid; second differences of core minus ref at 3.7e-14 of the largest.
  - **Effect on D-116:** the exclusion no longer applies. Parity's precision claim covers every period tested, and the separate materiality statement for zero- and all-default periods is withdrawn (`tolerances.md`, "The parity rule: what is claimed").
  - **Device code:** `sinh`, `asinh`, `cosh` and `expm1` are exercised in the nvcc compile check.
- **D-119 (Q17; owner; amends D-038 and D-095 on where SEs come from) — Core's SEs are the observed information at the refined estimate: the Hessian of the objective itself, by central differences in each axis's scaled coordinate, then the delta method to natural scale. The point estimate still comes from the grid stencil (D-095).**
  - **Step:** `kHessianStepFraction` = 0.15 × the grid-stencil SE on that axis (owner: a fixed fraction in 0.1–0.2), so the step scales with the surface's own curvature. The fraction is recorded in code and here. A sensitivity check compares 0.10 and 0.20 against 0.15: SEs move by at most 4.7e-5 relative (`TOL_SE_STEP_SENSITIVITY_REL` = 1e-4).
  - **Headline fits only:** 8 extra evaluations of the full log-likelihood. Resampling (W × L) never calls it, so the surface reuse of D-009 is unaffected.
  - **Near a bound:** an estimate within `kNearBoundSe` = 2 SEs of an axis bound (in scaled coordinates) sets the new flag `kFlagNearBound` (bit 5). Wald SEs are unreliable there in any engine. Profile-likelihood intervals are the right tool, and are scheduled for M2. SEs are still reported, but flagged.
  - **If the Hessian is not negative definite:** SEs are NaN and `kFlagFlatSurface` is set, as for the stencil.
  - **Measured:**
    - **Against the old stencil SEs:** −0.7% (PD) and −1.3% (ρ) on the M1.4 panel.
    - **Coarse- vs fine-grid SEs:** now 0.46% apart (was 0.97%).
    - **Against ref:** within 1.56% on the 4 panels not flagged near-bound (was 10.4% on 6). The 2 near-bound panels are within 4.2%; they are reported, not compared.
- **D-107 amended (owner) — a fourth separately named comparison tolerance: core's SEs across platforms against its regression values, `TOL_XREF_CROSS_PLATFORM_SE_REL` = 1e-10.** Second differences at step ≈ 0.15 SE amplify last-ulp libm differences in l_t. Estimates and log-likelihood stay at 1e-12. The regression values (`tests/golden/xref/core_estimates.csv`) were regenerated with MSVC 19.51, because SEs and some estimates changed with D-118 and D-119. Measured in CI prior to public release: 0 on Windows. On Linux (GCC 11, GCC 14 and Clang 18, identical), estimates and log-likelihood differ by up to 1.9e-14 (53× inside 1e-12) and SEs by up to 7.6e-12 (13× inside 1e-10). That confirms SEs needed their own tolerance.
- **D-120 — The per-run check flags a period when abs(l_t(rule) − l_t(doubled)) > max(1e-10, 64·ε·(abs(log C) + abs(log I))).** The cross-reference found periods at n = 10⁶ whose rule-vs-doubled difference is about 1.2e-10 purely from rounding in terms of size 6.5e5. A fixed 1e-10 would falsely flag them. Below the term-scaled floor no integration error can be detected anyway. `Objective::rounding_scale` supplies the term size.
  - **Detection:** with parity, no golden period exceeds the threshold, so the contract test also runs the same family two doublings coarser (GH 32 and K = 4, checked against GH 64 and K = 8). 42 periods are flagged there, and the check never reports less than 0.99998 of the true error. `TOL_QUAD_CHECK_MIN_DETECTION_RATIO` is tightened from 0.2 to 0.5 accordingly.

## M1.8 recovery harness, 2026-09-26

- **D-121 (owner's M1.8 notes) — The recovery harness reports three numbers per scenario and parameter, and does not let excluded replicates flatter coverage.**
  - **The three numbers:**
    1. bias (with its Monte Carlo SE) and RMSE over **all** replicates, grid-edge estimates included at their grid value;
    2. the fraction of replicates **flagged**, meaning they have no reliable Wald interval: on the grid edge, within 2 SEs of a bound (`kFlagNearBound`), or a flat surface;
    3. coverage of the 95% Wald interval **among unflagged replicates**, always labelled conditional. Dropping flagged replicates selects on the outcome, so this number alone overstates how well the intervals work.
  - **Verdict:**
    - **DEFERRED** when the flagged fraction is at least 5% (`TOL_RECOVERY_MAX_FLAGGED_FRACTION`, owner's figure). The coverage verdict then waits for M2's profile-likelihood intervals, the right tool near a bound. It is neither passed nor failed now.
    - Otherwise **PASS** if the conditional coverage lies in the Monte Carlo band 0.95 ± z·√(0.95·0.05/m), where m is the unflagged count, and **FAIL** if not. *(Amended by D-124: an out-of-band verdict is a reviewed KNOWN FINDING, not an engine failure.)*
    - The band uses z = 3.29, the two-sided 0.1% point (`TOL_RECOVERY_COVERAGE_BAND_Z`). Over 162 scenario-parameter verdicts that allows 0.16 expected false FAILs.
  - **The interval** is symmetric in the axis's logit coordinate, where the Hessian is taken (D-119): û ± z₀.₉₇₅·se_u. That keeps it inside (0, 1), and anyone can reproduce it from the reported SE and dv/du.
  - **"Flagged" is joint:** `kFlagNearBound` does not say which axis is near its bound, so a ρ̂ near the 0.5 upper bound also defers the PD verdict. That is conservative; per-axis reporting can come with M2's intervals.
  - **Other exit criteria:**
    - RMSE falls strictly with T for every (PD, ρ, n) and both parameters;
    - no replicate is flagged for quadrature or numerics.
    - The known small-T downward bias of ρ̂ is reported, not tested against zero.
  - **Design:**
    - Scenario matrix PD {0.1%, 1%, 5%} × ρ {0.02, 0.12, 0.24} × T {20, 40, 100} × n {100, 1000, 10000}, with n constant over time. Scenario id = ((i_PD·3 + i_ρ)·3 + i_T)·3 + i_n, which is also the DGP scenario number.
    - R = 1000 replicates. DGP seed `0x4D31385245434F56`.
    - The M1.7 estimation box: PD [1e-4, 0.2] × 61, ρ [1e-3, 0.5] × 41, logit axes.
    - Parity rule, D-119 SEs.
  - **Where it lives:**
    - The full run is a tool, `recovery_harness --write`, taking about 47 minutes on 24 threads. It writes `tests/golden/recovery/` (summary, replay values, manifest) and the generated `docs/methodology/recovery_results.md`.
    - CI replays replicates 0 and 1 of every scenario, a D-107 comparison of its own with `TOL_RECOVERY_REPLAY_REL` = 1e-12 and `TOL_RECOVERY_REPLAY_SE_REL` = 5e-10. The first CI run reused M1.7's SE tolerance of 1e-10. Linux measured 2.5e-10 on scenario 74 (n = 10⁴, where last-digit differences in terms of about 10⁵ are amplified by the Hessian's step), and `main` went red until this was fixed. M1.7's 1e-10 is unchanged.
    - CI also enforces the exit criteria on the committed summary, recomputing each verdict from its recorded counts under the current tolerances.
- **D-122 — `evaluate_surface` evaluates each distinct observation once per grid point and copies the row to periods that repeat it.** l_t depends only on the observation and θ, so the copy is bitwise what a second evaluation would give. A test compares against evaluating every period, for 1, 3 and all threads, and the M1.7 regression values are unchanged bit for bit.
  - **Why:** low-default panels repeat heavily. A T = 100 panel with n = 100 at PD 0.1% has 3 distinct d, and its zero-default periods are the costliest to integrate (D-118). Measured: 2.0–5.2 s → 0.12–1.9 s per serial T = 100 fit, which makes R = 1000 affordable.
  - **Contract:** objectives must provide `Obs ==`; `check_objective` asserts it with its own message. The duplicate search is O(T²) comparisons, negligible next to T × K integrals.
- **D-123 (owner) — Parity sign-off for the binomial-mixture MLE, the M1 exit, requires:**
  1. the methodology doc, the scipy script and the tolerance register complete for the binomial-mixture MLE (M1.9);
  2. the recovery tables committed as goldens with provenance (M1.8);
  3. the deferred coverage verdicts listed explicitly, each with its M2 dependency (profile-likelihood intervals);
  4. a tag, `v0.1.0-parity-binomial`, on a commit whose CI is green.

  The sign-off lists both the KNOWN FINDINGs and the DEFERRED verdicts (D-124), each re-assessed with M2's profile-likelihood intervals.

  **Process rule (owner):** `main`'s head is always a green commit. CI must pass on every pushed commit, docs-only included, before new work starts.
- **D-124 (owner; amends D-121) — The harness judges two kinds of criterion differently. Engine correctness must pass. Interval coverage is a property of the interval method, reported as PASS, DEFERRED or KNOWN FINDING, never as an engine FAIL.**
  - **Engine correctness, which must hold:**
    - RMSE falls with T;
    - no quadrature or numeric flags;
    - the replay agrees across platforms;
    - estimates and SEs agree with ref (M1.7).
  - **Coverage:** an out-of-band verdict is a KNOWN FINDING only when it is reviewed and listed, with its diagnosis, in `kKnownFindings`. Otherwise it is UNREVIEWED, which CI rejects, as it rejects a listed finding that returns inside the band. The set of findings is therefore pinned: any change needs a review.
  - **Why (the first R = 1000 run):** 19 of 162 verdicts fell outside the band, while every engine-correctness criterion held.
    - **14 small-T Wald undercoverage, at T = 20 and 40** (12 of PD, 2 of ρ). Coverage 0.910–0.927. The spread of û exceeds the root-mean-square Hessian SE by 1–10%. A diagnostic re-fit of nine scenarios gave 0.924–0.943 with a t(T−1) quantile and 0.93–0.955 with the actual spread.
    - **5 ρ overcoverage, 0.973–0.998.** SE ratio 0.89–1.00, with û skewed towards the lower bound in the logit coordinate (skewness −0.43 to −0.67, measured by a re-fit of all five). I first attributed these to selection by the flagged-replicate exclusion. That holds only partly: two of them have under 1% flagged.
  - **Evidence committed:** per scenario and parameter, the unflagged sd(û) and the root-mean-square SE in the logit coordinate, and their ratio (`summary.csv`, `recovery_results.md`).
  - **Platforms:** verdicts are pinned exactly on the platform that wrote the goldens (the `MANIFEST.json` compiler). `recovery_harness --check` elsewhere allows counts within 2 replicates (`TOL_RECOVERY_CROSS_PLATFORM_COUNT`), because a replicate within about 1e-12 of a threshold can fall either side under libm differences. CI replays 162 fits on every platform with flags compared exactly; the full check is too slow for CI.
  - **Not done, by decision:** switching to a t(T−1) quantile now. That would change the interval after seeing it fail. It is evaluated in M2 as a `native` option alongside profile likelihood.
  - **Harness fix:** the first `--write` silently wrote nothing, because it never created `tests/golden/recovery/`. Writers now create their directory and throw on any stream error. After writing, every output is checked to exist and be non-empty, and `recovery_goldens_meet_the_exit_criteria` checks the same in CI.
  - **D-122 accepted by the owner.**

## M1.9 methodology and scipy replication, 2026-09-26

- **D-125 (owner) — Tolerance values live in one neutral data file, `tests/tolerances.toml`. `tools/gen_tolerances.py` (standard library only) renders `tests/tolerances.hpp` from it, and CI checks the header is current.** Validation scripts read the TOML directly.
  - **Why:** a replication script should never have to parse the C++ header. One source of truth also keeps the C++ tests and the scripts from drifting apart.
  - **Unchanged:** the register page (`tolerances.md`) still gives each id's rationale, and `tolerance_register_sync` still checks that the page and the header list the same ids.
  - **Conversion:** the header was converted mechanically, comments included. The regenerated constants are textually identical.
- **D-126 (owner) — `validation/scipy/binomial_mixture_mle.py` replicates the method with ordinary tools, and deliberately does not port the engine's numerics.**
  - **Tools:**
    - `binom.logpmf` in log space;
    - `quad` over a range centred on the integrand's mode, found by Brent's method, extending until the integrand has dropped by e⁻⁵⁰; relative accuracy 1e-12;
    - Nelder–Mead in logit coordinates over the engine's box;
    - 125 lines, readable top to bottom.
  - **Three checks against the engine's committed results on the six M1.7 panels:**
    1. the panel log-likelihood at 702 surface points (`tests/golden/validation/`, kept current by `unit_xref: validation_goldens_are_current`);
    2. the log-likelihood at scipy's optimum minus that at the engine's estimate;
    3. the estimates, in units of the engine's SE.
  - **A fifth named kind of comparison (D-107):** `TOL_SCIPY_SURFACE_LL_REL` = 1e-11, `TOL_SCIPY_OPTIMUM_LL_ABS` = 1e-3, `TOL_SCIPY_ESTIMATE_SE` = 0.06. They were set from what the plain script achieves: 4.1e-13, 3.9e-4 and 0.0255 respectively. The ε-level core-vs-ref and cross-platform tolerances are not reused.
  - **CI:** the job `validation (scipy)` runs the script with scipy 1.18.1 and numpy 2.5.3 pinned, taking about five minutes.
- **D-127 — The methodology note `docs/methodology/binomial_mixture_mle.md` is self-contained for anyone replicating the fit.** It covers:
  - the model and likelihood;
  - both integration rules, and why fixed-node GH fails at large n (D-037);
  - the validated ranges;
  - the per-run check and its empirical factor;
  - the materiality statement;
  - estimation, SEs and every flag, including near-bound;
  - the known statistical findings of M1.8;
  - the scipy replication.

  **Parity sign-off (D-123):** `docs/methodology/parity_signoff.md` lists the 19 known findings and the 31 deferred scenarios (62 verdicts), each with its M2 dependency.

## M2 (first part): profile-likelihood intervals, 2026-09-26

M2's scope stays as STATE.md defines it: profile intervals first, because they clear the deferred verdicts, then the resampling engine and C ABI v0.

- **D-128 (owner's M2 notes) — Profile-likelihood intervals, `engine/profile.hpp`, for headline fits.** For parameter a with nuisance b, the profile is P_a(u) = max over u_b of ℓ, in logit coordinates. The steps:
  1. **Polished maximum:** ℓ_max is maximised off-grid. Brent's method runs on the PD profile around the refined estimate, and each profile value is itself a Brent maximisation over ρ. Using the grid's value would misplace the threshold by up to about 4e-4 of log-likelihood (the refinement error, M1.9).
  2. **Inner maximum:** bracketed by the surface's own argmax in the neighbouring grid columns, ±2 steps. The bracket widens if the maximum lands on an edge that is not a bound of the box.
  3. **Endpoints (owner: resolution at the crossing):** walk outwards along the grid profile. Its values are lower bounds on the true profile, so grid points above the threshold need no work. At the first grid point below, the objective is evaluated directly, and Brent's root finder solves P_a(u) = ℓ_max − c to 1e-9 in u.
  - **Accuracy:** each fit reports its worst endpoint residual in log-likelihood. The register bounds it (`TOL_PROFILE_ENDPOINT_RESIDUAL_LL` = 1e-8, observed 2.1e-9) against an independent dense profile.
  - **Cost:** about 870 panel evaluations per fit, serial and deterministic.
  - **Independent check:** `ref/` gained its own implementation, using golden-section over the whole range and bisection. Core agrees with it to 1.5e-10 in logit units, and on every truncation.
- **D-129 (owner) — The 95% interval uses the χ²₁ threshold, c = χ²₁(0.95)/2 = 1.9207294103470630, with truncation at the box, including near ρ = 0.**
  - **Why this matters:** when a parameter's true value is on the boundary of the parameter space (ρ = 0), the likelihood-ratio statistic for *testing* it is not χ²₁. It is the mixture ½χ²₀ + ½χ²₁ (Chernoff 1954; Self and Liang 1987), whose 5% critical value is 2.71, not 3.84.
  - **For confidence intervals:** standard practice is χ²₁ with the interval truncated at the boundary. That is what the engine does, and anyone replicating it should know the choice was deliberate.
  - **Testing ρ = 0:** a test of "no correlation" must use the mixture, not this interval. That is not offered in parity.
  - **The box:** its lower ρ bound is 10⁻³, not 0, so the interval is truncated there.
- **D-130 (owner) — An interval that does not reach the threshold before a bound of the box is truncated at the bound and flagged (`kIntervalLowerTruncated` / `kIntervalUpperTruncated`), never extrapolated.**
  - **Common case:** at low ρ and in low-default panels the lower ρ end is usually the bound. A panel with no defaults puts both the PD estimate and the lower PD end on the bound.
  - **No interval:** fits flagged flat or numeric get `kIntervalNotComputed`, and harness coverage counts them as not covering; they are never dropped.
  - **Owner's note on the tag:** the `v0.1.0-parity-binomial` annotation says "13 jobs"; the correct count is 12 (STATE.md records it). Tags stay immutable once pushed. In the public repository, `v0.1.0-parity-binomial` marks the initial release snapshot, which contains the M1 sign-off.

- **D-131 (owner) — The recovery harness reports the profile-likelihood interval, Wald and t(T−1) side by side. Profile coverage counts all replicates and nothing is deferred. The 26 out-of-band profile verdicts are pinned under two labels.**
  - **The three intervals:**
    - **Profile:** coverage over all 1,000 replicates, because a truncated interval is still an interval, and a not-computed one counts as not covering.
    - **Wald:** unchanged (conditional on unflagged, DEFERRED at 5% flagged).
    - **t(T−1):** a `native` comparison on the same replicates, reported without a verdict.
  - **The band:** 0.95 ± 3.29·√(0.95·0.05/m), a two-sided 0.1% band. About 0.16 of 162 verdicts would fall outside by chance, so out-of-band verdicts are real effects, not noise.
  - **Results (Wald → profile), the headline of M2a:**

    | Wald \ profile | PASS | CONSERVATIVE | KNOWN FINDING | total |
    |---|---|---|---|---|
    | PASS | 80 | 0 | 1 | 81 |
    | DEFERRED | 41 | 20 | 1 | 62 |
    | KNOWN FINDING | 15 | 0 | 4 | 19 |

    - Profile fixes 15 of the 19 Wald findings and 41 of the 62 deferred verdicts.
    - What remains is either conservative, or within a point of the band.
  - **Labels (owner): two, each valid only on its own side of the band.** A reviewed verdict that crosses to the other side becomes UNREVIEWED again.
    - **CONSERVATIVE (20):** above the band, 0.977–0.999. These are near-uninformative settings (n = 100 at PD 0.1%, or ρ = 0.02 with little data), where 25–100% of intervals are truncated at a bound of the box. Such intervals contain the truth more often than advertised: a safe-side property, not a problem.
    - **KNOWN FINDING (6):** below the band, 0.920–0.927 against 0.927, at T = 20–40. This is the genuine limitation: small-sample undercoverage of the likelihood-ratio interval.
      - **ρ cases:** most misses lie below the truth (61–65 against 9–16), which is ρ̂'s downward small-T bias.
      - **PD cases:** misses are symmetric.
  - **Endpoint accuracy (owner):** `TOL_PROFILE_ENDPOINT_RESIDUAL_LL` = 1e-7. The worst over all 81,000 fits was 1.7e-8, the margin chosen by the owner. The harness asserts it on every fit, and CI checks each scenario's recorded maximum.
  - **Not done, by decision:** changing the threshold or the estimator after seeing these results. Both remedies go to the backlog as `native` options, each to be evaluated against these same pinned verdicts:
    - a Bartlett-type correction of the 1.92 threshold;
    - a bias-corrected ρ̂, the small-sample bias benchmark already planned.

## M2b: resampling engine, 2026-09-26

- **D-132 (owner) — Resampling draws come from Philox4x32-10 in a key domain of their own: key = seed XOR `0x52534D50424F4F54` ("RSMPBOOT"), counter = (scheme, replicate, draw/2, 0).**
  - **Draws:** draw j uses words (0, 1) if j is even and (2, 3) if odd. The uniform is (k + ½)·2⁻⁵², exactly as in the DGP. An index in [0, m) is ⌊u·m⌋, clamped to m − 1, with bias below m·2⁻⁵².
  - **Key domain:** for any seed, the key differs from the DGP's, so resampling a simulated panel can never reuse or shift the stream that generated it.
  - **Addressing:** (seed, scheme, replicate, draw), recorded like the DGP's.
  - **Independent check:** a Python mirror in `tools/gen_dgp_tables.py` writes every draw of a reference iid and moving-block bootstrap, and the C++ matches it bit for bit.
- **D-133 (owner) — The block bootstrap is the moving-block bootstrap (Künsch 1989), non-circular, with a fixed block length ℓ from the caller.**
  - **Construction:** ⌈T/ℓ⌉ blocks with uniformly drawn starts in [0, T − ℓ], laid end to end and cut at T.
  - **Default:** ℓ = ⌈T^(1/3)⌉, the usual rate-optimal order.
  - **Later `native` options:** the stationary bootstrap (Politis–Romano) and data-driven block length (Politis–White).
  - **Scope (owner):** the DGP's periods are independent, so the moving-block bootstrap is the wrong tool for the recovery matrix. Its coverage study moves to M6 with the AR(1) DGP, where block should beat iid. It is tested for correctness now.
- **D-134 (owner) — The resampling engine never stores B × K surfaces.**
  - **Pass 1:** for each replicate, the fused, tiled `reduce_weighted` finds the grid argmax. Tile order is fixed, so results are bitwise identical for any thread count.
  - **Pass 2:** the 3×3 refinement stencil is recomputed from L and that replicate's weights.
  - **Memory:** O(B × tiles).
  - **Supplied index matrices:** alongside generated W, the engine accepts an externally supplied B × m period-index matrix and turns it into counts. Supplying another tool's draws reproduces its bootstrap replicate for replicate.
  - **Verified:**
    - bitwise equal to materialising every replicate's surface, at 1, 3 and all threads;
    - equal to refitting the explicitly resampled panel (observed difference 0, `TOL_RESAMPLE_VS_REFIT_REL` = 1e-12).
- **D-135 (owner) — Bootstrap intervals are percentile intervals.**
  - **Quantile rule:** Hyndman–Fan type 7, as in `numpy.percentile` and R.
  - **Scale:** percentile intervals are invariant under monotone reparametrisation, so the logit-scale interval recommended for bounded parameters is the same interval.
  - **Replicates used:** grid-edge replicates are kept at their grid value, and non-finite ones are excluded and counted.
  - **Not offered:** BCa is a later `native` option; the basic interval is not offered.
- **D-136 (owner) — The bootstrap coverage study covers the full recovery matrix (81 scenarios × 1,000 replicates) with the iid bootstrap, B = 999, and its expected verdicts were written down before the run.**
  - **Cost:** measured at 0.04–0.2 s per panel on one core, 15–30% of the fit, not the days first estimated. The study therefore runs now rather than waiting for M4.
  - **Stream per panel:** each recovery panel's bootstrap uses seed = recovery seed XOR (scenario << 32 | replicate).
  - **Predictions:** `docs/methodology/bootstrap_predictions.md`, committed before the harness ran.

- **D-137 (owner) — Bootstrap percentile intervals undercover: 125 of 162 verdicts fall below the band and none above. Profile-likelihood intervals are the recommended method for inference in this model. Bootstrap percentile intervals are provided for comparison and are not recommended for ρ.**
  - **Results against the predictions registered before the run** (`docs/methodology/bootstrap_predictions.md`, committed before the run; see D-146). The misses are recorded as well: showing where a prediction was wrong is what makes the pre-registration credible.

    | Prediction | Result | Held? |
    |---|---|---|
    | 1. Group A, ρ: at least 2/3 below the band, many below 0.90; not conservative | 31/31 below, 26 below 0.90 (0.431–0.922); 2,878 zero-width intervals | Yes, in direction and severity. The mechanism was only partly right: in n = 100, PD 0.1% (scenarios 0, 9, 18, 21) many intervals collapse at the *upper* bound, not the lower |
    | 2. Group A, PD: mixed; below mainly where PD̂ sits at its bound (n = 100, PD 0.1%) | 26/31 below (0.762–0.944); the 5 PASS are 4, 6, 7, 33, 60 | No: undercoverage is far more widespread than predicted |
    | 3. Group B (T = 20): about 0.88–0.93, most of 24 below; ρ misses below the truth | PD 11/12 below (0.875–0.936); ρ 12/12 below (0.803–0.887); ρ intervals below the truth in 1,883 of 1,940 misses | Direction yes; ρ worse than the predicted range |
    | 4. Group C (T = 40): PD mostly in the band; ρ mixed, about 0.91–0.93 | PD 7 PASS, 10 below; ρ 17/17 below (0.874–0.909) | No: worse than predicted for both |
    | 5. Group D (T = 100): both in the band | PD 18/21 PASS; ρ 6/21 PASS, 15 below (0.908–0.938) | PD yes; ρ no |
    | 6. At most 5 CONSERVATIVE | 0 | Yes |
    | Totals: below 50–90, in band 70–110, above ≤ 5 | below 125, in band 37, above 0 | Undercoverage underestimated |

  - **Diagnoses (two):**
    - **Boundary breakdown (group A, 57 verdicts).** The bootstrap is inconsistent near a boundary of the parameter space (Andrews 2000). Resampled estimates pile onto a bound of the box and the percentile interval collapses; 3,500 intervals had zero width.
    - **No bias or skewness correction (the other 68).** ρ̂ is biased downwards and skewed towards zero in logit coordinates. The bootstrap distribution reproduces that shape, and the percentile interval does not correct it, so its intervals sit below the truth. That is why ρ undercovers even at T = 100. At T = 20 the iid bootstrap's variance is also (T − 1)/T of the true one.
  - **Engine correctness is unaffected.** The fused reduction is bitwise equal to a materialised one, and a replicate equals the refit of its resampled panel.
  - **Recording (owner):**
    - Every verdict is pinned in `kKnownBootstrapFindings` for CI, which costs nothing to maintain and catches any change.
    - Documents for people describe the result at the method level, with a group table, and do not list the 125.
  - **BCa stays on the backlog** as a `native` option, to be evaluated against these pinned verdicts. It will not fix group A, since BCa also breaks down at a boundary.

## M2c: C ABI v0, 2026-09-27

Owner's notes: an ABI is hard to change once others build on it, so get versioning, ownership,
errors, the context, symbols and a pure-C test right from the start. Scope it to what exists
today, and leave anything not yet built out of v0 rather than stubbing it.

- **D-138 (owner's M2c notes) — Versioning: `vcal_abi_version()` returns (major << 16) | minor, and every public struct starts with a `uint32_t struct_size` set by the caller.**
  - **Version bumps:** an additive change (a field appended, a function added) bumps the minor. A change that could break an existing caller bumps the major. The first version is 0.1; the M8 freeze is 1.0.
  - **Struct sizes:** structs only grow, and the library honours what the caller's struct_size says:
    - a size below the 0.1 size is `VCAL_E_INVALID_ARGUMENT`;
    - a larger input struct is accepted only if every byte beyond the fields the library knows is zero; otherwise `VCAL_E_UNSUPPORTED`. A newer caller's nonzero setting is refused, never silently ignored;
    - a larger output struct has only its known fields written.
  - **Reserved fields:** fields named `reserved` must be zero on input. They are future flags an old caller leaves at zero.
  - **Layout:** no struct has implicit padding, so its size is the sum of its fields and a new field can never hide in trailing padding. Checked by `static_assert` in the library and `_Static_assert` in the C test. 64-bit platforms only.
- **D-139 (owner's M2c notes) — Ownership: the caller allocates every buffer, and the library never returns memory for the caller to free.** This matters on Windows, where the caller and the library may use different C runtimes.
  - **Pattern:** variable-size outputs take (buffer, capacity, required).
    - A NULL buffer is a size query: it stores the length and returns `VCAL_OK`.
    - A buffer that is too small stores the length and returns `VCAL_E_BUFFER_TOO_SMALL`, writing nothing.
  - **Multi-array outputs** (`vcal_resample`) use a struct of caller-owned array pointers with one capacity. Any pointer may be NULL, and adding an output later appends a pointer.
  - **Static strings:** `vcal_status_string` returns a static string, which is never freed.
  - **Output structs** are written only after everything has been computed.
- **D-140 (owner's M2c notes) — Errors: every function returns a status, and `vcal_last_error` copies out the calling thread's message. No C++ exception crosses the boundary.**
  - **Statuses:** OK, INVALID_ARGUMENT, BUFFER_TOO_SMALL, UNSUPPORTED, NUMERIC, OUT_OF_MEMORY, FP_ENVIRONMENT (the DGP refuses flush-to-zero, D-109) and INTERNAL.
  - **Catch-all:** every entry point runs inside a catch-all that converts exceptions to statuses. `std::bad_alloc` becomes OUT_OF_MEMORY, and anything unexpected becomes INTERNAL with its message.
  - **Message storage:** a fixed thread-local buffer, so recording a failure cannot itself fail.
  - **Clearing:** every status-returning call clears the message on entry, except `vcal_last_error` itself. `vcal_abi_version` and `vcal_status_string` leave it alone, so a caller can look up a status name before reading the message.
  - **Validation:** inputs are validated in the ABI with messages naming the field and period. The engine's own checks stay as a second line, and a disagreement between the two is reported as INTERNAL.
  - **OpenMP:** exceptions cannot leave an OpenMP region. The parallel loops in the engine do not allocate or throw, and every allocation happens outside them.
- **D-141 (owner's M2c notes; Q12, D-045) — An opaque context, from `vcal_context_create` / `vcal_context_destroy`, holds the profile and the thread count.**
  - **Profile:** `VCAL_PROFILE_PARITY` is the only value in v0.
  - **Threads:** n_threads = 0 means the OpenMP default. Results are bitwise identical for any count (§6), so the count is configuration, not a numerical input.
  - **Concurrency:** a context must not be used by two threads at once; distinct contexts are independent. v0's context is immutable, so this is stricter than needed. It keeps room for caches or GPU streams later, and relaxing a rule is compatible while tightening one is not.
  - **Every numerical call takes a context,** including `vcal_dgp_simulate`, which uses nothing from it today; the GPU DGP (R-1) will.
- **D-142 (owner's M2c notes) — Symbols: a shared library `vcal` whose only exports are the `vcal_*` functions of `include/vcal/vcal.h`.**
  - **Export and call macros:** `VCAL_API` is dllexport/dllimport on Windows and default visibility elsewhere, with `VCAL_STATIC` for a static consumer. `VCAL_CALL` is `__cdecl` on Windows, so a caller compiled with another default convention still links.
  - **Hidden by default:** everything else is hidden. On Linux a linker version script also hides the standard library's template instantiations, which would otherwise be exported with default visibility, and `--no-undefined` catches missing symbols at link time.
  - **Export check:** `abi_exports` compares `nm -D` (Linux) or `dumpbin /exports` (MSVC) with the header's declarations, and `abi_exports_fires` proves the comparison catches a difference.
  - **Provenance:** `vcal_build_info` reports the version, ABI, git commit and dirty flag (read at every build, never stale), compiler, C++ standard, build type, OpenMP, CUDA and profiles.
- **D-143 (owner's scope rule) — v0 covers what exists and nothing else: calibrate with profile intervals, surface, resample (with the weights it implies), DGP simulate, grid helpers, build info, errors and the context.**
  - **The single tuple served:** the binomial-mixture MLE, parity rule, ArgMax and CPU. So there is no model spec and no dispatch registry yet; both arrive with M3's second estimator, as appended fields or new functions.
  - **Left out:** backtests (M5) and a CUDA backend field (M4) are left out too, not stubbed.
  - **Resampling:** all six W sources are exposed, since the engine consumes any W. Supplied weights must be finite and non-negative, with at least one positive weight per row.
  - **Fixed at 95% in v0:** the profile-interval level. It is a fixed χ²₁ threshold in the engine, and a `level` field can be appended. The bootstrap percentile level is a field already, because the engine takes it.
  - **Default block length** (D-133): computed in integers as the smallest ℓ with ℓ³ ≥ T, so no floating-point cube root can round at a perfect cube.
  - **Grid bounds:** both axes must satisfy 0 < lo < hi < 1 whatever their scale, since both parameters are probabilities.
- **D-144 (owner's M2c notes) — The ABI is tested from plain C and against the engine.**
  - **`abi_c` (C11, linked like an outside consumer):**
    - versioning, struct_size handling in both directions, size queries, every status name, messages that name the field, and thread-count invariance;
    - calibrate, profile, surface, every resampling scheme's replicate count, and the iid bootstrap;
    - supplied identity indices reproduce the estimate bit for bit;
    - the DGP reproduces the Python mirror's reference panels bit for bit.
  - **`abi_c99_header_check`:** compiles the header alone as C99 with `-pedantic-errors` (C11 on MSVC, which has no C99 mode).
  - **`abi_engine_equivalence` (C++):** every number the ABI returns is bitwise what the engine gives when called directly. This covers calibrate, profile, surface, all six resampling schemes with their intervals and weights, the DGP, and the default grid. Together with the recovery replay, it ties the ABI to the pinned goldens without a tolerance of its own.

## CI economy, 2026-09-27

- **D-145 (owner) — CI spends runner minutes where they buy a check.** The aim is faster feedback and fewer redundant runs.
  - **Measured cost:** each Release job spent 7–8 minutes on the slow-labelled tests (`xref_estimates`, `xref_profile`, `recovery_replay`). The accuracy step then ran `recovery_replay` a second time, taking up to 4 minutes on GCC 11.
  - **What runs when** (`.github/workflows/ci.yml`; a `plan` job decides):

    | Trigger | Jobs | Slow tests |
    |---|---|---|
    | push to `main`, pull request | the full 12 | no |
    | push to any other branch | `linux-gcc11 / cpu-release`, `windows-msvc / cpu-debug`, generated files | no |
    | nightly (04:23 UTC), only if `main` moved in the last 25 hours | the full 12 | yes |
    | `workflow_dispatch` | the full 12 | yes |

  - **Prose-only commits start no run.** These are `*.md`, `LICENSE` and `NOTICE`, except the Markdown files the code reads, which the path filters re-include. This replaces the earlier practice of confirming CI on docs-only commits too.
    - **The re-include list is derived, not remembered (owner).** `ci_path_filter_sync` (`cmake/CiPathFilterCheck.cmake`) scans what CI executes: `tests/`, `cmake/`, `tools/`, `validation/` and the top-level CMakeLists.txt. It strips comments and Python docstrings, collects every reference that names an existing Markdown file, and fails if the workflow's list differs in either direction. `ci_path_filter_sync_fires` proves that it reports a file that is read but not listed, and a file that is listed but not read.
    - **Today's list:** `docs/methodology/tolerances.md` (the register sync), `docs/methodology/recovery_results.md` (the recovery tests) and `README.md`.
    - **The README's C example is tested (owner).** The owner's review assumed a test compiled it; it did not, since it had only been compiled once by hand. `readme_c_example` now compiles the first `c` block of README.md as C11 against the shared library, runs it, and requires its output to equal the `text` block that follows. A wrong documented number fails the test, and that was checked.
  - **Concurrency:** a newer push to a branch or pull request cancels its older run. Runs on `main` are never cancelled, so every commit that reaches `main` keeps its own result.
  - **Rules that follow:**
    - `main`'s head is green when the run for its latest non-prose commit is green. Prose-only commits after it inherit that result.
    - Work goes to a branch (reduced set). The full matrix runs through a pull request, or on `main` itself, and must be green before new work starts.
    - Before any tag or sign-off, a `workflow_dispatch` run on that commit, with the slow tests, must be green.
  - **Trade-off, accepted:** the cross-platform replay (`TOL_RECOVERY_REPLAY_*`, `TOL_PROFILE_CROSS_PLATFORM_REL`, `TOL_BOOTSTRAP_CROSS_PLATFORM_REL` on the replay) and the core-vs-ref fits now run nightly and on demand, not on every push. The cross-platform regression values (`xref_core_regression`, 1e-12) and the recovery exit criteria still run on every push.

## Pre-publication pass, 2026-09-27

- **D-146 (owner) — Before the first public release, the repository is described as what it is: an independent credit-risk research library for practitioners, researchers and students, designed so every result can be independently replicated, and a personal open-source project of its author.** Recorded as one pass, with the CI changes of D-145 merged in first. The public repository is published from a single squashed snapshot of the result.
  - **Audience.** CLAUDE.md, README.md and the docs describe that audience. Wording that tied the project to a particular kind of user, or to one user's review process, is replaced by what the property is for, in terms of independent replication: D-002, D-005, D-006, D-014, D-016, D-018, D-032, D-033, D-036, D-037, D-049, D-051, D-077, D-089, D-121, D-125, D-127, D-129 and the §3.7 profile paragraph of ARCHITECTURE.md. The README's guide for reproducing results is now headed "Replicating the results", and M8 delivers a "replication pack". What each of those decisions requires is unchanged.
  - **Toolchain floors are portability goals.** The project supports CUDA ≥ 11.8 and older host compilers (GCC 11, MSVC 14.39) for broad compatibility, with device-visible code in C++17; CMake stays at or below 3.26 so that it configures on the older CMake that long-term-support distributions ship. D-026, D-059, D-060, D-066, D-094, D-101, D-104 and D-105 are restated in those terms. Every version, flag, CI job and rationale is kept. The open item about a second toolchain target is withdrawn: the nvcc 12.8 CI job stays, as the newer end of the supported range and the first toolkit that builds sm_120.
  - **In-place rewording.** Like D-104, this pass edits earlier entries in place, the second exception to the append-only rule. Only rationale wording changed; no decision was reversed.
  - **Authorship.** NOTICE names Keith Yip as the copyright holder. LICENSE stays the unmodified Apache-2.0 text, as the Apache convention asks and as licence detection expects. The README says this is a personal open-source project. No organisation is named as a user, sponsor or context.
  - **ARCHITECTURE.md §3.0** states each layer (objective, integrator, surface engine, reducer, resampling, reference, DGP) as a language-agnostic contract: inputs, outputs and invariants, so it can be re-implemented without reading the C++. The C++ signatures stay, and those of resampling, the reference and the DGP are added.
  - **CI path filter (D-145).** `ci_path_filter_sync` now also scans CUDA sources, TOML and JSON under `tests/`, and `CMakePresets.json`. An independent grep of every non-Markdown file for `.md` references agrees with the derived list: `README.md` (its C example is compiled and run), `docs/methodology/tolerances.md` and `docs/methodology/recovery_results.md`. Every other `.md` mention in code is in a comment or a link written into generated Markdown.
  - **Local notes.** `.gitignore` keeps the single `*.local.md` rule for machine-specific notes (D-104).
  - **References that a squashed snapshot cannot resolve.** CI run numbers and commit hashes from before the public release are replaced by "in CI prior to public release" and similar neutral wording; the measured values themselves are unchanged. The tag `v0.1.0-parity-binomial` is placed on the initial release snapshot.
  - **The bootstrap predictions were pre-registered (D-136).** `docs/methodology/bootstrap_predictions.md` was committed on 2026-09-26 at 21:22 (UTC−4), about three hours before the bootstrap coverage results (2026-09-27, 00:21), and has not changed since: its SHA-256 is `e37a9887f2da4df3722541349071bcc4834d8a20cd029d8c6f134eac7d954666`. This ordering was verified prior to public release. The file itself is left as written, because it promises not to be edited after the run.

## Enforced green main, 2026-09-27

- **D-147 (owner; supersedes D-145's path filters, keeps its job plan) — `main` changes only through a pull request whose required check `ci-ok` has passed on a head that is up to date with `main`. The rule "main is always green" is enforced by a repository ruleset instead of being a convention.**
  - **The ruleset** is committed as `.github/rulesets/main.json` and imported under Settings → Rules → Rulesets. On the default branch it requires a pull request (no approving review needed, since this is a one-person project), requires `ci-ok` from GitHub Actions with the branch up to date, and forbids force-pushes and deletion. No one bypasses it, the owner included.
  - **One required check, `ci-ok`.** The CI job names come from the `plan` job's matrix and change with the trigger, so requiring them one by one would be fragile. `ci-ok` always runs (`if: always()`) and passes only if `plan` succeeded and every job the plan called for succeeded. A job the plan skipped counts as passing; a failed or cancelled one never does.
  - **No workflow-level path filter.** A run that never starts reports no check, so under D-145's `paths` filters a prose-only pull request would have waited forever for `ci-ok`. The filters are gone. `plan` now reads the diff (the pull request's base...head, or a push's before..after) and, when every changed file is prose, skips all other jobs; `ci-ok` still reports. Prose is defined as before: `*.md`, `LICENSE` and `NOTICE`, except the Markdown files the code reads, which `plan` lists in `code_md`. An unknown or unreachable base (a new branch, a force push, the nightly and `workflow_dispatch`) counts as code.
  - **Checks.** `ci_path_filter_sync` now compares the derived list with `plan`'s `code_md` list and also fails on any workflow-level `paths` or `paths-ignore` filter. `ci_path_filter_sync_fires` (a wrong list) and `ci_path_filter_sync_fires_paths` (a path filter) prove that both checks fire. The `plan` logic was simulated on prose-only, code-reading-Markdown, code and new-branch diffs, and the `ci-ok` decision on every combination of job results, before the first run.
  - **Rules that follow:** work goes to a branch and reaches `main` through a pull request, which runs the full matrix. The nightly and `workflow_dispatch` rules of D-145 are unchanged.

## Open

- **R-1 (revisit at M6) — GPU-side DGP.** Only matters for large parametric bootstraps (AR(1), D-043).
- **R-2 (revisit at M7, when the macro pipeline starts interacting with the estimator) — Z sign convention and macro sign filters.** Declare expected macro effects in economic terms ("worsens" or "improves credit conditions"). Map them to coefficient signs in one place, through the engine's Z convention, and test the mapping on a synthetic DGP.
  - **The engine's convention:** p(z) = Φ((Φ⁻¹(PD) − √ρ·z)/√(1 − ρ)), in `core/model/vasicek.hpp` and `dgp/`. A higher Z therefore means *better* credit conditions (a lower conditional PD).
  - **The risk:** a macro variable that worsens credit conditions has a negative coefficient on Z, but a positive one on a PD or default-rate scale. A sign filter written against the wrong scale silently keeps the wrong specifications.
