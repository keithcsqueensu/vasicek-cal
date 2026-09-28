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

## Research studies and variants: ground rules and addenda, 2026-09-27

Owner's notes: a roadmap of research studies and `native` variants (S-1 to S-20 here; the addenda D-149 to D-151 extend it to S-33), indexed in
`studies/README.md`. Planning only; nothing is implemented until the owner approves
the plan.

- **D-148 (owner) — Ground rules for research studies and variants.**
  - **Variants never change parity.** Every variant is a `native` option or a standalone study.
  - **Every study is pre-registered.** Before any run, `studies/<slug>/PREDICTION.md` is committed,
    stating the expected results, their direction and their rough size. Results are compared with
    it, and misses are reported, not explained away. This generalises D-136, whose bootstrap
    predictions were committed before the run and never edited afterwards.
  - **Explore on the study subset; pin on the full matrix.** Exploration uses a fixed subset of
    about 9 scenarios, spanning short, medium and long T and low and high ρ, that runs in about 10
    minutes on the CPU. The full 81 × 1,000 matrix is run only to pin a study's verdicts.
  - **Pinned verdicts follow the M1.8/M2a policy:** PASS, CONSERVATIVE, KNOWN FINDING or DEFERRED,
    each with a diagnosis (D-124, D-131).
  - **Every study has an entry in `studies/README.md`:** the question, the experiment,
    the prediction, the result and the status of any mitigation.
  - **The study subset (approved by the owner; recorded in D-150):** scenarios 29, 37, 72 (T = 20),
    4, 68, 49 (T = 40) and 7, 43, 51 (T = 100), one per (T, ρ) cell, at R = 1,000 with the recovery
    seeds, so every subset panel is the pinned run's panel and comparisons are pairwise. It holds
    3 of the 6 small-T profile findings, a truncation-conservative verdict, both bootstrap
    diagnoses, a Wald skewed-overcoverage finding and profile PASS verdicts, and touches all four
    D-136 groups. **Measured:** 16.7 s of single-thread CPU per replicate across the nine (fit
    10.5, profile 4.2, bootstrap 2.0), so 4.6 thread-hours at R = 1,000: about 12 minutes on 24
    threads at the measuring machine's speed, an estimated 7 on the development machine. The full
    matrix measured 50.6 thread-hours for the same steps.
  - **Placement:** S-1 to S-7 now or alongside M3; S-8 and S-9 with M3; S-10 and S-11 on the subset
    now and the full matrix after M4; S-12 to S-17 after M4; S-18 to S-20 deferred to M7 and
    recorded together (S-18 is R-2).
  - **Noted for the owner, not decided:** with a fixed n per period, the 1,000 replicates of a
    scenario contain far fewer distinct default counts than their periods do (measured on the
    subset: 10 to 2,074 distinct d across a scenario's 1,000 replicates, against 4,750 to 39,000
    surface rows evaluated today after D-122's deduplication within each panel).
    A surface cached by (n, d), with each panel a count vector over it, would make the parametric
    bootstrap (S-10, S-4b), the misspecification runs (S-11) and R = 10,000 (S-13) W × L jobs. It
    extends D-122 across panels, and results would be bitwise unchanged. A spike would confirm it
    before any study's placement changes.

- **D-149 (owner; extends D-148) — Addendum: studies S-21 to S-28, and two further ground rules.**
  - **Placement:**
    - **Now** (existing machinery): S-21 period influence, S-22 box sensitivity, S-23 intervals for derived quantities.
    - **Now, on synthetic data variants:** S-24 PD heterogeneity and S-25 PD trend (DGP variants), and S-26 varying portfolio size (a scenario definition only: the DGP already takes n per period, so no new mirror or hash).
    - **With M3:** S-27 the large-portfolio approximation and S-28 zero-default treatments for rate-based estimators.
  - **A variant is judged on every verdict it could change**, not only the findings it targets. A fix that moves the targeted findings into the band but pushes PASS verdicts out of it is reported as both.
  - **A new DGP variant is specified like the base DGP:** a prose description, a line-for-line Python mirror and a reference-panel hash (D-110–D-112), in its own key domain.
  - **Recorded while planning:** a diagnostic refit of replicates 0–39 of the 20 CONSERVATIVE profile verdicts (the pinned run's own panels, no variant) found that truncation is almost only at the lower bounds (ρ ≥ 1e-3, PD ≥ 1e-4). Upper-only truncation was at most 6 of 40 in any verdict. S-22 therefore gains a lower-bound arm beside the ρ ≤ 0.9 arm. This finding is recorded in its `studies/README.md` entry so that S-22's prediction is written knowing it.

- **D-150 (owner; extends D-148 and D-149) — Decisions on the D-149 points, and studies S-29 to S-33 under the theme "data and population instability".**
  - **Approved:**
    - S-22's lower-bound arm (ρ ≥ 1e-5, PD ≥ 1e-6), with the grid's spacing held fixed so it is not confounded with S-7;
    - S-21 shares the jackknife run with S-3 and S-5, with exact refits to measure the quadratic refinement's error;
    - S-24's separate stream key for the second segment;
    - S-25's oracle and two-stage detrending arms now, the joint fit in M7;
    - S-26 paired with the constant-n scenarios and repeated for the M3 estimators;
    - S-27 and S-28 in the S-8 estimator-comparison run, with S-28's "refuse" row reported as the refusal share.
  - **Misspecification studies report two targets** (S-11, S-24, S-25, S-29 to S-33):
    - **primary:** the pseudo-true value, the value the fit converges to on a very long panel. That is what the estimator estimates under misspecification, and coverage is judged against it;
    - **also reported:** the distance from the generating parameter, because that gap is the misspecification bias a practitioner cares about.
  - **Placement of S-29 to S-33:**
    - **Now** (single-segment engine and small DGP variants):
      - S-29: a two-grade version, with per-grade fits on the remapped history, a two-stage structural break and a portfolio-level fit with allocation;
      - S-30: per-bucket separate fits;
      - S-31: the bias, exclusion and a two-stage joining indicator;
      - S-32: in full;
      - S-33: at portfolio level.
    - **M6** (multi-grade calibration): S-29's probabilistic mapping, joint break and full scale; S-30's shared ρ and its cross with S-29; S-31's separate segment; S-33's late rating assignment.
  - **The study subset is approved** as proposed in D-148: the paired design of 9 scenarios (29, 37, 72, 4, 68, 49, 7, 43, 51) at R = 1,000 with the recovery seeds.
  - **Raised for the owner (all four decided in D-151):**
    - **Pseudo-true value for non-stationary designs.** A "very long panel" changes the design when the design itself changes over the sample (S-25, S-29, S-31, S-32, S-33). Proposed definition: the maximiser of the expected log-likelihood of the study's exact design, exact by quadrature where the recorded counts have a computable distribution, simulated otherwise. For a stationary design it equals the very-long-panel value.
    - **Multi-grade calibration is not in M6 as written.** M6 lists multi-segment correlated factors. K grades with their own PDs on one factor with a shared ρ (the ASRF model) is a different, (K + 1)-parameter model that needs its own optimiser. Proposed: add it to M6 explicitly.
    - **One DGP observation layer.** S-29, S-31, S-32 and S-33 all alter what is recorded from a true panel: grade relabelling, missed or spurious defaults, dropped defaulters, a segment added from a date. Proposed: build these once as a single layer with one prose description, one Python mirror and one hash.
    - **A grade scenario matrix** for S-29 and S-30, whose scales the single-PD subset and matrix cannot express.

- **D-151 (owner) — Decisions on the D-150 open points, and a priority order for the studies.**
  - **Multi-grade calibration becomes an explicit M6 item:** a single-factor model with K grades, each with its own PD, and a shared ρ.
    - All grades share the same Z_t, so each period's likelihood is still a one-dimensional integral, ∫ φ(z) ∏_k Binom(d_kt; n_kt, p_k(z)) dz. The quadrature stays cheap; only the optimiser grows.
    - Approach: a one-dimensional grid in ρ, with the K PDs maximised for each fixed ρ (Newton or coordinate ascent). This keeps the profile-likelihood machinery for ρ intact.
    - The adaptive rule's mode hint is revisited for the product integrand.
  - **The pseudo-true value** (for misspecification studies, D-150) is the maximiser of the expected log-likelihood under the study's exact design: by quadrature where the recorded counts have a computable distribution, by simulation otherwise. For a stationary design it equals the very-long-panel value.
  - **One DGP observation layer**, with one prose description, one Python mirror and one hash, shared by S-29 and S-31 to S-33.
  - **The grade scenario matrix** is defined explicitly in `studies/README.md`, with its own seed (`0x475241444553434E`, "GRADESCN") and description:
    - 3 scales (L, M, H) × T ∈ {20, 40, 100};
    - a 17-grade master scale, PD_k = 0.03% × 1.5^(k−1);
    - 10,000 obligors per period, spread over the grades by a discretised normal of standard deviation 4 grades;
    - ρ of 0.24, 0.12 and 0.02, and portfolio PDs of 0.61%, 1.54% and 3.42%;
    - R = 1,000.
  - **S-33's redesign is approved:** partial survivorship, and a misjudged recording-change date.
  - **S-32's misclassification-aware objective has two oracles:** a scipy script and `ref/`. The changed integrand is exactly where the adaptive rule's assumptions could fail, and `ref/` shares no quadrature code with core.
  - **Priority order.** 33 studies are a research programme, and running them in placement order would delay M3 indefinitely.
    1. **First batch, the highest value per hour:**
       1. S-23;
       2. S-13 targeted (R = 10,000 on the borderline scenarios, which settles whether the 6 small-T profile findings are real);
       3. the shared jackknife run (S-3, S-5, S-21);
       4. S-1 and S-2.
    2. **Then M3,** with S-8, S-9, S-15, S-27 and S-28 folded into it. S-15 moves from after M4 to M3.
    3. **The rest** keep their placement and are ordered at the end of M3.

- **D-152 (owner) — Monitoring implications, M5 features, an extension of S-14, and study S-34.**
  - **Every `studies/README.md` entry gets a monitoring implication**, filled in when the study finishes: the monitoring metric, threshold or data check the result supports, if any. M5's monitoring design draws its metrics from finished entries and is not frozen until the relevant studies are done.
  - **M5 features (artifacts, not studies):**
    - what-if recalibration with hypothetical future periods, reusing the per-period surfaces (a hypothetical period is one more row);
    - a threshold table mapping next period's default count to tiered statuses (within tolerance / warning / threshold exceeded), each tier tied to a stated rule;
    - parameter-shock propagation to the 99.9% conditional PD, sharing S-23's machinery;
    - conditional PD across factor levels, and reverse factor stress, expressed in adverse and benign terms rather than raw signs of Z;
    - an estimator comparison, a view of S-8 and S-11.

    Names are neutral and describe what each feature computes.
  - **S-14 is extended** to the threshold table's operating characteristics: the false-alarm rate for a correctly calibrated model, and the detection delay for a miscalibrated one.
  - **S-34, sensitivity to severe new periods:** how far PD̂, ρ̂ and the 99.9% conditional PD move after one or two severe periods, as a function of T.
    - **Placement:** the first batch, after S-1 and S-2. It needs S-23's interval for the conditional PD and none of M3's estimators.
    - **Cost:** it is cheap: each added period is one more surface row. About 12–20 minutes on the subset.

## Study index moved, 2026-09-27

- **D-153 (owner) — The study index moves from `docs/methodology/questions.md` to `studies/README.md`, retitled "Studies".**
  - **Why:** the index and the studies' own files (`studies/<slug>/PREDICTION.md`, later results) belong together, and a study's directory now sits beside the entry that indexes it. The index lists studies, so it is named for them.
  - **Wording:** where the old text called the index or its entries "questions", it now says "study" or "studies". Each study's own research question stays in its entry (the **Question:** line and the index table's short-question column).
  - **What changed elsewhere:** path references only, in README.md, STATE.md, S-23's `PREDICTION.md` and earlier D-entries (D-148 to D-152). No earlier entry was reworded beyond its paths, and S-23's pre-registered predictions are unchanged: its two edited lines name the index's new path and nothing else.
  - **No code reads the index,** so the CI prose filter's `code_md` list is unchanged (`ci_path_filter_sync` confirms it).

## Studies S-1 and S-2, 2026-09-27

Numbered from D-300 so that entries made in parallel sessions do not collide.

- **D-300 (owner) — S-1 is redesigned as a surface-level property check, not a recovery study.**
  - **Why:** z → −z is an exact invariance of the likelihood, not a statistical property, so recovery replicates and coverage verdicts add nothing that a direct comparison does not show.
  - **Design:** a test-only mirrored objective (+√ρ·z, with its own hint) against parity on 34 fixed panels chosen for the hard cases (zero- and all-default periods, high ρ, estimates at or near a bound, n up to 10⁶). The per-period surfaces are compared cell by cell, then estimates, SEs, profile and bootstrap intervals and flags. A wrong-hint control must be detected.
  - **Prediction:** ε-level agreement with bounds derived from the rules' exact node symmetry and reversed summation order, and identical flags; any larger difference is a finding. Registered in `studies/z-sign-invariance/PREDICTION.md` before any run.
  - **The recovery-matrix version** is kept in reserve, run only if the check finds a difference that needs statistical characterisation.

## S-13 targeted: scope and pre-registration, 2026-09-27

- **D-154 (owner) — S-13's targeted run covers the borderline profile verdicts, and its results never overwrite the pinned R = 1,000 verdicts.**
  - **Scope:** the 38 profile-likelihood verdicts for PD or ρ whose R = 1,000 coverage lies within 0.010 of a band edge, in 31 scenarios, at R = 10,000 (replicates 0–9,999; 0–999 are the pinned ones). The profile interval is the recommended one and S-4 is judged against its findings. Borderline Wald and bootstrap verdicts would take in 66 of the 81 scenarios, which is the whole-matrix run that stays after M4. The other families of the same scenarios are reported at R = 10,000 without predictions.
  - **Results are S-13's own,** in `studies/recovery-r10000/`, with their provenance. At R = 10,000 the band is 0.9428–0.9572, and group B's mean profile coverage (0.937 for ρ, 0.940 for PD, pinned) is below it. So many verdicts that PASS at R = 1,000 are expected to fail at that precision. That is a narrower band detecting a smaller shortfall, not a contradiction, and the pinned verdicts stay as they are. Whether any reviewed list changes is a separate decision.
  - **Pre-registered** in `studies/recovery-r10000/PREDICTION.md` (Q1–Q7), merged before the run. The predictions shrink each R = 1,000 coverage towards its group's mean (empirical Bayes), because the verdicts were selected for lying near an edge. The expected conclusion, if they hold: profile intervals at T = 20 undercover by about one point systematically, not only in the six flagged scenarios. That is the evidence S-4a is to be judged against.
## The shared jackknife run (S-3, S-5, S-21): pre-registration, 2026-09-27

- **D-155 (owner) — One pre-registration for the three studies that share the jackknife pass, with their definitions fixed before any run.**
  - **Where:** `studies/jackknife-bias-rho/PREDICTION.md`, predictions J1–J14 (S-3: J1–J6; S-5: J7–J10; S-21: J11–J14), merged before the run.
  - **S-3:** ρ̃ = T·ρ̂ − (T − 1)·mean ρ̂₍₋ₜ₎ on the natural scale, both from the grid refinement, so the correction measures bias and not a method difference. ρ̃ outside the box is set to the bound and counted. The profile interval is shifted in logit ρ; an end truncated at the box stays there.
  - **S-5:** BCa in the logit coordinate, with z₀ from the existing B = 999 replicates and the acceleration from the jackknife. An interval with infinite z₀ or no jackknife spread is not computed, and counts as not covering (D-131).
  - **S-21:** leave-one- and leave-two-period-out influence in Hessian SE units, against the realised Z_t. The resolution check compares the refined leave-one-out estimates with exact off-grid maxima on replicates 0–4 of every scenario.
  - **Pinning:** one fitting run writes the results, with the reviewed verdict labels for the new families kept in a separate file, so a review never requires a refit. S-3 and S-5 are `native` options reported beside the parity verdicts, which do not change.

## Study S-1 result, 2026-09-27

- **D-301 — S-1 result: the calibration is invariant under z → −z within rounding; one reviewed finding, ρ not identified at n = 1.**
  - **Result** (`studies/z-sign-invariance/RESULTS.md`): on 34 panels, the mirrored objective's hints are exact mirrors and its surfaces agree with parity at the ε level (worst cell 0.24 of its bound), with identical flags everywhere. The wrong-hint control is caught.
  - **Downstream, on 33 of 34 panels,** estimates, SEs, profile and bootstrap intervals agree within the registered bounds. Three expected sizes were missed, each inside its bound: C3 by 2.3×, C8 by 27× (the Hessian's second differences, the size seen across platforms), C10 by 2.4×.
  - **The finding: panel 26, n = 1 in every period.** A single obligor's likelihood is PD or 1 − PD, free of ρ, so the surface is flat along ρ to rounding (2.8e-14) and rounding picks ρ̂ (0.262 against 0.001). The flags and the profile interval (the whole box) agree. It is non-identification, not asymmetry, and any rounding perturbation (platform, compiler, summation order) moves ρ̂ the same way.
  - **Pinned:** `study_z_sign_invariance` (CTest, label slow) exits non-zero on any failure other than C5, C7 and C10 in panel 26.
  - **The registered statistical follow-up is not run (owner, confirmed).** The prediction's verdict rule says a finding that moves estimates is characterised on the study subset. That step is skipped, for these reasons:
    - The mechanism is exact: ρ is absent from the likelihood only when every period has n_t = 1. It is a structural property of the data, not a rate that sampling could estimate.
    - It cannot arise in any recovery scenario, since every one has n ≥ 100. A recovery run would therefore measure a rate that is zero by construction.
    - The nine subset panels checked here (panels 1–9, replicate 0 of scenarios 29, 37, 72, 4, 68, 49, 7, 43, 51) pass every check.
  - **Mitigation, decided by the owner:** a structural flag when no period has n_t ≥ 2, with the numbers still reported. It is its own parity decision, D-302. A "flat to rounding" test is not added: that condition depends on the platform and stays with the existing flat-surface flag.
  - **Monitoring implication:** no metric for the sign convention, since calibration outputs do not depend on it; factor-scale outputs state it in adverse and benign terms. A data check reports ρ as not identified when no period has n_t ≥ 2, or when ρ's profile interval spans the whole box.

## Unidentified ρ flagged, 2026-09-27

- **D-302 (owner) — A fit whose panel has no period with n ≥ 2 is flagged `kFlagRhoNotIdentified` (`VCAL_FLAG_RHO_NOT_IDENTIFIED`, bit 6); its numbers are still reported. A parity behaviour change, from S-1's finding (D-301).**
  - **The condition is exact and structural.** A period with n = 1 contributes E[p(Z)] = PD or 1 − PD, free of ρ, so with no period of n ≥ 2 the data carry no information about ρ at all. It is a check on the counts: it does not depend on rounding or the platform.
  - **Flag, not refusal.** PD, ρ, the SEs, the log-likelihood and the profile intervals are reported as before. The flag says ρ̂ is meaningless, and ρ's profile interval is the whole box, truncated at both ends. A consistent output shape is easier for consumers than NaN or an error.
  - **"Flat to rounding along an axis" gets no second flag.** That softer condition depends on the platform, and stays with the existing flat-surface flag.
  - **Where it lives:** `BinomialMixture::rho_identified(obs, periods)` states the condition in core. `engine::calibrate` sets the flag whenever an objective declares that member, so the engine stays generic, and objectives without it are never flagged. Headline fits only: resampled replicates do not carry it, since `replicate_estimates` sees W × L and not the counts.
  - **ABI:** a new flag value is additive and breaks no caller. The minor version goes to 0.2, so a caller can tell whether the library may set the bit. The struct sizes are unchanged.
  - **No golden changes:** no recovery, replay or cross-reference panel lacks a period with n ≥ 2. The recovery replay, which compares flags exactly, passes unchanged. S-1's panel 26 now carries the flag in both conventions, and the pinned S-1 finding is otherwise unchanged.
  - **Tests:** `unit_engine` (an n = 1 panel is flagged with finite numbers; one period with n = 2 lifts the flag; the reference panel is not flagged) and `abi_c` (the same through the C ABI, with ρ's profile interval the whole box).

## Study S-2: pre-registration, 2026-09-27

- **D-303 — S-2's method is fixed in its pre-registration, before the table is computed.**
  - **Where:** `studies/sample-size-table/PREDICTION.md` (Q1–Q10) and the stdlib-only script `studies/sample-size-table/sample_size_table.py`, both merged before the run. It needs no new fits: only the committed recovery summary is read.
  - **Targets:** ρ within ±0.05 absolute, and PD within ±25% relative, at 95% (1.96 × RMSE).
  - **Method:**
    - RMSE ≈ C/√T, with C fitted in logs over T ∈ {20, 40, 100}, gives T*;
    - a free-slope OLS fit is the scaling check: slope in [−0.75, −0.30] and every residual ≤ 0.15 in log RMSE;
    - a scenario with more than half its replicates on the grid edge is left out of the fit, and a cell with all three left out is "not estimable at any T studied";
    - T* outside [20, 100] is marked an extrapolation.
  - **Evidence for the predictions:** large-n theory only, SD(ρ̂) ≈ √2·ρ(1 − ρ)/√T and relative SD(PD̂) ≈ λ(Φ⁻¹PD)·√ρ/√T. No RMSE value from the summary was read in writing them.
  - **The optional refit at an implied T is not planned.** It would be added under an appended prediction if asked.

## Study S-2 result, 2026-09-27

- **D-304 — S-2 result: the sample-size planning table.** The table, the comparison and the diagnosis are in `studies/sample-size-table/RESULTS.md`.
  - **Scaling:** RMSE ∝ T^(−½) holds in all 24 cells where it could be checked. Two cells are unchecked and one is not estimable, all at PD 0.1% and n = 100.
  - **ρ within ±0.05:** about 35–40 years at ρ = 0.12 with n·PD ≥ 50, rising to 180–225 at n·PD = 1. More than 100 years at ρ = 0.24 everywhere studied.
  - **PD within ±25%:** 6–7, 38–39 and 82–86 years at 5% PD for ρ = 0.02, 0.12 and 0.24.
  - **The large-n formula** is close at n·PD ≥ 50. Finite n multiplies the years by 1.6–2 at n·PD = 10 and by 5–7 at n·PD = 1.
  - **Against the prediction:** Q1, Q2, Q3, Q5, Q7 and Q10 held; Q4, Q6, Q8 and Q9 missed. The common cause of Q4, Q6 and Q9: with few defaults, the estimators do better than predicted when ρ is low, since PD is then close to the pooled binomial estimate. Q8 missed in size (188 years against ≤ 170).
  - **Monitoring implication:** a data-sufficiency check printed beside every calibration: the expected half-width 1.96·C/√T and the years each target needs, with the shortfall stated when the history is shorter. No year-on-year threshold follows from it; that is S-34's and S-21's question.

- **D-305 (owner) — S-2's write-up marks its extrapolations, and adds a relative target for ρ as a post-hoc view.**
  - **Extrapolations:** the T^(−½) law was verified at T = 20, 40 and 100 only. In the grids, figures beyond 100 years are marked † and figures below 20 are marked ‡, each with a footnote; cells where the law could not be checked are marked §. D-304's figures above 100 years, namely 180–225 years for ρ at n·PD = 1, 110–120 at ρ = 0.24, and Q8's 188, are extrapolations. What is observed is that those targets are not met at T = 100.
  - **ρ within ±25% relative,** beside ±0.05 absolute. ±0.05 is about ±20% of ρ at 0.24 but ±40% at 0.12, which partly explains why high-ρ cells looked hard. Under the relative target the order reverses: 76–83 years at ρ = 0.24, 95–111 at 0.12 and 124–183† at 0.02, for large portfolios. That matches the large-n relative SD √2·(1 − ρ)/√T.
  - **Not pre-registered, and labelled post hoc.** With the slope fixed at −½, T\* scales as 1/target², so the relative grid follows from the already published C by arithmetic (T\*(relative) = T\*(absolute) × (0.05/(0.25·ρ))²). There was nothing left to predict blind. No new fits; the registered prediction is unchanged.

## Study S-23 result: intervals for the 99.9% conditional PD, 2026-09-27

- **D-156 (owner) — S-23 is pinned: q's profile-likelihood interval is the recommended interval for the 99.9% conditional PD; its delta-method Wald and bootstrap percentile intervals are reviewed findings. Seven predictions held and four did not.**
  - **What was added:** `engine/conditional_pd.hpp` computes q, its logit gradient, and a profile interval along q = c. The inner maximum is over ρ, bracketed by the grid surface. The ends are walked in logit(q) and solved against the objective. Two flags: truncated at q's limit in the box, and box-limited. The recovery harness adds three verdict families for q, with the same band and policy as PD and ρ: profile (over all replicates), delta-method Wald in logit(q) (on the unflagged replicates, DEFERRED as for PD and ρ), and the bootstrap percentile of q at the existing B = 999 replicates. `validation/scipy/conditional_pd_profile.py` cross-checks the end points on the 162 replay panels (`replay_panels.csv`, new).
  - **Nothing existing changed:** every pre-existing cell of `summary.csv` (66 columns × 81 rows) and `replay.csv` (18 × 162) is bitwise identical after the pinning run. The profile code's panel log-likelihood helpers were shared with the new header with the arithmetic unchanged.
  - **Result:** profile 64 PASS, 15 CONSERVATIVE (box-limited, near-uninformative data), 2 KNOWN FINDING (72 and 74, T = 20). Delta-method Wald 32 PASS, 18 KNOWN FINDING (17 below the band, biased low; 13 above it, SE overstated), 31 DEFERRED. Bootstrap 8 PASS, 73 KNOWN FINDING (each inheriting ρ's diagnosis in its scenario). Two Wald diagnoses are new, for q. One is low estimates with narrow intervals: q̂ carries ρ̂'s downward bias, and its delta-method SE moves with it (error and SE correlate at +0.64 to +0.88), so the lowest estimates get the narrowest intervals. The other is overcoverage from an SE that overstates the spread. The first diagnosis was written as "location, not width" when the review began. It was corrected before merging: removing the mean bias exactly (an oracle shift of the pinned estimates) lifts coverage only to 0.916–0.947. The reviewed labels are unchanged, and the results page was regenerated from the pinning run's saved fits (`--results-md`), whose summary reproduces the committed one byte for byte. The study's entry in `studies/README.md` gives the full result, the comparison and a mitigation list.
  - **Predictions:** 7 held (P1, P2, P6, P7, P8, P10, P11) and 4 did not (P3, P4, P5, P9); S1, on the subset, held (8 of 9). Every miss is recorded as a miss, including the favourable ones: the predictions underestimated how well profile intervals perform for q.
  - **Timing of the pre-run files (UTC−4).** Times are commit author dates, which the rebase onto main kept. The two commits were pushed to the study branch as they were made, as 1acd449 and 480e050; the rebase gave them new hashes.
    - `PREDICTION.md` was committed at 16:12 and revised after review at 16:27, both in #3. Its last change was D-153's path-only edit at 16:55. SHA-256 `5a8e31a6cfce3badc8ba401f5629ae67ed861e81cc7a85777a562f21ab03085b`.
    - The study-subset exploration started at 17:07.
    - The scipy tolerance `TOL_SCIPY_Q_PROFILE_ENDPOINT_S`, with the implementation, was committed at 17:30:53. The first full-matrix run started 9 s later, at 17:31:02.
    - `compare.py`, which scores P1–P11, was committed at 17:32:22, 80 s after that run started. The harness prints nothing until it exits (at 18:58), so no result had been seen. Where the prediction's wording admits two readings, the script takes the stricter one. SHA-256 `320fe2e121cae6e02a86fd3a069cab7e8d0563edf1d87f503187f40f35b800ab`.
  - **An acceptance failure, fixed before pinning:** the first full-matrix run failed "every profile interval contains q̂" in 11 of 81,000 fits. All were estimates on the box's corner (PD 10⁻⁴, ρ 10⁻³), whose truncated lower end was 1 ulp above q̂: it was formed through logit and back, with the corner's ρ recovered from the scan coordinate. A truncated end is now q at the limit point evaluated directly, and a limit at an end of the ρ axis uses the bound itself. A regression test covers two of the replicates. Coverage was unaffected, since the true q is never the corner's. The reviewed lists come from that run; the pinning run (`--write` at the fixed commit, 5,165 s) reproduced every verdict.
  - **Process change for later studies (owner):** a pinning run should be a `--write` run from the start, with reviewed verdict labels kept in a separate file, so that review never needs a refit (adopted in D-155 for the shared jackknife run). S-23 needed two full runs.
  - **scipy pitfalls, recorded for validators** (methodology note §8): Nelder–Mead stalls on a bound of the box when the maximum lies just inside it, and bounded Brent stops about √ε·|x| short of a bound. Before the script handled them they cost 10⁻⁷ to 10⁻⁶ in logit(q).

## S-1's bounds on the Windows libm, 2026-09-27

- **D-157 (owner) — Two of S-1's check bounds are set from the worst measured platform, and C6 joins the panel-26 review. Recorded as a post-registration change: S-1's prediction file is unchanged.**
  - **What failed:** S-1's test (`study_z_sign_invariance`, label `slow`) passed on Linux (glibc), where S-1 was pinned (D-301). It failed on Windows, where both MSVC and MinGW GCC use the UCRT libm, with the same figures on both. PR CI runs no slow tests, so the nightly and dispatch runs caught it; `main` was red there on Windows.
    - C8, panel 20 (n = 10⁵, ρ = 0.24): the correlation differs by 6.95·10⁻¹⁰ against a bound of 5·10⁻¹⁰.
    - C11, panel 20: the quadrature check's value differs by 3.64·10⁻¹² against an absolute 10⁻¹².
    - C6, panel 26: two bootstrap replicates' flags differ.
  - **Why these are rounding:**
    - **C8:** its bound was `TOL_RECOVERY_REPLAY_SE_REL`, measured on panels up to n = 10⁴. The Hessian's finite differences cancel more at n = 10⁵.
    - **C11:** it compares the check's own value, a difference of two rounded log-likelihoods. The engine treats such a difference as rounding below 64 ε times the term size (D-120), and an absolute 10⁻¹² ignores the term size at n = 10⁵.
    - **C6:** on panel 26 the bootstrap flags follow the argmax along ρ, which rounding decides because ρ is not identified. That is exactly D-301's reviewed finding, which already covers C5, C7 and C10 there.
  - **The change:** C8's bound becomes `TOL_ZSIGN_SE_REL` = 2·10⁻⁹ (observed, doubled and rounded up), and C11's becomes 10⁻¹² + 64 ε × the largest term size at the estimate. C6 is allowed to differ on panel 26 only. On MinGW GCC, C8 is now at 0.35 of its bound and C11 at 0.005. The verdict is again "FINDING, reviewed (D-301)" on every platform, and the test exits 0.
  - **Why not rewrite S-1:** by the prediction's rule, any excess is a finding. So these are recorded as platform findings, reviewed as rounding, and S-1's `PREDICTION.md` and D-301 stand as written.

## Study S-13 targeted result, 2026-09-28

- **D-158 (owner) — S-13's targeted run: profile-likelihood intervals at T = 20 undercover by about one point, systematically. The pinned R = 1,000 verdicts are unchanged (D-154), and the finding is S-13's own.**
  - **What was run:** R = 10,000 on the 38 borderline profile verdicts in 31 scenarios. The band at that precision is 0.9428–0.9572. Replicates 0–999 reproduce the pinned fits field for field. The run took 3 h 42 min on the development machine; `studies/recovery-r10000/MANIFEST.json` has the commit, command and hashes.
  - **Result:**
    - all 11 targeted group B (T = 20) verdicts are below the band (0.932–0.942), 7 of them PASS at R = 1,000;
    - the six small-T findings are real but not special;
    - q's profile interval shows the same shortfall at T = 20;
    - at T = 40 the shortfall is about half a point, for PD only;
    - near the bounds, the conservative verdicts stay above the band, and six PASS verdicts at the upper edge join them.
  - **Predictions:** 4 of 7 held (Q2, Q5, Q6, Q7). Q1 (by 0.0002), Q3 (41/PD inside the band) and Q4 (72/ρ at 0.9399, above its range) missed in size, not direction. The empirical-Bayes shrinkage was, if anything, slightly too pessimistic about the flagged six.
  - **Consequences:**
    - S-4a is judged against a systematic shortfall at T = 20, not against six findings;
    - STATE's description of the profile interval now says that at T = 20 its coverage is about 94%.
    - Whether the pinned reviewed lists should change is a separate decision, not taken here.
  - **Per-replicate data are committed** as Parquet (`fits_r10000.parquet`, 11.6 MB): each replicate's PD and ρ estimates with their profile interval ends and flags, exact doubles. So every S-13 figure can be recomputed without the 162 MB fits file (owner's suggestion).
  - **Noted, harmless:** `--save-fits` writes the Fit record's padding bytes, which are uninitialised, so two saved-fits files of identical fits can differ byte for byte. Field-by-field comparison is unaffected, and that is how the provenance check was made.

## The shared jackknife run (S-3, S-5, S-21) result, 2026-09-28

- **D-159 — The shared jackknife run is finished. BCa replaces the percentile interval for ρ; jackknife bias correction is a point-estimate adjustment that needs exact delete-one refits; bias correction plus Wald does not rescue q; period influence must be measured, not inferred. 11 of 23 predictions held.**
  - **What was run:** the registered arms on the full matrix (81 × 1,000, 3 h 42 min), then the polished arm on the study subset (9 × 1,000, 2 h 01 min), after its registration (J19–J23, #18) had merged. `studies/jackknife-bias-rho/MANIFEST.json` has the commits, commands and hashes; the run commits were rebased onto main afterwards with no code change. The scorer `compare.py` was committed after the full run started and before any of its results existed. The polished run's registered-arm columns equal the full run's, all 45, on its 9 scenarios.
  - **S-3:** the registered arm computes the delete-one estimates from the 3 × 3 grid refinement. It lost up to 0.143 of coverage and broke 24 ρ profile PASS verdicts, because ρ̃ multiplies the refinement's delete-one error by T − 1: badly at ρ = 0.02 (gaps up to 0.97 SE), and at T = 100 for every ρ. With exact delete-one maxima (the polished arm) the shifted interval is within 0.008 of parity in all 6 informative subset scenarios, the bias is removed, and the RMSE is 0–5% above ρ̂'s. Near a bound (group A) setting ρ̃ to the floor still costs about 0.07 even when exact.
  - **S-5:** BCa beats the percentile interval for ρ in all 50 informative scenarios, and is in the band in 20 of 21 at T = 100 (percentile: 6). For PD it changes nothing (26 of 50 in the band, as before): the PD percentile interval is short on width, which BCa does not correct; the normal-quantile and (T − 1)/T narrowness alone account for about half the shortfall at T = 20.
  - **S-21:** the largest leave-one-out change in ρ̂ is 0.4–0.7 SE (median), as predicted, and pairs add 1.65–1.87 ×. The most influential period is the one with the most extreme factor in only 25–80% of replicates, rising with the expected defaults per period (median 0.34 at n·PD ≤ 5, 0.75 at 500).
  - **q's Wald sub-arms (S-23's mitigation 3):** neither rescues the Wald interval. A symmetric interval with a steady SE of the right size would cover (0.944–0.957 outside ρ = 0.02), so skew is not the problem. The delta-method SE is slightly small, noisy and tied to the error; the jackknife SE is right on average but twice as noisy, which alone costs 2–4.5 points. S-23's mitigation now reads *tested and failed*, and profile-likelihood intervals are the only supported method for q, with the parametric bootstrap (S-10) the one alternative untested. The diagnosis is from `diagnose.py`, written after the results.
  - **Predictions: 11 of 23 held** (J1, J4, J7, J10, J11, J13, J17, J19, J21, J22, J23). Of the 12 misses, 7 are in size (J5, J6, J9, J12, J14, J18, and J20, which missed favourably) and 5 in direction or in an effect that did not appear (J2's T = 100 half, J3, J8, J15, J16). The direction misses in S-3 trace to the grid refinement: the polished arm behaves as J2 and J3's mechanism described, which explains the misses without rescoring them. Those for q come from the SE's noise and its link to the error, which the mechanism did not foresee. So the predictions were better at mechanism than at size, but not only at size.
  - **Review labels relabelled before commit, after the polished arm:** the 13 shifted-interval verdicts at ρ ≥ 0.12 first labelled "correction noise" are `refinement_resolution` (the polished arm brings 49 and 51 back to parity); the q sub-arms' labels follow the diagnosis above (`se_tracks_error`, `jackknife_se_noisy`, and `small_rho_resolution` at ρ = 0.02). No label was committed in its first form.
  - **Consequences:**
    - the pinned parity verdicts are unchanged: S-3 and S-5 are `native` options reported beside them;
    - a `native` jackknife bias correction, when built, must use exact delete-one refits, and is reported beside ρ̂, never as the centre of an interval;
    - a `native` BCa interval for ρ is recommended over the percentile interval; the percentile interval stays the parity bootstrap interval;
    - the full matrix was not polished, so the S-3 conclusion rests on the 6 informative subset scenarios.

## Study S-34: pre-registration, 2026-09-28

- **D-160 — S-34 is registered before any run, with its method fixed: added periods at the median default count of a fixed adverse factor level, exact off-grid estimates, and predictions from the model's large-n limit.**
  - **Method:** each recovery panel is extended by one or two periods at a 1-in-100 or 1-in-1,000 adverse factor level, each with the median default count at that level, so the added periods are the same for every replicate. PD̂, ρ̂ and q̂ are the exact off-grid maxima of the original and extended panels, since D-159 showed the grid refinement's error matters in differences. Shifts are reported in the original panel's SE units and in relative terms, with the extended panel's profile intervals, their coverage of the unchanged truth, and how often the new q̂ lies above the original interval.
  - **Mechanism:** in the large-n limit q̂'s shift in SE units depends only on T and the added periods, not on PD or ρ: about (3.6–3.8)/√T SE for one 1-in-100 period and (6.3–6.7)/√T for one 1-in-1,000. `large_n_reference.py` computes the reference from the model alone, reading no recovery panel. Finite n should dilute the shift where n·PD is small.
  - **Predictions K1–K10:** direction, scaling with T and with the number of periods, size against the reference where n·PD ≥ 100, dilution where n·PD ≤ 10, exceedance of the original q interval, coverage of the unchanged truth (predicted to rise after one 1-in-100 period, and to fall after two 1-in-1,000 periods at T = 20), relative size, and three scenarios (0, 3, 6) where the "severe" period has a median of zero defaults.
  - **Process:** the subset first, then one full-matrix run that pins the results, scored by a script committed before that run's results exist. S-34 is descriptive and adds no verdict family.

## Study S-34 result, 2026-09-28

- **D-162 — S-34 result: after a severe year q̂ jumps by an amount predictable from T alone (about 3.7/√T SE per 1-in-100 year, 6.5/√T per 1-in-1,000), and the refitted interval moves up rather than widening, so it covers the unchanged truth far less often. 8 of 10 predictions held.**
  - **What was run:** every recovery panel refitted after adding one or two periods at a 1-in-100 or 1-in-1,000 adverse factor level (median default count), with exact off-grid estimates and the PD, ρ and q profile intervals: the subset first (35 min), then the full matrix (81 × 1,000 × 4 variants, 6 h 12 min). `studies/severe-period-sensitivity/MANIFEST.json` has the commits, commands and hashes. The original panels' q coverage reproduces the pinned recovery summary in all 81 scenarios.
  - **Result:**
    - the jump matches the model's large-n limit: 0.93–1.02 × the reference for q where n·PD ≥ 100, and 0.61–1.33 × (median 0.98) where n·PD ≤ 10;
    - the refitted interval's lower end rises more than its upper end (for ρ in every B–D scenario and variant), so coverage of the unchanged truth falls: to 0.81–0.89 after one 1-in-1,000 year at T = 20, and to 0.26–0.51 after two;
    - the refitted q̂ exceeds the original interval in at most 6.3% of replicates after one severe year, and 28–45% after two 1-in-1,000 years at T = 20;
    - in 10 group A scenarios (n·PD ≤ 1) a "1-in-100" year lowers q̂.
  - **Predictions:** 8 of 10 held (K1–K4, K6, K7, K9, K10). K5 missed: the predicted dilution at small n·PD did not appear. K8 missed in size in its first part: after one 1-in-100 year coverage stayed level rather than rising. The subset had shown both before the full run; the predictions were not amended.
  - **Consequences:**
    - no estimator change;
    - M5's what-if recalibration reports the expected jump beside the refit;
    - the monitoring implications in the S-34 entry (the expected jump, the breach of the previous interval as a review trigger, and post-crisis refits reported with the number of severe years) go into M5's monitoring design (D-152).
  - **Per-replicate data are committed** as Parquet (`replicates.parquet`, 28 MB, exact in all 7,695,000 cells). That is larger than earlier studies' because the nine double columns per row (estimates and interval ends for five panels per replicate) do not compress.
## Phase 2 pre-registrations: S-8, S-9, S-15, S-27, S-28, 2026-09-28

- **D-161 — S-8, S-9, S-15, S-27 and S-28 are registered before M3's estimators exist, with the estimators defined in the registrations for M3 to implement.**
  - **Files:** `studies/mle-vs-mom/PREDICTION.md` for the estimator-comparison pass (S-8, S-27, S-28; D-150; E1–E13), with pointer files for S-27 and S-28; `studies/bayes-coverage/PREDICTION.md` (S-9; F1–F5); `studies/bayes-sbc/PREDICTION.md` (S-15; G1–G5).
  - **Definitions fixed there:** MoM in the joint-default-probability form, exact in finite n; the Vasicek-rate MLE with a profile-likelihood interval, and four zero-default treatments: refuse, drop, a true censored likelihood (the Vasicek rate CDF below the detection limit c_n = 1/(2n), instead of the density) and substitution (half a default, a continuity correction), the last two kept distinct at the owner's request; "indistinguishable" for S-27 (0.1 SE in the mean, 0.0227 in coverage); the grid posterior with flat and Jeffreys priors and cell-uniform marginals; SBC's draws and continuous rank statistic. If M3 has to define an estimator differently, the change is appended to the registration, dated, before the study runs.
  - **Evidence:** model-only reference scripts (`mle-vs-mom/model_reference.py`, `bayes-coverage/grid_spread.py`) and pinned results; no estimator of these studies has run on any panel.
  - **Design point found while registering:** on the parity grid (61 × 41) the estimator's spread is under half a grid spacing for PD in 8 scenarios of groups B–D, so a grid posterior there is set by the grid, not the data. At the owner's direction this becomes a rule of the estimator rather than a fixed finer grid: **at least 4 grid points per posterior SD on each axis** (logit coordinates), with local refinement, or a refusal after three failed refinements. 4 keeps a normal posterior's 95% interval within 0.001 of its mass at any grid offset. S-9 and S-15 run the estimator with the rule and, as a diagnostic, the parity grid without it, whose failure (overcoverage, and a hump in the SBC ranks) they predict in advance.
  - **E9 is an implementation check:** the refusal share follows exactly from the model, so the file says a miss there means a bug, not a surprise.
  - **Headline mechanisms:** MoM loses efficiency for ρ by a factor that grows with ρ, with lower PD and with T (1.04–2.68 in the large-n limit); the rate MLE's ρ̂ is biased up by binomial noise, about 1/n; the parity rate estimator refuses exactly 1 − (1 − P₀)^T of panels, which is most of the matrix at n ≤ 1,000.

## S-34 follow-up, per-replicate data policy, and the CI gate, 2026-09-28

- **D-163 (owner) — S-34's coverage drop is framed as selection; its group A oddity is a design artifact; per-replicate data files leave the repository; the required CI gate is reported by the pull_request run only.**
  - **Selection, not interval failure:** S-34 conditions on the added periods being severe, so the refit is pushed up by construction; lower coverage of the unchanged truth *given* a severe year is expected of a calibrated interval, whose coverage is an unconditional property (and reproduces as pinned on the original panels). The S-34 entry now says so, so it is not read as profile intervals breaking after a crisis.
  - **Group A oddity, checked:** where n·PD ≤ 1 the median default count at the "1-in-100" factor level is 0 (scenarios 0, 3, 6: PD̂ falls) or 1–2 (1, 4, 7, 15, 24, 30, 33: PD̂ rises but ρ̂ falls, because a 1–2 default year is milder than the periods that set ρ̂ at those counts). Both are artifacts of defining severity by the factor level with the median count; M5's what-if recalibration states its definition of severity.
  - **Per-replicate data are not committed from now on.** The runs are deterministic, so per-replicate files are regenerable from the committed code, seed and command. A study commits its summary tables and the per-replicate file's SHA-256 in its `MANIFEST.json`, and attaches the file to a GitHub release when it should be downloadable. S-34's 28 MB Parquet is removed from the tree and attached to the release `s34-data-2026-09-28` (https://github.com/keithcsqueensu/vasicek-cal/releases/tag/s34-data-2026-09-28); the manifest records its SHA-256 and how to regenerate it. It remains in the history of #21; the history is not rewritten. The earlier studies' committed Parquet files (S-13 11.6 MB, the jackknife run 5.8 MB and 0.8 MB) stay where they are.
  - **The CI gate:** a push to a pull request's branch starts a reduced run as well as the full pull_request run, and both reported a check named `ci-ok`. The reduced run's gate passed first and satisfied the required check while the full run's scipy job was still running, so #19 and #21 merged before it finished (it passed both times). Only the pull_request run now reports `ci-ok`; the other runs report `ci-ok (push)`, `ci-ok (schedule)` or `ci-ok (workflow_dispatch)`. The PR run's `ci-ok` already required every job, scipy included, so "`ci-ok` passed" now means every job of the full run passed.
## M3.1: the Vasicek-rate MLE, 2026-09-28

- **D-164 — The Vasicek-rate MLE is an objective of the existing engine, closed form per period, with its zero-rate treatment a type parameter; parity refuses rates of 0 or 1.**
  - **What:** `core/objectives/vasicek_rate.hpp`. The log density of an observed default rate, l = ½ log(1 − ρ) − ½ log ρ + x²/2 − (√(1 − ρ)x − c)²/(2ρ), x = Φ⁻¹(r), runs through the same surface, refinement, profile and q-interval code as the binomial objective. Its treatment of rates of 0 and 1 is part of its type, so `panel_error` can see it: `Refuse` (parity; `zero_rate_periods` names the periods, D-044), `Censor` (the Vasicek CDF below the detection limit), `Substitute` (a continuity correction, stated as a data edit), and `drop_zero_rate_periods` as an explicit edit before the fit. For counts the detection limit is 1/(2n). These are the definitions registered for S-28 (D-161).
  - **Why the grid path and not the closed form:** one mechanism for every objective (including M2b's W × L resampling), and the censored likelihood has no closed form. The closed form, the normal MLE of x, is kept as an independent check: the polished maximum agrees with it within 1.6·10⁻⁷ in logit units over 27,000 fits.
  - **Validation:** mpmath goldens at the binary64 inputs (`tools/gen_vasicek_rate_goldens.py`, in the CI's generated-files check), within 2.9 ε of the terms' size; a scipy replication (`validation/scipy/vasicek_rate_mle.py`, in the scipy CI job) on 63 panels including 9 censored count panels, within 1.8·10⁻⁷ in the estimates and 2.3·10⁻¹⁰ in the profile ends; four new tolerances, set from observed values.
  - **Recovery on its own model:** the n → ∞ limit of the recovery panels (the same factor draws), 27 cells × 1,000. Predictions V1–V5 were committed before the run; 4 held. The one miss (V5) is 1.8–2.2% near-bound flags at ρ = 0.24, T = 20, against ≤ 1% predicted: ρ̂ within two SEs of the 0.5 cap. 78 of 81 coverage verdicts are in the band; the three below (29/PD, 68/PD, 74/q) use the factor draws of pinned binomial small-T findings, so the dips belong to those draws at small T, not to the binomial likelihood. Reviewed in `unit_vasicek_rate`.
  - **Not in this step:** the C ABI (still v0, binomial only; exposing M3's estimators needs an ABI minor version, decided separately), and bootstrap intervals for the rate objective.

## M3.2: the method of moments, 2026-09-28

- **D-165 — The method of moments is an engine function over per-period sufficient statistics, joint-default-probability form, exact in finite n; its joint default probability is the n = 2, d = 2 binomial-mixture integral, so it shares the parity integrator. It is validated on panels of its own seed, not the recovery panels, so S-8's comparison is first computed in S-8's registered run.**
  - **What:** `engine/moments.hpp`. PD̂ = Σd/Σn and PD̂₂ = Σd(d − 1)/Σn(n − 1), weighted by a W row when resampling (D-042); ρ̂ solves E[p(Z)²](PD̂, ρ) = PD̂₂ by Brent in logit(ρ) over the box. No defaults: refused; PD̂₂ outside what the box's ρ range gives: ρ̂ at the floor or the cap, flagged; PD̂ outside the PD box: flagged. A rate-series form (mean rate, mean squared rate) for D-046's rate data.
  - **Why the n = 2 integral:** Φ₂(c, c; ρ) = E[p(Z)²] is exactly the binomial-mixture integral with n = 2, d = 2, so no bivariate-normal routine is added to core, and the moment uses the same validated integrator as the likelihood.
  - **Validation:** mpmath goldens (tanh-sinh at 50 digits; `tools/gen_moments_goldens.py` in the generated-files check), within 7.1·10⁻¹⁵ in log PD₂; inversion within 1.3·10⁻¹³; edge cases and weights; consistency on a 20,000-period panel; 108 reference panels (count and rate series, the recovery matrix's cells at n = 1,000, seed M3MOMVAL) reproduced on MSVC and GCC and replicated by `validation/scipy/method_of_moments.py` with Owen's T (exact Φ₂ at equal arguments), within 1.6·10⁻¹¹ in logit(ρ̂).
  - **Why its own seed:** S-8 compares MoM with the MLE on the recovery panels and is scored by a script committed before its run (D-161). Running MoM over the recovery matrix here would compute S-8's result first. Likewise M3.1's replay touched nine recovery count panels (replicate 0, PD 1%, n = 100) under the censored rate likelihood; that is disclosed for S-28, whose run covers 81,000 panels.
  - **Not in this step:** intervals (no likelihood; bootstrap via W, not validated here) and the C ABI.

## M3.3: the grid-Bayesian estimator, 2026-09-28

- **D-166 — The grid-Bayesian estimator carries the resolution rule (at least 4 grid points per posterior SD on each axis), refining on local grids whose extent follows the posterior's SD and tails; the registered extent (±8 of the smaller of the Hessian SE and the posterior SD) is amended before S-9 runs, because it dropped real mass for skewed posteriors.**
  - **What:** `engine/posterior.hpp`. The posterior on the grid in logit coordinates, flat (natural scale) or Jeffreys priors (√det I_u of one period, tabulated per n and interpolated), cell-uniform marginals, the resolution rule with local refinement or refusal, equal-tailed and HPD intervals. `docs/methodology/grid_bayesian.md`.
  - **The amendment, and why:** on the estimator's own validation panels (seed M3BAYVAL, T = 20, n = 100), the registered extent refused 3 of 16 fits: the ρ posterior is wide and skewed towards the box's floor, so ±8 Hessian SEs dropped posterior mass. The extent now follows the posterior SD and the marginal's 5·10⁻⁸ tails (snapped to cell edges, so the mass check and the quantiles agree), and the spacing follows the smaller of the SD and the SE, at 5 points per SD. Recorded as an addendum to `studies/bayes-coverage/PREDICTION.md` before any S-9 run; F1–F5 are unchanged. No recovery panel was used.
  - **Validation:** exact against a normal posterior (equal-tailed ends 0.009 SD, HPD 0.071 SD, both within bounds derived from the rule); the Jeffreys table against scipy's adaptive quad (1.45·10⁻⁹ in log); 16 reference fits replicated by `validation/scipy/grid_posterior.py` (3.0·10⁻⁵ in logit, the rule's decisions exactly); thread-count determinism.
  - **M3 is complete with this step:** the Vasicek-rate MLE (D-164), the method of moments (D-165) and the grid-Bayesian estimator (D-166). Still open, and decided separately: exposing them through the C ABI (a minor version), and the registered studies S-8, S-9, S-15, S-27 and S-28 that run on them.
## Spike: surface rows cached by (n, d) across replicates, 2026-09-28

- **D-167 (owner; approved after the spike) — Per-period surface rows depend only on the observation and the configuration that evaluates it, so they are computed once and shared across panels (replicates, scenarios with the same n, bootstrap panels), extending D-122. Bit-for-bit identical results; on the recovery subset the surface work falls 50× and the whole pass runs 2.1× faster.**
  - **The owner's two conditions, as built (`engine/surface_cache.hpp`):**
    - **The key is everything that determines a row,** not (n, d) alone: the objective's type (its precision policy and, for the rate objective, its zero-rate treatment) and its state if it has any (an objective with state must declare `cache_signature`, enforced at compile time); the integrator's type and every table and count, N included; each grid axis's bounds (exact bits), point count and scale; then the observation's bytes. `unit_surface_cache` changes each of these in turn (both bounds, spacing and scale of each axis; Gauss–Hermite N, the composite rule's panels and points, the doubled check rule; the precision policy; the rate objective's censored and substituted treatments; the detection limit) and asserts a miss whose row equals the direct evaluation.
    - **Rows are filled on demand,** the first time an observation appears, with a once-per-row guard so concurrent fits compute a row once: no precomputed range and no fixed bound, which serves S-27's large n.
  - **Adopted in:** `recovery::fit` (optional cache), the recovery harness (one cache per run: scenarios with the same n share rows), `study_jackknife_run` and `study_severe_period`. `engine::calibrate_cached` is `calibrate` with the surface from the cache.
  - **The full matrix through the cache:** `recovery_harness --check` (81 × 1,000, MSVC) reproduces the committed golden summary byte for byte (the same SHA-256; 0 mismatches), computing 4,412 surface rows and reusing 1,744,128, in 2,266 s against the pinned run's 5,165 s (2.3× faster).
  - **Mechanism:** `engine::calibrate_from_surface` (calibrate's steps after the surface, on a surface supplied by the caller; `calibrate` is now the surface followed by it, unchanged in behaviour). A caller that simulates a scenario's panels first computes each distinct count's row over the grid once (parallel over rows × points) and assembles each replicate's L from copies. The profile, bootstrap and q-interval code read L unchanged.
  - **Identity (the spike, `tests/studies/surface_cache_spike.cpp`):** on the study subset (9 scenarios × 1,000), every field of all 18,000 fits, cached and uncached, equals `recovery::fit` bit for bit.
  - **Measured (MSVC 19.51, 24 threads):**

    | | Surface rows | Surface CPU s | Other CPU s (calibrate, profile, bootstrap, q) | Total CPU s | Wall s |
    |---|---|---|---|---|---|
    | Uncached | 147,885 | 5,380 | 4,531 (106, 2,180, 1,552, 693) | 9,911 | 447.5 |
    | Cached | 2,922 | 88 (86 build + 3 assembly) | 4,710 (109, 2,245, 1,643, 714) | 4,799 | 217.3 |

    The surface was 54% of the work and becomes 2%; what remains is per replicate: profile solves (47%), the bootstrap (34%), the q interval (15%). Distinct counts per scenario: 10 (n·PD = 1) to 2,074 (n = 10⁴, PD 5%); the largest cache is 2,074 × 2,501 doubles, 41 MB.
  - **Projected, from the measured split (estimates, not runs):**
    - **Full recovery matrix (81 × 1,000):** about 2× faster (the surface share on the subset is representative of the matrix's mix of T and n).
    - **S-13 full (81 × 10,000):** about 2×: distinct counts grow far more slowly than R, so the surface becomes negligible and the per-replicate profile and bootstrap work sets the time; the index's 12–21 h becomes roughly 6–10 h on CPU.
    - **S-10 (parametric bootstrap) and S-4b (its Bartlett factor):** the largest gain. Rows depend on (n, d), not on the parameters a panel was drawn from, so one scenario's cache serves every replicate *and* every bootstrap panel. A bootstrap panel's surface drops from about 0.6 CPU s to a copy and a reduction (about a millisecond), so a full-matrix S-10 (81,000 replicates × 999 panels, fits without profiles) goes from months of CPU to roughly a day of CPU, an hour or two of wall time: feasible on CPU, where the index placed it after M4 (GPU). For S-4b the likelihood at θ̂, fixed for a replicate's bootstrap panels, is likewise one integral per distinct count.
  - **Limits:** the cache needs the grid and integrator fixed within a scenario (true of every study so far), and memory grows with distinct counts × grid points; at n ≥ 10⁵ (S-27) the counts' range widens and a bounded or on-demand cache would be needed. Resampling by W × L is unaffected.
  - **If adopted:** the recovery harness and the study tools build a scenario cache before fitting; a fast test asserts cached and uncached fits agree bit for bit on a few scenarios.

- **D-168 (owner) — The roadmap after M3: S-10 and S-4b move from "after M4" to the CPU, and M6 comes before M4.**
  - **Why:** with surface rows shared across panels (D-167), a parametric-bootstrap panel's surface is a copy, and what remains of S-10 and S-4b is per-replicate profile and bootstrap work, which parallelises on the CPU. M4 (GPU) is no longer on the critical path of any registered study; M6 unlocks the multi-grade, misspecification and data-instability studies.
  - **Order:** (1) merge the M3 PRs and the cache; (2) C ABI 0.3 exposing the three M3 estimators, the rate MLE's zero-default treatments as an enum of explicit options with no recommended default until S-28 reports; (3) the M3 folded studies (S-8/S-27/S-28, S-9, S-15) with the cache; (4) S-10 and S-4b on the CPU; (5) a short consolidation (tag v0.2; README updated with the study results); then M6; M4 becomes a later performance milestone.

## C ABI 0.3: the M3 estimators, 2026-09-28

- **D-169 (owner) — ABI 0.3 exposes the three M3 estimators in one minor version, adds the `native` profile, and offers the rate MLE's zero-default treatments as explicit options with no recommended default until S-28 reports.**
  - **Functions (additive; no existing struct or function changes):** `vcal_calibrate_rate` (count data, rates d/n with detection limits 1/(2n), or a `vcal_rate_series` with its own limits), `vcal_calibrate_moments` (counts or rates) and `vcal_calibrate_posterior` (flat or Jeffreys prior). New structs `vcal_rate_series` (32 bytes), `vcal_moments_estimate` (32) and `vcal_posterior` (112), padding-free and checked from C.
  - **Zero rates:** `VCAL_ZERO_RATES_REFUSE` (0) is the parity behaviour of D-044, and its refusal names the periods in the error message; `_CENSOR`, `_SUBSTITUTE` and `_DROP` are native enhancements and need a `VCAL_PROFILE_NATIVE` context (parity refuses them with a message saying so). The enum has no recommended value: which to prefer is S-28's question, and when it reports the recommendation is documentation plus a default value, with no ABI change.
  - **The native profile:** `VCAL_PROFILE_NATIVE` (1) is accepted by `vcal_context_create`; every existing function behaves as under parity. Build info reports `profiles=parity,native`.
  - **Limits, explicit errors rather than silent behaviour:** the Jeffreys prior needs the same n in every period (`VCAL_E_UNSUPPORTED` otherwise, since the table is per n); the posterior needs logit axes.
  - **Tests:** `abi_engine_equivalence` checks every new call against the engine bit for bit (the rate MLE under each treatment, on counts and on a rate series; MoM on counts and rates; the posterior under both priors, including the context's cached Jeffreys table) and the refusals; `abi_c_test` checks the layouts and calls each function from C.

## Open

- **R-1 (revisit at M6) — GPU-side DGP.** Only matters for large parametric bootstraps (AR(1), D-043).
- **R-2 (revisit at M7, when the macro pipeline starts interacting with the estimator) — Z sign convention and macro sign filters.** Declare expected macro effects in economic terms ("worsens" or "improves credit conditions"). Map them to coefficient signs in one place, through the engine's Z convention, and test the mapping on a synthetic DGP. Tracked as study S-18 (D-148).
  - **The engine's convention:** p(z) = Φ((Φ⁻¹(PD) − √ρ·z)/√(1 − ρ)), in `core/model/vasicek.hpp` and `dgp/`. A higher Z therefore means *better* credit conditions (a lower conditional PD).
  - **The risk:** a macro variable that worsens credit conditions has a negative coefficient on Z, but a positive one on a PD or default-rate scale. A sign filter written against the wrong scale silently keeps the wrong specifications.
