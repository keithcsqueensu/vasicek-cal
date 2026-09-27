# vasicek-cal — project instructions

## What this is
An independent credit-risk research library: an open-source C++/CUDA engine
for calibrating Vasicek/ASRF PD models, with benchmarking, backtesting, and
monitoring. It is for practitioners, researchers and students, and is designed
so every result can be independently replicated. Correctness, traceability,
and replicability therefore outrank cleverness. It is a personal open-source
project.

## Hard constraints
- Public or synthetic data only (synthetic DGP fixtures, FDIC, FRED, StatCan).
- Language: C++17/20 + CUDA for the core. No Rust in the core; an optional
  Rust CLI lives in a separate repo and links the C ABI.
- Dependencies: CUDA toolkit + CMake only. No other third-party libraries in
  the core without discussion.
- Stable, flat C ABI (extern "C", plain structs, pointers + scalars).
- License: Apache-2.0.

## Architecture
- core/: host+device math policies (logPhi, conditional PD, quadrature)
- Policies composed as templates: Objective × Integrator × Reducer × Precision
- Surface engine: evaluate scalar surface over parameter grid, then reduce
- Per-period surfaces L (T × grid) computed once; bootstrap, walk-forward,
  and jackknife are weight matrices W, so surfaces = W × L
- CPU backend (OpenMP) is first-class for correctness and CI; GPU must
  match CPU within tolerance
- ref/: independent FP64 reference that must NOT share quadrature code
- dgp/: counter-based RNG (Philox), deterministic across CPU/GPU
- Profiles: `parity` (textbook, replicable) vs `native` (enhancements)

## Methodology scope
- Estimators: binomial-mixture MLE, Vasicek-rate MLE, method of moments,
  grid Bayesian
- Benchmarks: multi-segment correlated factors (two-stage + pairwise),
  AR(1) factor, alternative mixing distributions, Pluto–Tasche, window and
  influence analysis
- Macro: transform library with lineage, pre-screen filters, batched
  specification fits, pluggable scorers, config-driven ranker
- Monitoring: backtests, run ledger, generic monitoring report

## Working rules
- Numerics in log space; document every tolerance.
- Every feature ships with synthetic recovery tests and golden values.
- Anything someone must reproduce to replicate a result must be describable
  in plain prose and implementable in a short scipy script.
- Build fat binaries for sm_89 and sm_120 plus PTX; don't tune only for one GPU.
- Portability goals: supports CUDA >= 11.8 and older host compilers (GCC 11,
  MSVC 14.39) for broad compatibility; device-visible code stays C++17.
- Keep STATE.md (current status, next steps) and DECISIONS.md (decision +
  rationale) up to date at the end of every session.