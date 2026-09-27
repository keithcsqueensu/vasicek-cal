# vasicek-cal

An independent credit-risk research library: a C++/CUDA engine for calibrating one-factor
Vasicek / ASRF probability-of-default models from default counts. It is for practitioners,
researchers and students, and is designed so every result can be independently replicated: every
method is described in plain prose, reproduced by a short scipy script, and cross-checked against
an independent reference implementation, with every tolerance documented.

vasicek-cal is a personal open-source project by Keith Yip, developed independently.

> **Status: pre-alpha research software, not for production use.**
> - **M1, the parity sign-off for the binomial-mixture MLE, is complete** and tagged
>   `v0.1.0-parity-binomial`.
> - **M2 is done:** profile-likelihood intervals, the resampling engine and the C ABI v0.1.
> - **Not yet built:** the other single-factor estimators (M3) and the CUDA backend (M4).

## What it does today

- **Estimates PD and asset correlation ρ** by maximum likelihood from a panel of periods, each
  with n obligors and d defaults. The model is the binomial mixture of the one-factor Vasicek
  model, in the **`parity` profile**: textbook, deterministic, and reproducible in scipy.
- **Accurate, log-space integration** for every period shape (D-118):
  - adaptive Gauss–Hermite (128 nodes, centred on the integrand's mode) when 0 < d < n;
  - composite Gauss–Legendre on a sinh map for zero- and all-default periods.

  Every fit re-checks each period with the rule doubled and flags any disagreement.
- **Grid search, then refinement.** The log-likelihood surface is evaluated on a (PD, ρ) grid, so
  there are no optimiser convergence failures. The maximum is refined with a 2-D quadratic,
  including the PD–ρ cross term. Standard errors come from the Hessian at the estimate.
- **Profile-likelihood intervals** (95%, χ²₁). They are truncated at the search box and flagged,
  never extrapolated, which makes them the right interval near ρ = 0 or with very few defaults.
- **Resampling as weights (W × L).** The costly integrals run once per panel, and each replicate is
  a weighted re-sum of them, so 999 bootstrap replicates cost 15–30% of one fit. Schemes:
  - the iid and moving-block bootstrap, with their own random-number key domain;
  - jackknife and walk-forward;
  - index matrices supplied by the caller, to reproduce another tool's replicates exactly.

  Bootstrap percentile intervals are also computed.
- **Flags instead of silent fixes:**
  - grid edge;
  - flat surface;
  - refinement rejected;
  - quadrature unconverged;
  - numeric;
  - near a bound;
  - interval truncated.
- **A synthetic data generator** (Philox4x32-10 with its own elementary functions). Its panels are
  bitwise identical on every platform and match a pure-Python mirror.
- **Deterministic results.** They are bitwise identical for any thread count, and agree across
  platforms within documented tolerances.
- **A flat C ABI** (`include/vcal/vcal.h`, plain C99): calibrate with profile intervals,
  resample, surface, simulate and build provenance. It is tested from a C11 program, and every
  number it returns is bitwise what the engine gives.

## Accuracy and replication evidence

Each figure comes from a CI-enforced test, and each tolerance is in the
[tolerance register](docs/methodology/tolerances.md).

| Check | Result |
|---|---|
| Special functions (log Φ, erfcx, Φ⁻¹, log C(n, d), …) vs mpmath at 50 digits | at most 6 ulp |
| Period log-likelihood, parity rule vs mpmath: PD 10⁻⁶–0.5, ρ up to 0.99, n up to 10⁶, including d = 0 and d = n | within 11.7 ε |
| Period log-likelihood vs the independent `ref/` implementation, 4,296 periods | within 35 ε (term-scaled); worst absolute 2.3·10⁻¹⁰ |
| Estimates vs `ref/` | within 0.026 SE; SEs within 1.6% |
| Profile-interval endpoints vs `ref/` | within 1.5·10⁻¹⁰ in logit units; truncation identical |
| Plain scipy replication (`binom.logpmf`, `quad`, Nelder–Mead) | surface within 4.1·10⁻¹³ relative; estimates within 0.026 SE |
| Across platforms (Linux GCC 11/14 and Clang 18, Windows MSVC) | estimates within 4.4·10⁻¹³, SEs within 2.5·10⁻¹⁰ relative |
| Recovery: 81 scenarios × 1,000 simulated panels | RMSE falls with T everywhere; no quadrature flags in 81,000 fits |
| Profile-interval coverage, same 162 scenario × parameter verdicts | 136 in the Monte Carlo band; 20 conservative (above it); 6 slightly below at T = 20–40 |
| Bootstrap percentile coverage (iid, B = 999), same 162 verdicts, with predictions committed before the run | 37 in the band, 125 below it: undercovers for ρ at every T and for PD at T ≤ 40 |

Known statistical behaviour, all reported rather than corrected:

- small-T downward bias in ρ̂;
- Wald undercoverage at T ≤ 40;
- conservative profile intervals in near-uninformative settings;
- bootstrap percentile intervals that undercover: they break down near a bound and correct neither bias nor skewness.

**Profile-likelihood intervals are the recommended method for inference in this model.**
Bootstrap percentile intervals are provided for comparison and are not recommended for ρ.

These are documented with evidence in [recovery.md](docs/methodology/recovery.md), and the M1
sign-off is in [parity_signoff.md](docs/methodology/parity_signoff.md).

## Replicating the results

Start with the methodology note. It is self-contained.

| Document | Contents |
|---|---|
| [binomial_mixture_mle.md](docs/methodology/binomial_mixture_mle.md) | Model, likelihood, both integration rules and why fixed-node Gauss–Hermite fails at large n, validated ranges, the per-run check, materiality, estimation and flags, profile intervals and boundary theory, known statistical behaviour |
| [tolerances.md](docs/methodology/tolerances.md) | Every tolerance: value, observed error, rationale and the test that enforces it. Values live in `tests/tolerances.toml` |
| [recovery.md](docs/methodology/recovery.md), [recovery_results.md](docs/methodology/recovery_results.md) | Recovery harness design and the full results tables (bias, RMSE, Wald / t / profile coverage) |
| [reference_implementation.md](docs/methodology/reference_implementation.md) | The independent `ref/` implementation |
| [dgp.md](docs/methodology/dgp.md) | The synthetic data generator, precise enough to re-implement |
| [parity_signoff.md](docs/methodology/parity_signoff.md) | The M1 exit: requirements, deferred verdicts and their M2 re-assessment |
| [studies/README.md](studies/README.md) | Index of research studies and `native` variants: each study's question, experiment, pre-registered prediction, result, mitigation (D-148, D-153) |
| `validation/scipy/binomial_mixture_mle.py` | A 125-line scipy replication, run in CI |

[DECISIONS.md](DECISIONS.md) records every design decision and its rationale (D-nnn).
[ARCHITECTURE.md](ARCHITECTURE.md) describes the design, and [STATE.md](STATE.md) the current
status and the full plan.

## Build and test

Requirements:

- CMake ≥ 3.25 and Ninja;
- a C++20 compiler. CI uses GCC 11, GCC 14, Clang 18 and MSVC; the device-visible code is C++17;
- OpenMP, optional. Without it the CPU backend runs serially and gives identical results.

On Windows, configure from a Visual Studio developer shell.

```sh
cmake --preset cpu-release
cmake --build --preset cpu-release
ctest --preset cpu-release        # the cpu-debug preset skips the slow cross-reference tests
```

The presets are `cpu-debug`, `cpu-release` and `cuda-release`. `cuda-release` needs nvcc ≥ 11.8:
the project supports CUDA ≥ 11.8 and older host compilers (GCC 11, MSVC 14.39) for broad
compatibility, which is why device-visible code is C++17.
For now it compile-checks all device-visible code for sm_89 (plus sm_120 with CUDA ≥ 12.8) and
compute_80 PTX, and runs the CPU suite.

### Using the engine from C (the ABI)

`include/vcal/vcal.h` is the one public header. It is plain C99, backed by the shared library
`vcal`. The caller owns all memory, every function returns a status, and `vcal_last_error` explains
a failure. Its conventions (versioning by `struct_size`, size queries) are in the header and in
[ARCHITECTURE.md §4](ARCHITECTURE.md#4-c-abi-includevcalvcalh).

```c
#include <stdio.h>
#include <string.h>
#include "vcal/vcal.h"

int main(void) {
    const int64_t n[] = {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000};
    const int64_t d[] = {2, 5, 1, 9, 3, 0, 4, 12, 6, 2};
    vcal_context* ctx = NULL;
    vcal_panel panel;
    vcal_grid grid;
    vcal_estimate est;
    vcal_profile_intervals ci;

    memset(&panel, 0, sizeof panel);
    panel.struct_size = (uint32_t)sizeof panel;  /* every struct: set struct_size, zero the rest */
    panel.n_periods = 10;
    panel.n_obligors = n;
    panel.n_defaults = d;
    memset(&grid, 0, sizeof grid);
    grid.struct_size = (uint32_t)sizeof grid;
    memset(&est, 0, sizeof est);
    est.struct_size = (uint32_t)sizeof est;
    memset(&ci, 0, sizeof ci);
    ci.struct_size = (uint32_t)sizeof ci;

    if (vcal_context_create(NULL, &ctx) != VCAL_OK || vcal_grid_default(&grid) != VCAL_OK ||
        vcal_calibrate(ctx, &panel, &grid, &est, &ci) != VCAL_OK) {
        char msg[1024];
        vcal_last_error(msg, sizeof msg, NULL);
        fprintf(stderr, "%s\n", msg);
        return 1;
    }
    printf("PD %.5f [%.5f, %.5f], rho %.4f [%.4f, %.4f], flags %u\n", est.pd, ci.pd_lo, ci.pd_hi,
           est.rho, ci.rho_lo, ci.rho_hi, est.flags);
    vcal_context_destroy(ctx);
    return 0;
}
```

It prints the same fit as the C++ example below (the `readme_c_example` test compiles this code
as C11 and checks the output):

```text
PD 0.00442 [0.00255, 0.00921], rho 0.0522 [0.0086, 0.1974], flags 0
```

`vcal_resample` runs the bootstrap, jackknife,
walk-forward or a supplied index or weight matrix, and `vcal_resample_weights` returns the W behind
it. `vcal_surface` exposes the surface for cell-by-cell replication, and `vcal_dgp_simulate`
generates synthetic panels. `vcal_build_info` reports the git commit, compiler and build settings.

### Using the engine (C++)

The engine is header-only (plus `dgp/` for simulation). A minimal calibration:

```cpp
#include "backends/cpu/cpu_backend.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"

using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;

// One observation per period: obligors n_t and defaults d_t.
const std::vector<Objective::Obs> panel = {{1000, 2}, {1000, 5}, {1000, 1}, {1000, 9}, {1000, 3},
                                           {1000, 0}, {1000, 4}, {1000, 12}, {1000, 6}, {1000, 2}};
const auto T = static_cast<std::int64_t>(panel.size());

// Search box on logit axes: PD in [1e-4, 0.2] (61 points), rho in [1e-3, 0.5] (41 points).
const vcal::Grid<2> grid{{{1e-4, 0.2, 61, vcal::AxisScale::Logit}, {1e-3, 0.5, 41, vcal::AxisScale::Logit}}};

const auto rule = vcal::quadrature::parity_rule();       // the parity integrator
const auto check = vcal::quadrature::parity_rule(true);  // the same rule doubled: the per-run check

std::vector<double> L;  // the T x K surface, kept for profile intervals
vcal::engine::Estimate2 est{};
vcal::engine::calibrate(vcal::backends::CpuBackend{}, Objective{}, rule, check, panel.data(), T, grid, L, est);
const auto ci = vcal::engine::profile_intervals(Objective{}, rule, panel.data(), T, grid, L, est);
// est.value[0], est.se[0]: PD and its SE;  est.value[1], est.se[1]: rho;  est.flags
// ci.lo[a], ci.hi[a]: 95% profile interval;  ci.flags[a]: truncated / not computed
```

With this panel it prints PD 0.00442 (SE 0.00122), with a profile interval of [0.00255, 0.00921],
and ρ 0.0522 (SE 0.0378), with a profile interval of [0.0086, 0.1974].

### Replication and generated files

Every table of coefficients, quadrature rule and golden value is generated by `tools/` at 50 digits
with mpmath, never typed by hand. Each set carries a `MANIFEST.json`, and CI requires exact
regeneration. Tolerance values live in `tests/tolerances.toml`, and the C++ header is generated from
them.

```sh
# the generated tables (mpmath 1.4.1) and the tolerance header (standard library only)
python tools/gen_special_tables.py --check
python tools/gen_quadrature_tables.py --check
python tools/gen_dgp_tables.py --check
python tools/gen_tolerances.py --check

# the plain scipy replication (about five minutes)
uv run --no-project --with scipy==1.18.1 --with numpy==2.5.3 python validation/scipy/binomial_mixture_mle.py

# the recovery harness: 81 scenarios x 1,000 replicates, about 70 minutes on 24 threads
build/cpu-release/tests/recovery_harness --replicates 50   # a quick look, writes nothing
build/cpu-release/tests/recovery_harness --write           # regenerate the recovery goldens
```

## Roadmap

| Milestone | Scope | Status |
|---|---|---|
| M1 | Core math, independent `ref/`, synthetic DGP, cross-reference, recovery harness, methodology and scipy replication, for the binomial-mixture MLE on the CPU | **Done**, tag `v0.1.0-parity-binomial` |
| M2a | Profile-likelihood intervals; recovery re-run with profile, Wald and t(T − 1) side by side | **Done** |
| M2b | Resampling engine: iid and moving-block bootstrap, jackknife, walk-forward and supplied index matrices as weight matrices W × L; percentile intervals; bootstrap coverage study | **Done** |
| M2c | C ABI v0.1: calibrate with profile intervals, resample, surface, DGP, build info; tested from C and bitwise against the engine | **Done** |
| M3 | Vasicek-rate MLE, method of moments, grid Bayesian; `native` profile options | Next |
| M4 | CUDA backend and CPU↔GPU parity suite. Today CI only compiles device code with nvcc 11.8 and 12.8 | Planned |
| M5 | Monitoring: backtests, run ledger, report data | Planned |
| M6 | Benchmark models: multi-segment, AR(1) factor, alternative mixing distributions, Pluto–Tasche, window and influence analysis | Planned |
| M7 | Macro pipeline: transforms with lineage, pre-screens, batched fits, scorers, ranker | Planned |
| M8 | Public-data demonstrations (FDIC, FRED, StatCan), the full scipy replication set, ABI 1.0 | Planned |

## Design principles

- **Replicability over cleverness.**
  - Anything needed to replicate a result is describable in prose and in a short scipy script.
  - Tolerances are named, documented and set from measured errors. Tests use the name, never a
    literal.
  - Each kind of comparison has its own tolerances: mpmath, `ref/`, cross-platform, recovery replay
    and scipy (D-107).
- **Independent checking.**
  - `ref/` shares no code or algorithm with the core, and a layering check enforces it.
  - The goldens are generated by mpmath or `ref/`, never by the engine certifying itself.
- **Log-space numerics.** Log Φ is the core primitive, integrals are computed as log ∫ exp, and
  sums are compensated, so tail PDs and periods with 10⁶ obligors neither underflow nor cancel.
- **Determinism.** Parallelism is only over independent outputs, with fixed summation order and
  fixed tiles. Periods that repeat an observation are evaluated once, with bit-identical results.
- **W × L.** Per-period surfaces are computed once. Every resampling scheme then becomes a weight
  matrix, so a bootstrap costs a weighted reduction, not a refit: 999 replicates cost 15–30% of
  one fit.
- **Say what is known.** Results that are a property of the statistics, not the engine, are
  reported with their diagnosis and pinned, so any change is reviewed. Examples: small-sample
  bias, and interval coverage at small T or near a bound.

## Repository layout

```
include/vcal/vcal.h   the public C ABI (plain C99)
abi/                  its implementation: the shared library vcal
core/                 header-only, host+device, C++17
  special/            log Φ, erfcx, Φ⁻¹, log C(n,d), log-add-exp, ...; generated/ from tools/
  quadrature/         Gauss–Hermite and composite Gauss–Legendre rules; the parity rule
  model/  objectives/ Vasicek conditional PD, binomial-mixture integrand and objective
  reducers/           ArgMax;  grid.hpp  precision.hpp
engine/               surface, weighted reduction, refinement, calibrate, profile intervals
backends/cpu/         OpenMP backend
resample/             weight matrices W (bootstrap, jackknife, walk-forward, supplied indices), replicate estimates
dgp/                  synthetic data generator (Philox4x32-10, own elementary functions)
ref/                  independent FP64 reference (standard library only)
tests/                unit/, xref/ (core vs ref), recovery/ (harness and tests), golden/ with
                      MANIFEST.json files, tolerances.toml -> tolerances.hpp, abi/ (C tests), cuda/, ...
tools/                generators: mpmath tables, DGP reference panels, tolerance header
validation/scipy/     scipy replication scripts (outside the build)
docs/methodology/     methodology notes, tolerance register, recovery results, sign-off
cmake/  .github/      build modules and checks; CI
```

## License

Apache-2.0, copyright 2026 Keith Yip. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

## Disclaimer

This is pre-alpha research software, provided "as is" without warranty of any kind (see LICENSE).
It has not been validated for production or regulatory use, and interfaces and results may change
without notice. Anyone relying on its output is responsible for validating it independently.
