# Correctness audit, November 2026

A fresh correctness audit of the library and its scripts (D-218), run because the earlier Copilot
audit's findings were never committed and could not be recovered. The audit looked for bugs only,
not style, naming or performance. It split the code five ways and reviewed each part independently:

1. `core/`, `common/`, `backends/`;
2. `engine/`, `resample/`;
3. `abi/`, `include/`, `ref/`, `dgp/`;
4. `benchmark_models/`, `validation/`, `tools/` (and the CI plan job);
5. every script under `studies/`.

Every finding below was checked against the code a second time before it was recorded. CONFIRMED
means it was reproduced or follows from the exact lines; PLAUSIBLE means the defect is real in the
code but no input reaching it in practice was found.

**Headline.** No finding changes a golden or a recorded verdict. Every study scorer that can run on
committed data was re-run and reproduced its recorded table exactly. The one verdict-changing defect
found is S-6's P5 slack, which D-217 already records.

## Findings

| # | Severity | Where | Defect | Status | Disposition |
|---|---|---|---|---|---|
| 1 | medium | `validation/scipy/{method_of_moments,vasicek_rate_mle,binomial_mixture_mle,grid_posterior,hac_se}.py` | Checks written `x > tol` and running `max(worst, x)` never fail on a NaN, so a NaN engine golden passed the independent scipy gate (reproduced: exit 0). | CONFIRMED | Fixed: `not x <= tol`, and NaN mapped to inf before a running max. |
| 2 | medium | `.github/workflows/ci.yml` (plan job) | `git diff --name-only` lists only a rename's destination, so renaming a code file or a code-relevant `.md` to a prose path skipped all CI and `ci-ok` still passed (shown on `43e9b87`). | CONFIRMED | Fixed: `--no-renames`. |
| 3 | medium | `engine/posterior.hpp` (`jeffreys_table`), `backends/cpu/cpu_backend.hpp` | The Jeffreys table allocates n + 1 doubles inside an OpenMP region; a `bad_alloc` there calls `std::terminate`, so a huge n killed the caller's process instead of returning `VCAL_E_OUT_OF_MEMORY`. A large n that fits in memory is effectively a hang. | CONFIRMED | Fixed: `CpuBackend::parallel_for` keeps the first exception and rethrows it after the region; the ABI refuses n > `VCAL_JEFFREYS_MAX_OBLIGORS` (10⁶) for the Jeffreys prior. |
| 4 | low | `engine/profile.hpp` (`panel_loglik`), `engine/derivs.hpp` (`panel_derivs_natural`) | The compensated sum turns a −inf term into NaN ((−inf) − (−inf)); `reducers::weighted_sum` already guards this. A NaN fails the conditional-PD walk where −inf would bracket. No trigger known with valid binomial inputs. | CONFIRMED (arithmetic) | Fixed: non-finite terms bypass the compensation, as in `weighted_sum`; bit for bit unchanged otherwise. |
| 5 | low | `core/model/binomial_mixture.hpp` | The natural-scale score and Hessian overflow to NaN for PD below about 1e-150 or ρ below about 1e-200; the log-likelihood stays correct. The ABI accepted such grids (lo > 0 only). | CONFIRMED | Fixed at the ABI: every axis needs lo ≥ `VCAL_GRID_MIN_PROBABILITY` (1e-100). The math is unchanged. |
| 6 | low | `abi/vcal.cpp` (`vcal_context_create`) | `n_threads` had no upper bound; libgomp may `exit()` on an absurd thread count. | PLAUSIBLE | Fixed: `n_threads` ≤ `VCAL_MAX_THREADS` (4096). |
| 7 | low | `engine/posterior.hpp` (`masses`, `jeffreys_log_prior`) | Grid points with a NaN log posterior silently got zero mass with no flag; the Jeffreys interpolation could make one with 0 × (−inf) next to a −inf table node. | PLAUSIBLE | Fixed: zero-weight corners are skipped (bit for bit the old value when finite), and a new flag, `kPosteriorNanCells` / `VCAL_POSTERIOR_NAN_CELLS`, reports dropped cells. |
| 8 | low | `abi/vcal.cpp` (`fit_rates`) | `vcal_calibrate_rate` wrote `*estimate` before the profile step could still fail, against the file's own all-or-nothing contract. | CONFIRMED | Fixed: both outputs computed, then written. |
| 9 | low | `engine/calibrate.hpp` | `out` was not reset on the early error returns. Current callers are safe. | CONFIRMED | Fixed: reset first. |
| 10 | latent | `studies/{jackknife-bias-rho,bayes-coverage,derived-quantity-intervals}/compare.py` | An inclusive "within 0.015" test of a difference of k/1,000 coverages fails at exactly 15 in 1,000 (`0.961 − 0.946` is `0.015000000000000013`). No committed result sits on a bound. | CONFIRMED | Not changed: these scorers belong to scored, registered studies. Rule for future scorers in D-218. |
| 11 | latent | `studies/parametric-bootstrap/compare.py` (consistency check) | Ends compared with `==`, so two NaNs counted as different. It could only stop the scoring, never change a verdict. | CONFIRMED | Fixed: NaN-aware equality. Its output on the committed fits is byte for byte unchanged. |
| 12 | docs | `include/vcal/vcal.h` | `VCAL_E_UNSUPPORTED` is documented only as "a field this library does not know", but is also returned for non-logit posterior axes, unequal n with Jeffreys and a missing DGP kernel; the posterior's logit-only requirement is not in the header. | CONFIRMED | Copilot-audit PR B (docs accuracy). |
| — | (known) | `studies/pluto-tasche/score.py` | P5's 1e-12 slack is narrower than the two solvers agree. | CONFIRMED | Already recorded in D-217; no change. |

## Reviewed and found clean

- **Philox and the DGP:** the constants, rounds and key bump match Random123 (the known-answer
  vectors were reproduced independently); the SSE2/AVX2 lanes; the counter layout (no stream
  collisions, no overflow); the uniform, which is never 0 or 1; the polar method; the AS241
  coefficients; the AR(1) panel at φ = 0 bit for bit equal to the base panel.
- **`ref/`:** includes only the standard library, so it shares no quadrature with `core/`; its
  math checks out.
- **Special functions and quadrature:** `log_phi`, `erfcx`, `probit` (under 1 ulp against mpmath),
  `lbinom`, Gauss–Hermite (fixed and adaptive), the composite Legendre rule, the SIMD kernels.
- **Engine and resampling:** weight matrices (iid and block bootstrap, jackknife, walk-forward), the
  bounded argmax, distinct-row compaction, quantiles, BCa, Bartlett, the surface engine and its
  row cache, refinement, the Brent and safeguarded Newton solvers, profile walks, conditional PD
  (gradients checked numerically), HAC, the method of moments, the posterior's moments and HPD.
- **Benchmarks:** the AR(1) filter, search and kernels; the Pluto–Tasche bounds against the paper
  and an independent scipy replication (which can fail: worst 6.3e-11 against a 2e-10 tolerance).
- **Tools and CI:** `check_manifests.py`, `check_data_files.py`, `gen_engine_flags.py`,
  `tablegen`/`gen_tolerances --check`, the `ci-ok` gate, `perf/project_runtime.py`.
- **Study scorers:** all re-run where the data are committed (S-3/5/21, S-8/27/28, S-9, S-13, S-15,
  S-23, S-34, the S-35 main arm), reproducing their recorded tables cell for cell; the rest read
  against their registrations.
