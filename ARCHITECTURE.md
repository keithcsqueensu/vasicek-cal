# vasicek-cal — Architecture

Status: **accepted**, 2026-09-26. Decision ids (D-nnn) refer to `DECISIONS.md`, which also lists every open item.
References like (D-0nn) point at the decision that fixed a design point; **[Qn]** marks one still open.

---

## 1. Model and notation

Single-factor ASRF. Periods `t = 1..T`, `n_t` obligors, `d_t` defaults, systematic factor
`Z_t ~ N(0,1)` (iid; the AR(1) variant is a benchmark model, §5.4).

```
x(z; PD, ρ) = (Φ⁻¹(PD) − √ρ · z) / √(1−ρ)
p(z)        = Φ(x(z))              conditional PD
1 − p(z)    = Φ(−x(z))             computed directly, never as 1 − p

ℓ_t(θ) = log C(n_t, d_t) + log ∫ exp( d_t·logΦ(x(z)) + (n_t−d_t)·logΦ(−x(z)) ) φ(z) dz
ℓ(θ)   = Σ_t w_t · ℓ_t(θ)          (w = 1 for a plain calibration)
```

Everything is in log space: `logΦ` is the core primitive, integrals are `log ∫ exp(g)`,
and the quadrature sum is a log-sum-exp. Without this, tail PDs and large-`n` periods underflow.

Useful closed forms used as tests (not as implementation):
`E[p(Z)] = PD` exactly, and `E[p(Z)²] = Φ₂(Φ⁻¹PD, Φ⁻¹PD; ρ)`.

---

## 2. Repository layout

```
vasicek-cal/
├── CMakeLists.txt  CMakePresets.json  LICENSE  NOTICE
├── cmake/                  arch lists, warning sets, FP-contraction flags per target
├── include/vcal/vcal.h     the public C ABI — the only installed header
├── core/                   header-only, host+device (VCAL_HD); no heap, no I/O, no threads, no exceptions
│   ├── precision.hpp       Precision policies, VCAL_HD macro
│   ├── special/            log_phi, erfcx, probit, log_add_exp, lbinom; generated/ from tools/
│   ├── model/              Vasicek1F (conditional PD), binomial-mixture log-integrand + adaptive hint
│   ├── quadrature/         Integrator contract, Gauss–Hermite policies; generated/ rules from tools/ (D-080)
│   ├── objectives/         per-period contribution policies
│   ├── reducers/           ArgMax, LogSumExp (posterior), moment inversion
│   └── grid.hpp            axes, flattening, point lookup
├── engine/                 backend-agnostic surface evaluation + weighted reduction templates
├── backends/
│   ├── cpu/                OpenMP
│   └── cuda/               .cu kernels; built only when VCAL_ENABLE_CUDA=ON (D-030)
├── resample/               builders for W: iid bootstrap, block bootstrap, jackknife, walk-forward, custom
├── dgp/                    Philox4x32-10, own deterministic log/exp/Φ/Φ⁻¹, polar factor, Bernoulli sums — CPU only,
│                           bitwise identical across platforms and to a Python mirror (D-052–D-057, D-109–D-111)
├── abi/                    extern "C" implementation, dispatch registry, status/error handling
├── ref/                    independent FP64 reference — links nothing above (D-031–D-033)
├── benchmark_models/       (M6) multi-segment, AR(1), alt. mixing, Pluto–Tasche, window/influence
├── macro/                  (M7) transforms + lineage, pre-screens, batched fits, scorers, ranker
├── monitoring/             (M5) backtests, run ledger, report data
├── tests/
│   ├── harness/            minimal in-house runner, registered with CTest (D-029)
│   ├── unit/               special functions, quadrature, grid, RNG known-answer tests
│   ├── crossref/           core vs ref on identical inputs
│   ├── recovery/           statistical recovery on DGP panels
│   ├── parity/             CPU vs GPU (M4)
│   ├── abi/                a pure-C translation unit that exercises vcal.h
│   └── golden/             checked-in golden values + provenance manifest (D-031–D-033)
├── validation/scipy/       short replication scripts matching the prose docs; outside the build (D-031–D-033)
├── tools/                  generators outside the build, e.g. gen_special_tables.py (mpmath → headers + goldens)
├── docs/methodology/       one plain-prose note per estimator/test + tolerances register
└── perf/                   performance benchmarks (distinct from benchmark_models/)
```

Naming: "benchmark" in the charter means *challenger models*, so that code lives in
`benchmark_models/`; performance measurement lives in `perf/`.

### Layering rules (enforced by CMake target graph + a CI include check)

| Target | May depend on | Must not |
|---|---|---|
| `core` | nothing | allocate, do I/O, throw, spawn threads |
| `engine` | core | know about a specific backend |
| `backends/*` | core, engine | — |
| `dgp` | nothing (own deterministic math, D-052) | use `core/special` — it's a third independent Φ |
| `ref` | C++ standard library only | include anything from core/engine/backends/dgp |
| `abi` | everything above except ref | let exceptions or C++ types cross the boundary |

---

## 3. Layer contracts and policy interfaces

### 3.0 Layer contracts (language-agnostic)

Each layer is specified here in plain terms, so it can be re-implemented in any language (a
short scipy script included) and checked against this engine without reading the C++. Numbers
are IEEE binary64 throughout, and every sum that the contract calls "fixed order" runs in
ascending index. The C++ interfaces that implement these contracts follow in §3.1–§3.7; the
tolerances that define "agrees" are in the register (§7).

**Objective** (§3.5; `core/objectives/`)

- *Inputs:* one period's observation (binomial mixture: obligors n_t ≥ 1 and defaults
  0 ≤ d_t ≤ n_t), one parameter point θ = (PD, ρ) with 0 < PD < 1 and 0 < ρ < 1, and an
  integrator.
- *Outputs:* the period's log-likelihood contribution
  ℓ_t(θ) = log C(n_t, d_t) + log ∫ p(z)^d_t (1 − p(z))^(n_t−d_t) φ(z) dz, with p(z) as in §1;
  a hint for the integrator (the mode and scale of the integrand, or, for d_t ∈ {0, n_t}, the
  centre and width of its one-sided drop); and the size of ℓ_t's terms, used by the per-run
  check's threshold.
- *Invariants:* ℓ_t depends on (observation, θ) only: no state, no other period. Equal
  observations give bitwise-equal values. Everything is computed in log space; 1 − p(z) is
  computed as Φ(−x), never as 1 − Φ(x). No allocation, I/O, threads or exceptions. Invalid
  panels are rejected before evaluation, with the offending period named.

**Integrator** (§3.4; `core/quadrature/`)

- *Inputs:* a log-integrand g (z ↦ log f(z) for a non-negative f) and the objective's hint.
- *Outputs:* log ∫ exp(g(z)) φ(z) dz, one binary64 number.
- *Invariants:* the result is a log-sum-exp over the rule's nodes, folded in node order. Nodes
  and weights come from fixed tables generated at 50 digits with mpmath (stored as
  log-weights), never computed at run time. The number of evaluations is fixed for a given
  rule: only where the nodes sit depends on the hint. Every rule has a doubled counterpart,
  used by the per-run convergence check. The parity rule is adaptive Gauss–Hermite with
  N = 128 for 0 < d < n, and composite Gauss–Legendre with 16 panels of 16 points for
  d ∈ {0, n}.

**Surface engine** (§3.7; `engine/`)

- *Inputs:* a panel of T observations; a grid of K parameter points (per axis: lower and upper
  bound, n ≥ 3 points, and a scale, uniform in the scaled coordinate); an objective, an
  integrator, and a backend that runs independent tasks.
- *Outputs:* stage 1, the surface L, a T × K row-major matrix with L[t][k] = ℓ_t(θ_k).
  Stage 2, for each row b of a B × T weight matrix W, the reducer's result over the weighted
  surface s_b[k] = Σ_t W[b][t]·L[t][k].
- *Invariants:* L is computed once per (panel, grid, model) and does not depend on W. Each
  Σ_t runs in fixed order with compensated (Neumaier) summation, and skips periods whose weight
  is zero. Parallel work is split only across independent outputs: grid points in stage 1,
  fixed 1024-point k-tiles in stage 2, merged in tile order. Results are therefore bitwise
  identical for any thread count. The B × K matrix is never stored. Periods with equal
  observations share one evaluated row, copied bit for bit.

**Reducer** (§3.6; `core/reducers/`)

- *Inputs:* a stream of (grid index k, value s[k]) pairs, in any order and split into any
  number of partial folds.
- *Outputs:* a summary of the stream. ArgMax returns the index and value of the maximum and the
  number of NaN values seen.
- *Invariants:* a reducer is an initial state, a fold of one pair, and a merge of two states.
  ArgMax's order is total (the larger value wins; on a tie the smaller k wins; NaN never wins
  and is counted), so its merge is exactly associative and commutative. A reducer whose merge
  is not exact is merged in a fixed order. Refinement, standard errors and flags are separate
  passes that read the neighbours of the maximum; they are not part of the fold. A flag is
  reported, never used to alter a result silently.

**Resampling** (§5.2; `resample/`)

- *Inputs:* the number of periods T and a scheme with its settings: the iid or moving-block
  bootstrap (seed, number of replicates B, block length), the jackknife, walk-forward windows
  (window length, first end, step), a caller-supplied B × m matrix of period indices, or a
  caller-supplied B × T weight matrix.
- *Outputs:* the B × T weight matrix W and, for the random schemes, the B × T period indices
  behind it. From W and the surface L: B replicate estimates, each with its flags, and type-7
  percentile intervals.
- *Invariants:* every random draw is a pure function of (seed, scheme, replicate, draw index):
  Philox4x32-10 with key = seed XOR the ASCII tag "RSMPBOOT" and counter
  (scheme, replicate, draw/2, 0). Replicates therefore do not depend on the thread count or on
  B, and never share a stream with the DGP for the same seed. The weight for a period is the
  number of times it is drawn: supplying another tool's indices reproduces its replicates one
  for one. The jackknife is W = 1·1ᵀ − I, and walk-forward rows are 0/1 windows. The
  integrals are never recomputed: each replicate is a re-weighted sum of L, which is exact for
  objectives that add over periods, and so excludes the AR(1) factor model (D-043).

**Reference** (`ref/`)

- *Inputs:* a panel of (n_t, d_t) and a search box for PD and ρ; for an interval, also the fit,
  which parameter, and the likelihood-ratio threshold.
- *Outputs:* log Φ, Φ⁻¹, log C(n, k), the log mixture integral, per-period and total
  log-likelihoods; the maximum-likelihood estimate with standard errors from a
  finite-difference Hessian, their correlation and a boundary flag; profile-likelihood
  interval endpoints with at-bound flags.
- *Invariants:* binary64, with the language's standard library only. It shares no code and no
  algorithm with the engine: a different Φ, a different quadrature (adaptive Gauss–Legendre
  with nodes computed at start-up) and a different optimiser (nested golden-section search).
  The layering check enforces the code part (§2). Its own accuracy is measured against the
  mpmath goldens. Speed is not a goal, and it is used only as an oracle in tests, never on the
  engine's path.

**DGP** (`dgp/`; full description in `docs/methodology/dgp.md`)

- *Inputs:* seed, scenario, replicate, PD, ρ, and the obligor count n_t of each period.
- *Outputs:* the default count d_t of each period and, optionally, its factor draw Z_t.
- *Invariants:* a panel is a pure function of its inputs and is bitwise identical on every
  platform and to the pure-Python mirror. Draws are Philox4x32-10 with counter
  (scenario, replicate, period, block). Adding scenarios, replicates or periods, or changing
  n_t in one period, leaves every other period unchanged. Z_t comes from the Marsaglia polar
  method; d_t counts obligors whose uniform falls below p(Z_t). The elementary functions (log,
  exp, Φ, Φ⁻¹) are the DGP's own, built from exactly or correctly rounded operations only. It
  runs with FP contraction off, refuses to run with flush-to-zero active, and runs on the CPU
  only.

The resampling, reference and DGP entry points in C++:

```cpp
namespace vcal::resample {   // resample/weights.hpp, resample/bootstrap.hpp
std::vector<std::int64_t> bootstrap_indices(std::uint64_t seed, Scheme, std::uint32_t replicates,
                                            std::int64_t periods, std::int64_t block_length = 1);
std::vector<double> weights_from_indices(const std::int64_t* indices, std::int64_t replicates,
                                         std::int64_t draws, std::int64_t periods);
std::vector<double> jackknife_weights(std::int64_t periods);
std::vector<double> walk_forward_weights(std::int64_t periods, std::int64_t window,
                                         std::int64_t first_end, std::int64_t step);
template <class Backend>
void replicate_estimates(const Backend&, const Grid<2>&, const double* L, std::int64_t periods,
                         const double* W, std::int64_t replicates, Replicate2* out);
PercentileInterval percentile_interval(const Replicate2* reps, std::int64_t replicates, int param,
                                       double level = 0.95);
}
namespace vcalref {          // ref/ref.hpp
double log_mixture(double pd, double rho, std::int64_t n, std::int64_t d);
double loglik(double pd, double rho, const std::vector<Period>& panel);
Fit fit(const std::vector<Period>& panel, double pd_lo, double pd_hi, double rho_lo, double rho_hi);
ProfileInterval profile_interval(const std::vector<Period>& panel, const Fit&, int param, double pd_lo,
                                 double pd_hi, double rho_lo, double rho_hi, double threshold);
}
namespace vcal::dgp {        // dgp/dgp.hpp
Status simulate_panel(const PanelSpec& spec, std::int64_t* defaults, double* factors);
}
```

### C++ policy interfaces

Device-visible code (`core/`, `engine/`, `backends/cuda/`, per D-058 and D-062) is C++17, so the
contracts below are written as concepts only for readability. In code, each one is a trait
(`is_precision_v<P>`, …) built from `std::void_t`/SFINAE, checked by a `static_assert` next
to each policy's definition, one message per clause (D-063). A contract break then fails at
the policy, with a readable message, not at a use site:

```cpp
template <class P, class = void> struct is_precision : std::false_type {};
template <class P>
struct is_precision<P, std::void_t<typename P::eval_t, typename P::accum_t>>
    : std::bool_constant<std::is_floating_point_v<typename P::eval_t> &&
                         std::is_same_v<typename P::accum_t, double>> {};
template <class P> inline constexpr bool is_precision_v = is_precision<P>::value;
```

### 3.1 Precision

```cpp
namespace vcal {

#if defined(__CUDACC__)
#  define VCAL_HD __host__ __device__ inline
#else
#  define VCAL_HD inline
#endif

// eval_t: arithmetic inside the integrand.  accum_t: log-sum-exp and Σ_t accumulation.
struct PrecisionF64   { using eval_t = double; using accum_t = double; };  // parity + native
struct PrecisionMixed { using eval_t = float;  using accum_t = double; };  // native only (D-035)(D-039)

template <class P>
concept Precision = std::floating_point<typename P::eval_t> &&
                    std::floating_point<typename P::accum_t>;
```

### 3.2 Special functions (core/special)

```cpp
namespace vcal::special {                        // binary64 only until M4 (D-076)
VCAL_HD double log_phi(double x);                // log Φ(x) on ℝ; upper tail via log1p(−Φ(−x)) (D-068)
VCAL_HD double erfcx(double x);                  // exp(x²)·erfc(x); mpmath-fitted Chebyshev (D-078)
VCAL_HD double probit(double p);                 // Φ⁻¹(p): Wichura AS241 (PPND16) (D-077)
VCAL_HD double probit_upper(double q);           // Φ⁻¹(1 − q), taking q directly (D-069)
VCAL_HD double log_add_exp(double a, double b);  // −∞-safe (D-071)
VCAL_HD double lbinom(int64_t n, int64_t k);     // log C(n,k): exact n ≤ 60, else Stirling differences (D-070)
VCAL_HD double stirling_error(int64_t n);        // Loader's δ(n): generated table n ≤ 30, series beyond
}
```

Coefficients, constants and the Stirling table live in `core/special/generated/`, written by
`tools/gen_special_tables.py` from mpmath (D-033, D-067); nothing there is typed by hand.
Accuracy per function is in the tolerance register (§7).

The same implementation runs on host and device so the CPU/GPU gap is only
compiler/FMA-level, not algorithmic. `ref/` uses different algorithms (std::erfc +
continued-fraction tail) on purpose.

### 3.3 Grid

```cpp
enum class AxisScale : int32_t { Linear, Log, Logit, Probit };  // uniform in the scaled coordinate u = S(v)
struct Axis {                                   // n >= 3 (refinement needs a 3-point stencil)
    double lo, hi; int32_t n; AxisScale scale;
    double step() const;                        // in u
    double scaled_at(int32_t i) const;          // u_lo + i·step
    double value_at(int32_t i) const;           // S^-1(u_i); endpoints exactly lo and hi
};
template <int D> struct Grid {
    Axis axis[D];
    int64_t size() const;                              // Π n_a
    void    unflatten(int64_t k, int32_t (&i)[D]) const;  // row-major, axis 0 slowest
    int64_t flatten(const int32_t (&i)[D]) const;
    void    values(int64_t k, double (&v)[D]) const;   // natural scale
};
inline constexpr double kDefaultRhoUpper = 0.5;        // configurable per grid (D-089)
```
`grid::to_scaled`, `grid::from_scaled` and `grid::dvalue_dscaled` (the delta method's
Jacobian) are defined per scale. `axis_error` and `grid_error` validate on the host.

### 3.4 Integrator

Computes `log ∫ exp(g(z)) φ(z) dz` for a log-integrand `g`.

```cpp
namespace vcal::quadrature {
// mode of h = g(z) − z²/2 and 1/√(−h''); fixed rules ignore it. For a one-sided integrand
// (d ∈ {0, n}) also the survival factor's half-point and width, used by the composite rule.
struct IntegrandHint { double mode; double scale; double center; double width; bool one_sided; };

// Standard-normal rule, non-owning: E[f(Z)] ≈ Σ exp(log_weight[i]) f(node[i]).
struct GaussHermiteRule { const double* node; const double* log_weight; int n; };

// Contract (C++17 trait, D-063): is_integrator_v<I> / check_integrator<I>() require
//   double log_integrate(const G& g, IntegrandHint) const   for g callable as double(double).

template <class P> struct GaussHermiteFixed    { GaussHermiteRule rule; /* log_integrate */ };
// Mode-centred, scaled (Liu & Pierce 1994): z = mode + scale·u. Required for large n (D-035, D-037).
template <class P> struct GaussHermiteAdaptive { GaussHermiteRule rule; /* log_integrate */ };

GaussHermiteRule gauss_hermite_rule(int n);  // host; n ∈ {8,16,24,32,48,64,96,128,256}, else {nullptr,nullptr,0}

// One-sided integrands (D-118): K panels of M Gauss–Legendre points under z = center + width·sinh(u),
// over the range where h has dropped by 50 from the mode. K and M are fixed at every parameter point.
struct GaussLegendreRule { const double* node; const double* weight; int n; };
template <class P> struct CompositeLegendre { GaussLegendreRule rule; int panels; /* log_integrate */ };
// hint.one_sided ? one_sided : interior
template <class P> struct SplitRule { GaussHermiteAdaptive<P> interior; CompositeLegendre<P> one_sided; };
GaussLegendreRule gauss_legendre_rule(int m);  // host; m ∈ {8,12,16,20,24,32}

// The parity rule: GH N = 128 and composite K = 16 × M = 16; doubled = true gives the per-run
// check (GH 256, K = 32), D-092, D-118.
SplitRule<PrecisionF64> parity_rule(bool doubled = false);
}
// For closed-form objectives (Vasicek-rate): no integral (M3).
```

The rules are mpmath-generated tables, not computed at runtime (D-080, D-081). Log-weights
keep the tiny outer weights at full relative precision. There is no eigensolver anywhere in
the engine. The tables are host arrays: the CPU backend passes pointers into them, and the
CUDA backend copies them to device memory once per context (M4).

Terms are folded in node order through a one-pass log-sum-exp, so results are deterministic.
The adaptive rule's exponent correction log φ(z) − log φ(u) is computed as (u − z)(u + z)/2,
so u² and z² never cancel.

The binomial-mixture log-integrand and its hint live in `core/model/binomial_mixture.hpp`
(D-084). The hint is the mode of the strictly concave h(z) = g(z) − z²/2, found by bracketed
Newton with a bisection fallback from the large-n start, plus scale = 1/√(−h''). Every step is
deterministic. Accuracy and convergence by N are in the tolerance register (§7).

### 3.5 Objective

Per-period contribution at one parameter point.

```cpp
template <class O, class P, class I>
concept Objective = requires(const O& o, const typename O::Obs& y,
                             const typename O::Theta& th, const I& in) {
    { O::n_params } -> std::convertible_to<int>;
    { o.log_contrib(y, th, in) } -> std::same_as<typename P::accum_t>;
    { O::hint(y, th) } -> std::same_as<IntegrandHint>;
    { O::theta_from_grid(std::declval<const double(&)[O::n_params]>()) } -> std::same_as<typename O::Theta>;
    { y == y } -> std::convertible_to<bool>;  // equal observations share one surface row (D-122)
};

struct Theta1F { double pd; double rho; };

template <Precision P>
struct BinomialMixture {
    struct Obs { int64_t n, d; friend bool operator==(const Obs&, const Obs&); };
    using Theta = Theta1F;
    static constexpr int n_params = 2;
    // mode: z* solving p(z*) = (d+½)/(n+1); scale from curvature of the log-integrand at z*.
    // d ∈ {0, n}: one_sided, center where (1−p)ⁿ (or pⁿ) = ½, width 1/(n·β·λ) there (D-118).
    static VCAL_HD IntegrandHint hint(const Obs&, const Theta&);
    // abs(log C) + abs(log I): the size of l_t's terms, for the check threshold (D-120).
    static VCAL_HD double rounding_scale(const Obs&, double l_t);
    template <class I> VCAL_HD typename P::accum_t log_contrib(const Obs&, const Theta&, const I&) const;
};

template <Precision P>
struct VasicekRate {            // log Vasicek density of an observed default rate; closed form
    struct Obs { double dr; };  // zero/one rates: (D-044)
    using Theta = Theta1F; static constexpr int n_params = 2;
    template <class I> VCAL_HD typename P::accum_t log_contrib(const Obs&, const Theta&, const I&) const;
};

// Method of moments does not produce a grid surface; it produces per-period
// sufficient statistics (d/n and d(d−1)/(n(n−1))) that W reweights. (D-042, D-043)
```

### 3.6 Reducer

Folds one weighted-surface row `s_b[k] = Σ_t W[b,t]·L[t,k]` into a result.

```cpp
template <class R, class P>
concept Reducer = requires(const R& r, typename R::State& s, const typename R::State& o,
                           int64_t k, typename P::accum_t v) {
    { r.init() } -> std::same_as<typename R::State>;
    r.push(s, k, v);   // fold value v at grid index k
    r.merge(s, o);     // combine partial states from different k-tiles
};

struct ArgMax {
    struct State { double best; int64_t k; };   // ties → smaller k, so merge is exactly
    /* init/push/merge */                        // associative & commutative (order-independent)
};

struct LogSumExpPosterior {                      // grid Bayesian: log-normaliser + moments
    struct State { double m; double s; double s1[2]; double s2[3]; };
    // floating merge is NOT associative → engine merges k-tiles in a fixed order (§6)
};
```

Post-reduction steps (sub-grid refinement (D-038), observed-information SEs, edge flags)
read the neighbours of `k*` in a second, cheap pass; they are not part of the fold.

### 3.7 Engine and backends

```cpp
// Stage 1: L[t][k] = O::log_contrib(obs_t, θ(grid, k), integrator). Computed once per (data, grid, model).
template <class O, class I, Precision P, class Backend>
Status evaluate_surface(Backend&, const O&, const I&,
                        std::span<const typename O::Obs> obs,
                        const Grid<O::n_params>&,
                        MatrixView<typename P::accum_t> L);          // T × K, row-major

// Stage 2: for each row b of W, reduce_k Σ_t W[b,t]·L[t,k].
// Fused and tiled over (b, k): the B × K matrix is never materialised.
template <class R, Precision P, class Backend>
Status reduce_weighted(Backend&, const R&,
                       ConstMatrixView<typename P::accum_t> L,       // T × K
                       ConstMatrixView<double> W,                    // B × T
                       std::span<typename R::Result> out);           // B

struct CpuBackend  { int32_t n_threads; };
struct CudaBackend { int32_t device; cudaStream_t stream; };
```

Templates are instantiated only for supported tuples, listed in one registry in `abi/`
that maps `(profile, objective, integrator, reducer, precision, backend)` enums to function
pointers. An unlisted tuple returns `VCAL_E_UNSUPPORTED`, it never falls back silently.
A **profile** is a whitelist over that registry plus fixed settings (D-035). The public `parity`
profile means textbook and scipy-replicable, not conformance to any particular organisation's
method; a user who needs one defines their own conformance profile as another whitelist over the
same registry (D-036).

---

## 4. C ABI (`include/vcal/vcal.h`)

The header is the reference: every struct, field and function is documented there. v0.1 (M2c,
D-138–D-144) is summarised below.

**Conventions**

- **Plain C99.** The header uses only `<stdint.h>` types, doubles and pointers. tests/abi compiles
  a consumer as C11 and the header alone as C99 with pedantic errors, so a C++-only construct
  fails the build.
- **Names.** Every function and type starts with `vcal_`, and every macro and enumerator with
  `VCAL_`. Exports go through `VCAL_API` (`__declspec(dllexport/dllimport)`, or default
  visibility). Calls use `VCAL_CALL`, which is `__cdecl` on Windows.
- **Versioning (D-138).** `vcal_abi_version()` returns `(major << 16) | minor`.
  - An additive change bumps the minor; a breaking one bumps the major (the Linux SOVERSION
    follows the major).
  - Every public struct starts with `uint32_t struct_size`, set by the caller, and structs only
    grow.
  - A struct_size below this version's size is invalid.
  - A larger input struct is accepted only if its unknown bytes are zero; otherwise the call
    returns `VCAL_E_UNSUPPORTED`, never a silent ignore.
  - A larger output struct gets its known fields written.
  - `reserved` fields must be zero.
  - Layouts have no implicit padding and are checked on both sides.
- **Memory (D-139).** The caller allocates everything; the library never returns memory to
  free and keeps no caller pointer.
  - Variable-size outputs take `(buffer, capacity, required)`. A NULL buffer is a size query.
  - A buffer that is too small gives `VCAL_E_BUFFER_TOO_SMALL`, with the length needed stored
    and nothing written.
- **Errors (D-140).** Every function returns a `vcal_status`; only `vcal_abi_version` and
  `vcal_status_string` cannot fail.
  - A failure stores a message for the calling thread, which `vcal_last_error` copies out. The
    message names the function and the field, e.g.
    `vcal_calibrate: panel.n_defaults[3] = 1001: need 0 <= d <= n_obligors[3] = 1000`.
  - Each entry point runs inside a catch-all, so no C++ exception crosses the boundary.
- **Context (D-141, D-045).** An opaque `vcal_context`, made by `vcal_context_create` and freed
  by `vcal_context_destroy`, holds the profile (`VCAL_PROFILE_PARITY`) and the thread count.
  - Results are bitwise identical for any thread count.
  - A context must not be used by two threads at once; distinct contexts are independent.

**Functions (v0.1; v0.3 adds the M3 estimators, D-169)**

| Function | What it does |
|---|---|
| `vcal_abi_version`, `vcal_build_info`, `vcal_last_error`, `vcal_status_string` | Version, build provenance (`key=value` lines: git commit and dirty flag, compiler, build type, OpenMP, CUDA, profiles), errors |
| `vcal_context_create`, `vcal_context_destroy` | The context |
| `vcal_grid_default`, `vcal_grid_values` | The recovery-study grid (D-115); an axis's exact natural-scale values |
| `vcal_calibrate` | Estimate, SEs, flags and the quadrature check (§5.1); with a non-NULL `vcal_profile_intervals`, the 95% profile-likelihood intervals too (D-128–D-130) |
| `vcal_surface` | The T × K per-period log-likelihood surface, for cell-by-cell replication |
| `vcal_resample` | Replicate estimates from W × L (§5.2) for the iid and moving-block bootstrap, jackknife, walk-forward, a supplied index matrix or a supplied weight matrix, with type-7 percentile intervals |
| `vcal_resample_weights` | The W a scheme implies, so any replicate can be rebuilt elsewhere |
| `vcal_dgp_simulate` | One synthetic panel, bitwise identical on every platform (D-052–D-057) |
| `vcal_calibrate_rate` (0.3) | The Vasicek-rate MLE on count data or a rate series, with profile intervals; zero rates by an explicit `VCAL_ZERO_RATES_*` option: refuse (parity, D-044) or, in a `VCAL_PROFILE_NATIVE` context, censor, substitute or drop. No default is recommended until S-28 reports |
| `vcal_calibrate_moments` (0.3) | The joint-default-probability method of moments, counts or rates |
| `vcal_calibrate_posterior` (0.3) | The grid posterior, flat or Jeffreys prior, with the resolution rule; equal-tailed and HPD intervals. Jeffreys tables are kept in the context per (n, grid) |

**Deliberately absent from v0.** These arrive as appended fields or new functions, under the
same conventions:

- model selection (objective, integrator, precision), the dispatch registry, and the `native`
  profile (M3);
- the CUDA backend and device selection (M4);
- backtests and the ledger (M5);
- multi-segment models (M6);
- macro fits (M7);
- a profile-interval level other than 95%.

---

## 5. Data flow

### 5.1 Calibrate

```
panel (n_t, d_t) ──► validate: n_t > 0, 0 ≤ d_t ≤ n_t, T ≥ 2 ──► Grid(spec)
                                                                  │
                     evaluate_surface(Objective, Integrator) ─────┴──► L [T × K]
                                                                         │
                     W = [1 … 1]  (1 × T) ──► reduce_weighted(ArgMax) ───┴──► k*, ℓ(k*)
                                                                                │
   2nd pass on the 3×3 stencil around k*, in scaled coordinates: quadratic vertex with cross
   term (D-095). SEs from the objective's own Hessian at the vertex, step 0.15 × the stencil
   SE, delta method (D-119). Flags: edge / flat surface / quadrature unconverged /
   near a bound (within 2 SEs)  (never silently clamp)                           │
                                                                                ▼
                                                                         vcal_estimate
```

`evaluate_surface` evaluates each distinct observation once per grid point and copies the row
to the periods that repeat it (D-122). The copy is bitwise what a second evaluation would give.

After calibration, `engine::profile_intervals` (D-128–D-130) gives each parameter's 95%
profile-likelihood interval. It is bracketed from the same L, solved against the objective at
off-grid points, and truncated at the box with a flag, never extrapolated. It is serial,
deterministic, and for headline fits only.

### 5.2 Bootstrap (and every other resampling scheme)

```
L [T × K]  ◄── computed once, exactly as in 5.1
seed ──► Philox4x32-10(key = seed ⊕ "RSMPBOOT", counter = (scheme, b, j/2, 0)) ──► W[b, ·]   (D-132)
          iid bootstrap:   multinomial(T; 1/T) counts
          block bootstrap: moving-block counts, block length ℓ (default ⌈T^(1/3)⌉, D-133)
          jackknife:       W = 1·1ᵀ − I            (B = T)
          walk-forward:    W[b, t] = 1{t ∈ window_b}
          custom:          caller-supplied W, or a B × m period-index matrix → counts (D-134)
W [B × T] ──► distinct rows (D-172): R [D × K], the panel's distinct rows in ascending observation
              order ((n, d) for the binomial objective); M [B × D], M[b,j] = Σ_{t in row j} W[b,t]
          ──► bounded argmax (D-172), per b: tile bounds Σ_j M[b,j]·max_{k∈tile} R[j,k] (weights ≥ 0),
              tiles in descending bound order, each evaluated in full, until the next bound is
              strictly below the best value: exactly the full grid's argmax. Falls back to, and is
              checked against, reduce_weighted, fused & tiled over (b, k):
                s = Σ_j M[b,j]·R[j,k]   (ascending j, compensated summation)
                ArgMax.push(k, s); tile states merged in fixed order
          ──► per-b refinement pass ──► B × vcal_replicates ──► percentile intervals (type 7, D-135)
```

Cost: the quadrature (expensive) runs `T·K` times total, independent of `B`. The fused
reduction is `O(B·D·K)` multiply-adds with `O(D·K + B·D)` memory, D ≤ T the number of distinct
observations (16 on average against T = 53 over the recovery subset); `B·K` is never stored.
The order over rows is defined by the observations, so a replicate's estimate does not depend on
how the panel's periods are laid out; against a sum over periods it differs only by rounding.
Validity: this is exact for objectives additive over periods; see (D-042, D-043) for the AR(1)
factor and method-of-moments cases where it isn't a pure `W × L`.

### 5.3 Backtest

```
calibration run (ledger id) ──► PD forecast per period / grade
observed panel (n_t, d_t) ──► per-period and pooled tests:
      exact binomial · Jeffreys (Beta(d+½, n−d+½)) · correlation-adjusted binomial
      (quantile of D under the ASRF mixture, via core integrator at given ρ) · traffic-light zones
  ──► vcal_backtest_result[] ──► ledger record {run ids, input hashes, spec, result hash}
  ──► monitoring report data ──► rendering (D-045)
```

### 5.4 Later milestones (shape only)

- **Multi-segment**: stage 1 = per-segment 2-D fits (as 5.1); stage 2 = pairwise
  composite likelihood over a 1-D grid in inter-segment correlation `r`, using 2-D quadrature
  per period. Still per-period additive → fits `W × L`.
- **AR(1) factor**: binomial likelihood needs a forward filter over a factor grid; not
  per-period separable → resampling via parametric bootstrap through `dgp/` (D-043; the DGP is CPU-only until R-1 is revisited at M6).
- **Pluto–Tasche**: per grade, root of the monotone 1-D surface
  `P(D ≤ d_obs | PD, ρ) = 1 − γ` (a root-finding reducer).
- **Macro**: many small specification fits batched as the `B` dimension.

---

## 6. Determinism contract

1. Same inputs + seed + build + backend + profile ⇒ **bitwise-identical** outputs,
   independent of thread count and launch configuration. Achieved by: each output cell owns
   its `Σ_t` (fixed ascending order, no cross-thread sums); exactly-associative merges
   (ArgMax with index tie-break) where possible; fixed tile decomposition + fixed merge tree
   otherwise. OpenMP is used only as `parallel for` over independent outputs — no
   `reduction` clauses (portable to MSVC's OpenMP 2.0, and deterministic).
2. Across platforms (glibc vs MSVC UCRT): within documented tolerance, not bitwise, about 3e-15
   relative on estimates (D-099). CPU vs GPU: within documented tolerance (D-041).
3. DGP: runs on the CPU only (D-034); the context's backend never changes a simulated panel.
   Panels are bitwise identical across OSes: own elementary functions, `/fp:strict` or
   `-ffp-contract=off`, and a CI hash check (D-052–D-057).

## 7. Numerics and tolerances

- `docs/methodology/tolerances.md` is the register: each entry has an id, value, why that
  value, and the test that enforces it. Tests reference tolerances by id, never by literal.
  Values live in `tests/tolerances.toml`; the C++ header is generated from it, and validation
  scripts read it directly (D-125).
- Quadrature convergence is checked, not assumed: at calibration, the headline point is
  re-evaluated with the doubled rule (GH 256, composite K = 32). A period differing by more
  than max(1e-10, 64·ε·(abs(log C) + abs(log I))) sets `VCAL_FLAG_QUADRATURE_UNCONVERGED` (D-092,
  D-120).
- Parity's integrator is split by period (D-118): adaptive Gauss–Hermite, N = 128, for
  0 < d < n; composite Gauss–Legendre, 16 × 16, for d ∈ {0, n}, whose integrand is a truncated
  prior that polynomial rules centred on a mode resolve slowly.

## 8. Test layers

| Layer | Question it answers | Oracle |
|---|---|---|
| unit | Is each primitive right? | high-precision tables, closed forms, polynomial exactness, Philox KATs |
| crossref | Does core compute the same surface/estimate as an independent implementation? | `ref/` |
| recovery | Does the estimator recover truth on synthetic data at the expected rate? | DGP ground truth. Engine correctness must hold; interval coverage is reported as PASS / DEFERRED / reviewed KNOWN FINDING (D-121, D-124) |
| parity | Does GPU match CPU? | CPU backend |
| abi | Is the C ABI usable from plain C, and does it give the engine's numbers? | a C11 consumer and a C99 header compile; the ABI bitwise against direct engine calls; the export table against the header (D-144) |
| golden | Did anything change that shouldn't have? | checked-in values with provenance (D-031–D-033) |

## 9. Build

- CMake ≥ 3.25 (D-047). C++20 (D-025) except device-visible code, which is C++17 and checked by
  a `-std=c++17` TU (D-058), Ninja via `CMakePresets.json`
  (`cpu-debug`, `cpu-release`). On Windows, configure from a VS developer shell.
- Every target goes through `vcal_target_defaults()` (`cmake/VcalTargetDefaults.cmake`):
  C++20, warnings (`/W4 /permissive-` or `-Wall -Wextra -Wpedantic -Wconversion -Wshadow`),
  optional warnings-as-errors (D-049), and FP-contraction off (D-048). `dgp/` additionally gets `/fp:strict` on MSVC (D-056).
- CUDA (M4): `VCAL_ENABLE_CUDA` defaults ON only if a CUDA compiler is found (D-030).
  Architectures: `89-real;80-virtual`, plus `120-real` when nvcc ≥ 12.8 (D-026, D-040).
  The resolved list is printed at configure time; `vcal_build_info()` reports `cuda=off` until M4 embeds it.
  CUDA TUs are C++17 (D-058). Floor: nvcc 11.8 with its host-compiler range (D-065).
  Compile-only CI jobs (nvcc 11.8 and 12.8, GCC 11) build the C++17 header TU and a kernel
  calling every device-visible function, with host calls from device code as errors (D-100).
- `cxx17_header_check` compiles a generated TU of every header in the C++17 zone (D-064);
  `cxx17_header_check_fires` and `contract_fail_precision_message` are compile-fail tests
  proving that the check and the contract messages fire.
- `layering_check` (CTest) enforces the §2 layering table (D-050).
- The C ABI is the shared library `vcal` (D-142). Symbols are hidden by default. On Linux a
  linker version script exports `vcal_*` only, and `--no-undefined` catches unresolved symbols.
  `abi_exports` compares the export table (`nm -D` or `dumpbin /exports`) with the header's
  declarations, and `abi_exports_fires` proves that the comparison catches a difference.
