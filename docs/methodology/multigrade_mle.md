# The multi-grade binomial MLE with a shared ρ (M6)

**Status: the likelihood is built (PR 1, D-412); the estimator (§3, §4) and the DGP layer (§6) are still the
plan of D-411.** Sections 1, 2 and 5 describe what the code does, with the measured numbers; the later
sections say what PR 2 and PR 3 will do. The plan was fixed before any code (D-411), as the P-item and
study rules do for their expectations (D-177). Decision numbers (D-nnn) point to `DECISIONS.md`.

The code: `core/model/multigrade_mixture.hpp` (the integrand, its batch evaluation, its hint and its
parameter derivatives) and `core/objectives/multigrade.hpp` (the period log-likelihood and its score and
Hessian). Its tests are `unit_multigrade`.

The design is D-151's: a single-factor model with K grades, each with its own PD, and one asset
correlation ρ shared by all grades. The one-dimensional grid in ρ, the K PDs maximised for each fixed
ρ, and the revisited mode hint are taken from D-151 as decided. Miu and Ozdemir (2008, SSRN 1365842)
estimate a model of the same structure, a common correlation with grade-specific default points, by
maximum likelihood on each grade's default rates; this note is the count-based version, with each grade
binomial given the factor.

## 1. Model and likelihood

A panel has T periods. In period t grade k has n_kt obligors and d_kt defaults, k = 1, …, K. Every
obligor in every grade loads on the same factor Z_t ~ N(0, 1), independent across periods, with the same
ρ. Given Z_t = z, grade k defaults independently with

    p_k(z) = Φ((Φ⁻¹(PD_k) − √ρ z) / √(1 − ρ)),

so d_kt given z is Binomial(n_kt, p_k(z)), and the grades are independent given z. The grades are
therefore *not* independent: they share z. Period t's log-likelihood is

    l_t(PD_1, …, PD_K, ρ) = Σ_k log C(n_kt, d_kt) + log ∫ ∏_k p_k(z)^d_kt (1 − p_k(z))^(n_kt − d_kt) φ(z) dz,

and the panel's is ℓ = Σ_t l_t. As in the single-grade note (D-070), the binomial coefficients are
kept, so log-likelihoods, likelihood ratios and information criteria are the true values. Still one
integral per period, so the quadrature stays cheap; only the optimiser grows (D-151).

A grade with n_kt = 0 in a period contributes nothing to that period. A grade may be empty in some
periods and not others.

**Three exact identities** anchor the build, and each is a test:

1. **K = 1 is the single-grade likelihood.** With one grade, l_t is `binomial_mixture`'s log_contrib.
   The K-grade integrand performs that integrand's operations grade by grade, in the same order, and
   with one non-empty grade the hint is the single-grade hint itself, so the identity holds **bit for
   bit**, value and derivatives, under the parity rule and its doubled check: 2,384 cases (every
   single-grade mixture golden and a sweep over PD, ρ, n and d), alone and with empty grades inserted
   around the grade.
2. **Two grades with equal PDs pool.** Binom(d₁; n₁, p) · Binom(d₂; n₂, p) =
   [C(n₁, d₁) C(n₂, d₂) / C(n₁ + n₂, d₁ + d₂)] · Binom(d₁ + d₂; n₁ + n₂, p) for any p. So with PD₁ = PD₂
   the two-grade l_t equals the single-grade l_t on the pooled counts plus that log-ratio of
   coefficients. To rounding, not bit for bit (different sums): worst 5.8 ε of the terms' size over 144
   pairs, within `TOL_MULTIGRADE_POOLING_EPS` (24).
3. **ρ → 0 factorises.** At ρ = 0, p_k(z) = PD_k for every z and l_t = Σ_k log Binom(d_kt; n_kt, PD_k),
   in closed form. At ρ = 10⁻¹⁴ the K-grade l_t is within 1.8e-13 of it (`TOL_MULTIGRADE_RHO_ZERO_ABS`,
   1e-12), and the gap grows 100-fold at ρ = 10⁻¹²: ρ's own first-order effect, not error.

## 2. The integral

The log-integrand is

    h(z) = Σ_k [d_kt log p_k(z) + (n_kt − d_kt) log(1 − p_k(z))] − z²/2.

Each grade's term is concave in z: log Φ is concave, and log p_k(z) and log(1 − p_k(z)) are log Φ of
affine functions of z. With the strictly concave −z²/2, **h is strictly concave**, so it has at most one
peak, and the single-grade rules carry over with the period classified by its counts:

- **Interior period** (some grade has 0 < d_kt < n_kt, or the period has both a grade with d = 0 and one
  with d = n): h has an interior peak. **Adaptive Gauss–Hermite, N = 128**, centred on the mode z* and
  scaled by σ = 1/√(−h″(z*)), as the single-grade rule (D-118).
- **One-sided period** (every non-empty grade has d_kt = 0, or every one has d_kt = n_kt): no interior
  peak; the normal prior is cut by the product of the grades' survival factors, ∏_k (1 − p_k)^n_k (or
  ∏_k p_k^n_k). **The composite Gauss–Legendre rule, 16 × 16**, with its transition point z_c where that
  product equals ½ and its width from the product's slope there, which reduces to the single-grade z_c
  and w at K = 1.

**The mode hint, revisited (D-151).** h′ is the sum of the grades' h′ terms, and h″ the sum of theirs,
so the existing Newton-with-bisection solve applies unchanged to the summed derivatives. Its start is
the precision-weighted mean of the grades' own single-grade starts, Σ τ_k z_k / (Σ τ_k + 1), with z_k
where p_k(z_k) = (d_k + ½)/(n_k + 1) and τ_k = −g_k″(z_k), shrunk toward the prior mean exactly as the
single-grade hint shrinks its start. (D-411 had proposed the z where the n-weighted mean of p_k equals the
pooled rate, one root solve; the weighted mean needs none and reduces to the single-grade start at K = 1.)

**The one-sided transition.** log S(z) = Σ_k n_k log Φ(∓x_k(z)) is monotone in z, so its root at −ln 2 is
bracketed from the mode by doubling steps and found by bisection to adjacent doubles; the width is
1/|d log S/dz| there. At K = 1 the closed-form single-grade z_c and w are used (through the single-grade
hint), so that path is unchanged.

**The risk the plan named, and how it was settled.** In a mixed interior period, a grade with d = 0 and
large n puts a steep cliff into the integrand next to the peak that other grades set. Its log is still
smooth and concave, and its curvature enters h″ and so the scale σ, but whether 128 Gauss–Hermite nodes
resolve it to the parity tolerance is what the single-grade D-116 finding warned about. It was measured
against 65 mpmath goldens (50 digits, `tools/gen_multigrade_goldens.py`): mixed periods with a d = 0
grade of n = 10⁴, 10⁵ and 10⁶ beside one or two interior grades at ρ from 0.02 to 0.9; two-grade interior
periods; d = 0 beside d = n; one-sided periods of up to five grades; D-151's 17-grade master scale in a
benign, an average and an adverse year at ρ of 0.02, 0.12 and 0.24; and five grades at ρ up to 0.9.
**The interior rule holds:** worst 3.6 ε, at exactly the named case (ρ = 0.9, a d = 0 grade of n = 10⁶
beside two interior grades), against `TOL_PARITY_MIXTURE_EPS` of 24; one-sided periods 5.6 ε. So the
classification stands as written, and no period moves to the composite rule. The period's rule depends
only on the counts, never on the parameters, so the likelihood is a smooth function of the parameters
(as the single-grade note requires).

**The derivatives** (D-170, extended): with E the expectation under the period's posterior for z, the
score is E[g_i] and the Hessian E[g_ij + g_i g_j] − E[g_i] E[g_j], on the value's own nodes. Grade k's term
depends on (PD_k, ρ) only, so g_{PD_k PD_l} = 0 for k ≠ l, and the Hessian's PD-PD cross terms are the
posterior covariances of the grades' scores: the matrix is dense, as §3 expects. Against
Richardson-extrapolated central differences in logit coordinates on every golden: 1.2e-10 (score) and
1.3e-9 (Hessian) for grades of n ≤ 10⁵, 1.2e-9 and 3.8e-8 above, within the single-grade tolerances
`TOL_SCORE_HESSIAN_FD_REL` and `..._LARGE_N`. The test's step is half the single-grade test's: at its 0.02
the 17-grade periods' ρ-ρ entry is limited by the differences' own O(h⁴) truncation, which falls 16-fold
per halving of the step.

**Cost.** Each node evaluates K log Φ pairs instead of one, through the existing erfcx kernels (one
batch per grade, P-7). At K = 17 a period is about 17 times a single-grade period's work, so a panel of
T = 100 periods on 17 grades costs about what 1,700 single-grade periods do. The derivatives' node
functions number (K + 1) + (K + 1)(K + 2)/2, 189 at K = 17. Not optimised yet: PR 2's estimator will be
timed, and the per-node erfcx sharing of P-21 is not used by the K-grade derivatives.

## 3. The estimate

Built in PR 2 (D-413) as `benchmark_models/multigrade/fit.hpp`, `multigrade::fit`. Parameters θ = (u_1, …, u_K, v),
with u_k = logit(PD_k) and v = logit(ρ).

- **Boxes.** Each PD_k in [1e-5, 0.5], ρ in [1e-3, ρ_hi] with ρ_hi = `kDefaultRhoUpper` (0.5), logit axes
  (`multigrade::Box`). The PD box is wider than the single-grade grid's [1e-4, 0.2], because a 17-grade master
  scale starts at 0.03% (D-151's scenario matrix), and above the floor of `VCAL_GRID_MIN_PROBABILITY` (D-218).
  One PD box for every grade: the plan's per-grade boxes are not built, since nothing yet needs them (D-413).
- **Outer: ρ on a grid (D-151).** 41 points uniform in v, endpoints exactly the box's. At each v_j the inner solve
  gives the profile ℓ_p(v_j) = max over the PDs of ℓ, and dℓ_p/dv, the v-component of the score at the inner
  maximum (the envelope theorem).
- **Inner: the K PDs at fixed v.** Safeguarded Newton in (u_1, …, u_K) with the analytic Hessian, which is
  dense: the grades share z, so the cross terms are the covariance, under each period's posterior for z, of the
  grades' scores. The step solves (−H + λI) s = g with λ raised from 0 until Cholesky succeeds; then a
  backtracking line search on ℓ, projected onto the box. A grade on a bound whose gradient, or whose step, points
  out of the box is held there (an active set). The solve stops when no coordinate moves by 10⁻¹⁰.
  - Start at the first grid point (ρ = 10⁻³): each grade's pooled rate, (Σ_t d_kt + ½)/(Σ_t n_kt + 1), clamped
    into the box.
  - Later grid points: warm-started from the previous one's solution.
- **Refinement of ρ.** From the grid argmax, safeguarded Newton on ℓ_p(v) inside the bracket of its two
  neighbours: dℓ_p/dv from the envelope theorem, d²ℓ_p/dv² the Schur complement H_vv − H_vu H_uu⁻¹ H_uv on the free
  grades, and a bisection step on the sign of dℓ_p/dv whenever Newton's step leaves the bracket or the curvature
  is not negative; to 10⁻¹⁰ in v. At a box end whose dℓ_p/dv points out of the box, ρ is that end. Then one joint
  Newton step in every free coordinate, kept only if it raises ℓ. (The plan's quadratic vertex before Newton is
  not used: Newton from the grid point, safeguarded by the bracket, converges in a few steps without it.)
- **Derivatives.** PR 1's analytic score and Hessian of l_t in θ (`multigrade_log_contrib_derivs`), summed over
  the periods and chained to (u, v): H_u,ij = J_i J_j H_θ,ij + δ_ij g_θ,i J_i (1 − 2θ_i), with J_i = θ_i(1 − θ_i).
- **Cost, measured.** A fit is 470 to 610 panel evaluations (value only, or value with score and Hessian), with
  about 40% of them taking derivatives. On one core of the cloud machine, a 2-grade fit at T = 20 takes about 0.1 s.
  A 17-grade fit at T = 100 takes 7.5 s, nearly all of it in the derivative evaluations: each node folds
  (K + 1) + (K + 1)(K + 2)/2 = 189 moments.

## 4. Standard errors, intervals and flags

- **SEs:** the observed information, −H from the analytic Hessian at the estimate, on the free coordinates (the
  grades off their bounds, and v unless ρ is on a bound), inverted by Cholesky in (u, v), then the delta method to
  the natural scale (as P-15, D-185). A coordinate on a bound has SE NaN. If −H is not positive definite, every
  SE is NaN and the fit is flagged `kFlagFlatSurface`.
- **Intervals:**
  - **ρ: the profile-likelihood interval** at `kProfileThreshold95`. Each end is bracketed by the grid's own
    profile values and solved by safeguarded Newton on ℓ_p − (ℓ̂ − χ²₁,₀.₉₅/2), to `kProfileRootTol` in v. If the
    profile stays above the threshold out to a box end, that end is the bound, flagged truncated (D-128–D-130).
  - **Each PD_k: the logit-scale Wald interval,** u_k ± 1.96 SE_u. NaN for a grade on a bound. Per-grade profile
    intervals stay a later item (§8). Each end needs a K-dimensional solve (the other PDs and ρ) at every Newton
    step on it. A 17-grade fit at T = 100 already takes 7.5 s, and 34 ends of such solves would cost several
    times that fit, so they wait until a study needs them.
- **Flags.** Fit-level, in `Fit::flags`:
  - `kFlagGradeEdge`: some grade's PD is on a box bound. Per grade, `kGradeAtLower`, `kGradeAtUpper` and
    `kGradeEmpty` (n = 0 in every period: not estimated, PD NaN) are in `Fit::grade_flags`. The expected edge
    case is a grade with no defaults over the whole panel, whose MLE is 0. Its PD is reported as the box's lower
    end exactly, its SE and interval NaN; nothing is clamped silently (owner's choice, 2026-11-10). ρ is still
    estimated from every grade.
  - engine `kFlagGridEdge`: ρ on a box end. `kFlagFlatSurface`: −H not negative definite.
    `kFlagQuadratureUnconverged`: the doubled rule at the estimate differs by more than
    max(`kQuadratureCheckFlagAbs` · T, 256 ε Σ_t |l_t|), as calibrate's check (D-120, D-178).
  - `kFlagNotConverged`: an inner solve or the refinement hit its iteration limit. `kFlagNumeric`: a profile
    value is not finite.

## 5. The independent reference and the replication

- **`ref/`:** `vcalref::log_mixture_grades` and `period_loglik_grades`, the K-grade integral by ref's own
  method: the maximum bracketed and found by golden-section search, then adaptive Gauss–Legendre on panels
  around the mode sized by the half-drop widths, out to e⁻⁶⁰. That is the single-grade `log_mixture`'s body,
  now shared by the two through one helper; it shares no quadrature code with core (D-031–D-033). (The plan
  had described ref's method as Brent's mode and a trapezoid window, which is not what ref does.) Core
  against ref on every K-grade golden: worst 3.5 ε, within `TOL_XREF_PERIOD_LL_EPS`; ref against mpmath:
  3.1 ε, within `TOL_REF_LOG_MIXTURE_EPS`. Both are the single-grade tolerances, unchanged. The comparison
  is in `unit_multigrade`.
- **Goldens from mpmath:** `tools/gen_multigrade_goldens.py` writes `tests/golden/multigrade/likelihood.csv`
  (65 cases, §2) with its manifest, and CI re-generates and compares it. Each value is `mp.quad` on h at 50
  digits, split at the mode ± k·σ and at every one-sided grade's own transition point, so a steep cliff
  beside the peak falls on a split point. Before writing, it checks that K = 1 with an empty grade and two
  grades with one PD give the single-grade integral to 10⁻³⁰. The two new tolerances
  (`TOL_MULTIGRADE_POOLING_EPS`, `TOL_MULTIGRADE_RHO_ZERO_ABS`) were set from their measurements, worst × 4
  rounded up, as the register's other entries; the plan's "stated before they are measured" was not kept,
  and every other check uses an existing tolerance.
- **`validation/scipy/multigrade_mle.py`:** the likelihood half. With `binom.logpmf` products and
  `scipy.integrate.quad` on the e⁻⁵⁰ window around the mode found by Brent's method, it reproduces every
  golden row plus its binomial coefficients to 1.5e-13 relative (`TOL_SCIPY_MULTIGRADE_LL_REL`, 1e-12). The
  estimate half (PR 2): on the two panels of `tests/golden/multigrade/estimates.csv` (written by
  `tests/studies/multigrade_goldens.cpp`, K = 2, T = 12 and K = 3, T = 10), `scipy.optimize.minimize`
  (L-BFGS-B in the logit coordinates inside the same box, from each grade's pooled rate and ρ = 0.1, with
  central-difference gradients of the `quad` likelihood) reaches the engine's estimate to 6.1e-10 in every
  logit coordinate (`TOL_SCIPY_MULTIGRADE_FIT_LOGIT`, 1e-4). It never beats the engine's log-likelihood
  (`TOL_SCIPY_MULTIGRADE_FIT_LL_ABS`, 1e-8), and its log-likelihood at the engine's estimate matches the
  committed value to 1.3e-13. `unit_multigrade_fit` refits the same panels on every platform against the
  committed estimates (`TOL_MULTIGRADE_FIT_REGRESSION_LOGIT`, 1e-7).

## 6. The data-generating process (PR 3)

Built in PR 3 (D-414). The design is D-150's, for S-24's segments; D-151 shares it with S-29 and S-31 to
S-33.

- **`dgp::simulate_panel_segments`** (`docs/methodology/dgp.md`, "Segments on one factor"). Each period
  has one Z_t, the base DGP's factor on the base key, and then each segment's defaults given Z_t by the
  base DGP's Bernoulli step.
  - Segment 0 draws on the base key. So one segment is `simulate_panel`'s panel bit for bit, and grade 1
    of a grade panel is the base panel at its PD.
  - Segment s ≥ 1 draws on the key of seed ⊕ `0x5345474D454E5453` ("SEGMENTS") ⊕ s, its own key domain
    (D-149). The segments are independent given Z_t.
  - A segment's draws do not depend on the other segments' PDs, sizes or number.
- **The kernels are unchanged.** The segment key goes through `count_defaults` like any seed, so every
  Bernoulli kernel already reproduces it.
- **Mirror and hash.** The Python mirror (`Mirror.panel_segments` in `tools/gen_dgp_tables.py`) writes six
  reference panels with their SHA-256, matched bit for bit under every kernel. The six span K = 1 (base
  panel 0), 2, 3, 5 and the 17-grade M scale.
- **The grade scenario matrix** (D-151, `studies/README.md`), in `tests/recovery/grade_matrix.hpp`: the
  master scale, the obligor rule (shares by `dgp::det_exp`, so the counts are the same on every
  platform), and the nine scenarios on seed "GRADESCN", with grade k as segment k − 1.
- **Not in PR 3: D-151's observation layer.** That layer relabels, misrecords and drops what is recorded,
  for S-29 and S-31 to S-33. Each of those studies defines its own thinning in its `PREDICTION.md`, so the
  layer is built with the first of them.

## 7. The checks

Every item ships with synthetic recovery tests and golden values (CLAUDE.md).

- **PR 1, the likelihood (done; `unit_multigrade`, six tests, and the scipy script):**
  - the three identities of §1: K = 1 bit for bit on both rules (2,384 cases), pooling (5.8 ε) and
    ρ → 0 (1.8e-13);
  - core against the 65 mpmath goldens within `TOL_PARITY_MIXTURE_EPS`: 3.6 ε interior, 5.6 ε one-sided,
    so the rule is not amended;
  - core against `ref/` (3.5 ε) and ref against mpmath (3.1 ε);
  - the derivatives against central differences (§2), and at K = 1 against `log_contrib_derivs` bit for
    bit;
  - the scipy replication of the likelihood (1.5e-13).
- **PR 2, the estimator (done; `unit_multigrade_fit`, seven tests and the slow recovery test, and the scipy
  estimate half).** Tolerances were stated before the tests ran (D-413):
  - **At K = 1, against `calibrate` polished by `profile_intervals`:** four panels. The estimate is within
    6.1e-10 in logit units (tolerance 1e-5), the log-likelihood within 1.1e-13 (1e-9), and the ρ interval's
    ends within 9.7e-11 (1e-6). The SEs match the engine's analytic SE at the same point to rounding (1e-8).
    The stated comparison against calibrate's own SE was the wrong comparand: that SE is central differences
    of the surface (D-119). It differs by 0.1% to 2.1% and is reported only.
  - **Pooling, as an exact identity:** a duplicated grade gives PD₁ = PD₂ (a split of 0) and the K = 1 fit on
    the doubled counts, to 4.5e-10 in logit units. The log-likelihoods agree to 9.1e-13 after the binomial
    coefficients. This replaces the plan's grades generated with one PD, which agree only to sampling error.
  - **A zero-default grade** (n = 100 or 2,000, beside two interior grades) has its PD exactly at the box's
    lower end, flagged, with SE and interval NaN. ρ moves by 0.005 and 0.17 SE from the fit without it
    (tolerance 0.25).
  - **An empty grade** inserted anywhere changes nothing, bit for bit; the fit is deterministic.
  - **The estimate is a maximum of ref's independent likelihood:** every 10⁻³ step lowers it, and ref and core
    agree at the estimate to 6e-13.
  - **The two committed panels** refit exactly on this machine (`TOL_MULTIGRADE_FIT_REGRESSION_LOGIT`, 1e-7).
  - **Recovery,** R = 100 per cell, on panels drawn in the test from `dgp::factor_draw` and
    `dgp::count_defaults` until PR 3. Three scales: 2 grades; 5 grades; and D-151's 17-grade master scale with
    10,000 obligors. Each at T ∈ {20, 100} and ρ ∈ {0.24, 0.12, 0.02}. All 18 cells pass:

    | Scale | T | ρ profile coverage | PD Wald coverage | Bias of logit ρ̂ (sd) |
    |---|---|---|---|---|
    | 2 grades | 20 | 0.94–0.95 | 0.915–0.95 | −0.06 to −0.21 (0.34–0.46) |
    | 2 grades | 100 | 0.93–0.95 | 0.92–0.96 | −0.02 to −0.03 (0.17–0.18) |
    | 5 grades | 20 | 0.90–0.98 | 0.91–0.946 | −0.07 to −0.12 (0.33–0.38) |
    | 5 grades | 100 | 0.97 | 0.93–0.958 | −0.02 to −0.03 (0.12–0.14) |
    | 17 grades | 20 | 0.89–0.92 | 0.919–0.948 | −0.10 to −0.19 (0.33–0.38) |
    | 17 grades | 100 | 0.95–0.97 | 0.953–0.977 | −0.00 to −0.01 (0.14–0.15) |

    The band's lower end is 0.878. No fit is flagged, none has ρ on the box, and none hit an iteration limit.
    At T = 20 the profile interval undercovers slightly, most on the 17-grade scale (0.89 at ρ = 0.24), as the
    single-grade fit does (0.937, `bartlett_profile.md`). ρ̂ is biased low by 0.06 to 0.21 in logit units at
    T = 20, and by at most 0.03 at T = 100. The run takes about 70 CPU-minutes, so it
    is labelled slow; `VCAL_MG_RECOVERY_R`, `_SCALE` and `_RHO` size and split it.
  - **The scipy replication of the estimate:** 6.1e-10 (§5).
- **PR 3, the DGP layer (done; D-414):**
  - the segment panels match the Python mirror bit for bit under every kernel (scalar, SSE2, AVX2), and
    their SHA-256 too;
  - one segment is the base panel bit for bit, and Z_t is the base factor;
  - each segment's count is `count_defaults` on its own key;
  - a segment's counts do not move when another segment's PD or size changes, or a segment is added;
  - two segments with one PD draw different uniforms;
  - invalid input is refused;
  - the grade matrix's counts equal the studies index's lists, summing to 10,000, with its portfolio PDs;
  - in CI, replicate 0 of the six T ∈ {20, 40} scenarios is fitted, converged, finite and not flat
    (11 s);
  - **Recovery on the nine grade scenarios** (slow, R = 100, D-413's assertions, expectation stated before the
    run). Eight pass. The table shows ρ's profile coverage, the PDs' Wald coverage, the bias of logit ρ̂ (sd),
    and the fits with a grade on the PD box:

    | Scenario | Scale, T | ρ profile | PD Wald | Bias (sd) | Grade on the box |
    |---|---|---|---|---|---|
    | 0 | L, 20 | 0.94 | 0.934 | −0.09 (0.35) | 11 |
    | 1 | L, 40 | 0.95 | 0.970 | −0.07 (0.23) | 2 |
    | 2 | L, 100 | 0.98 | 0.961 | −0.01 (0.13) | 0 |
    | 3 | M, 20 | 0.95 | 0.928 | −0.12 (0.32) | 33 |
    | 4 | M, 40 | 0.95 | 0.955 | −0.04 (0.24) | 10 |
    | 5 | M, 100 | 0.98 | 0.947 | −0.02 (0.13) | 0 |
    | **6** | **H, 20** | **0.87** | 0.951 | −0.09 (0.42) | 89 |
    | 7 | H, 40 | 0.94 | 0.948 | −0.05 (0.25) | 65 |
    | 8 | H, 100 | 0.93 | 0.952 | −0.01 (0.16) | 25 |

    No fit failed to converge, gave a non-finite value or a flat surface, and none put ρ on the box.

    **Scenario 6** (H, T = 20, ρ = 0.02) is below the band, 0.878. Rerun at R = 400 its coverage is 0.9125,
    and the band there is 0.914. Its misses lean low: 24 of the 35 intervals lie below the true ρ.
    - **The grade on the box is not the cause.** 358 of the 400 fits have the 0.03% grade (46 obligors,
      0.28 expected defaults) on the box, but only 28 of the 35 misses do, so it is not over-represented.
    - **It is the single-grade fit's known finding.** At T = 20 and ρ = 0.02, the single-grade profile
      interval undercovers too: 0.927 at PD 1% with 10,000 obligors, and 0.926 at PD 5% with 1,000. Both
      are reviewed KNOWN FINDINGs in the pinned recovery (D-154 expects about a point of undercoverage at
      T = 20). Scenario 6 is about one Monte Carlo standard error (0.014) below those.
    - **For comparison,** scenario 3 at R = 400 covers 0.9525.

    The cell is a finding of the method, not of the multi-grade fit. The owner made it a reviewed KNOWN
    FINDING (2026-11-11, D-414), as in the single-grade recovery. Scenario 6's profile verdict is listed in
    `kGradeMatrixKnownFindings`. Any other cell below the band fails the test, and so does scenario 6 if it
    returns to the band, so the list is revisited.

## 8. The PRs, and what is not in this work

1. **PR 1, the likelihood** (§1, §2, §5's ref, goldens and the likelihood half of the scipy script): done
   (D-412).
2. **PR 2, the estimator** (§3, §4, the estimator half of the scipy script, recovery on simple panels): done
   (D-413).
3. **PR 3, the DGP layer and the grade scenario matrix** (§6), with recovery on it: done (D-414).
4. **Later, each its own decision:** per-grade profile intervals (§4); a C ABI entry point
   (`vcal_calibrate_multigrade`), after the recovery study has run (owner's choice, 2026-11-10); the
   studies that need the model (S-29 full, S-30 shared ρ, S-40), each registered before its run; and the
   multi-grade Pluto–Tasche sequence (S-52), which needs PR 3's layer.

**Not in this model:** grade-specific ρ (S-40 compares it, as a separate estimator), rating migration
(Miu and Ozdemir model it; here the grade counts are taken as given), monotone PDs across grades (not
imposed; owner's choice, 2026-11-10), and serial correlation of the factor (the AR(1) model's; combining
the two is its own item).

## 9. Open points for the owner

PR 2 took the defaults proposed here, as no other answer had been given (D-413):

- **The PD box** [1e-5, 0.5] (§3), chosen for the 17-grade scale's 0.03% top grade: taken, one box for every
  grade.
- **Wald intervals for the PDs,** with profile intervals later (§4): taken.
- **Parity or native.** The model is textbook and scipy-replicable, so it fits the `parity` profile's
  definition (D-035); but parity's whitelist today holds only the single-grade estimators. `native` until the
  recovery study has run, then reconsidered with the ABI entry point: taken. The estimator has no ABI entry
  point, so nothing in the profiles changes yet.

Each can still be changed by the owner before PR 3.
