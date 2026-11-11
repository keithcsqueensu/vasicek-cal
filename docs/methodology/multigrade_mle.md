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

Parameters θ = (u_1, …, u_K, v), with u_k = logit(PD_k) and v the ρ axis's scaled coordinate.

- **Boxes.** Each PD_k in [1e-5, 0.5], ρ in [1e-3, ρ_hi] with ρ_hi = `kDefaultRhoUpper` (0.5), logit
  axes. The PD box is wider than the single-grade grid's [1e-4, 0.2], because a 17-grade master scale
  starts at 0.03% (D-151's scenario matrix), and above the floor of `VCAL_GRID_MIN_PROBABILITY` (D-218).
  Configurable per grade.
- **Outer: ρ on a grid (D-151).** 41 points on the logit ρ axis, as the recovery grid's ρ axis. At each
  ρ_j the inner solve gives the profile ℓ_p(ρ_j) = max over the PDs of ℓ.
- **Inner: the K PDs at fixed ρ.** Safeguarded Newton in (u_1, …, u_K) with the analytic K × K Hessian.
  The grades share z, so the Hessian is dense, not diagonal: the cross terms are the covariance, under
  each period's posterior for z, of the grades' scores. K ≤ 17 makes a dense Cholesky trivial.
  Backtracking line search; the box handled by an active set (a grade pinned at its bound leaves the
  free set while its gradient points outwards).
  - Start at the first ρ: each grade's pooled rate, (Σ_t d_kt + ½)/(Σ_t n_kt + 1), clamped into the box.
  - Later ρ: warm-started from the previous ρ's solution.
- **Refinement of ρ.** The grid argmax of ℓ_p, then a quadratic vertex through it and its neighbours,
  then Newton on ℓ_p. By the envelope theorem dℓ_p/dρ is ∂ℓ/∂ρ at the inner optimum; d²ℓ_p/dρ² is the
  Schur complement of the full Hessian, H_vv − H_vu H_uu⁻¹ H_uv. Then one final joint Newton step in all
  K + 1 coordinates from that point, accepted only if it raises ℓ.
- **Derivatives.** Analytic first and second derivatives of l_t in all K + 1 coordinates, from the
  same nodes as the value, extending D-170's `log_contrib_derivs` to the product integrand: each grade's
  derivative of log p_k and log(1 − p_k) at each node, combined under the node weights. K = 1 must
  reproduce `log_contrib_derivs` bit for bit.

## 4. Standard errors, intervals and flags

- **SEs:** the observed information, the full (K + 1) × (K + 1) analytic Hessian at the estimate,
  inverted in scaled coordinates, then the delta method to the natural scale (as P-15, D-185). Where it
  is not negative definite, the SEs are NaN and the fit is flagged, as calibrate does.
- **Intervals:**
  - **ρ: the profile-likelihood interval,** read off ℓ_p directly. The ends are solved off-grid by
    Newton on ℓ_p − (ℓ̂ − χ²₁,₀.₉₅/2) and truncated at the box with a flag, as `profile_intervals` does
    (D-128–D-130).
  - **Each PD_k: the logit-scale Wald interval** in PR 2. Per-grade profile intervals need a
    (K)-dimensional solve per point and per end; they are listed as their own later item (§8), decided
    on PR 2's measured cost.
- **Flags, per grade and for the fit:**
  - **edge**, per grade: PD_k at its box bound. The expected case is a grade with no defaults over the
    whole panel, whose MLE is 0. Its PD is the box's lower end, its SE NaN, and the flag says so;
    nothing is clamped silently (owner's choice, 2026-11-10). ρ is still estimated from all grades.
  - **ρ at the box edge**, **flat profile**, **quadrature check** (each period re-integrated by the
    doubled rule at the estimate, as calibrate's check, D-178), **not negative definite**, and **inner
    solve not converged** at any ρ.

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
  estimate half, by `scipy.optimize.minimize` (L-BFGS-B in the logit coordinates, from the same start) on
  two small committed panels, arrives with the estimator in PR 2.

## 6. The data-generating process (PR 3)

D-151 asks for **one DGP observation layer** shared by S-29 and S-31 to S-33, with one prose
description, one Python mirror and one reference-panel hash (D-149). `dgp::simulate_panel` draws one PD
per panel; the multi-grade draw is new: per period one Z_t, then each grade's defaults given Z_t, on a
counter-based Philox stream whose counter layout is fixed in `docs/methodology/dgp.md` before the code,
bitwise identical on every platform and to the Python mirror (D-052–D-057). The grade scenario matrix of
D-151 (seed "GRADESCN", 3 scales × T ∈ {20, 40, 100}, the 17-grade master scale PD_k = 0.03% × 1.5^(k−1),
10,000 obligors per period spread by a discretised normal of standard deviation 4 grades, ρ of 0.24, 0.12
and 0.02, R = 1,000) is written into `studies/README.md` with that layer.

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
- **PR 2, the estimator:**
  - at K = 1, the estimate against `calibrate`'s within solver tolerance (different searches: a 2-D grid
    against a ρ grid with an inner Newton), and the ρ profile interval against `profile_intervals`';
  - pooling: two grades generated with one PD fitted as two give PD₁ ≈ PD₂ and the pooled single-grade
    ρ̂, within solver tolerance;
  - a zero-default grade lands on the box edge with the flag, and ρ̂ is unchanged when it is removed
    beyond what its counts carry (stated tolerance);
  - recovery on synthetic panels (until PR 3, from `simulate_panel` per grade on a shared Z_t drawn in
    the test): K ∈ {2, 5, 17}, T ∈ {20, 100}, the three ρ of the scenario matrix; bias and Wald and
    profile coverage reported against tolerances stated before the run;
  - the scipy replication of the estimate on the two committed panels.
- **PR 3, the DGP layer:** the counter layout's known-answer tests, the Python mirror's bitwise match and
  the reference-panel hash; then recovery on the grade scenario matrix's subset in CI.

## 8. The PRs, and what is not in this work

1. **PR 1, the likelihood** (§1, §2, §5's ref, goldens and the likelihood half of the scipy script): done
   (D-412).
2. **PR 2, the estimator** (§3, §4, the estimator half of the scipy script, recovery on simple panels).
3. **PR 3, the DGP layer and the grade scenario matrix** (§6), with recovery on its subset.
4. **Later, each its own decision:** per-grade profile intervals (§4); a C ABI entry point
   (`vcal_calibrate_multigrade`), after the recovery study has run (owner's choice, 2026-11-10); the
   studies that need the model (S-29 full, S-30 shared ρ, S-40), each registered before its run; and the
   multi-grade Pluto–Tasche sequence (S-52), which needs PR 3's layer.

**Not in this model:** grade-specific ρ (S-40 compares it, as a separate estimator), rating migration
(Miu and Ozdemir model it; here the grade counts are taken as given), monotone PDs across grades (not
imposed; owner's choice, 2026-11-10), and serial correlation of the factor (the AR(1) model's; combining
the two is its own item).

## 9. Open points for the owner

None blocks PR 1. Recorded so they are decided before the PR that needs them:

- **The PD box** [1e-5, 0.5] (§3): chosen for the 17-grade scale's 0.03% top grade. PR 2.
- **Wald intervals for the PDs in PR 2,** with profile intervals later (§4). PR 2.
- **Parity or native.** The model is textbook and scipy-replicable, so it fits the `parity` profile's
  definition (D-035); but parity's whitelist today holds only the single-grade estimators. Proposed:
  `native` until the recovery study has run, then reconsidered with the ABI entry point.
