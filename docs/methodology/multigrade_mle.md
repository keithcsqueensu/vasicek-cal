# The multi-grade binomial MLE with a shared ρ (M6): the plan

**Status: plan, written before the build (D-411).** Nothing here is built yet. This note fixes the model,
the likelihood, the estimator and the checks before any code, as the P-item and study rules do for
their expectations (D-177). It becomes the methodology note once PR 1 lands; then each section is
rewritten from "will" to "does" and its measured numbers are filled in. Decision numbers (D-nnn) point
to `DECISIONS.md`.

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
   The K-grade integrand will be written so that at K = 1 it performs the same operations in the same
   order, so the identity holds **bit for bit**, on both rules.
2. **Two grades with equal PDs pool.** Binom(d₁; n₁, p) · Binom(d₂; n₂, p) =
   [C(n₁, d₁) C(n₂, d₂) / C(n₁ + n₂, d₁ + d₂)] · Binom(d₁ + d₂; n₁ + n₂, p) for any p. So with PD₁ = PD₂
   the two-grade l_t equals the single-grade l_t on the pooled counts plus that log-ratio of
   coefficients. To rounding, not bit for bit (different sums).
3. **ρ → 0 factorises.** At ρ = 0, p_k(z) = PD_k for every z and l_t = Σ_k log Binom(d_kt; n_kt, PD_k),
   in closed form. Checked at the axis's smallest ρ against the limit.

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
so the existing Newton-with-bisection solve applies unchanged to the summed derivatives. Its start
moves from "where p(z) = (d + ½)/(n + 1)" to the pooled analogue: the z at which the n-weighted mean of
p_k(z) equals (Σd + ½)/(Σn + 1). One root solve (bisection on a monotone function), then Newton on h′.

**The risk the plan names, and how it is settled.** In a mixed interior period, a grade with d = 0 and
large n puts a steep cliff into the integrand next to the peak that other grades set. Its log is still
smooth and concave, and its curvature enters h″ and so the scale σ, but whether 128 Gauss–Hermite nodes
resolve it to the parity tolerance is exactly what the single-grade D-116 finding warns about, and has
to be measured, not assumed. PR 1 measures the K-grade rule against mpmath on a stress matrix: mixed
periods with a d = 0 grade of n up to 10⁶ beside interior grades, ρ up to 0.9, K up to 17. **If the
interior rule misses `TOL_PARITY_MIXTURE_EPS` on any such case,** mixed periods move to the composite
rule on h's e⁻⁵⁰ window centred on the peak, and the classification above is amended in this note before
PR 1 merges. Either way, the period's rule depends only on the counts, never on the parameters, so the
likelihood stays a smooth function of the parameters (as the single-grade note requires).

**Cost.** Each node evaluates K log Φ pairs instead of one, through the existing erfcx kernels (P-7,
P-21). At K = 17 a period is about 17 times a single-grade period's work, so a panel of T = 100 periods
on 17 grades costs about what 1,700 single-grade periods do.

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

- **`ref/`** gains the K-grade period log-likelihood with its own rule (Brent's mode and a trapezoid
  window, as its single-grade code), sharing no quadrature code with core (D-031–D-033). Core and ref
  are compared on identical inputs in `tests/crossref/`.
- **Goldens from mpmath** (`tools/`): the K-grade period log-likelihood at 50 digits on the stress matrix
  of §2, and the three identities of §1. Tolerances are registered in `tolerances.md`, the new ones
  stated before they are measured.
- **`validation/scipy/multigrade_mle.py`:** the likelihood by `scipy.integrate.quad` on h, and the
  estimate by `scipy.optimize.minimize` (L-BFGS-B in the logit coordinates, from the same start) on two
  small committed panels, against the engine within stated tolerances. Short enough to read in a sitting.

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

- **PR 1, the likelihood:**
  - the three identities of §1 (K = 1 bit for bit on both rules; pooling and ρ → 0 to rounding);
  - core against mpmath goldens on the stress matrix, within `TOL_PARITY_MIXTURE_EPS`, or the rule
    amended as §2 says;
  - core against `ref/` on the same inputs;
  - the derivatives against central differences of the value, and at K = 1 against
    `log_contrib_derivs` bit for bit;
  - the scipy replication of the likelihood.
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

1. **PR 1, the likelihood** (§1, §2, §5's ref, goldens and the likelihood half of the scipy script).
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
