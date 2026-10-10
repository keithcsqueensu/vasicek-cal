# The AR(1)-aware binomial MLE (M6)

This note describes, in plain prose, how vasicek-cal estimates PD, the asset correlation ρ and the
factor's serial correlation φ from a series of default counts when the systematic factor follows an
AR(1) process. It is the M6 benchmark model of ARCHITECTURE §5.4 and the estimator of S-35's M6 arm.
Decision numbers (D-nnn) point to `DECISIONS.md`.

The accompanying script `validation/scipy/ar1_binomial_mle.py` reproduces the likelihood with ordinary
scipy tools, three ways and in two formulations of the factor's law, none sharing the engine's node
placement (section 7).

## 1. Model and likelihood

Period t has n_t obligors and d_t defaults. Given the factor Z_t, defaults are binomial with the
one-factor Vasicek conditional PD

    p(z) = Φ( (Φ⁻¹(PD) − √ρ · z) / √(1 − ρ) ),

as in the binomial-mixture MLE (`binomial_mixture_mle.md`). The factor is AR(1) with unit marginal
variance:

    Z_1 ~ N(0, 1),   Z_t | Z_{t−1} ~ N(φ Z_{t−1}, 1 − φ²),   |φ| < 1,

so every Z_t is N(0, 1) and φ = 0 is the independent model. This is the DGP's AR(1) variant
(`dgp.md`, "The AR(1) factor"; D-190). The likelihood is a T-dimensional integral,

    L(PD, ρ, φ) = ∫ Π_t C(n_t, d_t) p(z_t)^{d_t} (1 − p(z_t))^{n_t − d_t}  dP(z_1, …, z_T),

which no longer splits into a sum over periods. That is why the period-reweighting machinery (W × L)
does not apply, and why its intervals are resampled by parametric bootstrap through the DGP (D-043).

## 2. The forward filter

Write f_t(z) = C(n_t, d_t) p(z)^{d_t} (1 − p(z))^{n_t − d_t} and let π_t be the predictive density of
Z_t given d_1, …, d_{t−1} (π_1 = N(0, 1)). Then

    L = Π_t ∫ f_t(z) π_t(z) dz,

and the filtered density after period t is f_t π_t normalised; π_{t+1}(z) = ∫ filtered_t(u) K(u → z) du,
with K the transition density N(φu, 1 − φ²). The engine evaluates this recursion over a set of
quadrature nodes for each period, in log space throughout (log-sum-exp; no exponential of an
unnormalised sum).

**Each period's nodes are the binomial likelihood's own.** Period t is integrated by the parity rule:
adaptive Gauss–Hermite with 128 nodes, or the composite Gauss–Legendre rule (16 panels of 16 points)
when d_t is 0 or n_t (D-118). Those rules integrate against a normal density. The filter places them
for a Gaussian proxy q_t = N(μ_t, σ_t²) of where period t's mass lies: with z = μ_t + σ_t y,

    ∫ f_t(z) q_t(z) dz = ∫ f_t(μ_t + σ_t y) φ(y) dy,

and f_t(μ + σy) is again a Vasicek binomial term, with Φ⁻¹(PD) replaced by Φ⁻¹(PD) − √ρ μ and √ρ by
√ρ σ. So the existing rule and its mode-and-scale hint integrate it unchanged.

**The change of measure is explicit.** The rule integrates f_t q_t; the filter needs f_t π_t. At each
node z the integrand is multiplied by the exact ratio π_t(z)/q_t(z), computed over the previous period's
nodes z_j with normalised log-weights ℓ̂_j:

    log π_t(z) − log q_t(z) = LSE_j(ℓ̂_j + log K(z_j → z)) − LSE_j(ℓ̂_j) − log q_t(z).

In the first period the ratio is φ(z)/q_1(z). Any proxy gives the same integral up to quadrature error;
the proxy decides only where the nodes go.

## 3. Where the nodes go: forward, backward, forward

The placement is the one design choice that matters for accuracy, and it is fixed: exactly three
passes, never iterated (D-195).
1. **Forward, predictive placement:** q_t is the predictive moment-matched from the filtered weights,
   μ_t = φ m_{t−1}, σ_t² = φ² s²_{t−1} + 1 − φ², with m and s² the filtered mean and variance of Z_{t−1}.
2. **Backward:** the smoothed mean m̂_t and variance ŝ²_t of every Z_t over pass 1's nodes, by the
   standard backward recursion.
3. **Forward, neighbour placement:** q_t is the law of Z_t given its neighbours at their smoothed
   moments,

       μ_t = φ (m̂_{t−1} + m̂_{t+1}) / (1 + φ²),
       σ_t² = (1 − φ²)/(1 + φ²) + φ² (ŝ²_{t−1} + ŝ²_{t+1}) / (1 + φ²)²,

   with one neighbour at either end (μ = φ m̂, σ² = 1 − φ² + φ² ŝ²). Pass 3 gives the likelihood.

**Why not simpler placements.** Two simpler designs were tried and measured against an independent
dense filter (section 6) on adjacent periods with very different modes: a zero-default year next to a
crisis year, n up to 10,000, φ up to 0.95.
- **Nodes placed for the prior,** N(0, 1) in every period: off by up to 0.9 nats. The 128- and
  256-node versions agreed with each other to 4e-11 throughout, so **a doubled-node check cannot
  detect placement error**: both rules share the misplaced centre.
- **The predictive placement alone** (pass 1): off by up to 1e-4 nats when a later crisis selects an
  earlier period's tail, which a forward filter cannot see.

The neighbour placement is within 1.1e-13 (relative) of the dense filter on every such case. The
difference between passes 1 and 3 is therefore the per-run placement check (section 5).

**φ = 0 is the binomial likelihood bit for bit.** At φ = 0 every proxy is exactly N(0, 1) and every
log-ratio is exactly 0, so each period's term is the binomial-mixture MLE's and the periods are summed
the same way. There is no special case in the code: it is the recursion itself evaluating to the
independent model.

**Smooth in the parameters.** The placement moves with (PD, ρ, φ); a clamp, floor or branch in it would
put kinks into the likelihood. The variances are two-pass sums of squares (non-negative without a floor),
and the likelihood's second differences against the dense filter along lines in φ and ρ are at rounding
level (`TOL_AR1_SMOOTH_D2_REL`).

**Cost.** The period-to-period transition (128 × 128 node pairs, 256 for zero-default periods) dominates. Its
log-sum-exps run one lane per node (scalar, SSE2 or AVX2, chosen at run time; all three give the same numbers) with
the in-house exp the DGP also uses, which is identical on every platform (P-18, D-199). With AVX2, one evaluation
(three passes) takes about 2.4–5.2 ms at T = 20 and 12.5–28 ms at T = 100 on one thread, and its score about
four-fifths as much again.

## 4. The estimate

**The score** is analytic, by Fisher's identity: the expectation, under the smoothed law of the factor
given the data, of the gradient of the complete-data log-density. For PD and ρ that is the sum over
periods of the binomial term's parameter derivatives at each node, weighted by the smoothed node
weights. For φ it is the sum over adjacent pairs of ∂ log K/∂φ = a u/v − φ a²/v² + φ/v (a = z − φu,
v = 1 − φ²) under the smoothed pairwise weights. It matches Richardson-extrapolated differences to 7e-10
relative (`TOL_AR1_SCORE_REL`).

**The search** is fixed and deterministic, over the recovery box (PD in [1e-4, 0.2], ρ in [1e-3, 0.5])
and |φ| ≤ 0.95 (D-195: φ = 0 interior, so neither the intervals nor a test of φ = 0 meet a boundary
mixture), in coordinates (logit PD, logit ρ, φ):
1. **The scan:** φ ∈ {−0.6, −0.3, 0, 0.3, 0.6, 0.9}, maximising over (PD, ρ) at each, starting at the
   binomial fit at φ = 0 and warm-started outward.
2. **The polish:** a 3-D maximisation from each local maximum along the scan; the best is the estimate.

Each maximisation is a quasi-Newton search on the analytic score: the Hessian is central differences of
the score at the start, with BFGS updates in between; the step is damped until the Hessian is
negative definite, then a backtracking line search on the log-likelihood, projected onto the box. A
coordinate at a bound whose gradient or step points out of the box is held there. The differenced
Hessian is recomputed when the set of coordinates held at a bound changes, when the line search keeps
shortening the step, and after any step that needed damping: far from the maximum the differenced
Hessian can be indefinite, and BFGS keeps a definite matrix definite but cannot make an indefinite one
so. Carried forward, it would leave every later step damped towards gradient ascent, converging
linearly (D-197).

**Standard errors** come from the observed information (central differences of the score at the
estimate), by the delta method to the natural scale.

**Profile intervals** (95%) follow the binomial MLE's rule (D-129, D-130) in three dimensions: the
profile of a parameter maximises over the other two; each end solves profile = maximum − 1.9207 by
Newton on the profile, whose derivative is the score's component at the inner maximum (the envelope
theorem); an end that would lie beyond a bound is that bound, flagged truncated.

## 5. Checks and flags

- **The doubled rule:** the likelihood at the estimate with twice the nodes and panels, as the binomial
  MLE's per-run check (D-092, D-120).
- **The placement check:** |pass 3 − pass 1| at the estimate. It is the check a doubled rule cannot be.
- Both are flagged above max(10⁻¹⁰ T, 256 ε Σ_t |l_t|).
- **At a bound:** a coordinate of the estimate on the box's boundary (φ at ±0.95 included) is flagged,
  never clamped silently.
- **Diagnostics:** the number of local maxima along the scan and the log-likelihood gap to the second.

## 6. The independent reference

`vcalref::ar1_loglik` (ref/, standard library only, D-031) evaluates the same likelihood a different way:
one fixed uniform grid over [−10, 10] (step min(0.02, 0.4/√max n)), trapezoid weights, the transition
density built directly, and ref's own Φ and log C(n, d). Its grid depends only on the panel, never on the
parameters, so it is smooth by construction; it is slow by design. It is the oracle for the filter's
accuracy, smoothness and goldens.

## 7. Replicating it

`validation/scipy/ar1_binomial_mle.py` recomputes the log-likelihood of every row of
`tests/golden/ar1/likelihood.csv` (11 hand-made panels, including the adversarial ones, at 6
parameter points; written by `tests/studies/ar1_goldens.cpp` with the engine's and ref's values). It
uses `scipy.stats`, `scipy.special` and `scipy.integrate` only, in two formulations of the factor's
law, so that an error derived the same way in the engine and in one formulation (the transition
variance 1 − φ², say, or the prior-ratio term) shows up against the other:

- **The Markov formulation**, through the transition density N(φ z_{t−1}, 1 − φ²), as the engine:
  1. a forward filter on one uniform grid of composite-Simpson nodes over [−12, 12] (step
     min(0.02, 0.5/√max n)), with the transition matrix from `norm.logpdf` and the recursion in log
     space (`logsumexp`); every row;
  2. for T ≤ 3 and n ≤ 1000, nested one-dimensional `quad`: given Z₂, the other periods' factors are
     independent, each N(φ Z₂, 1 − φ²), so the likelihood is a 1-D integral over z₂ of products of 1-D
     integrals.
- **The joint formulation**, with no transition density at all:
  3. for T ≤ 3 and every n, adaptive cubature (`scipy.integrate.cubature`, Gauss–Kronrod 21) of the
     explicit joint density: (Z_1, …, Z_T) multivariate normal with covariance φ^|i−j| (a Toeplitz
     matrix, `multivariate_normal`) times the binomial terms. The region is the joint mode ± 20 Laplace
     scales, both found from the joint density itself (BFGS, then central differences): over the whole
     cube the rule's first evaluation can miss the narrow ridge of an n = 10⁴, φ = 0.95 panel and report
     convergence on nothing. The joint log-density is concave, so the mode is unique.

log C(n, d) is the log of the exact integer `math.comb(n, d)`; a difference of `gammaln` values loses
about 2·10⁻¹¹ to cancellation at n = 10⁴, which would dominate the comparison. The range is ±12
because the hardest row (PD 0.001, ρ 0.05, 90 defaults in 1,000) puts the factor near z = −8; at ±9
all three methods lost the same tail mass, 3·10⁻¹² relative.

Agreement with the engine, relative to max(1, |l|): 2.9·10⁻¹⁴ (grid filter, 66 rows), 7.5·10⁻¹⁵
(nested, 36 rows), 1.8·10⁻¹⁴ (joint cubature, 48 rows). Tolerances `TOL_SCIPY_AR1_LOGLIK_REL`
(10⁻¹²) and `TOL_SCIPY_AR1_JOINT_REL` (10⁻¹⁰). It runs in the CI job "validation (scipy)", about six
minutes.

### Recovery and the search's checks

`unit_ar1_search` (`tests/unit/ar1_search_test.cpp`) fits panels from the AR(1) DGP with known φ, on
their own seed, and checks:
- the estimate is a maximum of ref's independent likelihood: ref falls a step of 10⁻³ away in every
  coordinate and diagonal direction;
- no fixed-φ maximum on a grid of step 0.05 over [−0.95, 0.95], each cold-started, beats it (slow);
- at φ = 0 the inner search is the binomial-mixture MLE (logit gaps at most 1.2·10⁻⁸);
- on T = 200, n = 10,000 panels at φ = −0.5, 0.3 and 0.8 every parameter is within 4 SEs of the
  truth, unflagged (observed at most 1.4 SEs; slow);
- φ at ±0.95 is held exactly with its score pointing out of the box and flagged; an all-zero panel
  stays finite and flagged;
- fits are deterministic bit for bit;
- every solved profile end is on the threshold (to 10⁻¹³) for the best of four cold-started inner
  maxima. One start is not enough: at an end the inner problem can have two maxima.

## 8. What the estimator does and does not do

S-35's M6 arm (registered in `studies/lrpd-serial-correlation/PREDICTION.md`) tests whether modelling the
serial correlation, on top of a binomial MLE that already handles the nonlinearity of the default rate,
restores interval coverage and moves the median PD estimate. Its expected Fisher information
(`information_m6.csv`) shows φ weakly identified at short histories: at T = 20 its information-based SE is
at least 0.18 even at n = 10,000, and far larger at small n·PD, where about one default per period says
little about the factor path. There φ̂ is often near or at a bound, and the profile in φ is nearly flat.
