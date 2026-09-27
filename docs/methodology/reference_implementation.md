# The independent reference implementation (`ref/`)

`ref/` is a second, deliberately different implementation of the single-factor
binomial-mixture model. It exists so that the engine (`core/`, `engine/`) is checked against
something that shares none of its code or algorithms (D-012, D-031). If both agree, an error
would have to be made twice, independently, in the same way.

It is host-only C++ using only the standard library. A layering check fails the build if any
file in `ref/` includes anything else from the repository (D-050). Speed is not a goal.

Its own accuracy is established against the mpmath golden tables, which are data rather than
code (see [tolerances.md](tolerances.md), section "ref/ against mpmath").

## Model

For a period with n obligors and d defaults, and parameters PD and ρ (both in (0, 1)):

- threshold c = Φ⁻¹(PD);
- for a systematic factor value z, the conditional default probability is
  p(z) = Φ(x(z)) with x(z) = (c − √ρ z) / √(1 − ρ);
- the period log-likelihood is log C(n, d) + log E[p(Z)ᵈ (1 − p(Z))ⁿ⁻ᵈ], Z ~ N(0, 1);
- the panel log-likelihood is the sum over periods.

## log Φ(x)

Three regions:

1. **x ≤ −20.** log Φ(x) = −x²/2 − log √(2π) + log R(−x), where R(t) = Q(t)/φ(t) is the Mills
   ratio. R is computed from the continued fraction R(t) = 1/(t + 1/(t + 2/(t + 3/(t + …)))) by
   the modified Lentz algorithm, stopping when a step changes the value by less than half an
   ulp.
2. **−20 < x ≤ 2.** log Φ(x) = log(erfc(−x/√2)/2), using the standard library's `erfc`.
3. **x > 2.** log Φ(x) = log1p(−Q(x)) with Q(x) = φ(x) R(x), R from the same continued fraction.
   Here φ(x) is computed from the exact square x² = hi + lo, with hi = x·x and
   lo = fma(x, x, −hi). That is an error-free product, so rounding x² costs nothing. Using `erfc`
   in this region would amplify the rounding of x/√2 by about 2x².

Measured against mpmath: at most 11 ulp wherever the result is a normal number.

## Φ⁻¹(p)

- p = 0 gives −∞, p = 1 gives +∞; p outside [0, 1] gives NaN.
- **¼ ≤ p ≤ ¾.** Newton's method on erf(x/√2)/2 = p − ½, starting from x = (p − ½)√(2π).
  p − ½ is exact in this range, and erf keeps relative accuracy near 0.
- **p < ¼.** Newton's method on log Φ(x) = log p (log Φ from above, derivative φ(x)/Φ(x)),
  safeguarded by the bracket [−40, 0]. Any step that would leave the bracket is replaced by
  bisection.
- **p > ¾.** −Φ⁻¹(1 − p). 1 − p is exact for p ≥ ½.

Measured against mpmath: at most 4 ulp over p ∈ [5e-324, 1).

## log C(n, k)

With k replaced by min(k, n − k): the compensated (Neumaier) sum over i = 1 … k of
log((n − k + i)/i). Every term is non-negative, so nothing cancels. The direct sum is limited
to min(k, n − k) ≤ 10⁸ (returns NaN beyond). Measured against mpmath: at most 1 ulp.

## The mixture integral

log E[pᵈ(1 − p)ⁿ⁻ᵈ] = log ∫ exp(h(z)) dz − log √(2π), with
h(z) = d log Φ(x(z)) + (n − d) log Φ(−x(z)) − z²/2. A zero count's term is omitted.

1. **Mode.** h is concave, since log Φ is. Evaluate h at −1, 0, 1. If h rises to the right
   (or left), step outward with doubling steps until it falls, which brackets the maximum.
   Then golden-section search to a relative width of 10⁻¹². Call the mode z_m and h_m = h(z_m).
2. **Widths.** On each side, find by bisection the distance w at which h has fallen by ½ from
   h_m: one standard deviation for a Gaussian. The two sides get their own widths, w_L and w_R,
   because zero-default integrands are strongly asymmetric: a Gaussian on one side, a steep
   sigmoid cliff on the other.
3. **Panels.** Breakpoints at z_m, z_m ± w·2ʲ (j = 0, 1, 2, …), continuing on each side until
   h has fallen by more than 60 (e⁻⁶⁰ ≈ 10⁻²⁶ of the peak).
4. **Adaptive quadrature on each panel.** Integrate f(z) = exp(h(z) − h_m) with 10- and 20-point
   Gauss–Legendre rules. The nodes and weights are computed at start-up by Newton's method on
   the Legendre recurrence, not typed in. A panel is accepted if the two rules agree within
   either:
   - its share of an absolute tolerance, 10⁻¹⁶ (w_L + w_R); or
   - the integrand's own relative noise, 64 ε (1 + abs(h_m)). f is only as accurate as h, whose
     absolute error is about ε·abs(h_m).

   Otherwise the panel is halved, recursively, to a depth of at most 30. Panel results are
   added with compensated summation.
5. **Result.** h_m + log(total) − log √(2π).

Measured against mpmath: at most 7.5 ε·max(1, abs(log I)) on every golden case, including
ρ up to 0.99 and zero-default periods of 10 to 200 obligors. This is the regime where the
engine's adaptive Gauss–Hermite rule is not at full precision (D-086, D-092), which is why ref
takes a panel-based approach.

## Estimation

Maximum likelihood over a box [PD_lo, PD_hi] × [ρ_lo, ρ_hi], both inside (0, 1), in logit
coordinates u = log(v/(1 − v)):

- the profile log-likelihood in ρ is the maximum over PD found by golden-section search;
- that profile is maximised over ρ by a second golden-section search;
- both searches stop at a width of 10⁻⁸ logit units.

An estimate within 10 search tolerances of a bound is flagged `on_boundary`.

**Standard errors.** Take the Hessian of the log-likelihood in logit coordinates by central
differences with step 10⁻³, including the cross term. If it is negative definite, the
covariance is its negated inverse, converted to natural scale by the delta method: dv/du =
v(1 − v) for each parameter. Otherwise the SEs are NaN.

On a 20-period panel of 1,000 obligors each, a fit takes about 8 s (optimised build, one
thread).
