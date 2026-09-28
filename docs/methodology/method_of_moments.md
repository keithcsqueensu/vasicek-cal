# Method of moments (M3)

This note describes, in plain prose, the method-of-moments estimator of PD and ρ in the one-factor
Vasicek model, in its joint-default-probability form (`engine/moments.hpp`). Decision numbers
(D-nnn) point to `DECISIONS.md`. The accompanying script `validation/scipy/method_of_moments.py`
reproduces it with ordinary scipy tools (section 5).

## 1. The two moments

In period t, n_t obligors each default with probability p(Z_t), Z_t ~ N(0, 1), as in the binomial
model (`binomial_mixture_mle.md`). Two moments identify the parameters:

- the PD, the probability that one obligor defaults, E[p(Z)];
- the **joint default probability**, the probability that two given obligors of the same period
  both default,

      PD₂(PD, ρ) = E[p(Z)²] = Φ₂(c, c; ρ),   c = Φ⁻¹(PD),

  the bivariate normal distribution function at (c, c) with correlation ρ.

## 2. The estimator

From a panel of counts, both moments are estimated **exactly in finite n** by pooled ratios:

    PD̂ = Σ_t d_t / Σ_t n_t,        PD̂₂ = Σ_t d_t (d_t − 1) / Σ_t n_t (n_t − 1).

Given Z_t = z, d(d − 1)/(n(n − 1)) is unbiased for p(z)², so PD̂₂ is unbiased for PD₂ whatever n is:
no large-portfolio approximation is made. ρ̂ is then the ρ at which PD₂(PD̂, ρ) = PD̂₂. PD₂ increases
with ρ, so the root is unique; it is found by Brent's method in logit(ρ), to 10⁻¹² there, over the
estimation box's ρ range.

**Computing PD₂.** E[p(Z)²] is exactly the binomial-mixture integral of a period with n = 2
obligors and d = 2 defaults, so it is computed by the same log integrand and the same parity
integrator as the binomial MLE, in log space. Against mpmath (tanh-sinh quadrature at 50 digits,
an independent method) its logarithm agrees within 7.1·10⁻¹⁵ over the box's PD and ρ.

**Edge cases,** flagged and never extrapolated:

- no defaults at all: refused, no estimate;
- PD̂₂ at or below what the box's smallest ρ gives, including PD̂₂ ≤ PD̂² (no dependence in the
  data): ρ̂ is the box's floor, flagged;
- PD̂₂ at or above what its largest ρ gives: ρ̂ is the cap, flagged;
- PD̂ outside the PD box: flagged, and ρ̂ is still solved at PD̂.

**Zero-default periods need no treatment:** they enter the sums as zeros, which is what they are.

**Resampling (D-042).** The statistics are sums over periods, so a resampling weight row reweights
them: weight w_t multiplies period t's terms, and a weight of 2 gives exactly the panel with that
period twice. The estimator therefore uses the same W matrices as the likelihood-based estimators
(bootstrap, jackknife, walk-forward), without a surface.

**Rate series (D-046).** Without counts, the moment form uses PD̂ = the mean rate and PD̂₂ = the
mean squared rate, the same inversion. It has no binomial correction, so on count data it reads
binomial noise as dependence, as the rate MLE does (`vasicek_rate_mle.md`).

## 3. What is known about its statistical behaviour

MoM matches the first two moments of the rate on the natural scale, where the likelihood uses the
whole distribution. Rates are right-skewed, more so at low PD and high ρ, so MoM is less efficient
than the MLE for ρ, and its ρ̂ is biased low. How much, across the recovery matrix, is study S-8's
question (registered, D-161; its reference in the large-n limit is `studies/mle-vs-mom/model_reference.py`).
This step deliberately runs **no** comparison on the recovery panels: that is S-8's run, scored by
its own pre-registered script.

## 4. Validation in this step

- **PD₂ against mpmath** (`tests/golden/moments/joint_default.csv`, `tools/gen_moments_goldens.py`
  in the CI's generated-files check): within 7.1·10⁻¹⁵ in log PD₂ (`TOL_MOM_JOINT_DEFAULT_LOG_ABS`).
- **Inversion:** ρ from the golden PD₂ is recovered within 1.3·10⁻¹³ in logit (`TOL_MOM_INVERSION_U`).
- **Edge cases and weights,** as listed above.
- **Consistency:** on one panel of 20,000 periods of its own validation seed (not the recovery
  seed), PD̂ and ρ̂ lie within 0.97 and 1.74 of their large-T standard deviations from the truth
  (`TOL_MOM_CONSISTENCY_Z` = 4).
- **Reference results** (`tests/golden/moments/reference.csv`, written by `mom_reference`): 108
  panels, count panels and rate series for every (PD, ρ, T) of the recovery matrix at n = 1,000,
  under the validation seed M3MOMVAL. Recomputed by `unit_moments` on every platform
  (`TOL_MOM_REFERENCE_REL`).

## 5. Replicating with scipy

`validation/scipy/method_of_moments.py` recomputes every reference row with ordinary tools. It
evaluates Φ₂ at equal arguments exactly, by Owen's T function,

    Φ₂(h, h; ρ) = Φ(h) − 2·T(h, √((1 − ρ)/(1 + ρ))),

independent of the engine's quadrature, and solves for ρ̂ by `brentq` in logit(ρ). It agrees with
the engine within 6.7·10⁻¹⁶ in the moments and 1.6·10⁻¹¹ in logit(ρ̂) (`TOL_SCIPY_MOM_*`), and on
every flag. It runs in the CI job "validation (scipy)".

## 6. What is not here yet

- **Intervals:** MoM has no likelihood, so no profile interval. Bootstrap intervals come from the W
  machinery above; their coverage is not validated here.
- **The C ABI** exposes it from 0.3 as `vcal_calibrate_moments` (D-169).
