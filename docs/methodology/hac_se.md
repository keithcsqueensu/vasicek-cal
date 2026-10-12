# The HAC (Newey–West) sandwich standard error (native option)

This note describes, in plain prose, the HAC sandwich standard error of the binomial-mixture MLE
(`engine/hac.hpp`, D-209). It is a `native` option and changes no parity result. Study S-48
(`studies/hac-se/PREDICTION.md`) tests whether it is worth recommending. It is not in the C ABI. If S-48
recommends it, a later decision adds it there.

## 1. What it corrects

calibrate's standard errors come from the observed information: V = (−H)⁻¹, with H the Hessian of the
panel log-likelihood at the estimate (`binomial_mixture_mle.md`). That is the variance of the estimate
when the periods are independent. The estimate's error is, to first order, (−H)⁻¹ times the sum of the
per-period scores s_t = ∂l_t/∂θ. When periods are dependent, the scores are too, and the variance of
their sum exceeds the sum of their variances.

S-35 measured this under an AR(1) systematic factor. The iid intervals for PD are too narrow by
about √((1 + φ)/(1 − φ)): a factor of 2 at φ = 0.6, at every history length (D-192).

The sandwich replaces the information in the middle by the long-run variance S of the scores:

    V = H⁻¹ S H⁻¹.

With independent periods S estimates −H, and V reduces to the iid variance.

## 2. The estimator

Everything is computed at the reported estimate, in the axes' scaled coordinates u, the same
coordinates as the Wald interval. For the recovery grid these are logit PD and logit ρ.

1. **The scores.** For each period t, s_t = (∂l_t/∂u_PD, ∂l_t/∂u_ρ). These are the natural-scale
   derivatives the objective already provides (`log_contrib_derivs`, D-170), multiplied by dv/du per axis.
   Periods with the same observation (n_t, d_t) share one evaluation.
2. **Centring.** c_t = s_t − s̄, with s̄ the mean over the T periods. At an exact maximum s̄ = 0. The
   reported estimate is a refined grid point, whose score is small but not exactly zero, and centring
   keeps that from entering S.
3. **The long-run variance**, with Bartlett weights w_l = 1 − l/(L + 1):

       S = Σ_t c_t c_tᵀ + Σ_{l=1}^{L} w_l Σ_{t=l}^{T−1} (c_t c_{t−l}ᵀ + c_{t−l} c_tᵀ).

   - The Bartlett weights keep S positive semi-definite.
   - S is a plain sum in period order, with no small-sample factor.
   - Lags of T or more contribute nothing.
   - L = 0 gives the outer-product sandwich.
4. **The lag.** The fixed rule is L = ⌊4 (T/100)^(2/9)⌋, giving 2, 3 and 4 at T = 20, 40 and 100. A
   caller may pass another lag; S-48 reports L = 1 and twice the rule.
5. **The sandwich.** V = H⁻¹ S H⁻¹, using the same H as calibrate's own SE (the analytic Hessian of
   D-185). The HAC SE in u is √V_aa, and on the natural scale it is (dv/du)·√V_aa: the delta method, as
   for the iid SE.
6. **The Wald interval** is u_a ± z √V_aa, mapped back to the natural scale, with z = 1.959963984540054.
   It is symmetric in logit units, like the iid Wald interval.

**When there is none.** The HAC SE is not computed (flag `kIntervalNotComputed`; every value NaN) when:

- the fit carries `kFlagGridEdge`, `kFlagFlatSurface` or `kFlagNearBound`, the fits that get no iid
  Wald interval either. Both Wald intervals are therefore computed on the same fits, and S-48 counts the
  rest;
- H is not negative definite;
- V's diagonal is not positive and finite.

## 3. Why these choices

- **The scores, not the default rates.** The sandwich is the variance of the estimator itself. It needs
  no model of the dependence, only the assumption that the dependence dies out within a few lags.
- **The fixed lag rule.** A data-driven lag (Andrews, Newey–West 1994) would be one more estimated
  quantity at T = 20. The rule is the common default (Newey and West 1994's T^(2/9) rate with constant 4).
  Fixing it in advance keeps S-48 a test of one method, not of a tuning search.
- **Logit units.** These are the Wald interval's own coordinates in this engine (D-095). An interval
  symmetric there respects the parameter bounds.
- **The known weakness, not corrected here.** At small T the Bartlett estimator is biased down, from the
  kernel's down-weighting and from estimating the mean. Its own variability is large, with effective
  degrees of freedom of about 3T/(2(L + 1)). Under independence it therefore covers less than the iid
  Wald interval. S-48 measures that price beside the gain.

## 4. Tests and replication

- **`tests/unit/hac_test.cpp`:**
  - the lag rule;
  - S by hand on four score vectors (exact);
  - the SE against an independent assembly on the M1.4 panel at L = 0, 1, 2 and 4: scores by central
    differences of each period's log-likelihood, H by differences of their sum, and S by the textbook
    double sum Σ_{t,r} k(|t − r|) c_t c_rᵀ. Agreement is about 10⁻⁷, the differences' own accuracy;
  - invariance to reversing the periods, and dependence on any other reordering once L > 0;
  - no SE on a flagged fit;
  - regression values on the six M1.7 panels (`tests/golden/hac/validation.csv`).
- **`validation/scipy/hac_se.py`**, a plain scipy script (CI job "validation (scipy)"):
  - The scores are posterior means of ∂ log binom.pmf/∂u, integrated with `quad` around the
    integrand's mode, sharing nothing with the engine's moment accumulators. H is central differences
    of their sum.
  - It checks the lag rule; scipy's (−H)⁻¹ against the engine's iid SEs; and the HAC SEs and
    correlations at L = 1, 2 and 4 on the six M1.7 panels. It also checks that the fits with no HAC SE
    are exactly the flagged ones, two near-bound fits whose distance to the box is under 2 SEs by
    scipy's own count.
  - Agreement: 5.7·10⁻¹⁰ (`TOL_SCIPY_HAC_SE_REL`, 2·10⁻⁹).

**In a short script.** Fit the binomial MLE. At the estimate, take each period's score in logit units,
centre the scores, and form S with Bartlett weights as above. Take the Hessian H of the panel
log-likelihood in the same units. Then SE = √diag(H⁻¹ S H⁻¹). `validation/scipy/hac_se.py` does exactly
this in about 100 lines.

## References

- Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and
  autocorrelation consistent covariance matrix. *Econometrica* 55(3), 703–708.
- Newey, W. K. and West, K. D. (1994). Automatic lag selection in covariance matrix estimation. *Review
  of Economic Studies* 61(4), 631–653.
