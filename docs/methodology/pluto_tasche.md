# Pluto and Tasche's most-prudent upper bounds (native benchmark)

This note describes, in plain prose, the most-prudent upper bounds for the PD of a low-default portfolio
(`benchmark_models/pluto_tasche/bounds.hpp`; D-214, D-215, D-216), after Pluto and Tasche (2005); §6 relates them to
the paper section by section. They are a `native`
benchmark: they change no parity result and are not in the C ABI. Study S-6 (`studies/pluto-tasche/PREDICTION.md`)
measures their coverage and conservatism.

## 1. The question they answer

A portfolio has few or no defaults, so the PD cannot be estimated with any precision. The most-prudent bound asks
instead: what is the largest PD that is still consistent with what was seen? At confidence γ, it is the PD p at
which observing D or fewer defaults is just a (1 − γ) event:

    P(S ≤ D | p) = 1 − γ,

with S the number of defaults the model expects. Any larger PD would make the observed total less likely than
1 − γ. The bound is "most prudent" because it takes the largest such p.

Here a panel has T periods with n obligors in each, d_t defaults in period t, and D = Σ_t d_t in total. The
bounds depend on the panel only through D.

## 2. The two bounds

**(a) Independent.** Defaults are independent across obligors and periods, so S ~ Bin(N, p) with N = n T. The
bound is the one-sided Clopper–Pearson upper limit. It equals the γ quantile of a Beta(D + 1, N − D)
distribution. At D = 0 it has the closed form 1 − (1 − γ)^(1/N).

**(b) One factor.** Defaults within a period share a systematic factor, as in the Vasicek model. Each period's
count S_t has the mixed binomial law

    P(S_t = d) = ∫ C(n, d) p(z)^d (1 − p(z))^(n − d) φ(z) dz,
    p(z) = Φ((Φ⁻¹(p) − √ρ_a z) / √(1 − ρ_a)),

and the periods' factors are independent, so S = Σ_t S_t is the T-fold convolution of that law. The correlation
ρ_a is assumed, not estimated: S-6 runs it at the generating value (the oracle) and at a fixed 0.12.

As ρ_a → 0 the mixed binomial tends to Bin(n, p), and (b) tends to (a).

## 3. How they are computed

- **The per-period law.** P(S_t = d) for d = 0..min(D, n) is exp of the binomial-mixture objective's per-period
  log-likelihood, `log_contrib`: log C(n, d) from `lbinom`, plus the parity quadrature of the integral, all in log
  space (`binomial_mixture_mle.md`). Counts above n have probability 0.
- **The convolution.** Only counts 0..D matter, so every convolution is truncated at D. The T-fold power is taken
  by repeated squaring (about log₂ T squarings and as many products). Each term of a convolution is a
  log-sum-exp over its pairs, shifted by its own largest pair, so nothing underflows however far the trial p is
  from the root.
- **The CDF.** log P(S ≤ D | p) is the log-sum-exp of the convolved log PMF over 0..D. For (a) it is the
  log-sum-exp of log C(N, k) + k log p + (N − k) log(1 − p) over k = 0..D, with log p and log(1 − p) formed from
  logit p directly.
- **The solve.** P(S ≤ D | p) falls as p rises, because p(z) rises with p at every z. So g(u) = log P(S ≤ D | p(u))
  − log(1 − γ), with u = logit p, has one root. It is found by bisection from the bracket u ∈ [−30, 15] (p from
  about 1e-13 to 1 − 3e-7), keeping g ≥ 0 at the lower end and g < 0 at the upper, until the bracket is at most
  1e-10 wide (`TOL_PT_SOLVE_LOGIT`). The bound is the bracket's midpoint. If g does not change sign over the
  starting bracket, the bound is that end of the bracket, flagged; S-6's tool stops if any bound is flagged.

## 4. Coverage, exactly

Because the bound rises with D, it covers the true PD exactly when D ≥ D_min, the smallest total at which the
assumed law gives P_a(S ≤ D | PD) ≥ 1 − γ. So its coverage is P_true(S ≥ D_min), a finite sum, with no simulation.
When the assumed law is the true one (the oracle), coverage is at least γ by construction. S-6's `reference.py`
computes these sums with scipy, and the study's scorer checks every simulated panel against them.

## 5. Checks

- The independent bound at D = 0 against its closed form, to the solve's tolerance.
- The one-factor bound tends to the independent one as ρ_a → 0: the gap falls from ρ_a = 1e-2 to 1e-8, and is
  within `TOL_PT_RHO_LIMIT_LOGIT` at 1e-8.
- Monotonicity: the bounds rise strictly with D and with γ, and independent < (b, 0.02) < (b, 0.12) < (b, 0.24) at
  every total and γ tested.
- Goldens (`tests/golden/pluto_tasche/bounds.csv`), reproduced on every platform, and replicated by
  `validation/scipy/pluto_tasche.py` with `beta.ppf` and an independent trapezoid quadrature, convolution and
  brentq.

## 6. Against the paper

Read and checked against the working paper, arXiv:cond-mat/0411699v3 (4 April 2005), on 2026-10-13 (D-216).

- **What is the paper's.** The independent bound is its equations (2.2), (3.2) and (A.1): the binomial tail as a beta
  distribution function, the exact Clopper–Pearson bound. The one-factor bound for one period is its (4.1) and (4.3a),
  under the Basel one-factor model with ρ = 12% in its examples. Proposition A.1 (and A.2, for (4.1)) gives the
  uniqueness this note's §3 relies on.
- **What is not.** The paper's most-prudent principle is about grades: grade k's bound pools grade k with every worse
  grade, so the bounds are monotone across grades. S-6 has one grade, so its bound is the pooled bound for the whole
  portfolio; the sequence across grades is study S-52.
- **The multi-period bound is S-6's own construction, not the paper's.** The paper's §6 follows one cohort of N
  obligors for T years (each can default once, no newcomers), with the yearly factors jointly normal at correlation
  ϑ^|s−t|, and integrates over the T factors by Monte Carlo. S-6's bound (b) sums T periods, each with its own n
  obligors and an independent factor, so it is the paper's one-period (4.3a) per period, convolved. With independent
  defaults the paper pools the years exactly as S-6's (a) does (its §6, second paragraph). The cohort model is study
  S-50.
- **The published examples are reproduced** (`unit_pluto_tasche: pt_reproduces_the_published_tables`): Tables 1–8 and
  the two values in the text, within 0.01 percentage points (the printed rounding and the paper's numerics; Table 8
  agrees with its values rounded up).
- **Errata in the paper** (none changes its method):
  - Table 4 at γ = 75% prints 0.65% (and Table 11 repeats it). The bound is 0.638%: at 0.65%, P(X ≤ 3) = 0.237, below
    1 − γ = 0.25, so 0.65% does not satisfy the paper's own (3.2). scipy's `beta.ppf(0.75, 4, 797)` agrees with the
    engine.
  - Equation (4.3a) prints G(p, ρ, y)^i, with G the threshold (4.3b); the conditional PD Φ(G(p, ρ, y)) is meant, as in
    (4.1). Table 8 is reproduced only with Φ.
  - Table 8 gives grade A's 50% bound as 0.72% and Table 12 as 0.71%; the bound is 0.7106%.
  - The proof after Proposition A.2 is headed "Proof of Proposition A.1".

## Reference

Pluto, K. and Tasche, D. (2005). Estimating probabilities of default for low default portfolios. Working paper,
arXiv:cond-mat/0411699v3 (4 April 2005; v1 28 November 2004), https://arxiv.org/abs/cond-mat/0411699. Read in full.
Published as "Thinking positively", *Risk*, August 2005, pp. 72–78; the title is confirmed on Risk.net, and the month
and pages only by secondary listings (the article is paywalled). Listed in `docs/references.md`.
