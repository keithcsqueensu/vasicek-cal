# S-1 Z-sign invariance, calibration half: predictions made before the run

Pre-registration for study S-1 (`studies/README.md`, D-148). This file is committed **before**
any S-1 quantity has been computed: no surface, estimate or interval has been evaluated under the
mirrored sign convention. The results will be compared with it check by check, and misses will be
reported, not explained away.

**This file is not edited after the run.** A change of method goes in a new, dated section
appended at the end, committed before the changed run. Nothing above it changes.

**Design (owner, D-300).** S-1 tests an exact invariance, not a statistical property, so it is a
pass/fail property check on a small fixed set of panels, compared at the surface level first and
then downstream. It uses no recovery replicates and gives no coverage verdicts. The recovery-matrix
version described in the original entry is kept in reserve: it is run only if this check finds a
difference that needs statistical characterisation.

**Evidence used in writing it:** the code of the parity rule (`core/quadrature/`,
`core/model/binomial_mixture.hpp`), and one property of its tables, inspected while writing this
file: the stored Gauss–Hermite nodes and log-weights (N = 128 and 256) and the Gauss–Legendre
nodes and weights (M = 16) are exactly antisymmetric and symmetric as doubles,
node[i] = −node[N−1−i] and weight[i] = weight[N−1−i], bit for bit. Nothing else was run.

## The invariance

The engine writes p(z) = Φ((c − √ρ·z)/√(1 − ρ)), c = Φ⁻¹(PD), so a higher Z means better
conditions (R-2). The mirrored convention writes p⁻(z) = Φ((c + √ρ·z)/√(1 − ρ)). Since p⁻(z) = p(−z)
and φ is symmetric, every per-period likelihood ∫ p(z)^d (1 − p(z))^(n−d) φ(z) dz is the same under
both. So is everything computed from the per-period surfaces: estimates, SEs, profile intervals,
bootstrap intervals and flags.

## What is run

- **The mirrored objective (test-only, never parity):** the binomial-mixture objective with the
  threshold x⁻(z) = (c + √ρ·z)/√(1 − ρ), and its own quadrature hint derived for that convention,
  not obtained by negating the parity hint:
  - the interior mode: Newton with bisection on h⁻′(z) = +β(d·λ(x⁻) − (n − d)·λ(−x⁻)) − z, started
    from z_L = (√(1 − ρ)·x_L − c)/√ρ shrunk toward 0 as in parity;
  - for d ∈ {0, n}, the survival factor's half-point z_c = (√(1 − ρ)·x_c − c)/√ρ and its width.

  It is integrated with the parity rule itself (adaptive GH, N = 128, for 0 < d < n; composite
  Gauss–Legendre, 16 × 16, sinh map, for d ∈ {0, n}), unchanged, with the doubled rule as the check.
- **Both conventions** are fitted with the same engine calls as a user's fit and the recovery
  harness: `calibrate` on the recovery grid (61 × 41, logit, PD ∈ [1e-4, 0.2], ρ ∈ [1e-3, 0.5]),
  `profile_intervals`, and the iid bootstrap, B = 999, percentile intervals, with the same W for
  both conventions.
- **A positive control:** the mirrored integrand with a *wrong* hint, the parity hint unchanged
  (centred on the unmirrored mode and half-point). The comparison must detect it; if it does not,
  the check is not trusted.
- **Where it lives:** a test-only objective and a study executable under `tests/`, plus a fast
  CTest entry that re-runs the check, so the invariance stays guarded after the study.

## The panels (fixed here, 34)

The DGP panels use seed `0x53315A5349474E53` ("S1ZSIGNS"), scenario = the panel number, replicate 0,
with n the same in every period unless stated. The constructed panels are fixed count vectors.

| # | Source | PD | ρ | T | n | The hard case it holds |
|---|---|---|---|---|---|---|
| 1–9 | recovery, replicate 0 | | | | | the study subset's scenarios 29, 37, 72, 4, 68, 49, 7, 43, 51: a reference to the recovery panels |
| 10 | DGP | 0.1% | 0.02 | 40 | 100 | mostly zero-default periods |
| 11 | DGP | 0.05% | 0.30 | 20 | 50 | nearly all periods zero, PD near its lower bound |
| 12 | DGP | 15% | 0.40 | 20 | 5 | all-default periods (d = n), high ρ |
| 13 | DGP | 18% | 0.45 | 20 | 2 | many all-default periods, ρ near its cap |
| 14 | DGP | 1% | 0.45 | 20 | 1,000 | high ρ, near the cap |
| 15 | DGP | 5% | 0.49 | 40 | 10,000 | ρ at the cap, informative |
| 16 | DGP | 1% | 0.0015 | 40 | 10,000 | ρ near its lower bound |
| 17 | DGP | 0.02% | 0.10 | 40 | 10,000 | PD near its lower bound, many zeros |
| 18 | DGP | 1% | 0.12 | 20 | 1,000,000 | large n |
| 19 | DGP | 0.1% | 0.05 | 20 | 1,000,000 | large n, low PD |
| 20 | DGP | 5% | 0.24 | 20 | 100,000 | large n, high ρ |
| 21 | DGP | 2% | 0.15 | 20 | n_t = round(10^(1 + 4t/19)), t = 0 … 19 (10 to 100,000) | varying n |
| 22 | DGP | 19% | 0.20 | 20 | 1,000 | PD near its upper bound |
| 23 | constructed | | | 20 | 1,000 | every d = 0: estimate on the grid edge |
| 24 | constructed | | | 10 | 3 | every d = n: estimate on the grid edge |
| 25 | constructed | | | 20 | 1,000 | one default in period 0, the rest zero |
| 26 | constructed | | | 20 | 1 | d alternating 0, 1 (n = 1: ρ not identified, flat in ρ) |
| 27 | constructed | | | 20 | 1,000 | d = 10 every period: no dispersion, ρ at its lower bound |
| 28 | constructed | | | 20 | 1,000 | d alternating 0, 500: extreme dispersion, ρ at its cap |
| 29 | constructed | | | 20 | 4 | d cycling 0, 1, 2, 3, 4: zero, interior and all-default periods together |
| 30 | constructed | | | 20 | 1,000,000 | d cycling 0, 0, 1, 3, 0: large n with zeros |
| 31 | constructed | | | 20 | 10,000 | d = 2,000 every period (PD 20%, the upper bound) |
| 32 | constructed | | | 20 | 100 | d cycling 0, 100: zero and all-default only |
| 33 | constructed | | | 20 | 10,000 | d cycling 1, 1, 1, 1, 30: one bad period in five |
| 34 | constructed | | | 40 | 1,000 | d = 3 in period 39, zero elsewhere: a single late default |

The run prints a census (zero-default and all-default periods, estimates on or near a bound, flat
or flagged fits) as an acceptance check that the hard cases are present. The constructed panels
guarantee them whatever the DGP draws.

## Mechanism, and the bound it gives

In IEEE arithmetic the mirror is nearly exact:

- **The integrand.** With z⁻ = −z, √ρ·(−z) = −(√ρ·z) and c + (−a) = c − a exactly, so
  x⁻(−z) = x(z) and g⁻(−z) = g(z) bit for bit.
- **The hint.** Applied to h⁻(z) = h(−z), the same Newton-and-bisection algorithm meets every
  quantity negated exactly (the start, the slopes, the brackets, the midpoints), and z_c⁻ is the exact
  negation of z_c. The hint's mode and centre should be **exact negations**, and its scale and width
  **bitwise equal**, at every surface cell.
- **Gauss–Hermite (0 < d < n).** The nodes are exactly antisymmetric, so the mirrored rule visits
  the negated abscissae, and each term log wᵢ + (uᵢ − zᵢ)(uᵢ + zᵢ)/2 + g(zᵢ) is bitwise one of the
  parity terms. **The only difference is the summation order**, which is reversed. For 128 positive
  terms folded by the online log-sum-exp, the bound is about 2N·ε = 256ε ≈ 5.7e-14 in log I, plus the
  last roundings of log σ + log S and of log C(n, d) + log I, each ε times the terms' size.
- **Composite Gauss–Legendre (d ∈ {0, n}).** The domain's end points are found by the same bisection
  on mirrored arguments, so they are exact negations, and u_lo⁻ = −u_hi exactly (asinh is odd). But
  the panel midpoints are computed from the other end, u_lo + (k + ½)·du, so the nodes differ from
  the mirrored parity nodes by a few ulp of |u| (up to about 20). Each term moves by
  |dh/du|·δu, and h falls by at most 50 across the domain (kCompositeDrop), so the error is
  about 1e-13 and is bounded here by 1e-12, again plus the summation-order and final roundings.

Downstream quantities are functions of these surfaces, so their differences are rounding
perturbations, of the same kind as a cross-platform replay (D-107). The bounds below reuse the
registered cross-platform tolerances where one exists.

## Predictions: acceptance checks (all must hold)

ε = 2.22e-16. "Scale" is the size of a period's terms, |log C(n, d)| + |log I| (the D-120 rounding
scale).

| # | Quantity | Predicted: bound in every case | Expected size |
|---|---|---|---|
| C1 | The hint: mode and centre negated, scale and width equal | bitwise, every cell | bitwise |
| C2 | Per-period surface cells ℓ_t(PD, ρ), 0 < d < n: \|Δ\| | ≤ 256ε + 4ε·scale | a few ulp of the scale |
| C3 | Per-period surface cells, d ∈ {0, n}: \|Δ\| | ≤ 1e-12 + 4ε·scale | ≤ 1e-13 |
| C4 | The panel surface Σ_t ℓ_t at every grid point | ≤ the sum of the per-period bounds | ≤ 1e-12 |
| C5 | Grid argmax k*, and the bootstrap replicates' grid argmax | identical in every fit and every replicate | identical |
| C6 | Flags: the fit's flags, the quadrature-check count, the profile flags (truncation, not computed), the bootstrap edge and excluded counts | identical in every fit | identical |
| C7 | PD̂, ρ̂ and the log-likelihood at the estimate, relative | ≤ 1e-12 (`TOL_RECOVERY_REPLAY_REL`) | ≤ 1e-13 |
| C8 | SEs and their correlation, relative | ≤ 5e-10 (`TOL_RECOVERY_REPLAY_SE_REL`) | ≤ 1e-11 |
| C9 | Profile end points, in logit units; the polished maximum, relative | ≤ 2e-9 (the root solver's tolerance, `TOL_XREF_PROFILE_ENDPOINT_U`); ≤ 1e-12 | ≤ 1e-11; ≤ 1e-13 |
| C10 | Bootstrap percentile end points, relative | ≤ 2e-11 (`TOL_BOOTSTRAP_CROSS_PLATFORM_REL`) | ≤ 1e-12 |
| C11 | The quadrature check's max \|ℓ_t(rule) − ℓ_t(doubled)\| at the estimate, absolute difference | ≤ 1e-12 | ≤ 1e-13 |
| C12 | The positive control (wrong hint) | fails C2 or C3 by at least 1e-6 in some large-n panel (18, 19, 20, 30), and sets the quadrature-check flag in at least one fit | as stated |

Estimates on a grid edge and intervals truncated at the box are compared like any other: an edge
estimate should be the same grid value in both conventions, and a truncated end point the same
bound.

**Reported, not predicted:** the fraction of surface cells that agree bit for bit, and the largest
difference for each check, per panel.

## Verdict rules

- **PASS:** C1 to C11 hold on all 34 panels and the control fires (C12). The calibration half of
  the Z convention is then shown to be invariant within rounding, and the recovery-matrix version
  is not run.
- **FINDING:** any check C1 to C11 fails, however small the excess. It is reported with the
  panel, the cell or output, and a diagnosis (for example, a hint that is not an exact mirror, or a
  rule that is asymmetric beyond rounding). A finding that affects estimates or intervals by more
  than rounding is then characterised statistically on the study subset, under a prediction
  appended here before that run.
- **NOT TRUSTED:** the control does not fire (C12 fails). The comparison is repaired and re-run
  under an appended prediction before any verdict.

## Cost

Two conventions and the control on 34 panels, each a full fit with profile intervals and B = 999
bootstrap: about 100 fits at 1–2 s of single-thread CPU each, so **a few minutes on 4 cores**. The
original entry's estimate ("under a minute") was for the replay alone.

## How the comparison will be reported

Each check C1–C12 gets a row: predicted bound, largest observed value and where, held or not. The
result and the comparison go in S-1's entry in `studies/README.md`, with its monitoring
implication, and in a D-entry.
