# S-1 Z-sign invariance, calibration half: results

Run on 2026-09-27, after the prediction ([`PREDICTION.md`](PREDICTION.md)) was merged, with
`study_z_sign_invariance` (`tests/studies/z_sign_invariance.cpp`), cpu-release, GCC 13.3 on Linux,
4 threads: 28 s wall, 88 s CPU for the three arms (parity, mirrored, wrong-hint control) on the 34
panels.

**Verdict: FINDING, reviewed (D-301).** The calibration is invariant under z → −z within rounding
on 33 of the 34 panels, and at the surface level on all 34: every surface cell, every quadrature
hint and every flag agrees, and the control is caught. The one exception, panel 26 (n = 1 in
every period), is not an asymmetry. There ρ is not identified: the likelihood of a single obligor
is PD or 1 − PD whatever ρ is. So the panel surface is flat along ρ to rounding (a spread of 2.8e-14
across the 41 ρ values, against at least 0.025 on every other panel), and the grid argmax along ρ is
decided by rounding noise. Parity reports ρ̂ = 0.262, the mirrored convention 0.001. Both flag the
fit (the PD estimate is on the grid edge), and both give the same correct profile interval for ρ,
the whole box [0.001, 0.5], truncated at both ends.

## Prediction against result

| # | Check | Predicted bound | Expected size | Largest observed (outside panel 26) | Bound held | Expected size held |
|---|---|---|---|---|---|---|
| C1 | hint mirrored | bitwise | bitwise | bitwise, 912,865 cells | yes | yes |
| C2 | cells, 0 < d < n | 256ε + 4ε·scale | a few ulp of the scale | 5.7e-14 (0.24 of the bound); at most 1.9 ε × scale on any panel | yes | yes |
| C3 | cells, d ∈ {0, n} | 1e-12 + 4ε·scale | ≤ 1e-13 | 2.3e-13 (0.10 of the bound) | yes | **no**, by 2.3× |
| C4 | panel surface | sum of the per-period bounds | ≤ 1e-12 | 9.1e-13 (0.13) | yes | yes |
| C5 | grid argmax, headline and bootstrap | identical | identical | identical on 33 panels; **768 differ, all in panel 26** (1 headline, 767 of its 999 bootstrap replicates) | **no** (panel 26) | — |
| C6 | flags | identical | identical | identical, 34,102 comparisons, including panel 26 | yes | yes |
| C7 | PD̂, ρ̂, log-likelihood | 1e-12 relative | ≤ 1e-13 | 8.5e-14; **panel 26's ρ̂: 0.262 against 0.001** | **no** (panel 26) | yes, elsewhere |
| C8 | SEs, correlation | 5e-10 relative | ≤ 1e-11 | 2.7e-10 (0.53 of the bound; panel 20, n = 10⁵) | yes | **no**, by 27× |
| C9 | profile end points; l_max | 2e-9 logit; 1e-12 relative | ≤ 1e-11; ≤ 1e-13 | 2.7e-12; 2.3e-14 | yes | yes |
| C10 | bootstrap percentile ends | 2e-11 relative | ≤ 1e-12 | 2.4e-12; **panel 26's upper ρ end differs by 0.996 relative** | **no** (panel 26) | **no**, by 2.4×, elsewhere |
| C11 | quadrature-check difference | 1e-12 | ≤ 1e-13 | 5.7e-14 | yes | yes |
| C12 | wrong-hint control fires | ≥ 1e-6 excess; ≥ 1 flagged fit | as stated | 2.7e6 excess; the check flagged 17 of 34 fits | yes | yes |

**Misses, as the prediction requires, not explained away:**

- **Panel 26 fails C5, C7 and C10.** The mechanism above was not foreseen when the panel was chosen,
  even though its row in the prediction says "ρ not identified, flat in ρ". The prediction expected
  the flatness to show up as identical flags and intervals, which it did, but not that a rounding
  perturbation would move the point estimate across the whole axis.
- **Three sizes were underestimated,** each still inside its registered bound:
  - C3 by 2.3×: the composite rule's node offsets matter more than estimated at n = 10⁶ (panel 30);
  - C8 by 27×: the Hessian's central second differences, at a step of 0.15 SE, amplify
    surface differences of about 1e-13 on a log-likelihood near 10⁶ (panel 20). That is the size
    already seen across platforms (2.5e-10, D-107);
  - C10 by 2.4×.

**Not predicted, reported:** the share of surface cells that agree bit for bit is 27–100% by panel
(the census). Mirrored GH differs from parity only in the order of summation, so many cells agree
exactly.

## Diagnosis of the finding

- **Where the likelihood depends on ρ.** The per-period likelihood is E[p(Z)^d (1 − p(Z))^(n − d)].
  At n = 1 that is E[p(Z)] = PD or 1 − PD, with no ρ in it. ρ enters only through periods with
  n_t ≥ 2. A panel whose periods all have n_t = 1 therefore cannot identify ρ, and the engine's
  ρ̂ is then an arbitrary grid value chosen by rounding: a different sign convention, platform,
  compiler or summation order can choose another.
- **What protects the user today.** The fit carries `kFlagGridEdge` because PD̂ is on the grid
  edge, and the profile interval for ρ is the whole box, truncated at both ends. Nothing flags ρ̂
  itself as unidentified. With PD inside the box (for example 20 single-obligor periods with 3
  defaults), the edge flag would depend on where rounding puts ρ̂: on an interior ρ nothing
  guarantees a flag, since the refinement flags a flat surface only if its stencil's Hessian is
  not negative definite, and a rounding-level surface does not decide that. This case was not run;
  it is an inference from the code.
- **Why the recovery-matrix version is not run,** although the prediction's verdict rules say a
  finding that moves estimates "is then characterised statistically on the study subset". The
  mechanism is exact and fully explained, and it cannot arise in the recovery panels: every
  recovery scenario has n ≥ 100, where ρ is identified, and the nine subset panels here (1–9) pass
  every check. A statistical run would measure nothing. This departure from the registered rule is
  put to the owner in D-301.
- **Proposed mitigation (not implemented; parity is unchanged).** Either a data check that refuses
  or flags ρ when no period has n_t ≥ 2, or a flat-direction flag: if the panel surface's spread
  along an axis, at the argmax, is below a rounding threshold (such as the D-120 threshold summed
  over periods), set `kFlagFlatSurface` and report no point estimate for that parameter. Either
  one changes parity's flags, so it is the owner's decision.

## Generated output

`study_z_sign_invariance`, verbatim. "rho spread" is the range of the panel surface Σ_t ℓ_t over the
41 ρ values in the PD column of the parity argmax.

### Census

| # | Panel | T | n | d = 0 | d = n | PD hat | rho hat | fit flags | profile truncation (PD, rho) | rho spread of S at the PD column of k* | rho at k* (parity, mirrored) | bitwise cells | worst C2 (eps x scale) | worst C3 (eps x scale) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | recovery 29 r0 | 20 | 100000 | 0 | 0 | 0.0091994 | 0.0174345 | -| -, - | 40.5 | 0.0185, 0.0185 | 0.996 | 0.52 | 0.00 |
| 2 | recovery 37 r0 | 20 | 10000 | 2 | 0 | 0.010855 | 0.148364 | -| -, - | 63.1 | 0.151, 0.151 | 0.889 | 1.22 | 10.00 |
| 3 | recovery 72 r0 | 20 | 1000 | 4 | 0 | 0.0544876 | 0.282674 | -| -, - | 58.6 | 0.2966, 0.2966 | 0.762 | 1.22 | 12.00 |
| 4 | recovery 4 r0 | 40 | 10000 | 14 | 0 | 0.00107266 | 0.0135892 | near-bound | -, lo | 33.5 | 0.01317, 0.01317 | 0.702 | 1.06 | 10.00 |
| 5 | recovery 68 r0 | 40 | 100000 | 0 | 0 | 0.0445686 | 0.107441 | -| -, - | 1.35e+03 | 0.1118, 0.1118 | 0.998 | 0.61 | 0.00 |
| 6 | recovery 49 r0 | 40 | 10000 | 6 | 0 | 0.0103551 | 0.226029 | -| -, - | 244 | 0.2299, 0.2299 | 0.850 | 1.22 | 10.00 |
| 7 | recovery 7 r0 | 100 | 10000 | 43 | 0 | 0.000998942 | 0.0255829 | -| -, lo | 64.7 | 0.02593, 0.02593 | 0.670 | 1.06 | 10.00 |
| 8 | recovery 43 r0 | 100 | 10000 | 4 | 0 | 0.00943672 | 0.0881129 | -| -, - | 209 | 0.08186, 0.08186 | 0.921 | 1.22 | 10.00 |
| 9 | recovery 51 r0 | 100 | 1000 | 73 | 0 | 0.00789686 | 0.366672 | near-bound | -, hi | 67.6 | 0.3733, 0.3733 | 0.439 | 1.22 | 12.00 |
| 10 | DGP PD 0.001 rho 0.02 T 40 n 100 | 40 | 1000 | 36 | 0 | 0.00104462 | 0.001 | edge | -, lo | 2.75 | 0.001, 0.001 | 0.345 | 1.22 | 12.00 |
| 11 | DGP PD 0.0005 rho 0.3 T 20 n 50 | 20 | 500 | 20 | 0 | 0.0001 | 0.5 | edge | lo, lo+hi | 0.0254 | 0.5, 0.5 | 0.268 | 0.00 | 12.00 |
| 12 | DGP PD 0.15 rho 0.4 T 20 n 5 | 20 | 50 | 8 | 0 | 0.161502 | 0.001 | edge | hi, lo | 4.14 | 0.001, 0.001 | 0.500 | 1.40 | 12.00 |
| 13 | DGP PD 0.18 rho 0.45 T 20 n 2 | 20 | 20 | 15 | 2 | 0.161502 | 0.5 | edge | hi, lo+hi | 1.6 | 0.5, 0.5 | 0.313 | 1.92 | 14.00 |
| 14 | DGP PD 0.01 rho 0.45 T 20 n 1000 | 20 | 10000 | 12 | 0 | 0.00439969 | 0.416544 | near-bound | -, hi | 71.7 | 0.4145, 0.4145 | 0.622 | 1.06 | 10.00 |
| 15 | DGP PD 0.05 rho 0.49 T 40 n 10000 | 40 | 100000 | 3 | 0 | 0.0319273 | 0.422905 | near-bound | -, hi | 4.57e+03 | 0.4145, 0.4145 | 0.956 | 0.97 | 10.00 |
| 16 | DGP PD 0.01 rho 0.0015 T 40 n 10000 | 40 | 100000 | 0 | 0 | 0.00950608 | 0.00334122 | rejected | -, - | 107 | 0.003341, 0.003341 | 0.996 | 0.55 | 0.00 |
| 17 | DGP PD 0.0002 rho 0.1 T 40 n 10000 | 40 | 100000 | 16 | 0 | 0.000258889 | 0.116654 | -| -, - | 34.4 | 0.1118, 0.1118 | 0.784 | 1.26 | 10.00 |
| 18 | DGP PD 0.01 rho 0.12 T 20 n 1000000 | 20 | 10000000 | 0 | 0 | 0.00929432 | 0.132159 | -| -, - | 1.6e+03 | 0.1302, 0.1302 | 1.000 | 0.32 | 0.00 |
| 19 | DGP PD 0.001 rho 0.05 T 20 n 1000000 | 20 | 10000000 | 0 | 0 | 0.00131959 | 0.0374765 | -| -, - | 322 | 0.03624, 0.03624 | 1.000 | 0.49 | 0.00 |
| 20 | DGP PD 0.05 rho 0.24 T 20 n 100000 | 20 | 1000000 | 0 | 0 | 0.0408727 | 0.172079 | -| -, - | 1.93e+03 | 0.1745, 0.1745 | 1.000 | 0.43 | 0.00 |
| 21 | DGP PD 0.02 rho 0.15 T 20 n 10..1e5 | 20 | 10..100000 | 3 | 0 | 0.022783 | 0.133068 | -| -, - | 248 | 0.1302, 0.1302 | 0.830 | 1.48 | 16.00 |
| 22 | DGP PD 0.19 rho 0.2 T 20 n 1000 | 20 | 10000 | 0 | 0 | 0.177213 | 0.202816 | near-bound | hi, - | 695 | 0.2008, 0.2008 | 0.995 | 0.83 | 0.00 |
| 23 | every d = 0 (n 1000) | 20 | 10000 | 20 | 0 | 0.0001 | 0.5 | edge | lo, lo+hi | 1.4 | 0.5, 0.5 | 0.423 | 0.00 | 10.00 |
| 24 | every d = n (n 3) | 10 | 30 | 0 | 10 | 0.2 | 0.5 | edge | hi, hi | 18.3 | 0.5, 0.5 | 0.666 | 0.00 | 2.50 |
| 25 | one default, period 0 | 20 | 10000 | 19 | 0 | 0.0001 | 0.151012 | edge | lo, lo+hi | 0.366 | 0.151, 0.151 | 0.444 | 1.01 | 10.00 |
| 26 | n = 1, d 0,1 | 20 | 10 | 10 | 10 | 0.2 | 0.26192 | edge | hi, lo+hi | 2.84e-14 | 0.2619, 0.001 | 0.315 | 0.00 | 16.00 |
| 27 | d = 10 every period | 20 | 10000 | 0 | 0 | 0.00950608 | 0.001 | edge | -, lo | 51.3 | 0.001, 0.001 | 0.960 | 0.54 | 0.00 |
| 28 | d 0,500 | 20 | 10000 | 10 | 0 | 0.102598 | 0.5 | edge | -, hi | 3.96e+03 | 0.5, 0.5 | 0.711 | 0.37 | 10.00 |
| 29 | n = 4, d 0..4 | 20 | 40 | 4 | 4 | 0.2 | 0.5 | edge | hi, hi | 13.8 | 0.5, 0.5 | 0.635 | 1.84 | 14.00 |
| 30 | n = 1e6, d 0,0,1,3,0 | 20 | 10000000 | 12 | 0 | 0.0001 | 0.5 | edge | lo, hi | 1.2e+03 | 0.5, 0.5 | 0.842 | 1.27 | 8.00 |
| 31 | d = 2000 of 10000 | 20 | 100000 | 0 | 0 | 0.2 | 0.001 | edge | hi, lo | 68.4 | 0.001, 0.001 | 1.000 | 0.00 | 0.00 |
| 32 | n = 100, d 0,100 | 20 | 1000 | 10 | 10 | 0.2 | 0.5 | edge | hi, hi | 1.63e+03 | 0.5, 0.5 | 0.618 | 0.00 | 12.00 |
| 33 | n = 1e4, d 1,1,1,1,30 | 20 | 100000 | 0 | 0 | 0.000581329 | 0.15542 | -| -, - | 109 | 0.151, 0.151 | 0.899 | 1.26 | 0.00 |
| 34 | one late period, d = 3 | 40 | 10000 | 39 | 0 | 0.0001 | 0.5 | edge | lo, hi | 5.4 | 0.5, 0.5 | 0.435 | 0.88 | 10.00 |

### Checks

| # | Quantity | Bound | Compared | Over the bound | Largest | Largest / bound | Where | Largest outside panel 26 | Held |
|---|---|---|---|---|---|---|---|---|---|
| C1 | hint: mode and centre negated, scale and width equal | bitwise | 912865 | 0 | 0 | 0 | - | 0 | yes |
| C2 | surface cells, 0 < d < n: |diff| | 256 eps + 4 eps scale | 1628151 | 0 | 5.68e-14 | 0.241 | panel 30 t 3 k 995 | 7.28e-12 | yes |
| C3 | surface cells, d in {0, n}: |diff| | 1e-12 + 4 eps scale | 1047919 | 0 | 2.27e-13 | 0.103 | panel 30 t 0 k 1929 | 2.27e-13 | yes |
| C4 | panel surface sum_t l_t: |diff| | sum of per-period bounds | 85034 | 0 | 9.09e-13 | 0.133 | panel 33 k 2347 | 7.28e-12 | yes |
| C5 | grid argmax k*, and every bootstrap replicate's | identical | 34000 | 768 in panel 26; | 1 | inf | panel 26 k* | 0 | NO |
| C6 | flags: fit, quad-check count, profile, bootstrap edge/excluded/per-replicate | identical | 34102 | 0 | 0 | 0 | - | 0 | yes |
| C7 | PD, rho, loglik at the estimate, relative | 1e-12 | 102 | 1 in panel 26; | 0.996 | 9.96e+11 | panel 26 rho | 8.48e-14 | NO |
| C8 | SEs and correlation, relative | 5e-10 | 102 | 0 | 2.67e-10 | 0.535 | panel 20 se PD | 2.67e-10 | yes |
| C9 | profile end points, logit units | 2e-9 | 136 | 0 | 2.74e-12 | 0.00137 | panel 18 hi rho | 2.74e-12 | yes |
| C9 | polished maximum l_max, relative | 1e-12 | 34 | 0 | 2.25e-14 | 0.0225 | panel 18 l_max | 2.25e-14 | yes |
| C10 | bootstrap percentile end points, relative | 2e-11 | 136 | 1 in panel 26; | 0.996 | 4.98e+10 | panel 26 boot hi rho | 2.37e-12 | NO |
| C11 | quadrature check max |l(rule) - l(doubled)|, absolute diff | 1e-12 | 34 | 0 | 5.68e-14 | 0.0568 | panel 22 quad-check max | 5.68e-14 | yes |
| C12 | control (wrong hint): excess over the C2/C3 bound in panels 18, 19, 20, 30; fits flagged by the quadrature check | >= 1e-6; >= 1 | 34 | | 2.69e+06; 17 (panels 2, 3, 5, 6, 8, 9, 14, 15, 17, 18, 19, 20, 21, 22, 28, 33, 34) | | | | yes |

Verdict: FINDING, reviewed (D-301): panel 26 only, rho not identified at n = 1
