# S-53 Block bootstraps under an AR(1) factor: predictions made before the run

Pre-registration for study S-53 (`studies/README.md`, D-133, D-415). This file is committed **before** the
stationary bootstrap or the data-driven block length exists in the engine, and before any block bootstrap
has been run on a study panel. The results will be compared with it prediction by prediction, and misses
will be reported, not explained away.

**This file is not edited after the run.** A change of method or a sharper prediction goes in a new, dated
section appended at the end, committed before the run it concerns.

## Disclosure: what was seen before this was written

- **S-35's results on these panels** (D-192, `lrpd-serial-correlation/summary.csv`): the iid bootstrap's
  coverage at φ = 0, 0.3 and 0.6. Pooled PD coverage at φ = 0.6 was B 0.633, C 0.680, D 0.706; ρ 0.659,
  0.756, 0.819.
- **S-48's results on these panels** (D-211, `hac-se/`): the Bartlett-kernel HAC Wald interval, lag
  ⌊4(T/100)^(2/9)⌋ (2, 3, 4), gave PD coverage at φ = 0.6 of B 0.721, C 0.805, D 0.852, with a cost at
  φ = 0 of B 0.043, C 0.030, D 0.017. For ρ it covered less than the iid Wald interval.
- **S-35's M6 arm** (D-208): the AR(1)-aware MLE's PD coverage at φ = 0.6 was B 0.910, C 0.929, D 0.942.
- **Nothing else.** No block bootstrap has been run on any study panel, and no stationary bootstrap or block
  length selector exists in the engine yet.

**Evidence used in writing it:** the model and pinned results only, through [`reference.py`](reference.py)
(output `reference.csv`). It uses the pinned φ = 0 iid-bootstrap coverage of the recovery summary and S-35's
reference machinery, imported unchanged. The selector's behaviour comes from model draws on a Gaussian AR(1)
proxy series, not from study panels. As a check of the mapping, the reference's iid-bootstrap PD coverage at
φ = 0.6 (B 0.648, C 0.684, D 0.706) lies within 0.015 of what S-35 measured. The mapping was not tuned to it.

## The question

When the systematic factor is serially correlated, how much of the iid bootstrap's lost coverage do block
bootstraps restore, and does a data-driven block length (Politis and White 2004; Patton, Politis and White
2009) beat the fixed rate-optimal length? D-133 deferred the moving-block bootstrap's coverage study to M6
with AR(1) data, and named the stationary bootstrap (Politis and Romano 1994) and the data-driven length as
later `native` options. This study covers all three.

## What is run

- **The panels:** S-35's, bit for bit. These are the recovery matrix, 81 scenarios × 1,000 replicates, from
  `dgp::simulate_panel_ar1` on the recovery seed, scenario ids and replicates, with φ ∈ {0, 0.3, 0.6}.
- **The fit:** `recovery::fit_counts`, as the recovery harness fits, with B = 999 bootstrap replicates and
  95% percentile intervals (type-7 quantiles), grid-edge replicates kept at the bound (D-135, D-136).
- **Four bootstrap arms per panel, each on that panel's bootstrap stream** (`bootstrap_seed(scenario,
  replicate)` in the resampling key domain, D-132):
  1. **iid** (scheme 1): S-35's arm, rerun.
  2. **Moving block, fixed length** (scheme 2, D-133): ℓ = ⌈T^(1/3)⌉, so 3, 4 and 5 at T = 20, 40 and 100.
  3. **Moving block, data-driven length** (scheme 2): ℓ = max(1, ⌊b_CB + ½⌋), with b_CB from the selector
     below.
  4. **Stationary, data-driven mean length** (a new scheme 3): mean block length b = max(1, b_SB).
- **The stationary bootstrap** (to be specified in `docs/methodology/bootstrap.md` before its code):
  - positions i = 0 … T − 1 of a replicate, circular;
  - position 0 starts a block at index ⌊u·T⌋ from draw 0;
  - at position i ≥ 1, draw 2i − 1 is a coin: if u < 1/b a new block starts at ⌊u′·T⌋ from draw 2i, otherwise
    the index is the previous one plus 1, modulo T;
  - every draw's address depends on the position only, never on earlier coins, so the stream is fixed for a
    given (seed, scheme, replicate).
- **The block length** (Politis and White 2004, with Patton, Politis and White 2009's corrected constants),
  from the panel's series x_t = Φ⁻¹((d_t + ½)/(n + 1)), the probit of each period's continuity-corrected
  default rate (owner's choice, 2026-11-11). It is `politis_white` in `reference.py`, which the engine's
  version must match:
  - x is demeaned, and the autocovariances use the 1/T divisor;
  - K_N = max(5, ⌈√log₁₀ T⌉), m_max = ⌈√T⌉ + K_N, b_max = ⌈min(3√T, T/3)⌉, and the threshold is
    c = 2√(log₁₀ T / T);
  - m̂ is the smallest m whose next K_N autocorrelations all lie within ±c (m_max if none), and
    M = min(2·max(m̂, 1), m_max);
  - with the flat-top window λ(s) = 1 for |s| ≤ ½, 2(1 − |s|) for ½ < |s| ≤ 1 and 0 otherwise, at s = k/M:
    ĝ = Σ_{|k|≤M} λ R̂(k) and Ĝ = Σ λ |k| R̂(k);
  - b_SB = (2Ĝ²/(2ĝ²))^(1/3) T^(1/3) and b_CB = (2Ĝ²/((4/3)ĝ²))^(1/3) T^(1/3), each capped at b_max, and
    1 if R̂(0) ≤ 0 or ĝ ≤ 0.

  The selected lengths are recorded per panel.
- **Checks before scoring:**
  - at every φ, the iid arm's per-scenario covered counts for PD and ρ must equal S-35's exactly;
  - the moving-block and stationary resampling indices must match the Python mirror's reference draws bit
    for bit;
  - the selector must match `reference.py`'s `politis_white` on its committed goldens.

  Any check failing stops the scoring.
- **Order:** the study subset (D-150, 9 scenarios) first, scored and reported as a subset; then the full
  matrix. The runtime projection (D-187) is recorded before the full run.

## Mechanism behind the predictions

- **What a block bootstrap captures.** To first order, a bootstrap of periods reproduces the per-period
  score's long-run covariance I + 2 Σ w_k C_k, with weights w_k:
  - iid: w_k = 0;
  - moving block: the Bartlett window 1 − k/ℓ for k < ℓ;
  - stationary: Politis–Romano's (1 − k/T) q^k + (k/T) q^(T−k), q = 1 − 1/b.

  The truth is I + 2 Σ_{k<T} (1 − k/T) C_k, with C_k = Σ_j a_j a_j′ φ^(jk) from S-35's Hermite expansion. The
  remaining shortfall r maps the pinned φ = 0 iid-bootstrap coverage through S-35's coverage mapping, with
  S-35's extra downward bias δ for ρ.
- **Short blocks capture only the short lags.** At φ = 0.6, PD's score correlation decays like 0.6^k. So
  ℓ = 3 to 5 recovers part of the lost variance, not all of it. The reference puts the fixed rule's PD
  coverage at φ = 0.6 at B 0.80, C 0.85, D 0.88 (iid 0.65, 0.68, 0.71). ρ's correlation decays like 0.36^k,
  so short blocks recover more of it: 0.77, 0.85, 0.90.
- **The selector undershoots at these T.** On the Gaussian proxy, the flat-top window and the small-sample
  shrinkage of R̂(k) give medians well below the selector's own population target:

  | | Oracle b_CB | Selected median | ℓ = 1 |
  |---|---|---|---|
  | φ = 0.6, T = 100 (D) | 7.7 | 4.0 | 1% |
  | φ = 0.6, T = 20 (B) | 4.6 | 1.8 | 27% |
  | φ = 0 | 1 | 1.4 | 55% |

  So the data-driven arms sit between iid and the fixed rule. Their PD coverage at φ = 0.6 is B 0.73,
  C 0.79–0.80, D 0.86–0.87.
- **The stationary and moving-block bootstraps at the same length nearly agree.** Their weights differ
  little when b ≪ T. The reference differs by at most 0.01.
- **What the first-order model leaves out:** the extra variance of using fewer, longer blocks. With T = 20
  and ℓ = 3 there are about 7 blocks per replicate. That costs coverage the model does not see (S-48's HAC
  interval cost 0.043 at T = 20 and φ = 0). The model therefore favours long blocks, and is held loosely for
  the fixed rule at T = 20 (B4 gives a range rather than a point).

## Predictions

Coverage is pooled over all replicates of a group's scenarios (groups A–D as in
`docs/methodology/bootstrap_predictions.md`); a non-finite interval misses. "Reference" is
`reference.csv`'s group mean. Group A is reported, not predicted.

| # | Prediction |
|---|---|
| B1 | **The fixed-length moving block follows the model where T ≥ 40:** in groups C and D at φ = 0.3 and 0.6, its pooled PD coverage is within ±0.04 of the reference (C 0.904 / 0.851, D 0.920 / 0.880), and its ρ coverage within ±0.05 (C 0.887 / 0.853, D 0.918 / 0.898) |
| B2 | **The data-driven arms follow the model where T ≥ 40:** in groups C and D at φ = 0.3 and 0.6, each data-driven arm's pooled PD coverage is within ±0.04 of its reference (moving block C 0.867 / 0.794, D 0.898 / 0.859; stationary C 0.866 / 0.804, D 0.897 / 0.866), and ρ within ±0.05 (moving block C 0.877 / 0.826, D 0.913 / 0.890; stationary C 0.876 / 0.827, D 0.912 / 0.890) |
| B3 | **The order at φ = 0.6, for PD:** in each of groups B, C and D, each data-driven arm covers more than iid by at least 0.02, and the fixed-length moving block covers more than both data-driven arms: by at least 0.02 in B and C, and by any margin in D (reference 0.014–0.021) |
| B4 | **T = 20 is held loosely:** in group B at φ = 0.6, the fixed-length moving block's PD coverage is 0.70–0.81 (reference 0.799, which omits the cost of few blocks), and each data-driven arm's is 0.66–0.75 (reference 0.725) |
| B5 | **Blocking costs little under independence:** at φ = 0, every block arm's pooled PD and ρ coverage is at most 0.06 below iid's in group B, at most 0.04 in C and at most 0.03 in D. The data-driven arms cost at most 0.02 in every group, since they mostly pick ℓ near 1 |
| B6 | **The two data-driven schemes agree:** in every group (B, C, D) and φ, the stationary and the data-driven moving block's pooled PD and ρ coverages are within 0.03 of each other |
| B7 | **The selector undershoots its target:** at φ = 0.6, the median selected b_CB over group D's panels is 2.5–5.5 (proxy 4.0; oracle 7.7), and over group B's 1.2–2.8 (proxy 1.8). At φ = 0, the share of panels with ℓ = 1 is 0.35–0.75 in every group (proxy 0.55) |
| B8 | **No block bootstrap reaches the band for PD at φ = 0.6:** in every group, every block arm's pooled PD coverage is below 0.927 |
| B9 | **The fixed-length moving block does at least as well as the HAC Wald for PD:** at φ = 0.6, in each of groups B, C and D, its pooled PD coverage is at least S-48's HAC coverage minus 0.02 (B 0.721, C 0.805, D 0.852) |

Reported without a prediction:
- group A;
- interval widths in logit units against iid's;
- the distribution of the selected lengths by scenario and φ, and their correlation with the panel's own lag-1
  autocorrelation;
- the AR(1) MLE (S-35's M6 arm) beside every arm, as the model-based benchmark;
- grid-edge replicate counts, and the runtime.

## How the comparison will be reported

- **Per prediction:** HELD or NOT HELD, with the numbers, in `results.txt` and in the studies index entry.
- **A miss is reported as a miss.** If its cause is found afterwards, it is recorded as a diagnosis, not a
  revision.
- **The recommendation** (a `native` default, or none) is the owner's decision after the run, with its own
  D-entry. If the fixed rule beats the data-driven lengths, as predicted, the data-driven lengths are still
  kept as options, since longer histories than this matrix's would favour them.
