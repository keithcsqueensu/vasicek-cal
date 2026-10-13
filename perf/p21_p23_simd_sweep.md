# P-21 to P-23: a SIMD sweep of the CPU paths, and a scope note for P-6 (D-213)

**Proposals only.** A sweep of the CPU paths for vectorisation not covered by P-6, P-7 (D-212), P-17 or P-18.
D-213 records it; each item still needs its own decision entry, benchmark and identity test under the
[P-item rules](README.md#rules-for-every-p-item-d-177) before it merges, and its place in the order is the
owner's call. All three keep results bit for bit: one lane per independent evaluation, each lane performing
the scalar operations in the scalar order, FMA contraction off, every kernel path tested as P-17.

**Relation to P-7 (D-212, #67).** The measurements below were taken before P-7. P-7 computes every node's
`erfcx(|x|/√2)` in one vector batch for the *value* (`BinomialLogIntegrand::batch`). It measured
`log_contrib` 33% faster but `log_contrib_derivs` only 0–11%, because the derivative functor's two
`inverse_mills` calls recompute that same `erfcx` per node, which is exactly P-21's target. After P-7, P-21's
natural form is to hand the batch's per-node `erfcx` values (and exp pair) to the derivative functor, rather
than caching them in the scalar integrand.

**Measurements** (`simd_scout`; parity rule, PD 1%, ρ 0.12; one thread on a cloud VM, Intel Xeon 2.1 GHz,
GCC 13.3, Release; indicative only, to be repeated on the owner's machine before any decision):

| (n, d) | `log_contrib` | `log_contrib_derivs` | derivs / value | hint alone | hint / value | integrand evaluations (nodes + search) |
|---|---|---|---|---|---|---|
| 100, 1 | 14.8 µs | 24.9 µs | 1.68 | 1.7 µs | 12% | 128 + 0 |
| 1,000, 5 | 14.6 | 25.7 | 1.76 | 1.5 | 10% | 128 + 0 |
| 10,000, 20 | 17.0 | 26.7 | 1.57 | 5.1 | 30% | 128 + 0 |
| 100, 0 | 50.2 | 61.7 | 1.23 | 1.3 | 2.5% | 256 + 113 |
| 1,000, 0 | 52.0 | 71.6 | 1.38 | 1.2 | 2.3% | 256 + 115 |
| 10,000, 0 | 53.8 | 64.6 | 1.20 | 1.6 | 2.9% | 256 + 116 |

Per call: `log_phi_pair` 68–99 ns (two runs), `inverse_mills(x)` + `inverse_mills(-x)` 90 ns, `erfcx` 27–30 ns,
the composite rule's node transform (`sinh`, `cosh`, `log`) 28–32 ns. Two things stand out. The derivative
path pays about as much for its inverse-Mills pair as for the value itself. And a one-sided period
(d ∈ {0, n}, common in low-PD panels) costs about 3.4× an interior one, a third of it in the extent search.

## P-21: one `erfcx` per node on the derivative path

**What.** At each node `log_contrib_derivs` calls `log_phi_pair(x)` (one `erfcx(|x|/√2)`, one exp pair
`exp(−s.hi/2)·exp(−s.lo/2)`), then `BinomialLogIntegrandDerivs` calls `inverse_mills(x)` and
`inverse_mills(−x)`, which recompute the same `erfcx(|x|/√2)` twice and the same exp pair once. The node's x
is the same expression in both (`(c − √ρ z)/√(1 − ρ)`). So λ(±x) can be finished from the pair's own
e = `erfcx(|x|/√2)` and E = `exp(−s.hi/2)·exp(−s.lo/2)`, each by `inverse_mills`'s own expression for its
sign: λ(−|x|) = √(2/π)/e and λ(|x|) = (1/√(2π))·E/(1 − 0.5·e·E), as `log_phi_pair` does for log Φ. Three
`erfcx` and four `exp` per node become one and two. `detail::mixture_slopes` (the hint's Newton) calls the same
inverse-Mills pair at 0 < d < n and shares the change.

**Changes results:** no. Each value is the expression it replaces, on the same operands; the test is
`inverse_mills_pair_is_two_inverse_mills_calls_bit_for_bit` on the `log_phi_pair` test's arguments, plus
`log_contrib_derivs` (value, score, Hessian) bit for bit on the recovery grid's distinct rows.

**Cost of the change.** The accumulator calls the derivative functor only for nodes whose weight does not
underflow, after the value. The tails must be handed from the value to the derivatives: either the integrand
caches them per node, or one fused integrand returns the value and the five node functions together (it
then computes them for underflowed nodes too, discarded as before). Scalar code; no kernel paths.

**Expected:** about a third of a derivative evaluation at 0 < d < n (about 65 of 195 ns per node, from the
per-call costs above), less at d ∈ {0, n} where only one inverse-Mills call is made. That path is the
profile solves: W* of S-4b (39.8% of the S-10 run, D-199), the analytic SE (13.5%) and the refinement.
The run-level saving should be projected with the D-187 method before it is built.

## P-22: SIMD derivative node functions

**What.** After P-7 and P-21, everything `BinomialLogIntegrandDerivs` does per node is arithmetic on
(z, x, λ(x), λ(−x)): x_ρ, x_ρρ, g_x, g_xx and the five node functions, then `m_[j] += w·f[j]` and the
rescale `m_[j] *= r`. Lanes over nodes (four in AVX2, two in SSE2), the node functions written to scratch
and folded in node order by the existing accumulator, as P-18 does with its pair terms. Divisions stay
exact divisions, no reciprocals. The remaining libm calls (the exp pair, the fold's exp) stay scalar.

**Changes results:** no; scalar, SSE2 and AVX2 kernels, each bit for bit against the scalar functor, the
`VCAL_*_KERNEL` override and build-info record as P-17 and P-18.

**Depends on** P-7's node batching (the same block of nodes, the same scratch) and on P-21 (without it the
lanes wait on three scalar `erfcx`).

## P-23: hint and extent search across grid points

**Rejected (owner, 2026-10-13; D-405), before it was built.** The row fill it targets runs inside the surface-row
cache's fill (P-1), once per distinct (n, d) per scenario. In S-10's tool that fill is 2.4 s once per scenario,
against about 1.5 s of per-replicate work; at R = 1,000 it is under 0.2% of the run. The proposal below is kept as
written.

**What.** Two serial searches sit in front of every integral: the hint's bracket and Newton
(`binomial_mixture_hint`, 10–30% of an interior period's cost, 30% at n = 10,000) and the composite rule's
extent bisections (113–116 of about 370 integrand evaluations at d ∈ {0, n}). Neither vectorises within
one integral: each step depends on the last. Across grid points they do. A surface row is one observation
at every grid point, so neighbouring grid points share n and d, take the same branches (one-sided or not)
and differ only in (PD, ρ). The row fill (`evaluate_surface`, `SurfaceRowCache::row`) would evaluate four
grid points per call, each lane running its own search in the scalar order, finished lanes masked until
the slowest lane stops (Newton's and the bisection's iteration counts differ per lane). The predicate
`g(z) − z²/2 < target` is one `log_phi` per lane, which P-7's vector `erfcx` serves.

**Changes results:** no: each lane's iterates are the scalar iterates. Test: every row of the recovery
subset's distinct observations, bit for bit against the scalar row fill, under every kernel, including
lanes that stop at the iteration cap and the masked tail when K is not a multiple of the lane count.

**Considered and not proposed:** a four-way section search inside one integral (4 probes per step, about
2.3× fewer steps). It returns the same end point only if `below` is monotone along the ray in floating
point, which rounding near the target does not guarantee; it would be a within-rounding change with its own
flip check, while lanes over grid points give the same gain bit for bit.

**Cost of the change.** The row fill becomes a blocked loop over grid points with a scalar tail; the
`parallel_for` stays over blocks, so results still do not depend on the thread count. The profile solves
evaluate off-grid points one at a time and do not benefit.

## Scope note for P-6

P-6's entry describes a SIMD reduction across grid points, as in `reduce_weighted`, but since D-172 the engine's default is the bounded argmax
(`ArgmaxSearch::Bounded`); `reduce_weighted` runs only when its guard fails. The 5–7% of the run that P-6
targets is therefore mostly `weighted_sum` inside visited tiles. Proposed for P-6's own decision entry:

- **Lanes = one tile row.** Tiles are 4 × 4 grid points (`kArgmaxTile`), and k is contiguous along the ρ
  axis, so a tile row is four adjacent k: one AVX2 vector, two SSE2. L is read at stride K per period with
  no gather; each lane keeps its own Neumaier sum and non-finite tally in the scalar order. On the recovery
  grid (41 ρ points) the last tile column is one point wide and is masked.
- **The tile bounds,** `sum_j w_j U[j][tile]` and its magnitude for every tile of every replicate, with
  lanes over tiles (`U` is contiguous in the tile index).
- **The `U` table,** a per-tile `fmax`, exact in any order.
- **The summed-surface calls** with unit weights in `posterior.hpp`, `profile.hpp` and
  `conditional_pd.hpp`, which call the same `weighted_sum`.

Not SIMD, noted while reading: each replicate sorts all tile bounds (176 on the recovery grid) and usually
visits a few; a heap pops them in the same order (bound descending, then index) for less work. Measure
before deciding.

## Unlocked only by the in-house `log`/`exp` decision (D-199)

These remain scalar under bit-for-bit rules because they call libm; they belong in that decision's
cost-benefit, not in new P-items: `log` and `log1p` in `log_phi`, its exp pair, the online log-sum-exp's
`exp` per node, the composite rule's `sinh`, `cosh` and `log` per node (28–32 ns, about a third of a
one-sided node), the posterior's mass exponentials and the AR(1) filter's `weighted_moments`. With an
in-house `log`/`exp`, a whole `log_contrib` could run with lanes over grid points and the online fold kept
per lane, with no scratch buffer: P-23's layout extended to the integral itself.

## Checked and not proposed

- **Bootstrap index draws** (`resample/weights.hpp`): scalar Philox, one block per two draws. P-17's
  lanes would apply bit for bit, but S-10 draws about 5 × 10⁴ blocks per panel against thousands of
  quadratures per distinct row: well under 1%.
- **Grid-point setup** (`make_vasicek1f` per grid point: `probit(PD)` and two square roots, recomputed for
  every ρ): tabulating per axis value is exact, but it is tens of ns against 15–50 µs per evaluation.
- **The DGP's factor draws** (one polar Box–Muller pair per period): T per panel, negligible beside the
  Bernoulli count P-17 already vectorised.

