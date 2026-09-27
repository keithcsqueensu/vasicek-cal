# The synthetic data-generating process (`dgp/`)

The DGP simulates default panels from the one-factor Vasicek/ASRF model, for recovery tests
and golden fixtures. Its contract is stronger than accuracy: **a panel is a pure function of
its inputs and is bitwise identical on every supported platform** (D-052 to D-057, D-109,
D-110). Everything below is part of that contract. Anyone can reproduce any panel exactly
from this description. `tools/gen_dgp_tables.py` contains a line-for-line Python
implementation that does exactly that.

## Inputs

A 64-bit seed, a scenario number, a replicate number (both 32-bit), PD and ρ (both strictly
between 0 and 1), and the number of obligors n_t for each period t = 0, 1, …, T−1.

## Random numbers: Philox4x32-10

The generator is Philox4x32 with 10 rounds (Salmon et al., 2011), exactly as in the Random123
library, and checked against Random123's published known answers.
- **Key:** the seed (low 32 bits, then high 32 bits).
- **128-bit counter:** (scenario, replicate, period, block).

Each evaluation yields four 32-bit words. Because every draw is addressed by
(scenario, replicate, period, block), adding scenarios, replicates, periods or obligors never
changes any existing draw.

**Uniforms.** Two words a, b give U = (k + ½)·2⁻⁵², with k = (a >> 6)·2²⁶ + (b >> 6), the top
26 bits of each. U is exactly representable in binary64 and lies strictly inside (0, 1).

## One period

1. **Factor Z_t (blocks 0 … 65535).** For block = 0, 1, 2, …: take U₁ from words (0, 1) and U₂
   from words (2, 3); set V₁ = 2U₁ − 1, V₂ = 2U₂ − 1, S = V₁² + V₂². If 0 < S < 1, then
   Z_t = V₁·√(−2 log S / S), and stop (Marsaglia's polar method, first normal of the pair).
   About 21% of blocks are rejected. If all 65,536 blocks were rejected, generation fails
   instead of borrowing from the default draws; this never happens in practice.
2. **Conditional PD.** p_t = Φ((Φ⁻¹(PD) − √ρ·Z_t) / √(1 − ρ)).
3. **Defaults.** d_t counts obligors i = 0 … n_t − 1 with U_i < p_t. Obligor i uses block
   65536 + ⌊i/2⌋: words (0, 1) if i is even, (2, 3) if odd. These blocks never overlap the
   factor's, so the number of rejections in step 1 cannot shift them.

## Arithmetic

Only operations that IEEE 754 makes exact or correctly rounded are used: + − × ÷, square root,
floor, multiplication by a power of two (ldexp), and reading or writing a double's bits.
No platform maths library is called, since `log`, `exp` and `erfc` differ between platforms
in the last bit. The code is compiled without fused multiply-adds (`-ffp-contract=off`,
`/fp:strict`). The DGP refuses to run if flush-to-zero or denormals-are-zero is active: it
checks at start-up that DBL_MIN/2 is non-zero and doubles back to DBL_MIN.

**log x.** Write x = m·2ᵉ with m ∈ (√½, √2]; subnormal x are first scaled by 2⁵⁴. With
s = (m − 1)/(m + 1) and z = s², log m = 2s + 2s·z·P(z), where
P(z) = Σ_{k=1..10} z^(k−1)/(2k + 1), evaluated by Horner from k = 10 down with coefficients
computed as 1/(2k + 1). Then log x = (e·ln2_lo + log m) + e·ln2_hi. Here ln2_hi is ln 2 with
32 significant bits and ln2_lo = ln 2 − ln2_hi, both from mpmath.

**exp x.** (+∞ above 800, 0 below −800.) k = ⌊x/ln 2 + ½⌋ (x multiplied by the constant
1/ln 2); r = (x − k·ln2_hi) − k·ln2_lo; e^r = Σ_{n=0..16} rⁿ/n! by Horner, with coefficients
1/n! computed by successive division; the result is ldexp(e^r, k).

**Φ(x).** φ(x) = e^(−x²/2)/√(2π).
- x < −2: φ(x)·R(−x).
- x > 2: 1 − φ(x)·R(x).
- Otherwise the series ½ + φ(x)·(x + x³/3 + x⁵/(3·5) + …), each term the previous times
  x²/(2n + 1), for exactly n = 1 … 32 (the worst case, |x| = 2, converges in 23).

R(t) is the Mills ratio, from its continued fraction 1/(t + 1/(t + 2/(t + …))) by the modified
Lentz algorithm, for exactly 128 steps (the worst case, t = 2, converges in 106). Both counts
are fixed, with no data-dependent exits (D-113).

**Φ⁻¹(p).** Wichura's AS241 (PPND16) on the log above.

## Checking a reproduction

`tests/golden/dgp/reference_panels.csv` lists six reference panels and their exact output (d_t
and Z_t as hex floats). `reference_panels.sha256` holds the SHA-256 of their canonical bytes:
for each panel in order and each period in order, d_t as a little-endian int64 followed by Z_t
as a little-endian IEEE binary64. Every CI platform must reproduce both exactly.
