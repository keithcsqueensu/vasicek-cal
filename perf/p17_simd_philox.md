# P-17: SIMD Philox for the DGP's Bernoulli draws (D-194)

Changes no result: every panel is bit for bit the same, under every kernel, on every platform. Machine:
Intel Core Ultra 9 275HX (8 P-cores, 16 E-cores, 24 threads); MSVC 19.51 and GCC 16.2 (WinLibs), Release.

## What changed

The DGP's cost is the Bernoulli count d_t: one Philox4x32-10 call per two obligors, then a uniform and a
comparison per obligor. `count_defaults` (`dgp/bernoulli.cpp`) now has three kernels:

- **scalar:** the specification's loop, operation for operation; any CPU.
- **SSE2:** 4 lanes × 4 registers interleaved; every x86-64 CPU, no compiler flag.
- **AVX2:** 8 lanes × 2 registers, in `dgp/bernoulli_avx2.cpp`, the only file compiled for AVX2. It runs
  only after a run-time check that the CPU and the OS support AVX2. It contains nothing with external
  linkage but its entry point, so no AVX2-compiled inline function can be merged into code that runs on
  other CPUs.

The vector kernels (`dgp/philox_lanes.hpp`) put one Philox block in each 32-bit lane. The 32 × 32 → 64-bit
products come from `mul_epu32` on the even lanes and on the odd lanes shifted down, and the halves are put
back in their lanes, so each lane computes exactly `philox4x32`. They also never form a uniform. Obligor i
defaults iff U_i < p, with U_i = (k + ½)·2⁻⁵². That holds iff k < K, where K = ⌈p·2⁵² − ½⌉ clipped to
[0, 2⁵²] (`bernoulli_threshold`; both floating-point steps are exact). With k = hi·2²⁶ + lo the comparison
is hi < K_hi or (hi = K_hi and lo < K_lo). All operands are below 2³¹, so SSE2's signed 32-bit compares
are exact. `docs/methodology/dgp.md` states the integer form for a reproduction.

**The kernel in use.** The default is the fastest available. The environment variable
`VCAL_DGP_KERNEL=scalar|sse2|avx2` forces one; a value naming a kernel the machine lacks makes simulation
refuse (`Status::KernelUnavailable`, `VCAL_E_UNSUPPORTED` through the ABI) rather than use another. The
kernel is recorded in `vcal_build_info` (`dgp_kernel=`), in every timing sample (`# dgp_kernel,`), in
S-10's checkpoint configuration, and on stderr by every tool that simulates panels.

## Identity

- **Every path on one machine.** `unit_dgp` runs both reference-panel SHA-256 checks (base and AR(1),
  against the Python mirror) under scalar, SSE2 and AVX2 in turn. It also compares 4,630 counts of each
  vector kernel with the scalar kernel:
  - every n from 0 to 80 (each tail length after a vector step, odd and even), up to n = 2²⁰ + 3;
  - p at the threshold's edges and beyond either end (0, −0, 1, 2, ∞, NaN, the smallest subnormal);
  - p exactly on the first, middle and last obligors' uniforms, which are threshold steps, and one ulp
    either side.

  Every CI job runs the full suite, and the runners have AVX2. So GCC 11 and 14, Clang 18, MSVC and both
  nvcc containers each run all three paths against the reference hashes, and their report steps print each
  kernel's hashes. `unit_dgp_kernel_{scalar,sse2,avx2,none-such}` rerun the panel checks with the kernel
  forced through `VCAL_DGP_KERNEL`, and check that the start-up selection follows it, including refusal for
  an invalid value. On a machine without a kernel, its entry exits 77 and CTest reports it as skipped, with
  the reason, never as passed. This was checked on a build with `-DVCAL_DGP_AVX2=OFF`.
- **The threshold at its edges.** `bernoulli_threshold_at_its_edges` checks:
  - K = 0 at and below U(0) = 2⁻⁵³;
  - K = 1 one ulp above, the smallest p that changes K;
  - K = 2⁵² − 1 at U(2⁵² − 1) = 1 − 2⁻⁵³, and K = 2⁵² from the next double, which is 1;
  - for 20,011 keys, K(U(k)) = k, K one ulp below = k, and K one ulp above = k + 1. The keys include both
    sides of the 2²⁶ word boundary and pseudo-random ones, with junk in the 6 discarded bits.
  - and U < p ⇔ k < K directly on the words.

  A deliberate mutant with "≤" semantics (K = ⌊x⌋ + 1) fails it and the count comparison at once.
- **Against simulate_panel.** The benchmark below checks 64 panels per cell and kernel against
  `simulate_panel` (CTest `p17_simd_philox_identity` runs it briefly).
- **On a study's rows.** In the A/B/B/A runs below, all four `study_parametric_bootstrap` row files are
  byte-identical.

## Speed: the kernels

`p17_simd_philox --seconds 2`: one thread, ms per panel at PD 1%, ρ 0.12, kernels timed in alternating order.

| T | n | scalar (MSVC) | AVX2 (MSVC) | SSE2 / AVX2 speed-up, MSVC | SSE2 / AVX2 speed-up, GCC |
|---|---|---|---|---|---|
| 20 | 100 | 0.0179 | 0.0111 | 1.41× / 1.61× | 1.20× / 1.33× |
| 20 | 1,000 | 0.0961 | 0.0300 | 2.14× / 3.20× | 1.80× / 2.60× |
| 20 | 10,000 | 0.888 | 0.216 | 2.43× / 4.12× | 2.08× / 3.53× |
| 100 | 100 | 0.0897 | 0.0557 | 1.41× / 1.61× | 1.20× / 1.33× |
| 100 | 1,000 | 0.480 | 0.150 | 2.14× / 3.20× | 1.81× / 2.60× |
| 100 | 10,000 | 4.43 | 1.08 | 2.41× / 4.10× | 2.08× / 3.56× |

(T = 40 matches T = 20 per period.) The time per panel is linear in T. At n = 100 the factor draw and Φ,
which are unchanged, are a large share of the period. GCC's scalar loop is faster than MSVC's, so its
ratios are smaller, but the AVX2 times are about the same with either compiler (0.97 and 1.08 ms at
T = 100, n = 10,000). A scratch prototype of the kernels built with clang 23 (not the committed build) gave
AVX2 2.4–3.4× its scalar loop at n = 10,000.

## Speed: a parametric-bootstrap run

`study_parametric_bootstrap --replicates 2` over all 81 scenarios, `main` at 1d2a868 (A) and this branch
(B). The order was A, B, B, A, after a 120 s warm-up (D-193). Thread-seconds per sample, and the projected
S-10/S-4 run (`perf/project_runtime.py`, 24 threads, efficiency 0.997):

| run | wall s | total thread-s | simulate | projected run |
|---|---|---|---|---|
| A | 120 | 1,286 | 243 | 4.47 h (4.27–4.66) |
| B | 111 | 1,115 | 77 | 3.51 h (3.30–3.71) |
| B | 111 | 1,111 | 75 | 3.45 h (3.25–3.63) |
| A | 119 | 1,262 | 241 | 4.40 h (4.19–4.59) |

The simulate phase falls 3.2×, and the projected run 21% (381,839 → 299,703 thread-s on the means). The
21% is a projection, calibrated and not validated (D-187): the next parametric-bootstrap full run records
its projection before it starts and confirms or corrects it against the actual time. The
other phases are unchanged within noise. W* (S-4b) also draws panels, but its cost is the fits. These are
projections from the sampled tool, not a full run. Both binaries predate D-193's fill split, so these
absolute hours are not comparable with D-193's; the A/B ratio is the result.
