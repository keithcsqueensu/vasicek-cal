// SPDX-License-Identifier: Apache-2.0
//
// Synthetic data-generating process for the one-factor Vasicek/ASRF model (M1.6; D-034,
// D-052..D-057, D-061, D-109, D-110). Runs on the CPU only.
//
// Contract: a panel is a pure function of (seed, scenario, replicate, PD, rho, n_t) and is
// bitwise identical on every supported platform. For that, the DGP:
//   - draws from Philox4x32-10 keyed by the seed, with counter (scenario, replicate, period,
//     block). Adding scenarios, replicates or periods never changes an existing panel;
//   - uses its own log, exp, Phi and Phi^-1 (dgp/det_math.cpp), built only from + - * /, sqrt,
//     floor and ldexp, which IEEE makes exact or correctly rounded, instead of the platform libm;
//   - is compiled with FP contraction off and strict FP semantics (D-048, D-056);
//   - refuses to run if flush-to-zero or denormals-are-zero is active (D-109).
// A pure-Python mirror (tools/gen_dgp_tables.py) produces reference panels that the C++
// must reproduce bit for bit on every CI platform (D-057).
//
// Per period t (block ranges D-109):
//   Z_t: Marsaglia polar method on Philox blocks 0 .. kFactorBlocks-1, one candidate pair per
//        block, the first accepted pair's first normal. Exhausting the range is an error.
//   p_t = Phi((Phi^-1(PD) - sqrt(rho) Z_t) / sqrt(1 - rho)).
//   d_t = number of obligors i < n_t with U_i < p_t, where obligor i uses block
//         kFactorBlocks + i/2, words (0,1) if i is even and (2,3) if odd.
//   A uniform from two words (a, b) is (k + 1/2) * 2^-52, k = (a >> 6) << 26 | (b >> 6):
//   exactly representable, never 0 or 1.
// Full description: docs/methodology/dgp.md.
#pragma once

#include <cstdint>

#include "dgp/philox.hpp"

namespace vcal::dgp {

inline constexpr std::uint32_t kFactorBlocks = 1u << 16;
// Bernoulli blocks run from kFactorBlocks to 2^32 - 1, two obligors per block.
inline constexpr std::int64_t kMaxObligors = 2 * (std::int64_t{1} << 32) - 2 * std::int64_t{kFactorBlocks};

enum class Status : std::int32_t {
    Ok = 0,
    InvalidSpec = 1,
    FpEnvironment = 2,
    FactorStreamExhausted = 3,
    KernelUnavailable = 4  // VCAL_DGP_KERNEL names a kernel this CPU or build lacks (count_defaults)
};

struct PanelSpec {
    std::uint64_t seed;
    std::uint32_t scenario;
    std::uint32_t replicate;
    double pd;   // 0 < pd < 1
    double rho;  // 0 < rho < 1
    const std::int64_t* obligors;  // n_t, 0 <= n_t <= kMaxObligors
    std::int64_t periods;          // 1 .. 2^32
};

// defaults[t] = d_t; factors[t] = Z_t if factors is non-null.
Status simulate_panel(const PanelSpec& spec, std::int64_t* defaults, double* factors);

// The AR(1)-factor variant (S-35, S-11; docs/methodology/dgp.md, "The AR(1) factor"; D-190):
//   Z_0 = e_0,  Z_t = phi Z_{t-1} + sqrt(1 - phi^2) e_t  (t >= 1),
// with e_t the base DGP's factor draw for period t (same key and counters) and every Bernoulli draw the
// base DGP's, so Z_t ~ N(0, 1) in every period and the variant draws no stream of its own: panels at
// different phi share their random numbers (common random numbers), and phi = 0 gives simulate_panel's
// panel bit for bit. Requires -1 < phi < 1. factors[t] = Z_t, the AR(1) factor.
Status simulate_panel_ar1(const PanelSpec& spec, double phi, std::int64_t* defaults, double* factors);

// --- building blocks, exported for tests and the methodology's worked examples -----------------

double det_log(double x);     // natural log
double det_exp(double x);     // e^x
double det_ncdf(double x);    // Phi(x)
double det_probit(double p);  // Phi^-1(p), Wichura AS241 on det_log
double uniform_from_words(std::uint32_t a, std::uint32_t b);

// True if subnormals are neither flushed on output nor treated as zero on input.
bool fp_environment_ok();

// The Bernoulli count d_t of one period: the number of obligors i < n with U_i < p, U_i drawn as above
// (P-17, D-194). Kernels: Scalar is the definition above, operation for operation; Sse2 (x86-64) and
// Avx2 (x86-64, when the CPU and the OS support it, checked at run time) compare integer keys against a
// threshold, 8 or 16 Philox blocks at a time, and return the same count bit for bit. count_defaults
// returns -1 if n is out of range or the kernel is not available on this CPU or build.
//
// simulate_panel and simulate_panel_ar1 use bernoulli_kernel_selected(): the fastest available, unless the
// environment variable VCAL_DGP_KERNEL (scalar, sse2 or avx2) names one, read once when first needed. A
// value that is not one of those, or names a kernel this CPU or build lacks, selects Unavailable, and the
// simulations then return Status::KernelUnavailable rather than quietly use another kernel.
// select_bernoulli_kernel changes the selection (the tests run every path on one machine); it is
// process-wide and must not be called while another thread is simulating. False, and no change, if the
// kernel is not available. Selecting Unavailable is allowed: the simulations then refuse, as above.
enum class BernoulliKernel : std::int32_t { Unavailable = -1, Scalar = 0, Sse2 = 1, Avx2 = 2 };
bool bernoulli_kernel_available(BernoulliKernel kernel);
BernoulliKernel bernoulli_kernel_selected();
bool select_bernoulli_kernel(BernoulliKernel kernel);
const char* bernoulli_kernel_name(BernoulliKernel kernel);  // "scalar", "sse2", "avx2", "unavailable"
std::int64_t count_defaults(std::uint64_t seed, std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                            std::int64_t n, double p, BernoulliKernel kernel);
// K, the number of integer keys k in [0, 2^52) whose uniform (k + 1/2) 2^-52 is below p: U < p iff k < K
// (docs/methodology/dgp.md, "An equivalent integer form of step 3"). 0 for p <= 2^-53 or NaN, 2^52 for
// p > 1 - 2^-53.
std::uint64_t bernoulli_threshold(double p);

// The factor draw for one period, with an explicit block budget (tests use a tiny one).
Status factor_draw(std::uint64_t seed, std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                   std::uint32_t max_blocks, double* z);

}  // namespace vcal::dgp
