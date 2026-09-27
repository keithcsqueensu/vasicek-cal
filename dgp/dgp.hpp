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

enum class Status : std::int32_t { Ok = 0, InvalidSpec = 1, FpEnvironment = 2, FactorStreamExhausted = 3 };

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

// --- building blocks, exported for tests and the methodology's worked examples -----------------

double det_log(double x);     // natural log
double det_exp(double x);     // e^x
double det_ncdf(double x);    // Phi(x)
double det_probit(double p);  // Phi^-1(p), Wichura AS241 on det_log
double uniform_from_words(std::uint32_t a, std::uint32_t b);

// True if subnormals are neither flushed on output nor treated as zero on input.
bool fp_environment_ok();

// The factor draw for one period, with an explicit block budget (tests use a tiny one).
Status factor_draw(std::uint64_t seed, std::uint32_t scenario, std::uint32_t replicate, std::uint32_t period,
                   std::uint32_t max_blocks, double* z);

}  // namespace vcal::dgp
