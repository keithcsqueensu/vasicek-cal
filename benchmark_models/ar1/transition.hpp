// SPDX-License-Identifier: Apache-2.0
//
// The AR(1) filter's transition sums (P-18, D-199): the log-sum-exps over node pairs that dominate an evaluation and
// its score, with the in-house exp (common/detmath/det_exp.hpp, shared with the DGP), so they can run one lane per
// output node. Three kernels, chosen at run time like the DGP's Bernoulli count (P-17): scalar, SSE2 (two lanes) and
// AVX2 (four lanes, only on a CPU that has it). Every lane does the scalar kernel's operations in the scalar
// kernel's order, so the three give the same results bit for bit (tests/unit/ar1_transition_test.cpp).
// VCAL_AR1_KERNEL=scalar|sse2|avx2 forces one; a value this CPU or build lacks selects none, and the kernels then
// refuse to run (the tests skip visibly, exit code 77).
//
// The terms, with A_k = 0.5 y_k^2 + log sigma (the proxy's terms at destination node k), c = 1 / (2 (1 - phi^2)),
// hlv = 0.5 log(1 - phi^2) and pu_j = phi u_j, in this operation order (filter.hpp's log_ratio):
//     R(j, k) = A_k - ((z_k - pu_j)^2 c + hlv)
//   forward:   for each destination k, LSE over sources j = 0, 1, ... of  lhat_j + R(j, k)
//   backward:  for each source j, LSE over destinations k = 0, 1, ... of  (R(j, k) + v_k) + b_k,
//              optionally writing each term to e[j * nn + k]
// Each LSE is the online one of core's OnlineLogSumExp (quadrature/integrator.hpp), term by term in index order:
// NaN and +inf are flagged, -inf is skipped, and otherwise, with d = x - max,
//     d > 0:   sum = sum * exp(-d) + 1,  max = x;     else:  sum = sum + exp(d)
// (exp(-|d|) in both branches). A term with d < kExpZeroBelow is skipped: exp is exactly +0 there and sum >= 1, so
// sum + 0 = sum and skipping changes nothing.
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace vcal::benchmark_models::ar1::kernel {

enum class Kernel { Scalar, Sse2, Avx2, Unavailable };

bool available(Kernel k);
const char* name(Kernel k);
Kernel selected();       // VCAL_AR1_KERNEL if set, else the fastest available; select() overrides
bool select(Kernel k);   // false if k is not available (Unavailable selects none: every call then refuses)

// exp(x) is exactly +0 for x < kExpZeroBelow, for det_exp and for the platform's std::exp alike (tested).
inline constexpr double kExpZeroBelow = -750.0;

// Optional counts (scalar kernel only; the measurement of P-18 item 1): terms seen, and those skipped because
// their exp is exactly 0.
struct Counts {
    std::uint64_t terms = 0;
    std::uint64_t below = 0;
};

// One LSE per output lane: its running max and sum, and flags (1: a NaN term, 2: a +inf term).
enum : unsigned char { kFlagNan = 1, kFlagInf = 2 };

// All return false (and compute nothing) if no kernel is selected.
bool forward(const double* pu, const double* lhat, std::size_t np, const double* z, const double* a, std::size_t nd,
             double c, double hlv, double* mx, double* sum, unsigned char* flags, Counts* counts = nullptr);
bool backward(const double* pu, std::size_t nc, const double* z, const double* a, const double* v, const double* b,
              std::size_t nn, double c, double hlv, double* mx, double* sum, unsigned char* flags, double* e,
              Counts* counts = nullptr);
// out[i] = det_exp(in[i] - shift), as the score's pair weights.
bool exp_shifted(const double* in, std::size_t n, double shift, double* out, Counts* counts = nullptr);
// det_exp(x) for each x, through the selected kernel (for the tests).
bool exp_lanes(const double* in, std::size_t n, double* out);

// The LSE's value, as OnlineLogSumExp::result.
inline double lse_result(double mx, double sum, unsigned char flags) {
    if (flags & kFlagNan) return static_cast<double>(NAN);
    if (flags & kFlagInf) return HUGE_VAL;
    if (mx == -HUGE_VAL) return -HUGE_VAL;
    return mx + std::log(sum);
}

}  // namespace vcal::benchmark_models::ar1::kernel
