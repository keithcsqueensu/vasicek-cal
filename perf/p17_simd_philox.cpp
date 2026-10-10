// SPDX-License-Identifier: Apache-2.0
//
// P-17's benchmark (perf/p17_simd_philox.md, D-194): panel simulation with each Bernoulli kernel, on
// the recovery grid's (T, n) cells (PD 1%, rho 0.12; the work per period depends on n alone).
//
// "Before" is the scalar kernel, which is the replaced loop operation for operation. Each panel is
// simulated through the DGP's public building blocks (factor_draw, det_probit, det_ncdf, count_defaults)
// with the kernel chosen explicitly, and every panel is checked against simulate_panel bit for bit.
// Single-threaded; the kernels are timed in alternating order (A B C, C B A, ...) so that a drift in
// the machine's clock affects them alike.
//
//   p17_simd_philox [--seconds S]   S per (cell, kernel), default 2, split over 4 rounds
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "dgp/dgp.hpp"

namespace {

namespace dgp = vcal::dgp;

constexpr std::uint64_t kSeed = 0x5EED17ull;
constexpr double kPd = 0.01;
constexpr double kRho = 0.12;

// simulate_panel's loop (dgp/dgp.cpp) with the kernel chosen by the caller.
void simulate_with(dgp::BernoulliKernel kernel, std::uint32_t scenario, std::uint32_t replicate, std::int64_t T,
                   std::int64_t n, std::int64_t* d) {
    const double c = dgp::det_probit(kPd);
    const double sr = std::sqrt(kRho);
    const double s1 = std::sqrt(1.0 - kRho);
    for (std::int64_t t = 0; t < T; ++t) {
        const auto period = static_cast<std::uint32_t>(t);
        double z = 0.0;
        if (dgp::factor_draw(kSeed, scenario, replicate, period, dgp::kFactorBlocks, &z) != dgp::Status::Ok) std::abort();
        const double p = dgp::det_ncdf((c - sr * z) / s1);
        d[t] = dgp::count_defaults(kSeed, scenario, replicate, period, n, p, kernel);
    }
}

const char* name(dgp::BernoulliKernel k) {
    return k == dgp::BernoulliKernel::Scalar ? "scalar" : k == dgp::BernoulliKernel::Sse2 ? "sse2" : "avx2";
}

}  // namespace

int main(int argc, char** argv) {
    double seconds = 2.0;
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::strcmp(argv[i], "--seconds") == 0) seconds = std::atof(argv[i + 1]);
    }
    std::vector<dgp::BernoulliKernel> kernels;
    for (auto k : {dgp::BernoulliKernel::Scalar, dgp::BernoulliKernel::Sse2, dgp::BernoulliKernel::Avx2}) {
        if (dgp::bernoulli_kernel_available(k)) kernels.push_back(k);
    }
    std::printf("p17_simd_philox: selected kernel %s; %.1f s per (cell, kernel), single thread\n",
                name(dgp::bernoulli_kernel_selected()), seconds);
    std::printf("%-5s %-6s", "T", "n");
    for (auto k : kernels) std::printf(" %12s", name(k));
    for (std::size_t j = 1; j < kernels.size(); ++j) std::printf(" %9s", name(kernels[j]));
    std::printf("   (ms per panel; speed-up against scalar)\n");

    constexpr int kRounds = 4;
    long long mismatches = 0;
    for (std::int64_t T : {20, 40, 100}) {
        for (std::int64_t n : {100, 1000, 10000}) {
            const auto scenario = static_cast<std::uint32_t>(T * 100000 + n);
            std::vector<std::int64_t> d(static_cast<std::size_t>(T)), want(d.size());
            const std::vector<std::int64_t> ns(d.size(), n);
            std::vector<double> sec(kernels.size(), 0.0);
            std::vector<long long> panels(kernels.size(), 0);
            std::uint32_t replicate = 0;
            for (int round = 0; round < kRounds; ++round) {
                for (std::size_t q = 0; q < kernels.size(); ++q) {
                    const std::size_t j = round % 2 == 0 ? q : kernels.size() - 1 - q;
                    const auto t0 = std::chrono::steady_clock::now();
                    double elapsed = 0.0;
                    while (elapsed < seconds / kRounds) {
                        simulate_with(kernels[j], scenario, replicate, T, n, d.data());
                        ++panels[j];
                        ++replicate;
                        elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                    }
                    sec[j] += elapsed;
                }
            }
            // The identity check, untimed: the first 64 replicates of the cell with every kernel.
            for (std::uint32_t r = 0; r < 64; ++r) {
                const dgp::PanelSpec spec{kSeed, scenario, r, kPd, kRho, ns.data(), T};
                if (dgp::simulate_panel(spec, want.data(), nullptr) != dgp::Status::Ok) std::abort();
                for (auto k : kernels) {
                    simulate_with(k, scenario, r, T, n, d.data());
                    if (d != want) ++mismatches;
                }
            }
            std::printf("%-5lld %-6lld", static_cast<long long>(T), static_cast<long long>(n));
            std::vector<double> ms(kernels.size());
            for (std::size_t j = 0; j < kernels.size(); ++j) {
                ms[j] = 1e3 * sec[j] / static_cast<double>(panels[j]);
                std::printf(" %12.4f", ms[j]);
            }
            for (std::size_t j = 1; j < kernels.size(); ++j) std::printf(" %8.2fx", ms[0] / ms[j]);
            std::printf("\n");
        }
    }
    std::printf("panels differing from simulate_panel: %lld\n", mismatches);
    return mismatches == 0 ? 0 : 1;
}
