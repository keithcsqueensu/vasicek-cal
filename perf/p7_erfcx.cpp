// SPDX-License-Identifier: Apache-2.0
//
// P-7's benchmark (perf/p7_erfcx.md, D-212): the per-period log-likelihood through the parity rule
// (log_contrib) and with its score and Hessian (log_contrib_derivs), on the recovery grid's cells, with each
// erfcx kernel; and erfcx alone on the arguments those periods' nodes produce.
//
// The same source builds on main, where there are no kernels: it then times main's code, the "before". On this
// branch "none" is the path with no kernel selected (every node through the scalar erfcx), and every kernel's
// values are checked against it bit for bit before timing. Single-threaded; kernels are timed in alternating order
// (A B C D, D C B A, ...) so a drift in the machine's clock affects them alike.
//
//   p7_erfcx [--seconds S]   S per (measure, kernel), default 2, split over 4 rounds
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "core/model/binomial_mixture.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "core/special/erfcx.hpp"

#if __has_include("core/special/erfcx_batch.hpp")
#include "core/special/erfcx_batch.hpp"
#define P7_HAS_KERNELS 1
#endif

namespace {

using Obj = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;

struct Period {
    Obj::Obs y;
    Obj::Theta th;
};

// The recovery grid's PD, rho and n, each period with d near n PD (0 < d < n, the Gauss-Hermite path).
std::vector<Period> periods() {
    std::vector<Period> out;
    for (const double pd : {0.001, 0.01, 0.05}) {
        for (const double rho : {0.02, 0.12, 0.24}) {
            for (const std::int64_t n : {100, 1000, 10000}) {
                for (const double f : {0.25, 1.0, 3.0}) {
                    auto d = static_cast<std::int64_t>(std::lround(f * pd * static_cast<double>(n)));
                    if (d < 1) d = 1;
                    if (d >= n) d = n - 1;
                    out.push_back({{n, d}, {pd, rho}});
                }
            }
        }
    }
    return out;
}

double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

struct Measure {
    const char* label;
    double (*run)(const std::vector<Period>&, const vcal::quadrature::SplitRule<vcal::PrecisionF64>&);
};

double run_value(const std::vector<Period>& p, const vcal::quadrature::SplitRule<vcal::PrecisionF64>& rule) {
    const Obj obj{};
    double s = 0.0;
    for (const auto& q : p) s += obj.log_contrib(q.y, q.th, rule);
    return s;
}

double run_derivs(const std::vector<Period>& p, const vcal::quadrature::SplitRule<vcal::PrecisionF64>& rule) {
    const Obj obj{};
    double s = 0.0, g[2], h[3];
    for (const auto& q : p) s += obj.log_contrib_derivs(q.y, q.th, rule, g, h) + g[0] + h[1];
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    double seconds = 2.0;
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::strcmp(argv[i], "--seconds") == 0) seconds = std::atof(argv[i + 1]);
    }
    const auto ps = periods();
    const auto rule = vcal::quadrature::parity_rule(false);
    const Measure measures[] = {{"log_contrib", run_value}, {"log_contrib_derivs", run_derivs}};

#if defined(P7_HAS_KERNELS)
    namespace k = vcal::special::erfcx_batch;
    struct K {
        const char* name;
        k::Kernel kernel;
    };
    std::vector<K> ks{{"none", k::Kernel::Unavailable}};
    for (const auto kk : {k::Kernel::Scalar, k::Kernel::Sse2, k::Kernel::Avx2}) {
        if (k::available(kk)) ks.push_back({k::name(kk), kk});
    }
    // Identity first: every kernel's sums equal "none"'s bit for bit.
    for (const auto& m : measures) {
        k::select(k::Kernel::Unavailable);
        const double want = m.run(ps, rule);
        for (const auto& kk : ks) {
            k::select(kk.kernel);
            const double got = m.run(ps, rule);
            if (std::memcmp(&got, &want, sizeof got) != 0) {
                std::printf("MISMATCH %s under %s\n", m.label, kk.name);
                return 1;
            }
        }
    }
#else
    struct K {
        const char* name;
    };
    std::vector<K> ks{{"main"}};
#endif

    std::printf("periods per pass: %zu (recovery grid cells, 0 < d < n)\n", ps.size());
    for (const auto& m : measures) {
        std::vector<double> total(ks.size(), 0.0);
        std::vector<long> passes(ks.size(), 0);
        for (int round = 0; round < 4; ++round) {
            for (std::size_t j = 0; j < ks.size(); ++j) {
                const std::size_t idx = (round % 2 == 0) ? j : ks.size() - 1 - j;
#if defined(P7_HAS_KERNELS)
                k::select(ks[idx].kernel);
#endif
                const double t0 = now();
                volatile double sink = 0.0;
                long n = 0;
                while (now() - t0 < seconds / 4.0) {
                    sink = sink + m.run(ps, rule);
                    ++n;
                }
                total[idx] += now() - t0;
                passes[idx] += n;
            }
        }
        for (std::size_t j = 0; j < ks.size(); ++j) {
            const double us = 1e6 * total[j] / (static_cast<double>(passes[j]) * static_cast<double>(ps.size()));
            std::printf("%-20s %-8s %8.3f us per period\n", m.label, ks[j].name, us);
        }
    }
    return 0;
}
