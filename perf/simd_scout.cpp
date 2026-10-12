// SPDX-License-Identifier: Apache-2.0
//
// The measurements behind perf/README.md's proposed SIMD items (P-21 to P-25). Single-threaded, one
// observation at a time, PD 1% and rho 0.12, the parity rule (128-node adaptive Gauss-Hermite; 16 x 16
// composite Gauss-Legendre for d in {0, n}). It times, per period:
//   value    BinomialMixture::log_contrib
//   derivs   log_contrib_derivs (the profile solves, W* of S-4b and the analytic SE)
//   hint     binomial_mixture_hint alone
// and counts the integrand's evaluations, split into the rule's nodes and the composite rule's extent search.
// It also times the per-node pieces (log_phi_pair, the inverse-Mills pair, erfcx, the composite rule's node
// transform). Changes nothing; run by hand: `simd_scout [--reps N]`.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "core/model/binomial_mixture.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "core/special/erfcx.hpp"
#include "core/special/inverse_mills.hpp"
#include "core/special/log_phi.hpp"

namespace {

using Clock = std::chrono::steady_clock;
volatile double g_sink = 0.0;

double seconds_since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

struct Counted {
    vcal::model::BinomialLogIntegrand f;
    std::int64_t* count;
    double operator()(double z) const {
        ++*count;
        return f(z);
    }
};

}  // namespace

int main(int argc, char** argv) {
    int reps = 2000;
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::strcmp(argv[i], "--reps") == 0) reps = std::atoi(argv[i + 1]);
    }
    using Obj = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
    const auto rule = vcal::quadrature::parity_rule();
    const Obj obj{};
    const double pd = 0.01, rho = 0.12;
    const Obj::Theta th{pd, rho};
    const std::int64_t cases[][2] = {{100, 1}, {1000, 5}, {10000, 20}, {100, 0}, {1000, 0}, {10000, 0}};

    std::printf("n,d,value_us,derivs_us,derivs_over_value,hint_us,hint_share_of_value,g_evals,nodes,search_evals\n");
    double s = 0.0;
    for (const auto& c : cases) {
        const Obj::Obs y{c[0], c[1]};
        auto t0 = Clock::now();
        for (int r = 0; r < reps; ++r) s += obj.log_contrib(y, th, rule);
        const double value_us = seconds_since(t0) / reps * 1e6;
        double grad[2], hess[3];
        t0 = Clock::now();
        for (int r = 0; r < reps; ++r) s += obj.log_contrib_derivs(y, th, rule, grad, hess);
        const double derivs_us = seconds_since(t0) / reps * 1e6;
        const auto f = vcal::model::make_binomial_log_integrand(vcal::model::make_vasicek1f(pd, rho), c[0], c[1]);
        t0 = Clock::now();
        for (int r = 0; r < reps; ++r) s += vcal::model::binomial_mixture_hint(f).mode;
        const double hint_us = seconds_since(t0) / reps * 1e6;
        std::int64_t count = 0;
        const Counted counted{f, &count};
        const auto hint = vcal::model::binomial_mixture_hint(f);
        s += rule.log_integrate(counted, hint);
        const std::int64_t nodes = hint.one_sided ? std::int64_t{rule.one_sided.rule.n} * rule.one_sided.panels
                                                  : std::int64_t{rule.interior.rule.n};
        std::printf("%lld,%lld,%.2f,%.2f,%.2f,%.2f,%.3f,%lld,%lld,%lld\n", static_cast<long long>(c[0]),
                    static_cast<long long>(c[1]), value_us, derivs_us, derivs_us / value_us, hint_us, hint_us / value_us,
                    static_cast<long long>(count), static_cast<long long>(nodes), static_cast<long long>(count - nodes));
    }

    const int m = reps * 1000;
    std::vector<double> xs(1024);
    for (std::size_t i = 0; i < xs.size(); ++i) xs[i] = -6.0 + 12.0 * static_cast<double>(i) / 1023.0;
    const auto x_at = [&](int i) { return xs[static_cast<std::size_t>(i & 1023)]; };
    auto t0 = Clock::now();
    for (int i = 0; i < m; ++i) {
        const auto p = vcal::special::log_phi_pair(x_at(i));
        s += p.at_x + p.at_neg_x;
    }
    const double pair_ns = seconds_since(t0) / m * 1e9;
    t0 = Clock::now();
    for (int i = 0; i < m; ++i) s += vcal::special::inverse_mills(x_at(i)) + vcal::special::inverse_mills(-x_at(i));
    const double mills_ns = seconds_since(t0) / m * 1e9;
    t0 = Clock::now();
    for (int i = 0; i < m; ++i) s += vcal::special::erfcx(std::fabs(x_at(i)) * 0.7071067811865476);
    const double erfcx_ns = seconds_since(t0) / m * 1e9;
    t0 = Clock::now();
    for (int i = 0; i < m; ++i) {
        const double u = 0.5 * x_at(i);
        s += std::log(0.01 * 0.3 * std::cosh(u)) + (0.1 + 0.3 * std::sinh(u));
    }
    const double transform_ns = seconds_since(t0) / m * 1e9;
    g_sink = s;
    std::printf("\nper call (ns): log_phi_pair %.1f, inverse_mills(x)+inverse_mills(-x) %.1f, erfcx %.1f, "
                "composite node transform (sinh, cosh, log) %.1f\n",
                pair_ns, mills_ns, erfcx_ns, transform_ns);
    return 0;
}
