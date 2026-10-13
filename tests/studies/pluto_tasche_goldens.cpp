// SPDX-License-Identifier: Apache-2.0
//
// Writes tests/golden/pluto_tasche/bounds.csv (S-6; D-215): Pluto and Tasche's most-prudent upper bounds at fixed
// inputs, from benchmark_models/pluto_tasche/bounds.hpp with the parity rule. Bounds only; no panels, no fits. The
// cases span the study's cells (n 100 and 1,000; T 20 and 100) at totals from 0 to well above the mean, at three
// confidences, for the independent bound and the one-factor bound at the study's three rho_a.
// unit_pluto_tasche checks the engine against them; validation/scipy/pluto_tasche.py replicates them.
//
//   pluto_tasche_goldens --out FILE
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

#include "benchmark_models/pluto_tasche/bounds.hpp"
#include "core/quadrature/parity.hpp"

namespace {

namespace pt = vcal::benchmark_models::pluto_tasche;

std::string hex(double v) {
    char b[64];
    std::snprintf(b, sizeof b, "%a", v);
    return b;
}

}  // namespace

int main(int argc, char** argv) {
    std::string out;
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::strcmp(argv[i], "--out") == 0) out = argv[i + 1];
    }
    if (out.empty()) {
        std::fprintf(stderr, "usage: pluto_tasche_goldens --out FILE\n");
        return 2;
    }
    const auto rule = vcal::quadrature::parity_rule();
    pt::FactorWork w;
    std::string text =
        "# Pluto-Tasche most-prudent upper bounds (tests/studies/pluto_tasche_goldens.cpp; D-215): "
        "benchmark_models/pluto_tasche/bounds.hpp, parity rule\n"
        "# rho_a is 0 for the independent bound; logit_bound is the bisection's midpoint, bound = logit^-1\n"
        "method,n,periods,d_total,rho_a_hex,gamma_hex,logit_bound_hex,bound\n";
    char buf[512];
    for (const std::int64_t n : {100, 1000}) {
        for (const std::int64_t T : {20, 100}) {
            for (const std::int64_t D : {0, 1, 5, 25}) {
                for (const double g : {0.5, 0.95, 0.99}) {
                    for (const double rho : {0.0, 0.02, 0.12, 0.24}) {
                        const pt::Bound b =
                            rho == 0.0 ? pt::independent_bound(n, T, D, g) : pt::factor_bound(n, T, D, rho, g, rule, w);
                        if (b.flags != 0) {
                            std::fprintf(stderr, "flagged bound: n %lld T %lld D %lld rho %g gamma %g\n",
                                         static_cast<long long>(n), static_cast<long long>(T), static_cast<long long>(D),
                                         rho, g);
                            return 1;
                        }
                        std::snprintf(buf, sizeof buf, "%s,%lld,%lld,%lld,%s,%s,%s,%.17g\n",
                                      rho == 0.0 ? "independent" : "factor", static_cast<long long>(n),
                                      static_cast<long long>(T), static_cast<long long>(D), hex(rho).c_str(),
                                      hex(g).c_str(), hex(b.u).c_str(), b.pd);
                        text += buf;
                    }
                }
            }
        }
    }
    std::ofstream f(out, std::ios::binary);
    f << text;
    f.close();
    if (!f) {
        std::fprintf(stderr, "cannot write %s\n", out.c_str());
        return 1;
    }
    return 0;
}
