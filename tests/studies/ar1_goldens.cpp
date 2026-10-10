// SPDX-License-Identifier: Apache-2.0
//
// Writes tests/golden/ar1/likelihood.csv (M6; D-195): the AR(1)-factor log-likelihood at fixed parameters on
// hand-made panels, from the engine (benchmark_models/ar1/filter.hpp, parity rule) and from ref's independent
// dense filter. Likelihood evaluations only; no fit. The panels are the hard cases for node placement (adjacent
// periods with very different modes, n up to 10^4), ordinary ones, and one whose n varies by period; T = 2 and 3
// rows are what validation/scipy/ar1_binomial_mle.py also integrates directly (nested, and over the joint density).
//
//   ar1_goldens --out FILE
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "benchmark_models/ar1/filter.hpp"
#include "core/quadrature/parity.hpp"
#include "ref/ref.hpp"

namespace {

namespace ar1 = vcal::benchmark_models::ar1;

struct Panel {
    std::vector<std::int64_t> n, d;
};

std::string hex(double v) {
    char b[64];
    std::snprintf(b, sizeof b, "%a", v);
    return b;
}

std::string join(const std::vector<std::int64_t>& v) {
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i) s += (i ? " " : "") + std::to_string(v[i]);
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    std::string out;
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::strcmp(argv[i], "--out") == 0) out = argv[i + 1];
    }
    if (out.empty()) {
        std::fprintf(stderr, "usage: ar1_goldens --out FILE\n");
        return 2;
    }
    const std::vector<Panel> panels = {
        {{20, 20}, {0, 3}},
        {{20, 20}, {3, 0}},
        {{100, 100, 100}, {0, 9, 0}},
        {{100, 100, 100}, {9, 0, 9}},
        {{1000, 1000}, {0, 90}},
        {{1000, 1000, 1000}, {90, 0, 90}},
        {{10000, 10000}, {0, 900}},
        {{10000, 10000, 10000}, {900, 0, 900}},
        {{100, 100, 100, 100, 100, 100}, {0, 3, 1, 7, 0, 2}},
        {{10000, 10000, 10000, 10000, 10000}, {50, 120, 80, 0, 300}},
        {{500, 800, 1000, 1200, 900, 700}, {3, 10, 25, 5, 0, 8}},
    };
    struct Theta {
        double pd, rho, phi;
    };
    const std::vector<Theta> thetas = {{0.01, 0.12, -0.6}, {0.01, 0.12, 0.0},  {0.01, 0.12, 0.6},
                                       {0.01, 0.12, 0.95}, {0.05, 0.3, 0.5},   {0.001, 0.05, 0.3}};
    std::ofstream f(out, std::ios::binary);
    if (!f) return 2;
    f << "# AR(1)-factor log-likelihood goldens (tests/studies/ar1_goldens.cpp; D-195): engine and ref, no fits\n";
    f << "panel,n,d,pd_hex,rho_hex,phi_hex,loglik_hex,ref_loglik_hex,placement_check_hex\n";
    const auto rule = vcal::quadrature::parity_rule();
    ar1::FilterWork w;
    for (std::size_t p = 0; p < panels.size(); ++p) {
        std::vector<ar1::Obs> obs;
        std::vector<vcalref::Period> q;
        for (std::size_t t = 0; t < panels[p].n.size(); ++t) {
            obs.push_back({panels[p].n[t], panels[p].d[t]});
            q.push_back({panels[p].n[t], panels[p].d[t]});
        }
        for (const auto& th : thetas) {
            const auto l = ar1::log_likelihood(obs.data(), static_cast<std::int64_t>(obs.size()), {th.pd, th.rho, th.phi},
                                               rule, w);
            const double r = vcalref::ar1_loglik(th.pd, th.rho, th.phi, q);
            f << p << ',' << join(panels[p].n) << ',' << join(panels[p].d) << ',' << hex(th.pd) << ',' << hex(th.rho)
              << ',' << hex(th.phi) << ',' << hex(l.value) << ',' << hex(r) << ',' << hex(l.placement) << '\n';
        }
    }
    std::fprintf(stderr, "wrote %s\n", out.c_str());
    return 0;
}
