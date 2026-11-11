// SPDX-License-Identifier: Apache-2.0
//
// Writes tests/golden/multigrade/estimates.csv (M6; D-413): two small shared-factor panels and the multi-grade
// MLE's estimate on each (benchmark_models/multigrade/fit.hpp, parity rule, the default box). The counts are what
// validation/scipy/multigrade_mle.py maximises on its own (the estimator half of the replication); the estimates are
// what unit_multigrade_fit checks the engine against on every platform.
//
// Panels: per period one Z_t from dgp::factor_draw on scenario S << 6, then each grade's defaults given Z_t from
// dgp::count_defaults on scenario (S << 6) + k, seed 0x4D474749444E ("MGGOLD").
//
//   multigrade_goldens --out FILE
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "benchmark_models/multigrade/fit.hpp"
#include "dgp/dgp.hpp"

namespace {

namespace mg = vcal::benchmark_models::multigrade;

constexpr std::uint64_t kSeed = 0x4D474749444Eull;

std::string hex(double v) {
    char b[64];
    std::snprintf(b, sizeof b, "%a", v);
    return b;
}

std::string dec(double v) {
    char b[64];
    std::snprintf(b, sizeof b, "%.10g", v);
    return b;
}

struct Spec {
    const char* name;
    std::uint32_t scenario;
    std::vector<double> pd;
    std::vector<std::int64_t> n;
    double rho;
    std::int64_t T;
};

bool panel(const Spec& s, std::vector<mg::GradeCount>& y) {
    const auto kernel = vcal::dgp::BernoulliKernel::Scalar;  // the definition, on every platform
    for (std::int64_t t = 0; t < s.T; ++t) {
        const auto period = static_cast<std::uint32_t>(t);
        double z = 0.0;
        if (vcal::dgp::factor_draw(kSeed, s.scenario << 6, 0, period, vcal::dgp::kFactorBlocks, &z) !=
            vcal::dgp::Status::Ok) {
            return false;
        }
        for (std::size_t k = 0; k < s.pd.size(); ++k) {
            const double c = vcal::dgp::det_probit(s.pd[k]);
            const double p = vcal::dgp::det_ncdf((c - std::sqrt(s.rho) * z) / std::sqrt(1.0 - s.rho));
            const std::int64_t d = vcal::dgp::count_defaults(kSeed, (s.scenario << 6) + static_cast<std::uint32_t>(k),
                                                             0, period, s.n[k], p, kernel);
            if (d < 0) return false;
            y.push_back({s.n[k], d});
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3 || std::strcmp(argv[1], "--out") != 0) {
        std::fprintf(stderr, "usage: multigrade_goldens --out FILE\n");
        return 2;
    }
    const std::vector<Spec> specs = {
        {"two grades", 1, {0.01, 0.04}, {800, 400}, 0.12, 12},
        {"three grades", 2, {0.003, 0.02, 0.08}, {1500, 600, 200}, 0.2, 10},
    };
    std::ofstream out(argv[2], std::ios::binary);
    out << "# The multi-grade shared-rho MLE on two shared-factor panels (M6, D-413), written by\n"
           "# tests/studies/multigrade_goldens.cpp: counts as \"n:d;n:d\" per period, periods separated by \"|\"; the\n"
           "# estimate (PDs \";\"-separated), its log-likelihood and the 95% rho profile interval. Default box\n"
           "# PD [1e-5, 0.5], rho [1e-3, 0.5]. Doubles exact in *_hex; decimal columns for readers only.\n"
           "name,K,T,counts,pd_hex,rho_hex,loglik_hex,rho_lo_hex,rho_hi_hex,flags,pd,rho,loglik\n";
    for (const auto& s : specs) {
        std::vector<mg::GradeCount> y;
        if (!panel(s, y)) {
            std::fprintf(stderr, "panel %s: the DGP refused\n", s.name);
            return 1;
        }
        const int K = static_cast<int>(s.pd.size());
        mg::Fit f;
        if (mg::fit(y.data(), s.T, K, mg::Box{}, f) != mg::Status::Ok) {
            std::fprintf(stderr, "panel %s: the fit failed\n", s.name);
            return 1;
        }
        std::string counts, pd_hex, pd_dec;
        for (std::int64_t t = 0; t < s.T; ++t) {
            if (t > 0) counts += "|";
            for (int k = 0; k < K; ++k) {
                const auto& o = y[static_cast<std::size_t>(t * K + k)];
                if (k > 0) counts += ";";
                counts += std::to_string(o.n) + ":" + std::to_string(o.d);
            }
        }
        for (int k = 0; k < K; ++k) {
            pd_hex += (k > 0 ? ";" : "") + hex(f.pd[static_cast<std::size_t>(k)]);
            pd_dec += (k > 0 ? ";" : "") + dec(f.pd[static_cast<std::size_t>(k)]);
        }
        out << s.name << "," << K << "," << s.T << "," << counts << "," << pd_hex << "," << hex(f.rho) << ","
            << hex(f.loglik) << "," << hex(f.rho_lo) << "," << hex(f.rho_hi) << "," << f.flags << "," << pd_dec << ","
            << dec(f.rho) << "," << dec(f.loglik) << "\n";
    }
    return out.good() ? 0 : 1;
}
