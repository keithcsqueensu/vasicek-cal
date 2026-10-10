// SPDX-License-Identifier: Apache-2.0
//
// S-35's M6 arm, model-based reference numbers (owner, 2026-10-10; D-196): the expected Fisher information of
// the AR(1)-aware binomial likelihood (benchmark_models/ar1/) at the TRUE parameters of each recovery scenario,
// for phi in {0, 0.3, 0.6}. Like S-35's reference.py it uses the model only: no estimate is ever computed. Per
// (scenario, phi), M panels are drawn from the AR(1) DGP and the score s(theta_true) is evaluated on each; the
// information is the mean of s s' (E[s] = 0 at the truth; its sample mean is reported as a check). Only scores
// and information are computed and stored: no fit, no maximisation, no phi-hat.
//
// Seed 0x1F020261010 (written 0x1F0_2026_1010): fresh, disjoint from S-35's recovery seed (kSeed), from the
// validation run's reserved seed 0x7A11DA7E and from the scratch seed 0xA11CE5EED (scenarios 900-901).
//
// What it is for (owner's condition): marking where asymptotics can be trusted, not predicting the spread
// where they cannot. From the inverse information, per-panel asymptotic SEs of PD, rho and phi. Stated before
// this tool was first run:
//   - phi: credible where the information-based SE is at most 0.15 (so the Wald half-width, 0.29, fits inside
//     the box from phi = 0.6, 0.35 from the bound);
//   - PD and rho: credible where 1.96 times the information-based SE in logit units is at most the distance, in
//     logit units, from the true value to the nearer bound of the estimation box (PD [1e-4, 0.2], rho
//     [1e-3, 0.5]).
// Where a parameter is not credible, the registration's predictions for it are qualitative.
//
//   study_ar1_information --out FILE [--panels M] [--threads N]
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "benchmark_models/ar1/filter.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "tests/harness/dgp_provenance.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

namespace ar1 = vcal::benchmark_models::ar1;
namespace rc = vcal::recovery;

constexpr std::uint64_t kInformationSeed = 0x1F020261010ull;
constexpr double kPhis[] = {0.0, 0.3, 0.6};
constexpr double kPhiSeCredible = 0.15;
constexpr double kBox[2][2] = {{1e-4, 0.2}, {1e-3, 0.5}};

struct Row {
    double mean_s[3];
    double info[3][3];   // per panel
    double se[3];        // natural scale, from the inverse information; NaN if not positive definite
    double se_logit[2];  // PD and rho, logit units
    double mc_rel_phi;   // Monte Carlo relative SE of the information's phi-phi entry
    bool credible[3];
    std::int64_t nonfinite;  // panels whose likelihood or score was not finite (excluded)
};

double logit(double p) { return std::log(p / (1.0 - p)); }

}  // namespace

int main(int argc, char** argv) {
    if (!vcal::test::announce_dgp_kernel("study_ar1_information")) return 2;
    std::string out_path;
    std::int64_t M = 500;
    int threads = 1;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (std::strcmp(argv[i], "--panels") == 0 && i + 1 < argc) {
            M = std::atoll(argv[++i]);
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = std::atoi(argv[++i]);
        } else {
            std::fprintf(stderr, "usage: study_ar1_information --out FILE [--panels M] [--threads N]\n");
            return 2;
        }
    }
    if (out_path.empty() || M < 2) {
        std::fprintf(stderr, "study_ar1_information: --out is required and --panels must be at least 2\n");
        return 2;
    }
    constexpr std::size_t kPhiCount = sizeof(kPhis) / sizeof(kPhis[0]);
    const std::int64_t jobs = static_cast<std::int64_t>(rc::kScenarios * kPhiCount);
    std::vector<Row> rows(static_cast<std::size_t>(jobs));
    vcal::backends::CpuBackend{threads}.parallel_for(jobs, [&](std::int64_t job) {
        const auto id = static_cast<std::uint32_t>(job / static_cast<std::int64_t>(kPhiCount));
        const double phi = kPhis[static_cast<std::size_t>(job) % kPhiCount];
        const auto s = rc::scenario(id);
        const auto rule = vcal::quadrature::parity_rule();
        ar1::FilterWork work;
        std::vector<std::int64_t> n(static_cast<std::size_t>(s.periods), s.obligors), d(n.size());
        std::vector<ar1::Obs> obs(n.size());
        double sum[3] = {0.0, 0.0, 0.0}, sum2[3][3] = {}, sum_pp2 = 0.0;
        std::int64_t used = 0, nonfinite = 0;
        for (std::int64_t m = 0; m < M; ++m) {
            const vcal::dgp::PanelSpec spec{kInformationSeed, id, static_cast<std::uint32_t>(m), s.pd, s.rho, n.data(),
                                            s.periods};
            if (vcal::dgp::simulate_panel_ar1(spec, phi, d.data(), nullptr) != vcal::dgp::Status::Ok) std::abort();
            for (std::size_t t = 0; t < n.size(); ++t) obs[t] = {s.obligors, d[t]};
            const ar1::Theta th{s.pd, s.rho, phi};
            const double l = ar1::log_likelihood(obs.data(), s.periods, th, rule, work).value;
            double g[3];
            if (std::isfinite(l)) ar1::score(obs.data(), s.periods, th, work, g);
            if (!std::isfinite(l) || !std::isfinite(g[0]) || !std::isfinite(g[1]) || !std::isfinite(g[2])) {
                ++nonfinite;
                continue;
            }
            ++used;
            for (int i = 0; i < 3; ++i) {
                sum[i] += g[i];
                for (int j = 0; j < 3; ++j) sum2[i][j] += g[i] * g[j];
            }
            sum_pp2 += (g[2] * g[2]) * (g[2] * g[2]);
        }
        Row& r = rows[static_cast<std::size_t>(job)];
        r.nonfinite = nonfinite;
        const double u = static_cast<double>(used);
        for (int i = 0; i < 3; ++i) {
            r.mean_s[i] = sum[i] / u;
            for (int j = 0; j < 3; ++j) r.info[i][j] = sum2[i][j] / u;
        }
        const double ipp = r.info[2][2];
        r.mc_rel_phi = std::sqrt(std::fmax(sum_pp2 / u - ipp * ipp, 0.0) / u) / ipp;
        // Inverse of the 3 x 3 information by cofactors; SEs where it is positive definite.
        const auto& A = r.info;
        const double c00 = A[1][1] * A[2][2] - A[1][2] * A[2][1];
        const double c11 = A[0][0] * A[2][2] - A[0][2] * A[2][0];
        const double c22 = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        const double det = A[0][0] * c00 - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                           A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
        const bool pos_def = A[0][0] > 0.0 && c22 > 0.0 && det > 0.0;
        const double inv_diag[3] = {c00 / det, c11 / det, c22 / det};
        const double truth[2] = {s.pd, s.rho};
        for (int i = 0; i < 3; ++i) r.se[i] = pos_def ? std::sqrt(inv_diag[i]) : static_cast<double>(NAN);
        for (int i = 0; i < 2; ++i) {
            r.se_logit[i] = r.se[i] / (truth[i] * (1.0 - truth[i]));
            const double dist = std::fmin(logit(truth[i]) - logit(kBox[i][0]), logit(kBox[i][1]) - logit(truth[i]));
            r.credible[i] = pos_def && 1.96 * r.se_logit[i] <= dist;
        }
        r.credible[2] = pos_def && r.se[2] <= kPhiSeCredible;
    });

    std::ofstream f(out_path, std::ios::binary);
    if (!f) {
        std::fprintf(stderr, "study_ar1_information: cannot write %s\n", out_path.c_str());
        return 2;
    }
    f << "# S-35 M6 arm: expected Fisher information at the true parameters (study_ar1_information; D-196)\n"
      << "# seed 0x1F020261010, panels per (scenario, phi) " << M << ", dgp_kernel " << vcal::test::dgp_kernel() << "\n"
      << "scenario,pd,rho,periods,obligors,phi,panels_used,nonfinite,mean_s_pd,mean_s_rho,mean_s_phi,"
         "i_pd_pd,i_rho_rho,i_phi_phi,i_pd_rho,i_pd_phi,i_rho_phi,se_pd,se_rho,se_phi,se_logit_pd,se_logit_rho,"
         "mc_rel_i_phi_phi,credible_pd,credible_rho,credible_phi\n";
    char buf[1024];
    for (std::int64_t job = 0; job < jobs; ++job) {
        const auto id = static_cast<std::uint32_t>(job / static_cast<std::int64_t>(kPhiCount));
        const double phi = kPhis[static_cast<std::size_t>(job) % kPhiCount];
        const auto s = rc::scenario(id);
        const Row& r = rows[static_cast<std::size_t>(job)];
        std::snprintf(buf, sizeof buf,
                      "%u,%.17g,%.17g,%lld,%lld,%.1f,%lld,%lld,%.6e,%.6e,%.6e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.6e,"
                      "%.6e,%.6e,%.6e,%.6e,%.4f,%d,%d,%d\n",
                      id, s.pd, s.rho, static_cast<long long>(s.periods), static_cast<long long>(s.obligors), phi,
                      static_cast<long long>(M - r.nonfinite), static_cast<long long>(r.nonfinite), r.mean_s[0],
                      r.mean_s[1], r.mean_s[2], r.info[0][0], r.info[1][1], r.info[2][2], r.info[0][1], r.info[0][2],
                      r.info[1][2], r.se[0], r.se[1], r.se[2], r.se_logit[0], r.se_logit[1], r.mc_rel_phi,
                      r.credible[0] ? 1 : 0, r.credible[1] ? 1 : 0, r.credible[2] ? 1 : 0);
        f << buf;
    }
    std::fprintf(stderr, "wrote %s\n", out_path.c_str());
    return 0;
}
