// SPDX-License-Identifier: Apache-2.0
//
// S-35's M6 arm (studies/lrpd-serial-correlation/PREDICTION.md, the 2026-10-10 addendum, L16-L28; D-195,
// D-196): the AR(1)-aware binomial MLE (benchmark_models/ar1/) on S-35's panels: the recovery matrix under the
// AR(1)-factor DGP (dgp::simulate_panel_ar1) at each phi, on the recovery seed, scenario ids and replicates, the
// panels study_lrpd_serial fits. Per panel:
//   - the binomial-mixture MLE (calibrate on the recovery grid, with the shared row cache): the search's start at
//     phi = 0 and the reference of the first check before scoring;
//   - the M6 estimate (ar1::fit): PD, rho, phi, SEs, log-likelihood, the phi = 0 scan point, the doubled-rule and
//     placement checks, the scan's local maxima and the gap to the second, flags;
//   - its 3-D profile intervals for PD and rho, and for phi on the study subset (D-150) only;
//   - on the study subset, replicates 0-99, phi = 0.6 only: the parametric-bootstrap comparator, B = 199 panels
//     from the AR(1) DGP at the panel's own estimate, panel b of (s, r) on seed 0xB0071A61, scenario s,
//     replicate 256 r + b, each refitted for the point estimate; percentile intervals (type 7) for PD and rho.
// Nothing is scored here.
//
// --seed S runs the same design on another seed (the validation run, D-196/D-197: 0x7A11DA7E); the bootstrap
// comparator is S-35's only and is not run then, unless --comparator (with --boot-panels B, for smoke tests).
//
//   study_lrpd_serial_m6 [--phis 0,0.3,0.6] [--replicates R] [--scenarios ID,...] [--seed S] [--out FILE]
//                        [--threads N] [--timing-out FILE]  (a runtime timing sample, perf/README.md, D-187)
//                        [--comparator] [--boot-panels B]
//
// Doubles are written to 17 significant digits. Results do not depend on the thread count.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

#include "benchmark_models/ar1/search.hpp"
#include "resample/bootstrap.hpp"
#include "tests/harness/dgp_provenance.hpp"
#include "tests/harness/timing_sample.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace e = vcal::engine;
namespace ar1 = vcal::benchmark_models::ar1;
using Clock = std::chrono::steady_clock;

constexpr std::uint64_t kBootSeed = 0xB0071A61ull;          // the comparator's panels (PREDICTION.md)
constexpr std::uint32_t kBootPanels = 199;                  // B (--boot-panels: smoke tests only)
constexpr std::uint32_t kBootReplicates = 100;              // replicates 0-99
constexpr double kBootPhi = 0.6;                            // phi = 0.6 only
constexpr std::uint32_t kSubset[] = {29, 37, 72, 4, 68, 49, 7, 43, 51};  // the study subset (D-150)
const ar1::Box kBox{1e-4, 0.2, 1e-3, kDefaultRhoUpper};    // the recovery box (D-089, D-115)

enum Phase { kSimulate, kBinomial, kFit, kProfile, kProfilePhi, kBootstrap, kCacheFill, kPhases };
const char* const kPhaseNames[kPhases] = {"simulate", "binomial MLE", "M6 fit", "M6 profile PD rho",
                                          "M6 profile phi", "M6 bootstrap", "cache fill"};

struct Timer {
    double* acc;
    Clock::time_point t0 = Clock::now();
    ~Timer() { *acc += std::chrono::duration<double>(Clock::now() - t0).count(); }
};

bool in_subset(std::uint32_t id) { return std::find(std::begin(kSubset), std::end(kSubset), id) != std::end(kSubset); }

struct Row {
    double binom_pd, binom_rho, binom_loglik;
    std::uint32_t binom_flags;
    ar1::Estimate3 est;
    bool fitted;
    ar1::ProfileIntervals3 prof;
    double boot_lo[2], boot_hi[2];
    std::uint32_t boot_used, boot_excluded;  // comparator panels refitted with finite estimates / not
};

double nan() { return std::nan(""); }

std::vector<ar1::Obs> as_obs(const std::vector<std::int64_t>& d, std::int64_t n) {
    std::vector<ar1::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {n, d[t]};
    return obs;
}

// The binomial-mixture MLE on the recovery grid, its surface from the shared row cache (fills timed apart).
e::Estimate2 binomial(const rc::Scenario& s, const std::vector<std::int64_t>& d, e::SurfaceRowCache& cache, double* sec,
                      Phase phase) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    std::vector<rc::Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {s.obligors, d[t]};
    std::vector<double> L;
    e::Estimate2 est{};
    const vcal::test::FillMeter fills;
    {
        Timer t{&sec[phase]};
        if (e::calibrate_cached(backends::CpuBackend{1}, cache, rc::Objective{}, primary, check, obs.data(), s.periods,
                                rc::grid(), L, est) != e::Status::Ok) {
            std::abort();
        }
    }
    sec[phase] -= fills.seconds();  // a fill is paid once per run, not per replicate (D-193)
    sec[kCacheFill] += fills.seconds();
    return est;
}

Row fit_one(const rc::Scenario& s, std::uint32_t r, double phi, const std::vector<std::int64_t>& d, bool boot,
            std::uint32_t boot_panels, e::SurfaceRowCache& cache, double* sec) {
    Row row{};
    row.boot_lo[0] = row.boot_lo[1] = row.boot_hi[0] = row.boot_hi[1] = nan();
    for (int a = 0; a < 3; ++a) row.prof.lo[a] = row.prof.hi[a] = nan();
    for (auto& f : row.prof.flags) f = e::kIntervalNotComputed;
    const auto T = s.periods;
    {
        const e::Estimate2 est = binomial(s, d, cache, sec, kBinomial);
        row.binom_pd = est.value[0];
        row.binom_rho = est.value[1];
        row.binom_loglik = est.loglik;
        row.binom_flags = est.flags;
    }
    const auto obs = as_obs(d, s.obligors);
    {
        Timer t{&sec[kFit]};
        row.fitted = ar1::fit(obs.data(), T, kBox, row.binom_pd, row.binom_rho, row.est);
    }
    if (!row.fitted) return row;
    {
        Timer t{&sec[kProfile]};
        const bool which[3] = {true, true, false};
        ar1::profile_intervals(obs.data(), T, kBox, row.est, which, row.prof);
    }
    if (in_subset(s.id)) {
        Timer t{&sec[kProfilePhi]};
        const bool which[3] = {false, false, true};
        ar1::ProfileIntervals3 p{};
        ar1::profile_intervals(obs.data(), T, kBox, row.est, which, p);
        row.prof.lo[2] = p.lo[2];
        row.prof.hi[2] = p.hi[2];
        row.prof.flags[2] = p.flags[2];
        row.prof.residual_max = std::fmax(row.prof.residual_max, p.residual_max);
    }
    if (boot && in_subset(s.id) && r < kBootReplicates && phi == kBootPhi) {
        std::vector<double> pds, rhos;
        std::vector<std::int64_t> n(static_cast<std::size_t>(T), s.obligors), db(n.size());
        for (std::uint32_t b = 0; b < boot_panels; ++b) {
            {
                Timer t{&sec[kBootstrap]};
                const dgp::PanelSpec spec{kBootSeed, s.id, 256u * r + b, row.est.value[0], row.est.value[1], n.data(), T};
                if (dgp::simulate_panel_ar1(spec, row.est.value[2], db.data(), nullptr) != dgp::Status::Ok) std::abort();
            }
            // Each refit is the estimator exactly as on the panel itself: the scan starts at the bootstrap panel's
            // own binomial MLE. Point estimate only (PREDICTION.md).
            const e::Estimate2 bin = binomial(s, db, cache, sec, kBootstrap);
            const auto ob = as_obs(db, s.obligors);
            ar1::Estimate3 eb{};
            bool ok;
            {
                Timer t{&sec[kBootstrap]};
                ok = ar1::fit(ob.data(), T, kBox, bin.value[0], bin.value[1], eb);
            }
            if (ok && std::isfinite(eb.value[0]) && std::isfinite(eb.value[1])) {
                pds.push_back(eb.value[0]);
                rhos.push_back(eb.value[1]);
            } else {
                ++row.boot_excluded;
            }
        }
        row.boot_used = static_cast<std::uint32_t>(pds.size());
        Timer t{&sec[kBootstrap]};
        if (!pds.empty()) {
            std::sort(pds.begin(), pds.end());
            std::sort(rhos.begin(), rhos.end());
            row.boot_lo[0] = resample::quantile_type7(pds, 0.025);
            row.boot_hi[0] = resample::quantile_type7(pds, 0.975);
            row.boot_lo[1] = resample::quantile_type7(rhos, 0.025);
            row.boot_hi[1] = resample::quantile_type7(rhos, 0.975);
        }
    }
    return row;
}

std::vector<double> parse_list(const char* p) {
    std::vector<double> out;
    while (*p) {
        char* end = nullptr;
        out.push_back(std::strtod(p, &end));
        if (end == p) std::abort();
        p = *end == ',' ? end + 1 : end;
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (!vcal::test::announce_dgp_kernel("study_lrpd_serial_m6")) return 2;
    std::vector<double> phis = {0.0, 0.3, 0.6};
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids;
    std::uint64_t seed = rc::kSeed;
    std::string out = "lrpd_serial_m6.csv", timing_out;
    unsigned threads = std::thread::hardware_concurrency();
    bool force_comparator = false;
    std::uint32_t boot_panels = kBootPanels;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--phis") == 0 && i + 1 < argc) {
            phis = parse_list(argv[++i]);
        } else if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            for (const double v : parse_list(argv[++i])) ids.push_back(static_cast<std::uint32_t>(v));
        } else if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            seed = std::strtoull(argv[++i], nullptr, 0);
        } else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out = argv[++i];
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--timing-out") == 0 && i + 1 < argc) {
            timing_out = argv[++i];
        } else if (std::strcmp(argv[i], "--comparator") == 0) {
            force_comparator = true;
        } else if (std::strcmp(argv[i], "--boot-panels") == 0 && i + 1 < argc) {
            boot_panels = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else {
            std::fprintf(stderr, "usage: study_lrpd_serial_m6 [--phis 0,0.3,0.6] [--replicates R] [--scenarios ID,...] "
                                 "[--seed S] [--out FILE] [--threads N] [--timing-out FILE] [--comparator] [--boot-panels B]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    if (threads == 0) threads = 1;
    const bool boot = seed == rc::kSeed || force_comparator;  // the comparator belongs to S-35's run
    const std::size_t P = phis.size(), S = ids.size();
    std::vector<Row> rows(S * R * P);
    std::vector<double> seconds(S * R * kPhases, 0.0);
    e::SurfaceRowCache cache;
    std::atomic<std::size_t> next{0};
    const auto t0 = Clock::now();
    const auto worker = [&] {
        for (std::size_t job; (job = next.fetch_add(1)) < S * R;) {
            const std::size_t k = job / R;
            const auto r = static_cast<std::uint32_t>(job % R);
            const rc::Scenario s = rc::scenario(ids[k]);
            double* sec = &seconds[job * kPhases];
            for (std::size_t p = 0; p < P; ++p) {
                std::vector<std::int64_t> n(static_cast<std::size_t>(s.periods), s.obligors), d(n.size());
                {
                    Timer t{&sec[kSimulate]};
                    const dgp::PanelSpec spec{seed, s.id, r, s.pd, s.rho, n.data(), s.periods};
                    if (dgp::simulate_panel_ar1(spec, phis[p], d.data(), nullptr) != dgp::Status::Ok) std::abort();
                }
                rows[job * P + p] = fit_one(s, r, phis[p], d, boot, boot_panels, cache, sec);
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    const double wall = std::chrono::duration<double>(Clock::now() - t0).count();

    char buf[2048];
    std::snprintf(buf, sizeof buf, "# S-35 M6 arm (D-196): one row per (phi, scenario, replicate); seed 0x%llX; dgp_kernel %s\n",
                  static_cast<unsigned long long>(seed), vcal::test::dgp_kernel().c_str());
    std::string text = buf;
    text += "# profile flags: engine kInterval*; phi's interval on the study subset only; boot_*: the comparator "
            "(subset, replicates 0-99, phi = 0.6, S-35's seed only)\n"
            "phi,scenario,replicate,pd_true,rho_true,periods,obligors,binom_pd,binom_rho,binom_loglik,binom_flags,"
            "fitted,pd,rho,phi_hat,se_pd,se_rho,se_phi,loglik,phi0_pd,phi0_rho,phi0_loglik,doubled_check,"
            "placement_check,scan_local_maxima,scan_unconverged,second_maximum_gap,evaluations,scores,flags,"
            "prof_lo_pd,prof_hi_pd,prof_flags_pd,prof_lo_rho,prof_hi_rho,prof_flags_rho,prof_lo_phi,prof_hi_phi,"
            "prof_flags_phi,prof_residual,boot_lo_pd,boot_hi_pd,boot_lo_rho,boot_hi_rho,boot_used,boot_excluded\n";
    for (std::size_t p = 0; p < P; ++p) {
        for (std::size_t k = 0; k < S; ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const Row& w = rows[(k * R + r) * P + p];
                const auto& x = w.est;
                std::snprintf(buf, sizeof buf,
                              "%.17g,%u,%u,%.17g,%.17g,%lld,%lld,%.17g,%.17g,%.17g,%u,%d,%.17g,%.17g,%.17g,%.17g,%.17g,"
                              "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%lld,%lld,%u,%.17g,%.17g,%u,%.17g,"
                              "%.17g,%u,%.17g,%.17g,%u,%.17g,%.17g,%.17g,%.17g,%.17g,%u,%u\n",
                              phis[p], s.id, r, s.pd, s.rho, static_cast<long long>(s.periods),
                              static_cast<long long>(s.obligors), w.binom_pd, w.binom_rho, w.binom_loglik,
                              w.binom_flags, w.fitted ? 1 : 0, x.value[0], x.value[1], x.value[2], x.se[0], x.se[1],
                              x.se[2], x.loglik, x.phi0_value[0], x.phi0_value[1], x.phi0_loglik, x.doubled_check,
                              x.placement_check, x.scan_local_maxima, x.scan_unconverged, x.second_maximum_gap,
                              static_cast<long long>(x.evaluations), static_cast<long long>(x.scores), x.flags,
                              w.prof.lo[0], w.prof.hi[0], w.prof.flags[0], w.prof.lo[1], w.prof.hi[1], w.prof.flags[1],
                              w.prof.lo[2], w.prof.hi[2], w.prof.flags[2], w.prof.residual_max, w.boot_lo[0],
                              w.boot_hi[0], w.boot_lo[1], w.boot_hi[1], w.boot_used, w.boot_excluded);
                text += buf;
            }
        }
    }
    {
        std::ofstream f(out, std::ios::binary);
        f << text;
        f.close();
        if (!f) {
            std::fprintf(stderr, "cannot write %s\n", out.c_str());
            return 1;
        }
    }
    double total[kPhases] = {};
    for (std::size_t j = 0; j < S * R; ++j) {
        for (int ph = 0; ph < kPhases; ++ph) total[ph] += seconds[j * kPhases + static_cast<std::size_t>(ph)];
    }
    double all = 0.0;
    for (int ph = 0; ph < kPhases; ++ph) all += total[ph];
    std::printf("%u replicates x %zu scenarios x %zu phi in %.0f s (%u threads); CPU s %.0f:", R, S, P, wall, threads, all);
    for (int ph = 0; ph < kPhases; ++ph) std::printf(" %s %.0f (%.0f%%);", kPhaseNames[ph], total[ph], 100.0 * total[ph] / all);
    std::printf("\n");
    if (!timing_out.empty()) {
        // Shares (perf/README.md): every phase runs on every replicate, except the comparator, on replicates 0-99
        // of 1,000 (0.1). The phi interval runs on every replicate of the subset's scenarios only, which the
        // sample's per-cell means already carry (other scenarios record 0). The row cache's fills run once.
        std::vector<std::string> phases(kPhaseNames, kPhaseNames + kPhases), shares(kPhases, "1");
        shares[kBootstrap] = "0.1";
        shares[kCacheFill] = "once";
        vcal::test::TimingSample sample("study_lrpd_serial_m6", phases, shares);
        for (std::size_t k = 0; k < S; ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const double* sec = &seconds[(k * R + r) * kPhases];
                sample.record(s.id, s.periods, s.obligors, r, std::vector<double>(sec, sec + kPhases));
            }
        }
        if (!sample.write(timing_out)) {
            std::fprintf(stderr, "cannot write %s\n", timing_out.c_str());
            return 1;
        }
    }
    return 0;
}
