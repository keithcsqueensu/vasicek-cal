// SPDX-License-Identifier: Apache-2.0
//
// The S-48 run (studies/hac-se/PREDICTION.md, H1-H6; D-209): S-35's panels, bit for bit (the recovery
// matrix under the AR(1)-factor DGP, dgp::simulate_panel_ar1, on the recovery seed, scenario ids and
// replicates; phi = 0 is the pinned recovery run's panels), each fitted exactly as S-35 fitted it
// (recovery::fit_counts with the profile arm), and at the estimate the HAC sandwich SE (engine/hac.hpp)
// at three lags: the rule's L = floor(4 (T/100)^(2/9)), L = 1 and 2 x the rule. Per fit and lag, the
// HAC Wald interval's coverage of PD and rho in logit units (u_hat +/- z se_u), and at the rule's lag the
// same interval with Kiefer and Vogelsang's fixed-b Bartlett critical value at b = L / T (reported
// without a prediction).
// Nothing is scored here: studies/hac-se/score.py scores the rows, after checking that the fits are S-35's.
//
//   study_hac_se [--phis 0,0.3,0.6] [--replicates R] [--scenarios ID,...] [--out FILE] [--threads N]
//                [--timing-out FILE]  (a runtime timing sample, perf/README.md, D-187)
//
// Doubles are written to 17 significant digits. Results do not depend on the thread count.
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

#include "engine/hac.hpp"
#include "tests/harness/dgp_provenance.hpp"
#include "tests/harness/timing_sample.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace e = vcal::engine;
using Clock = std::chrono::steady_clock;

enum Phase { kSimulate, kBinomial, kHac, kCacheFill, kPhases };
const char* const kPhaseNames[kPhases] = {"simulate", "binomial MLE", "HAC SE", "cache fill"};

struct Timer {
    double* acc;
    Clock::time_point t0 = Clock::now();
    ~Timer() { *acc += std::chrono::duration<double>(Clock::now() - t0).count(); }
};

// -1 not computed (or flagged, for Wald), 0 misses, 1 covers.
std::int8_t ind(bool computed, bool covers) { return static_cast<std::int8_t>(!computed ? -1 : covers ? 1 : 0); }

constexpr int kLags = 3;  // the rule, L = 1, 2 x the rule
const char* const kLagTags[kLags] = {"", "_l1", "_l2x"};

// Kiefer and Vogelsang (2005): the two-sided 95% fixed-b critical value for the Bartlett kernel.
double fixed_b_critical(double b) { return 1.9600 + 2.9694 * b + 0.4160 * b * b - 0.5324 * b * b * b; }

struct Hac {
    std::int32_t lag;
    std::uint32_t flags;
    double se_u[2];
    std::int8_t cover[2];  // -1 flagged fit (no iid Wald interval); a fit without a HAC SE misses
};

struct Row {
    double pd, rho, se_u[2];
    std::uint32_t flags;
    std::int8_t wald[2], prof[2];
    std::uint32_t prof_flags[2];
    Hac hac[kLags];
    double fixed_b_c;
    std::int8_t fixed_b[2];
};

Row fit_row(const rc::Scenario& s, std::uint32_t r, const std::vector<std::int64_t>& d, e::SurfaceRowCache& cache,
            double* sec) {
    Row row{};
    rc::Fit f{};
    {
        const vcal::test::FillMeter fills;
        {
            Timer t{&sec[kBinomial]};
            f = rc::fit_counts(s, r, d, &cache, rc::kArmProfile);
        }
        sec[kBinomial] -= fills.seconds();  // a fill is paid once per run, not per replicate (D-193)
        sec[kCacheFill] += fills.seconds();
    }
    Timer t{&sec[kHac]};
    const Grid<2> g = rc::grid();
    const double truth[2] = {s.pd, s.rho};
    const bool flagged = (f.flags & rc::kNoReliableInterval) != 0;
    row.pd = f.value[0];
    row.rho = f.value[1];
    row.flags = f.flags;
    double u_hat[2], u_true[2];
    for (int a = 0; a < 2; ++a) {
        const AxisScale sc = g.axis[a].scale;
        u_hat[a] = grid::to_scaled(sc, f.value[a]);
        u_true[a] = grid::to_scaled(sc, truth[a]);
        row.se_u[a] = f.se[a] / grid::dvalue_dscaled(sc, u_hat[a]);
        row.wald[a] = ind(!flagged, rc::covers(g, a, f, truth[a]));
        row.prof[a] = ind(!(f.prof_flags[a] & e::kIntervalNotComputed), rc::profile_covers(a, f, truth[a]));
        row.prof_flags[a] = f.prof_flags[a];
    }

    // The HAC SE needs only the estimate and its flags; the rest of Estimate2 is unused by hac_se.
    e::Estimate2 est{};
    est.value[0] = f.value[0];
    est.value[1] = f.value[1];
    est.flags = f.flags;
    std::vector<rc::Objective::Obs> obs(d.size());
    for (std::size_t k = 0; k < d.size(); ++k) obs[k] = {s.obligors, d[k]};
    static const auto primary = quadrature::parity_rule();
    const std::int32_t rule = e::newey_west_lag(s.periods);
    const std::int32_t lags[kLags] = {rule, 1, 2 * rule};
    for (int l = 0; l < kLags; ++l) {
        const e::HacSe2 h = e::hac_se(rc::Objective{}, primary, obs.data(), s.periods, g, est, lags[l]);
        Hac& out = row.hac[l];
        out.lag = h.lag;
        out.flags = h.flags;
        const bool have = !(h.flags & e::kIntervalNotComputed);
        for (int a = 0; a < 2; ++a) {
            out.se_u[a] = h.se_scaled[a];
            out.cover[a] = ind(!flagged, have && std::fabs(u_hat[a] - u_true[a]) <= rc::kZ975 * h.se_scaled[a]);
        }
        if (l == 0) {
            row.fixed_b_c = fixed_b_critical(static_cast<double>(rule) / static_cast<double>(s.periods));
            for (int a = 0; a < 2; ++a) {
                row.fixed_b[a] =
                    ind(!flagged, have && std::fabs(u_hat[a] - u_true[a]) <= row.fixed_b_c * h.se_scaled[a]);
            }
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
    if (!vcal::test::announce_dgp_kernel("study_hac_se")) return 2;
    std::vector<double> phis = {0.0, 0.3, 0.6};
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids;
    std::string out = "hac_se.csv", timing_out;
    unsigned threads = std::thread::hardware_concurrency();
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--phis") == 0 && i + 1 < argc) {
            phis = parse_list(argv[++i]);
        } else if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            for (const double v : parse_list(argv[++i])) ids.push_back(static_cast<std::uint32_t>(v));
        } else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out = argv[++i];
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--timing-out") == 0 && i + 1 < argc) {
            timing_out = argv[++i];
        } else {
            std::fprintf(stderr, "usage: study_hac_se [--phis 0,0.3,0.6] [--replicates R] [--scenarios ID,...] "
                                 "[--out FILE] [--threads N] [--timing-out FILE]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    if (threads == 0) threads = 1;
    const std::size_t P = phis.size();

    // Jobs: (scenario, replicate). A job fits its panel at every phi, as S-35's tool does.
    const std::size_t S = ids.size();
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
                    const dgp::PanelSpec spec{rc::kSeed, s.id, r, s.pd, s.rho, n.data(), s.periods};
                    if (dgp::simulate_panel_ar1(spec, phis[p], d.data(), nullptr) != dgp::Status::Ok) std::abort();
                }
                rows[job * P + p] = fit_row(s, r, d, cache, sec);
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    const double wall = std::chrono::duration<double>(Clock::now() - t0).count();

    char buf[1024];
    std::snprintf(buf, sizeof buf, "# S-48 (D-209): one row per (phi, scenario, replicate); seed 0x%016llX; dgp_kernel %s\n",
                  static_cast<unsigned long long>(rc::kSeed), vcal::test::dgp_kernel().c_str());
    std::string text = std::string(buf) +
                       "# coverage indicators: 1 covers, 0 misses, -1 not computed (Wald and HAC: flagged fit; a fit "
                       "without a HAC SE misses)\n"
                       "# se_*_u: SEs in logit units; hac_flags: engine::kIntervalNotComputed when no HAC SE\n"
                       "phi,scenario,replicate,pd_true,rho_true,periods,obligors,pd,rho,flags,se_pd_u,se_rho_u,"
                       "wald_pd,wald_rho,prof_pd,prof_rho,prof_flags_pd,prof_flags_rho";
    for (const char* tag : kLagTags) {
        text += std::string(",lag") + tag + ",hac_flags" + tag + ",hac_se_pd_u" + tag + ",hac_se_rho_u" + tag +
                ",hac_pd" + tag + ",hac_rho" + tag;
    }
    text += ",fixed_b_c,fixed_b_pd,fixed_b_rho\n";
    for (std::size_t p = 0; p < P; ++p) {
        for (std::size_t k = 0; k < S; ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const Row& w = rows[(k * R + r) * P + p];
                std::snprintf(buf, sizeof buf, "%.17g,%u,%u,%.17g,%.17g,%lld,%lld,%.17g,%.17g,%u,%.17g,%.17g,%d,%d,%d,%d,%u,%u",
                              phis[p], s.id, r, s.pd, s.rho, static_cast<long long>(s.periods),
                              static_cast<long long>(s.obligors), w.pd, w.rho, w.flags, w.se_u[0], w.se_u[1], w.wald[0],
                              w.wald[1], w.prof[0], w.prof[1], w.prof_flags[0], w.prof_flags[1]);
                text += buf;
                for (const Hac& h : w.hac) {
                    std::snprintf(buf, sizeof buf, ",%d,%u,%.17g,%.17g,%d,%d", h.lag, h.flags, h.se_u[0], h.se_u[1],
                                  h.cover[0], h.cover[1]);
                    text += buf;
                }
                std::snprintf(buf, sizeof buf, ",%.17g,%d,%d\n", w.fixed_b_c, w.fixed_b[0], w.fixed_b[1]);
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
    std::printf("%u replicates x %zu scenarios x %zu phi in %.0f s (%u threads; %zu cached surface rows); CPU s %.0f:",
                R, S, P, wall, threads, cache.size(), all);
    for (int ph = 0; ph < kPhases; ++ph) std::printf(" %s %.0f (%.0f%%);", kPhaseNames[ph], total[ph], 100.0 * total[ph] / all);
    std::printf("\n");
    if (!timing_out.empty()) {
        std::vector<std::string> phases(kPhaseNames, kPhaseNames + kPhases), shares(kPhases, "1");
        shares[kCacheFill] = "once";
        vcal::test::TimingSample sample("study_hac_se", phases, shares);
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
