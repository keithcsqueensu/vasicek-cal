// SPDX-License-Identifier: Apache-2.0
//
// The S-6 run (studies/pluto-tasche/PREDICTION.md, P1-P6; D-214, D-215): Pluto and Tasche's most-prudent upper
// bounds on the recovery grid's low-default cells (n PD <= 1: 27 scenarios), against the profile and Jeffreys upper
// limits on the same panels.
//
// Replicate r of scenario id is the base DGP panel (kSeed, id, r, PD, rho, n, T) of dgp/ (dgp::simulate_panel), n
// constant over the T periods; kSeed is S-6's own ("S06PLUTO", D-214), not the recovery seed. Per panel, with
// D = sum_t d_t, at each gamma in {0.5, 0.75, 0.9, 0.95, 0.99}:
//   - the bounds (benchmark_models/pluto_tasche/bounds.hpp, parity rule): independent; one factor at the
//     generating rho (oracle) and at rho_a = 0.12 (fixed). They depend on the panel only through D, so each
//     scenario's bounds are solved once per distinct total and reused;
//   - the profile upper limit: the binomial-mixture MLE's PD profile upper end at the threshold z_gamma^2 / 2
//     (engine::profile_intervals; D-131's conventions: a truncated end is still an end, a not-computed interval
//     does not cover). At gamma = 0.5 it is the polished MLE, the profile's own maximum;
//   - the Jeffreys upper limit: the gamma quantile of the Jeffreys posterior's PD marginal, resolution rule on
//     (engine::grid_posterior, as S-9's and S-10's tools; D-176). A refused or numeric posterior has none.
// The fit is the recovery fit's (evaluate_surface_cached, then calibrate_from_surface with the quadrature check).
// Nothing is scored here: studies/pluto-tasche/score.py scores the rows against reference.csv.
//
//   study_pluto_tasche [--replicates R] [--scenarios ID,...] [--out FILE] [--threads N] [--timing-out FILE]
//
// One row per (scenario, replicate, gamma). Doubles are written to 17 significant digits. Results do not depend on
// the thread count.
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "benchmark_models/pluto_tasche/bounds.hpp"
#include "engine/posterior.hpp"
#include "engine/profile.hpp"
#include "tests/harness/dgp_provenance.hpp"
#include "tests/harness/timing_sample.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace e = vcal::engine;
namespace pt = vcal::benchmark_models::pluto_tasche;
using Clock = std::chrono::steady_clock;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

inline constexpr std::uint64_t kSeed = 0x533036504C55544Full;  // "S06PLUTO" (D-214)
inline constexpr double kRhoFixed = 0.12;
constexpr int kGammas = 5;
inline constexpr double kGamma[kGammas] = {0.5, 0.75, 0.9, 0.95, 0.99};
// Phi^-1(gamma) (scipy.special.ndtri, 17 digits); the one-sided profile threshold is z^2 / 2. gamma = 0.5 has none.
inline constexpr double kZ[kGammas] = {0.0, 0.67448975019608171, 1.2815515655446004, 1.6448536269514722,
                                       2.3263478740408408};
constexpr int kMethods = 3;  // independent, oracle, fixed
const char* const kMethodNames[kMethods] = {"independent", "oracle", "fixed"};

enum Phase { kSimulate, kFit, kProfile, kJeffreys, kBounds, kCacheFill, kPhases };
const char* const kPhaseNames[kPhases] = {"simulate", "fit", "profile", "Jeffreys", "bounds", "cache fill"};

struct Timer {
    double* acc;
    Clock::time_point t0 = Clock::now();
    ~Timer() { *acc += std::chrono::duration<double>(Clock::now() - t0).count(); }
};

// The study's cells: n PD <= 1 (the slack covers rounding in n PD).
bool low_default(const rc::Scenario& s) { return static_cast<double>(s.obligors) * s.pd <= 1.0 + 1e-9; }

struct Bounds {
    double u[kGammas][kMethods];  // logit
    double pd[kGammas][kMethods];
};

// Each scenario's bounds, per distinct total, solved on first use. A bound is a pure function of (n, T, D, rho_a,
// gamma), so which job solves it first cannot change it.
class BoundCache {
public:
    const Bounds& get(const rc::Scenario& s, std::int64_t D, double* sec) {
        {
            const std::lock_guard<std::mutex> lock(mutex_);
            const auto it = cache_.find({s.id, D});
            if (it != cache_.end()) return it->second;
        }
        Timer t{sec};
        static const auto rule = quadrature::parity_rule();
        pt::FactorWork w;
        Bounds b{};
        for (int g = 0; g < kGammas; ++g) {
            const pt::Bound got[kMethods] = {pt::independent_bound(s.obligors, s.periods, D, kGamma[g]),
                                             pt::factor_bound(s.obligors, s.periods, D, s.rho, kGamma[g], rule, w),
                                             pt::factor_bound(s.obligors, s.periods, D, kRhoFixed, kGamma[g], rule, w)};
            for (int m = 0; m < kMethods; ++m) {
                if (got[m].flags != 0) {
                    std::fprintf(stderr, "flagged bound (%u): scenario %u, D %lld, gamma %g, %s\n", got[m].flags, s.id,
                                 static_cast<long long>(D), kGamma[g], kMethodNames[m]);
                    std::abort();
                }
                b.u[g][m] = got[m].u;
                b.pd[g][m] = got[m].pd;
            }
        }
        const std::lock_guard<std::mutex> lock(mutex_);
        return cache_.emplace(std::make_pair(s.id, D), b).first->second;
    }
    std::size_t size() const {
        const std::lock_guard<std::mutex> lock(mutex_);
        return cache_.size();
    }

private:
    mutable std::mutex mutex_;
    std::map<std::pair<std::uint32_t, std::int64_t>, Bounds> cache_;  // node-based: references stay valid
};

struct Row {
    std::int64_t d_total;
    std::uint32_t fit_flags;
    double bound[kGammas][kMethods];
    double profile[kGammas];
    std::uint32_t profile_flags[kGammas];
    double jeffreys[kGammas];
    std::uint32_t jeffreys_flags;
};

struct Context {
    e::SurfaceRowCache cache;
    std::map<std::int64_t, e::JeffreysTable> jeffreys;
    BoundCache bounds;
};

Row run_panel(Context& ctx, const rc::Scenario& s, std::uint32_t r, double* sec) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    const std::int64_t T = s.periods;
    Row row{};
    std::vector<std::int64_t> n(static_cast<std::size_t>(T), s.obligors), d(n.size());
    {
        Timer t{&sec[kSimulate]};
        const dgp::PanelSpec spec{kSeed, s.id, r, s.pd, s.rho, n.data(), T};
        if (dgp::simulate_panel(spec, d.data(), nullptr) != dgp::Status::Ok) std::abort();
    }
    for (const std::int64_t x : d) row.d_total += x;

    const Bounds& b = ctx.bounds.get(s, row.d_total, &sec[kBounds]);
    for (int k = 0; k < kGammas; ++k) {
        for (int m = 0; m < kMethods; ++m) row.bound[k][m] = b.pd[k][m];
    }

    // The fit, as the recovery fit (calibrate_cached in two calls; identical).
    std::vector<rc::Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {s.obligors, d[t]};
    std::vector<double> L;
    e::Estimate2 est{};
    {
        const vcal::test::FillMeter fills;
        {
            Timer t{&sec[kFit]};
            L.assign(static_cast<std::size_t>(T * g.size()), 0.0);
            e::evaluate_surface_cached(backends::CpuBackend{1}, ctx.cache, rc::Objective{}, primary, obs.data(), T, g,
                                       L.data());
            if (e::calibrate_from_surface(backends::CpuBackend{1}, rc::Objective{}, primary, check, obs.data(), T, g, L,
                                          est, e::kHessianStepFraction, true) != e::Status::Ok) {
                std::abort();
            }
        }
        sec[kFit] -= fills.seconds();  // a fill is paid once per run, not per replicate (D-193)
        sec[kCacheFill] += fills.seconds();
    }
    row.fit_flags = est.flags;

    {
        Timer t{&sec[kProfile]};
        double polished = kNaN;
        std::uint32_t polished_flags = e::kIntervalNotComputed;
        for (int k = 1; k < kGammas; ++k) {
            const auto prof =
                e::profile_intervals(rc::Objective{}, primary, obs.data(), T, g, L, est, 0.5 * kZ[k] * kZ[k]);
            const bool computed = !(prof.flags[0] & e::kIntervalNotComputed);
            row.profile[k] = computed ? prof.hi[0] : kNaN;
            row.profile_flags[k] = prof.flags[0];
            if (k == 1) {
                polished = computed ? prof.max_at[0] : kNaN;
                polished_flags = prof.flags[0] & e::kIntervalNotComputed;
            }
        }
        row.profile[0] = polished;
        row.profile_flags[0] = polished_flags;
    }

    {
        Timer t{&sec[kJeffreys]};
        double se_u[2];
        for (int a = 0; a < 2; ++a) {
            se_u[a] = est.se[a] / grid::dvalue_dscaled(g.axis[a].scale, grid::to_scaled(g.axis[a].scale, est.value[a]));
        }
        const auto p = e::grid_posterior(backends::CpuBackend{1}, rc::Objective{}, primary, obs.data(), T, g, L,
                                         e::Prior::Jeffreys, &ctx.jeffreys.at(s.obligors), se_u, true);
        row.jeffreys_flags = p.flags;
        const bool have = !(p.flags & (e::kPosteriorRefused | e::kPosteriorNumeric));
        for (int k = 0; k < kGammas; ++k) {
            row.jeffreys[k] =
                have ? grid::from_scaled(p.grid.axis[0].scale, e::posterior_detail::quantile(p.marginal[0], kGamma[k]))
                     : kNaN;
        }
    }
    return row;
}

std::vector<std::uint32_t> parse_ids(const char* p) {
    std::vector<std::uint32_t> out;
    while (*p) {
        char* end = nullptr;
        out.push_back(static_cast<std::uint32_t>(std::strtoul(p, &end, 10)));
        if (end == p) std::abort();
        p = *end == ',' ? end + 1 : end;
    }
    return out;
}

// 1 covers (the limit is at or above the true PD), 0 does not; a missing limit does not cover (D-131).
int covers(double limit, double truth) { return limit >= truth ? 1 : 0; }

// 17 significant digits; NaN always as "nan" (MSVC prints a computed NaN as "-nan(ind)").
std::string num(double v) {
    if (std::isnan(v)) return "nan";
    char b[64];
    std::snprintf(b, sizeof b, "%.17g", v);
    return b;
}

}  // namespace

int main(int argc, char** argv) {
    if (!vcal::test::announce_dgp_kernel("study_pluto_tasche")) return 2;
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids;
    std::string out = "pluto_tasche.csv", timing_out;
    unsigned threads = std::thread::hardware_concurrency();
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            ids = parse_ids(argv[++i]);
        } else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out = argv[++i];
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--timing-out") == 0 && i + 1 < argc) {
            timing_out = argv[++i];
        } else {
            std::fprintf(stderr, "usage: study_pluto_tasche [--replicates R] [--scenarios ID,...] [--out FILE] "
                                 "[--threads N] [--timing-out FILE]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
            if (low_default(rc::scenario(id))) ids.push_back(id);
        }
    }
    for (const auto id : ids) {
        if (id >= rc::kScenarios || !low_default(rc::scenario(id))) {
            std::fprintf(stderr, "scenario %u is not one of the study's low-default cells (n PD <= 1)\n", id);
            return 2;
        }
    }
    if (threads == 0) threads = 1;
    const auto t0 = Clock::now();

    Context ctx;
    for (const auto id : ids) {
        const auto n = rc::scenario(id).obligors;
        if (ctx.jeffreys.count(n) == 0) {
            ctx.jeffreys.emplace(n, e::jeffreys_table(backends::CpuBackend{static_cast<int>(threads)}, rc::Objective{},
                                                      quadrature::parity_rule(), n, rc::grid()));
        }
    }
    const double tables = std::chrono::duration<double>(Clock::now() - t0).count();

    const std::size_t S = ids.size();
    std::vector<Row> rows(S * R);
    std::vector<double> seconds(S * R * kPhases, 0.0);
    std::atomic<std::size_t> next{0};
    const auto worker = [&] {
        for (std::size_t job; (job = next.fetch_add(1)) < S * R;) {
            const rc::Scenario s = rc::scenario(ids[job / R]);
            const auto r = static_cast<std::uint32_t>(job % R);
            rows[job] = run_panel(ctx, s, r, &seconds[job * kPhases]);
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    const double wall = std::chrono::duration<double>(Clock::now() - t0).count();

    char buf[2048];
    std::snprintf(buf, sizeof buf,
                  "# S-6 (D-214, D-215): one row per (scenario, replicate, gamma); seed 0x%016llX; dgp_kernel %s\n",
                  static_cast<unsigned long long>(kSeed), vcal::test::dgp_kernel().c_str());
    std::string text = std::string(buf) +
                       "# bounds: independent, oracle (rho_a = rho), fixed (rho_a = 0.12); profile: the PD profile's "
                       "upper end at z_gamma^2 / 2 (the polished MLE at gamma = 0.5); jeffreys: the posterior PD "
                       "marginal's gamma quantile\n"
                       "# cov_*: 1 if the limit is >= pd_true, else 0 (a missing limit, NaN, does not cover); flags: "
                       "engine bits (fit: calibrate; profile: kInterval*; jeffreys: kPosterior*)\n"
                       "scenario,replicate,pd_true,rho_true,periods,obligors,gamma,d_total,independent,oracle,fixed,"
                       "profile,jeffreys,cov_independent,cov_oracle,cov_fixed,cov_profile,cov_jeffreys,fit_flags,"
                       "profile_flags,jeffreys_flags\n";
    for (std::size_t k = 0; k < S; ++k) {
        const rc::Scenario s = rc::scenario(ids[k]);
        for (std::uint32_t r = 0; r < R; ++r) {
            const Row& w = rows[k * R + r];
            for (int g = 0; g < kGammas; ++g) {
                const double* bd = w.bound[g];
                std::snprintf(buf, sizeof buf, "%u,%u,%s,%s,%lld,%lld,%s,%lld,%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%u,%u,%u\n",
                              s.id, r, num(s.pd).c_str(), num(s.rho).c_str(), static_cast<long long>(s.periods),
                              static_cast<long long>(s.obligors), num(kGamma[g]).c_str(),
                              static_cast<long long>(w.d_total), num(bd[0]).c_str(), num(bd[1]).c_str(),
                              num(bd[2]).c_str(), num(w.profile[g]).c_str(), num(w.jeffreys[g]).c_str(),
                              covers(bd[0], s.pd), covers(bd[1], s.pd), covers(bd[2], s.pd), covers(w.profile[g], s.pd),
                              covers(w.jeffreys[g], s.pd), w.fit_flags, w.profile_flags[g], w.jeffreys_flags);
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
    std::printf("%u replicates x %zu scenarios in %.0f s (%u threads; Jeffreys tables %.0f s; %zu cached surface rows; "
                "%zu distinct totals); CPU s %.0f:",
                R, S, wall, threads, tables, ctx.cache.size(), ctx.bounds.size(), all);
    for (int ph = 0; ph < kPhases; ++ph) std::printf(" %s %.0f (%.0f%%);", kPhaseNames[ph], total[ph], 100.0 * total[ph] / all);
    std::printf("\n");
    if (!timing_out.empty()) {
        // The bounds are solved once per distinct total, so, like the cache fills, they are a one-time cost.
        std::vector<std::string> phases(kPhaseNames, kPhaseNames + kPhases), shares(kPhases, "1");
        shares[kBounds] = "once";
        shares[kCacheFill] = "once";
        vcal::test::TimingSample sample("study_pluto_tasche", phases, shares);
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
