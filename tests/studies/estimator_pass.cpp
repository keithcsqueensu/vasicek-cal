// SPDX-License-Identifier: Apache-2.0
//
// The estimator-comparison pass: studies S-8, S-27 and S-28 (studies/mle-vs-mom/PREDICTION.md,
// E1-E13; D-150, D-161). Every panel is fitted by every estimator of the pass, and one row per
// (scenario, replicate) is written with each estimator's raw result. Nothing is scored here: the
// predictions are scored by studies/mle-vs-mom/compare.py from these rows.
//
// Panels: the recovery panels (scenarios 0-80, n <= 10^4), and for S-27 the 54 scenarios that extend
// n to 10^5 and 10^6 with the same PD, rho and T: id 81 + 2 (9 i_PD + 3 i_rho + i_T) + (i_n - 3),
// under the recovery seed and DGP (D-109). A panel is simulated once and shared by the estimators.
//
// Estimators, as registered:
//   - binomial MLE (the benchmark): recovery::fit with the profile arm, as pinned (no bootstrap);
//   - method of moments on counts (S-8): engine::mom_from_counts on the recovery box;
//   - Vasicek-rate MLE (S-27, S-28) on the rates d/n with detection limit 1/(2n), four treatments of
//     periods with d = 0 or d = n: refuse (parity), drop (the periods removed; refused if fewer than 3
//     remain), censor and substitute. The estimate is the polished maximum (the profile code's
//     max_at, "maximised over the box off the grid"); the 95% interval is the profile-likelihood
//     interval. A panel without such periods gets one fit, shared by all four treatments (they are
//     the same likelihood there);
//   - method of moments on rates (D-046), with zero rates as zeros, and with the drop treatment.
// Coverage of the true PD and rho by a profile interval: 1 covers, 0 misses, -1 not estimated;
// an interval that was not computed does not cover (D-131).
//
//   study_estimator_pass [--replicates R] [--scenarios ID,...] [--out FILE] [--threads N]
//
// --out writes the rows (CSV, doubles to 17 significant digits, an exact round trip) for conversion
// to a committed Parquet file. Results do not depend on the thread count.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "core/objectives/vasicek_rate.hpp"
#include "engine/moments.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace o = vcal::objectives;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

inline constexpr std::uint32_t kPassScenarios = 135;
inline constexpr std::int64_t kLargeN[2] = {100000, 1000000};

// Scenario ids 0-80 are the recovery matrix; 81-134 the registered n = 10^5 and 10^6 extension.
rc::Scenario pass_scenario(std::uint32_t id) {
    if (id < rc::kScenarios) return rc::scenario(id);
    const std::uint32_t j = id - rc::kScenarios, i_n = j % 2, cell = j / 2;
    const std::uint32_t i_t = cell % 3, i_rho = (cell / 3) % 3, i_pd = cell / 9;
    return {id, rc::kPds[i_pd], rc::kRhos[i_rho], rc::kPeriods[i_t], kLargeN[i_n]};
}

enum RateStatus : std::int8_t { kEstimated = 0, kRefused = 1, kEngineError = 2 };

struct RateResult {
    std::int8_t status = kRefused;
    double pd = kNaN, rho = kNaN;
    std::uint32_t flags = 0, prof_flags[2] = {0, 0};
    std::int8_t cover[2] = {-1, -1};
};

struct MomResult {
    double pd = kNaN, rho = kNaN;
    std::uint32_t flags = 0;
};

struct Row {
    std::int32_t zero_periods = 0, full_periods = 0;  // periods with d = 0, with d = n
    // binomial MLE
    double b_pd, b_rho, b_se_pd, b_se_rho;
    std::uint32_t b_flags, b_prof_flags[2];
    std::int8_t b_cover[2];
    MomResult mom;
    RateResult rate[4];  // refuse, drop, censor, substitute
    MomResult mom_rate, mom_rate_drop;
};

const char* const kTreatments[4] = {"refuse", "drop", "censor", "substitute"};

std::int8_t profile_cover(const engine::ProfileIntervals2& p, int a, double truth) {
    if (p.flags[a] & engine::kIntervalNotComputed) return 0;
    return static_cast<std::int8_t>(p.lo[a] <= truth && truth <= p.hi[a]);
}

template <class O>
RateResult rate_fit(const rc::Scenario& s, const std::vector<o::RateObs>& obs) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    const auto T = static_cast<std::int64_t>(obs.size());
    std::vector<double> L;
    engine::Estimate2 est{};
    RateResult out;
    if (engine::calibrate(backends::CpuBackend{1}, O{}, primary, check, obs.data(), T, g, L, est) != engine::Status::Ok) {
        out.status = kEngineError;
        return out;
    }
    const auto prof = engine::profile_intervals(O{}, primary, obs.data(), T, g, L, est);
    out.status = kEstimated;
    out.pd = prof.max_at[0];
    out.rho = prof.max_at[1];
    out.flags = est.flags;
    out.prof_flags[0] = prof.flags[0];
    out.prof_flags[1] = prof.flags[1];
    out.cover[0] = profile_cover(prof, 0, s.pd);
    out.cover[1] = profile_cover(prof, 1, s.rho);
    return out;
}

MomResult mom(const engine::MomEstimate& m) { return {m.pd, m.rho, m.flags}; }

Row pass_fit(const rc::Scenario& s, std::uint32_t replicate, engine::SurfaceRowCache* cache) {
    using Refuse = o::VasicekRate<PrecisionF64, o::ZeroRates::Refuse>;
    using Censor = o::VasicekRate<PrecisionF64, o::ZeroRates::Censor>;
    using Substitute = o::VasicekRate<PrecisionF64, o::ZeroRates::Substitute>;
    static const auto primary = quadrature::parity_rule();
    const Grid<2> g = rc::grid();
    Row row{};

    // The binomial MLE exactly as the recovery harness fits it (it simulates the same panel).
    const rc::Fit f = rc::fit(s, replicate, cache, rc::kArmProfile);
    row.b_pd = f.value[0];
    row.b_rho = f.value[1];
    row.b_se_pd = f.se[0];
    row.b_se_rho = f.se[1];
    row.b_flags = f.flags;
    for (int a = 0; a < 2; ++a) {
        row.b_prof_flags[a] = f.prof_flags[a];
        row.b_cover[a] = static_cast<std::int8_t>(rc::profile_covers(a, f, a == 0 ? s.pd : s.rho));
    }

    const std::vector<std::int64_t> d = rc::panel(s, replicate);
    const auto T = static_cast<std::int64_t>(d.size());
    std::vector<rc::Objective::Obs> counts(d.size());
    std::vector<o::RateObs> rates(d.size());
    std::vector<double> r(d.size()), r_kept;
    for (std::size_t t = 0; t < d.size(); ++t) {
        counts[t] = {s.obligors, d[t]};
        rates[t] = o::rate_obs(s.obligors, d[t]);
        r[t] = rates[t].rate;
        row.zero_periods += d[t] == 0;
        row.full_periods += d[t] == s.obligors;
        if (!o::is_boundary_rate(rates[t])) r_kept.push_back(r[t]);
    }

    row.mom = mom(engine::mom_from_counts(primary, counts.data(), T, g.axis[0].lo, g.axis[0].hi, g.axis[1].lo,
                                          g.axis[1].hi));
    row.mom_rate = mom(engine::mom_from_rates(primary, r.data(), T, g.axis[0].lo, g.axis[0].hi, g.axis[1].lo,
                                              g.axis[1].hi));
    if (r_kept.size() >= 3) {
        row.mom_rate_drop = mom(engine::mom_from_rates(primary, r_kept.data(), static_cast<std::int64_t>(r_kept.size()),
                                                       g.axis[0].lo, g.axis[0].hi, g.axis[1].lo, g.axis[1].hi));
    } else {
        row.mom_rate_drop.flags = engine::kMomRefused;
    }

    if (row.zero_periods + row.full_periods == 0) {
        const RateResult one = rate_fit<Refuse>(s, rates);
        for (auto& x : row.rate) x = one;
    } else {
        row.rate[0] = RateResult{};  // refuse
        const auto kept = o::drop_zero_rate_periods(rates.data(), T);
        row.rate[1] = kept.size() >= 3 ? rate_fit<Refuse>(s, kept) : RateResult{};
        row.rate[2] = rate_fit<Censor>(s, rates);
        row.rate[3] = rate_fit<Substitute>(s, rates);
    }
    return row;
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids;
    std::string out_path;
    unsigned threads = std::thread::hardware_concurrency();
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            for (const char* c = argv[++i]; *c != '\0';) {
                char* end = nullptr;
                const unsigned long id = std::strtoul(c, &end, 10);
                if (end == c || id >= kPassScenarios) {
                    std::fprintf(stderr, "--scenarios: ids 0..%u, comma-separated\n", kPassScenarios - 1);
                    return 2;
                }
                ids.push_back(static_cast<std::uint32_t>(id));
                c = *end == ',' ? end + 1 : end;
            }
        } else {
            std::fprintf(stderr,
                         "usage: study_estimator_pass [--replicates R] [--scenarios ID,...] [--out FILE] [--threads N]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < kPassScenarios; ++id) ids.push_back(id);
    }
    if (threads == 0) threads = 1;

    // Jobs, heaviest first (n, then T), so the long n = 10^6 fits do not finish last.
    std::vector<std::pair<std::size_t, std::uint32_t>> jobs;  // (scenario index, replicate)
    for (std::size_t k = 0; k < ids.size(); ++k) {
        for (std::uint32_t r = 0; r < R; ++r) jobs.emplace_back(k, r);
    }
    std::stable_sort(jobs.begin(), jobs.end(), [&](const auto& a, const auto& b) {
        const rc::Scenario sa = pass_scenario(ids[a.first]), sb = pass_scenario(ids[b.first]);
        return sa.obligors * sa.periods > sb.obligors * sb.periods;
    });

    // Surface rows are shared across panels for n <= 10^4 (D-167); at n >= 10^5 rows rarely repeat,
    // so those surfaces are computed directly (identical bit for bit either way).
    engine::SurfaceRowCache cache;
    std::vector<Row> rows(ids.size() * R);
    std::atomic<std::size_t> next{0};
    std::mutex progress;
    std::size_t done = 0;
    const auto start = std::chrono::steady_clock::now();
    auto worker = [&] {
        for (;;) {
            const std::size_t j = next.fetch_add(1);
            if (j >= jobs.size()) return;
            const auto [k, r] = jobs[j];
            const rc::Scenario s = pass_scenario(ids[k]);
            rows[k * R + r] = pass_fit(s, r, s.obligors <= 10000 ? &cache : nullptr);
            const std::lock_guard<std::mutex> lock(progress);
            if (++done % 1000 == 0 || done == jobs.size()) {
                std::fprintf(stderr, "fits %zu/%zu, %.0f s\n", done, jobs.size(),
                             std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    std::printf("%u replicates x %zu scenarios in %.0f s (%u threads; %zu cached surface rows)\n", R, ids.size(),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(), threads, cache.size());

    // A one-line digest per scenario, for reading the run; the scoring is compare.py's.
    for (std::size_t k = 0; k < ids.size(); ++k) {
        const rc::Scenario s = pass_scenario(ids[k]);
        std::int64_t est[4] = {}, mom_refused = 0;
        for (std::uint32_t r = 0; r < R; ++r) {
            const Row& w = rows[k * R + r];
            for (int m = 0; m < 4; ++m) est[m] += w.rate[m].status == kEstimated;
            mom_refused += (w.mom.flags & engine::kMomRefused) != 0;
        }
        std::printf("%3u PD %-5g rho %-4g T %-3lld n %-7lld | rate MLE estimates: refuse %4lld drop %4lld censor %4lld "
                    "substitute %4lld | MoM refused %lld\n",
                    s.id, s.pd, s.rho, static_cast<long long>(s.periods), static_cast<long long>(s.obligors),
                    static_cast<long long>(est[0]), static_cast<long long>(est[1]), static_cast<long long>(est[2]),
                    static_cast<long long>(est[3]), static_cast<long long>(mom_refused));
    }

    if (!out_path.empty()) {
        std::ofstream out(out_path, std::ios::binary);
        out << "scenario,replicate,pd_true,rho_true,periods,obligors,zero_periods,full_periods,"
               "bin_pd,bin_rho,bin_se_pd,bin_se_rho,bin_flags,bin_prof_flags_pd,bin_prof_flags_rho,bin_cover_pd,"
               "bin_cover_rho,mom_pd,mom_rho,mom_flags";
        for (const char* t : kTreatments) {
            out << ',' << t << "_status," << t << "_pd," << t << "_rho," << t << "_flags," << t << "_prof_flags_pd,"
                << t << "_prof_flags_rho," << t << "_cover_pd," << t << "_cover_rho";
        }
        out << ",momr_pd,momr_rho,momr_flags,momrd_pd,momrd_rho,momrd_flags\n";
        char line[512];
        for (std::size_t k = 0; k < ids.size(); ++k) {
            const rc::Scenario s = pass_scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const Row& w = rows[k * R + r];
                std::snprintf(line, sizeof line,
                              "%u,%u,%.17g,%.17g,%lld,%lld,%d,%d,%.17g,%.17g,%.17g,%.17g,%u,%u,%u,%d,%d,%.17g,%.17g,%u",
                              s.id, r, s.pd, s.rho, static_cast<long long>(s.periods),
                              static_cast<long long>(s.obligors), w.zero_periods, w.full_periods, w.b_pd, w.b_rho,
                              w.b_se_pd, w.b_se_rho, w.b_flags, w.b_prof_flags[0], w.b_prof_flags[1], w.b_cover[0],
                              w.b_cover[1], w.mom.pd, w.mom.rho, w.mom.flags);
                out << line;
                for (const RateResult& x : w.rate) {
                    std::snprintf(line, sizeof line, ",%d,%.17g,%.17g,%u,%u,%u,%d,%d", x.status, x.pd, x.rho, x.flags,
                                  x.prof_flags[0], x.prof_flags[1], x.cover[0], x.cover[1]);
                    out << line;
                }
                std::snprintf(line, sizeof line, ",%.17g,%.17g,%u,%.17g,%.17g,%u\n", w.mom_rate.pd, w.mom_rate.rho,
                              w.mom_rate.flags, w.mom_rate_drop.pd, w.mom_rate_drop.rho, w.mom_rate_drop.flags);
                out << line;
            }
        }
        if (!out) {
            std::fprintf(stderr, "cannot write %s\n", out_path.c_str());
            return 1;
        }
        std::printf("wrote %s\n", out_path.c_str());
    }
    return 0;
}
