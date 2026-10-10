// SPDX-License-Identifier: Apache-2.0
//
// The S-35 run (studies/lrpd-serial-correlation/PREDICTION.md, L1-L11; D-190): the recovery matrix
// under the AR(1)-factor DGP (dgp::simulate_panel_ar1) at each phi, on the recovery seed, scenario ids and
// replicates, so every phi shares its random numbers and phi = 0 is the pinned run's panels. Per panel:
//   - the simple average a = (1/T) sum d_t / n and its t interval a +/- t_{T-1,0.975} s / sqrt(T)
//     (s with the T - 1 denominator);
//   - the parity binomial-mixture MLE, fitted exactly as the recovery harness fits it
//     (recovery::fit_counts, every arm), with its coverage indicators computed by the harness's own
//     functions;
// and per (PD, rho, T) cell and replicate, the Vasicek-rate MLE on the n -> infinity rates
// r_t = p(Z_t) of the AR(1) factor of the cell's n = 10,000 scenario, fitted as rate_recovery_harness's
// rate_fit fits them (refuse; calibrate, profile intervals, q's interval); and per panel the observed-rate
// arm, the rate MLE on d_t / n with S-28's censored treatment, as study_estimator_pass fits it.
// Nothing is scored here: studies/lrpd-serial-correlation/compare.py scores the rows, after checking that
// the phi = 0 rows reproduce the pinned summaries' counts.
//
//   study_lrpd_serial [--phis 0,0.3,0.6] [--replicates R] [--scenarios ID,...] [--out FILE]
//                     [--rate-out FILE] [--threads N]
//                     [--timing-out FILE]  (a runtime timing sample, perf/README.md, D-187)
//                     [--count-repeats]    (P-9 on the AR(1) DGP: count repeated panels, fit nothing)
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

#include "core/objectives/vasicek_rate.hpp"
#include "engine/conditional_pd.hpp"
#include "tests/harness/timing_sample.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace o = vcal::objectives;
namespace e = vcal::engine;
using Clock = std::chrono::steady_clock;
using RateObjective = o::VasicekRate<PrecisionF64, o::ZeroRates::Refuse>;
using CensorObjective = o::VasicekRate<PrecisionF64, o::ZeroRates::Censor>;

enum Phase { kSimulate, kBinomial, kRate, kObserved, kCacheFill, kPhases };
const char* const kPhaseNames[kPhases] = {"simulate", "binomial MLE", "rate MLE", "observed-rate MLE",
                                          "cache fill"};

struct Timer {
    double* acc;
    Clock::time_point t0 = Clock::now();
    ~Timer() { *acc += std::chrono::duration<double>(Clock::now() - t0).count(); }
};

// -1 not computed (or flagged, for Wald and t), 0 misses, 1 covers.
std::int8_t ind(bool computed, bool covers) { return static_cast<std::int8_t>(!computed ? -1 : covers ? 1 : 0); }

struct Row {
    double avg, avg_lo, avg_hi;
    double pd, rho;
    std::uint32_t flags;
    std::int8_t wald[2], t[2], prof[2], boot[2], q_prof;
    std::uint32_t prof_flags[2];
};

// The observed-rate arm (PREDICTION.md, addendum: L12-L15): the rate MLE on d_t / n, fitted exactly as
// study_estimator_pass fits its censored treatment, so that phi = 0 reproduces the pinned S-8 rows.
struct ObservedRow {
    double pd = std::nan(""), rho = std::nan("");
    std::uint32_t flags = 0;
    std::int8_t cover[2] = {-1, -1};  // -1 not estimated (an engine error); a not-computed interval misses
};

template <class O>
ObservedRow observed_fit(const rc::Scenario& s, const std::vector<o::RateObs>& obs) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    const auto T = static_cast<std::int64_t>(obs.size());
    std::vector<double> L;
    e::Estimate2 est{};
    ObservedRow out;
    if (e::calibrate(backends::CpuBackend{1}, O{}, primary, check, obs.data(), T, g, L, est) != e::Status::Ok) return out;
    const auto prof = e::profile_intervals(O{}, primary, obs.data(), T, g, L, est);
    out.pd = prof.max_at[0];
    out.rho = prof.max_at[1];
    out.flags = est.flags;
    const double truth[2] = {s.pd, s.rho};
    for (int a = 0; a < 2; ++a) {
        out.cover[a] = static_cast<std::int8_t>(!(prof.flags[a] & e::kIntervalNotComputed) && prof.lo[a] <= truth[a] &&
                                                truth[a] <= prof.hi[a]);
    }
    return out;
}

ObservedRow observed_row(const rc::Scenario& s, const std::vector<std::int64_t>& d) {
    std::vector<o::RateObs> rates(d.size());
    bool boundary = false;
    for (std::size_t t = 0; t < d.size(); ++t) {
        rates[t] = o::rate_obs(s.obligors, d[t]);
        boundary = boundary || o::is_boundary_rate(rates[t]);
    }
    return boundary ? observed_fit<CensorObjective>(s, rates) : observed_fit<RateObjective>(s, rates);
}

struct RateRow {
    double pd, rho;
    std::uint32_t flags;
    std::int8_t cover[3];  // PD, rho, q profile intervals
};

Row binomial_row(const rc::Scenario& s, std::uint32_t r, const std::vector<std::int64_t>& d, e::SurfaceRowCache& cache) {
    Row row{};
    const auto T = static_cast<double>(s.periods);
    const auto n = static_cast<double>(s.obligors);
    double sum = 0.0;
    for (const auto x : d) sum += static_cast<double>(x) / n;
    row.avg = sum / T;
    double ss = 0.0;
    for (const auto x : d) {
        const double dev = static_cast<double>(x) / n - row.avg;
        ss += dev * dev;
    }
    const double half = rc::t_quantile(s.periods) * std::sqrt(ss / (T - 1.0)) / std::sqrt(T);
    row.avg_lo = row.avg - half;
    row.avg_hi = row.avg + half;

    const rc::Fit f = rc::fit_counts(s, r, d, &cache, rc::kArmAll);
    const Grid<2> g = rc::grid();
    const double truth[2] = {s.pd, s.rho};
    const bool flagged = (f.flags & rc::kNoReliableInterval) != 0;
    row.pd = f.value[0];
    row.rho = f.value[1];
    row.flags = f.flags;
    for (int a = 0; a < 2; ++a) {
        row.wald[a] = ind(!flagged, rc::covers(g, a, f, truth[a]));
        row.t[a] = ind(!flagged, rc::covers(g, a, f, truth[a], rc::t_quantile(s.periods)));
        row.prof[a] = ind(!(f.prof_flags[a] & e::kIntervalNotComputed), rc::profile_covers(a, f, truth[a]));
        row.boot[a] = ind(true, rc::bootstrap_covers(a, f, truth[a]));
        row.prof_flags[a] = f.prof_flags[a];
    }
    row.q_prof = ind(!(f.q_prof_flags & e::kIntervalNotComputed), rc::q_profile_covers(f, rc::q_truth(s)));
    return row;
}

// The rate panel: the AR(1) factor of the cell's n = 10,000 scenario (zero obligors: only Z_t is used),
// mapped through the model with the DGP's own functions, as rate_recovery_harness builds it at phi = 0.
RateRow rate_row(const rc::Scenario& s, std::uint32_t r, double phi) {
    std::vector<std::int64_t> zero(static_cast<std::size_t>(s.periods), 0), d(zero.size());
    std::vector<double> z(zero.size());
    const dgp::PanelSpec spec{rc::kSeed, s.id, r, s.pd, s.rho, zero.data(), s.periods};
    if (dgp::simulate_panel_ar1(spec, phi, d.data(), z.data()) != dgp::Status::Ok) std::abort();
    const double c = dgp::det_probit(s.pd), sr = std::sqrt(s.rho), s1 = std::sqrt(1.0 - s.rho);
    std::vector<o::RateObs> obs(z.size());
    for (std::size_t t = 0; t < obs.size(); ++t) obs[t] = {dgp::det_ncdf((c - sr * z[t]) / s1), 0.25};
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    std::vector<double> L;
    e::Estimate2 est{};
    if (e::calibrate(backends::CpuBackend{1}, RateObjective{}, primary, check, obs.data(), s.periods, g, L, est) !=
        e::Status::Ok) {
        std::abort();
    }
    const auto prof = e::profile_intervals(RateObjective{}, primary, obs.data(), s.periods, g, L, est);
    const auto qi = e::conditional_pd_interval(RateObjective{}, primary, obs.data(), s.periods, g, L, est, prof);
    RateRow out{};
    out.pd = prof.max_at[0];
    out.rho = prof.max_at[1];
    out.flags = est.flags;
    const double truth[3] = {s.pd, s.rho, e::conditional_pd(s.pd, s.rho)};
    for (int a = 0; a < 3; ++a) {
        const double lo = a < 2 ? prof.lo[a] : qi.lo, hi = a < 2 ? prof.hi[a] : qi.hi;
        const std::uint32_t f = a < 2 ? prof.flags[a] : qi.flags;
        out.cover[a] = ind(!(f & e::kIntervalNotComputed), lo <= truth[a] && truth[a] <= hi);
    }
    return out;
}

// P-9's measurement on the AR(1) DGP (perf/p9_repeat_rate.md, D-189): the panels exactly as the run draws
// them, no fits. Per scenario and phi, and per (T, n) cell and phi, the panels that repeat an earlier one
// under the ordered key, the only key under which a memoised fit is bit for bit the fit it replaces.
int count_repeats(const std::vector<std::uint32_t>& ids, std::uint32_t R, const std::vector<double>& phis) {
    std::printf("phi,scope,scenario,periods,obligors,panels,repeats_ordered\n");
    for (const double phi : phis) {
        std::vector<std::pair<std::int64_t, std::int64_t>> cells;
        std::vector<std::vector<std::vector<std::int64_t>>> by_cell;
        for (const std::uint32_t id : ids) {
            const rc::Scenario s = rc::scenario(id);
            std::vector<std::vector<std::int64_t>> panels(R);
            for (std::uint32_t r = 0; r < R; ++r) {
                std::vector<std::int64_t> n(static_cast<std::size_t>(s.periods), s.obligors), d(n.size());
                const dgp::PanelSpec spec{rc::kSeed, s.id, r, s.pd, s.rho, n.data(), s.periods};
                if (dgp::simulate_panel_ar1(spec, phi, d.data(), nullptr) != dgp::Status::Ok) std::abort();
                panels[r] = d;
            }
            const auto cell = std::make_pair(s.periods, s.obligors);
            std::size_t c = 0;
            while (c < cells.size() && cells[c] != cell) ++c;
            if (c == cells.size()) {
                cells.push_back(cell);
                by_cell.emplace_back();
            }
            by_cell[c].insert(by_cell[c].end(), panels.begin(), panels.end());
            std::sort(panels.begin(), panels.end());
            const auto distinct = static_cast<std::size_t>(std::unique(panels.begin(), panels.end()) - panels.begin());
            std::printf("%.17g,scenario,%u,%lld,%lld,%u,%zu\n", phi, id, static_cast<long long>(s.periods),
                        static_cast<long long>(s.obligors), R, R - distinct);
        }
        for (std::size_t c = 0; c < cells.size(); ++c) {
            auto& v = by_cell[c];
            const std::size_t total = v.size();
            std::sort(v.begin(), v.end());
            const auto distinct = static_cast<std::size_t>(std::unique(v.begin(), v.end()) - v.begin());
            std::printf("%.17g,cell,,%lld,%lld,%zu,%zu\n", phi, static_cast<long long>(cells[c].first),
                        static_cast<long long>(cells[c].second), total, total - distinct);
        }
    }
    return 0;
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
    std::vector<double> phis = {0.0, 0.3, 0.6};
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids;
    std::string out = "lrpd_serial.csv", rate_out = "lrpd_serial_rate.csv", timing_out;
    unsigned threads = std::thread::hardware_concurrency();
    bool count_only = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--phis") == 0 && i + 1 < argc) {
            phis = parse_list(argv[++i]);
        } else if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            for (const double v : parse_list(argv[++i])) ids.push_back(static_cast<std::uint32_t>(v));
        } else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out = argv[++i];
        } else if (std::strcmp(argv[i], "--rate-out") == 0 && i + 1 < argc) {
            rate_out = argv[++i];
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--timing-out") == 0 && i + 1 < argc) {
            timing_out = argv[++i];
        } else if (std::strcmp(argv[i], "--count-repeats") == 0) {
            count_only = true;
        } else {
            std::fprintf(stderr, "usage: study_lrpd_serial [--phis 0,0.3,0.6] [--replicates R] [--scenarios ID,...] "
                                 "[--out FILE] [--rate-out FILE] [--threads N] [--timing-out FILE]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    if (threads == 0) threads = 1;
    const std::size_t P = phis.size();
    if (count_only) return count_repeats(ids, R, phis);

    // Jobs: (scenario, replicate). A job fits its panel at every phi; the rate MLE is fitted in the jobs of
    // n = 10,000 scenarios, whose factor draws the rate panels use (one per (PD, rho, T) cell).
    const std::size_t S = ids.size();
    std::vector<Row> rows(S * R * P);
    std::vector<RateRow> rate_rows(S * R * P);
    std::vector<ObservedRow> observed_rows(S * R * P);
    std::vector<char> has_rate(S, 0);
    for (std::size_t k = 0; k < S; ++k) has_rate[k] = rc::scenario(ids[k]).obligors == 10000 ? 1 : 0;
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
                {
                    const vcal::test::FillMeter fills;
                    {
                        Timer t{&sec[kBinomial]};
                        rows[job * P + p] = binomial_row(s, r, d, cache);
                    }
                    sec[kBinomial] -= fills.seconds();  // a fill is paid once per run, not per replicate (D-193)
                    sec[kCacheFill] += fills.seconds();
                }
                {
                    Timer t{&sec[kObserved]};
                    observed_rows[job * P + p] = observed_row(s, d);
                }
                if (has_rate[k]) {
                    Timer t{&sec[kRate]};
                    rate_rows[job * P + p] = rate_row(s, r, phis[p]);
                }
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    const double wall = std::chrono::duration<double>(Clock::now() - t0).count();

    char buf[512];
    std::string text = "# S-35 (D-190): one row per (phi, scenario, replicate)\n"
                       "# coverage indicators: 1 covers, 0 misses, -1 not computed (Wald and t: flagged fit)\n"
                       "phi,scenario,replicate,pd_true,rho_true,periods,obligors,avg,avg_lo,avg_hi,pd,rho,flags,"
                       "wald_pd,wald_rho,t_pd,t_rho,prof_pd,prof_rho,boot_pd,boot_rho,q_prof,prof_flags_pd,"
                       "prof_flags_rho,obs_pd,obs_rho,obs_flags,obs_cover_pd,obs_cover_rho\n";
    std::string rtext = "# S-35 rate MLE on the n -> infinity rates: one row per (phi, PD, rho, T, replicate)\n"
                        "phi,z_source_scenario,replicate,pd_true,rho_true,periods,pd,rho,flags,cover_pd,cover_rho,"
                        "cover_q\n";
    for (std::size_t p = 0; p < P; ++p) {
        for (std::size_t k = 0; k < S; ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const std::size_t job = k * R + r;
                const Row& w = rows[job * P + p];
                const ObservedRow& ob = observed_rows[job * P + p];
                std::snprintf(buf, sizeof buf,
                              "%.17g,%u,%u,%.17g,%.17g,%lld,%lld,%.17g,%.17g,%.17g,%.17g,%.17g,%u,%d,%d,%d,%d,%d,%d,%d,%d,"
                              "%d,%u,%u,%.17g,%.17g,%u,%d,%d\n",
                              phis[p], s.id, r, s.pd, s.rho, static_cast<long long>(s.periods),
                              static_cast<long long>(s.obligors), w.avg, w.avg_lo, w.avg_hi, w.pd, w.rho, w.flags,
                              w.wald[0], w.wald[1], w.t[0], w.t[1], w.prof[0], w.prof[1], w.boot[0], w.boot[1],
                              w.q_prof, w.prof_flags[0], w.prof_flags[1], ob.pd, ob.rho, ob.flags, ob.cover[0],
                              ob.cover[1]);
                text += buf;
                if (!has_rate[k]) continue;
                const RateRow& q = rate_rows[job * P + p];
                std::snprintf(buf, sizeof buf, "%.17g,%u,%u,%.17g,%.17g,%lld,%.17g,%.17g,%u,%d,%d,%d\n", phis[p], s.id, r,
                              s.pd, s.rho, static_cast<long long>(s.periods), q.pd, q.rho, q.flags, q.cover[0],
                              q.cover[1], q.cover[2]);
                rtext += buf;
            }
        }
    }
    for (const auto& [path, body] : {std::make_pair(out, &text), std::make_pair(rate_out, &rtext)}) {
        std::ofstream f(path, std::ios::binary);
        f << *body;
        f.close();
        if (!f) {
            std::fprintf(stderr, "cannot write %s\n", path.c_str());
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
        vcal::test::TimingSample sample("study_lrpd_serial", phases, shares);
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
