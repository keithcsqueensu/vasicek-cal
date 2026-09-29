// SPDX-License-Identifier: Apache-2.0
//
// The S-10 / S-4 run (studies/parametric-bootstrap/PREDICTION.md, K1-K5; studies/bartlett-profile/
// PREDICTION.md, H1-H6; D-176). Every recovery panel is fitted as the recovery harness fits it, then:
//   - S-10: B = 999 parametric-bootstrap panels simulated at theta-hat with the recovery DGP (D-109)
//     under seed S10PBOOT, scenario id as the panel's, replicate 1,000 r + b; each fitted as the
//     recovery panels are (the parity grid through the row cache). Per parameter, in logit: the
//     percentile interval (type-7 2.5% and 97.5% points of the finite bootstrap estimates) and the
//     studentised interval (t*_b = (u*_b - u-hat)/se*_b, se from analytic_se_scaled, the analytic
//     observed information; bootstrap panels without it are left out, and with more than 10% left out
//     the interval is not computed);
//   - the comparators on the same rows: the profile interval (recomputed) and, for replicates 0-49,
//     the Jeffreys equal-tailed interval with the resolution rule, computed as S-9 computed it, for a
//     bit-for-bit check against S-9's rows (the other replicates are joined from them);
//   - S-4: W0 = 2 (l_max - P_a(truth)) per parameter (profile_lr_statistic); S-4b's factor k_b = the
//     mean of W*_b = 2 (l*_max - P*_a(theta-hat)) over the first 199 of the panel's own S-10 bootstrap
//     panels (b = 0..198, asserted: the same panels, identified by index and content hash), 1 where
//     theta-hat_a is on the box's bound or more than 10% of the W*_b are not finite. For replicates
//     0-199 a second pass solves the corrected intervals (threshold c k / 2) with S-4a's in-sample factor
//     (the scenario's mean W0) and with k_b, for widths.
// Nothing is scored here: studies/parametric-bootstrap/compare.py scores K1-K5 and H1-H6 from the rows.
//
//   study_parametric_bootstrap [--replicates R] [--scenarios ID,...] [--out FILE] [--threads N]
//
// Doubles are written to 17 significant digits (an exact round trip). Results do not depend on the
// thread count. The phase timings (CPU seconds summed over fits) are printed per scenario.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "engine/posterior.hpp"
#include "resample/bootstrap.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace e = vcal::engine;
using Clock = std::chrono::steady_clock;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

inline constexpr std::uint64_t kBootSeed = 0x53313050424F4F54ull;  // "S10PBOOT"
inline constexpr std::uint32_t kB = 999;                           // S-10's bootstrap panels
inline constexpr std::uint32_t kBartlettPanels = 199;              // S-4b: the first 199 of them
inline constexpr std::uint32_t kJeffreysCheck = 50;                // replicates 0-49
inline constexpr std::uint32_t kWidthReplicates = 200;             // S-4 corrected intervals solved
inline constexpr double kLeftOutMax = 0.10;
inline constexpr double kC = 2.0 * e::kProfileThreshold95;         // chi2_1 0.95 quantile

double logit(double v) { return std::log(v) - std::log1p(-v); }
double inv_logit(double u) { return grid::from_scaled(AxisScale::Logit, u); }

enum Phase { kSimulate, kSurface, kFit, kSe, kProfile, kLr0, kLrBoot, kJeffreys, kWidths, kPhases };
const char* const kPhaseNames[kPhases] = {"simulate", "surface", "argmax+refine+SE", "analytic SE", "profile",
                                          "W0", "W* (S-4b)", "Jeffreys", "S-4 widths"};

struct Timer {
    double* acc;
    Clock::time_point t0 = Clock::now();
    ~Timer() { *acc += std::chrono::duration<double>(Clock::now() - t0).count(); }
};

// FNV-1a over a panel's default counts: identifies the panel S-4b uses as S-10's own.
std::uint64_t fnv(std::uint64_t h, std::int64_t x) {
    for (int k = 0; k < 8; ++k) {
        h ^= static_cast<std::uint64_t>(x >> (8 * k)) & 0xFFu;
        h *= 1099511628211ull;
    }
    return h;
}
std::uint64_t panel_hash(const std::vector<std::int64_t>& d) {
    std::uint64_t h = 1469598103934665603ull;
    for (const auto x : d) h = fnv(h, x);
    return h;
}
std::uint64_t panel_hash(const std::vector<rc::Objective::Obs>& obs) {
    std::uint64_t h = 1469598103934665603ull;
    for (const auto& y : obs) h = fnv(h, y.d);
    return h;
}

struct Row {
    double pd, rho;
    std::uint32_t flags;
    double se_calib_u[2], se_analytic_u[2];
    double prof_lo[2], prof_hi[2];
    std::uint32_t prof_flags[2];
    double w0[2];
    double jef_lo[2], jef_hi[2];
    std::uint32_t jef_flags;
    double pct_lo[2], pct_hi[2];
    std::int32_t pct_left_out[2];
    double stud_lo[2], stud_hi[2];
    std::int32_t stud_left_out[2];
    double kb[2];
    std::int32_t kb_fallback[2], wstar_finite[2];
    double s4a_lo[2], s4a_hi[2], s4b_lo[2], s4b_hi[2];
    std::uint32_t s4a_flags[2], s4b_flags[2];
    double seconds[kPhases];
};

struct Context {
    e::SurfaceRowCache cache;
    std::map<std::int64_t, e::JeffreysTable> jeffreys;
};

// The recovery fit's calibrate_cached, split in its two calls so the phases can be timed; identical.
e::Estimate2 fit_panel(Context& ctx, const std::vector<rc::Objective::Obs>& obs, std::vector<double>& L,
                       double* sec) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    const auto T = static_cast<std::int64_t>(obs.size());
    e::Estimate2 est{};
    if (vcal::grid_error(g) != nullptr || rc::Objective::panel_error(obs.data(), T) != nullptr) std::abort();
    {
        Timer t{&sec[kSurface]};
        L.assign(static_cast<std::size_t>(T * g.size()), 0.0);
        e::evaluate_surface_cached(backends::CpuBackend{1}, ctx.cache, rc::Objective{}, primary, obs.data(), T, g,
                                   L.data());
    }
    Timer t{&sec[kFit]};
    if (e::calibrate_from_surface(backends::CpuBackend{1}, rc::Objective{}, primary, check, obs.data(), T, g, L, est) !=
        e::Status::Ok) {
        std::abort();
    }
    return est;
}

std::vector<rc::Objective::Obs> as_obs(const std::vector<std::int64_t>& d, std::int64_t n) {
    std::vector<rc::Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {n, d[t]};
    return obs;
}

Row fit_one(Context& ctx, const rc::Scenario& s, std::uint32_t r) {
    static const auto primary = quadrature::parity_rule();
    const Grid<2> g = rc::grid();
    Row row{};
    double* sec = row.seconds;
    std::vector<std::int64_t> d;
    {
        Timer t{&sec[kSimulate]};
        d = rc::panel(s, r);
    }
    const auto obs = as_obs(d, s.obligors);
    const auto T = static_cast<std::int64_t>(obs.size());
    std::vector<double> L;
    const e::Estimate2 est = fit_panel(ctx, obs, L, sec);
    row.pd = est.value[0];
    row.rho = est.value[1];
    row.flags = est.flags;
    const double v_hat[2] = {est.value[0], est.value[1]};
    double u_hat[2], se_calib[2];
    for (int a = 0; a < 2; ++a) {
        u_hat[a] = logit(v_hat[a]);
        se_calib[a] = est.se[a] / grid::dvalue_dscaled(g.axis[a].scale, grid::to_scaled(g.axis[a].scale, est.value[a]));
        row.se_calib_u[a] = se_calib[a];
    }
    e::AnalyticSe se_hat{};
    {
        Timer t{&sec[kSe]};
        se_hat = e::analytic_se_scaled(rc::Objective{}, primary, obs.data(), T, g, v_hat);
    }
    for (int a = 0; a < 2; ++a) row.se_analytic_u[a] = se_hat.ok ? se_hat.se_u[a] : kNaN;
    {
        Timer t{&sec[kProfile]};
        const auto prof = e::profile_intervals(rc::Objective{}, primary, obs.data(), T, g, L, est);
        for (int a = 0; a < 2; ++a) {
            row.prof_lo[a] = prof.lo[a];
            row.prof_hi[a] = prof.hi[a];
            row.prof_flags[a] = prof.flags[a];
        }
    }
    {
        Timer t{&sec[kLr0]};
        const double truth[2] = {s.pd, s.rho};
        const auto lr = e::profile_lr_statistic(rc::Objective{}, primary, obs.data(), T, g, L, est, truth);
        row.w0[0] = lr.w[0];
        row.w0[1] = lr.w[1];
    }
    row.jef_lo[0] = row.jef_lo[1] = row.jef_hi[0] = row.jef_hi[1] = kNaN;
    if (r < kJeffreysCheck) {
        Timer t{&sec[kJeffreys]};
        // As S-9's tool: the posterior's local grid seeded by calibrate's SE in u.
        const auto p = e::grid_posterior(backends::CpuBackend{1}, rc::Objective{}, primary, obs.data(), T, g, L,
                                         e::Prior::Jeffreys, &ctx.jeffreys.at(s.obligors), se_calib, true);
        row.jef_flags = p.flags;
        for (int a = 0; a < 2; ++a) {
            row.jef_lo[a] = p.et_lo[a];
            row.jef_hi[a] = p.et_hi[a];
        }
    }

    // S-10's bootstrap panels, and S-4b's W* on the first 199 of them.
    std::vector<double> ub[2], tb[2];
    std::vector<double> wstar[2];
    // (index, content hash) of S-10's first 199 panels as simulated, and of the panels S-4b fits.
    std::vector<std::pair<std::uint32_t, std::uint64_t>> s10_panels, s4b_panels;
    std::vector<std::int64_t> n(static_cast<std::size_t>(T), s.obligors), db(n.size());
    std::vector<double> Lb;
    for (std::uint32_t b = 0; b < kB; ++b) {
        {
            Timer t{&sec[kSimulate]};
            const dgp::PanelSpec spec{kBootSeed, s.id, 1000u * r + b, v_hat[0], v_hat[1], n.data(), T};
            if (dgp::simulate_panel(spec, db.data(), nullptr) != dgp::Status::Ok) std::abort();
        }
        if (b < kBartlettPanels) s10_panels.emplace_back(b, panel_hash(db));
        const auto obs_b = as_obs(db, s.obligors);
        const e::Estimate2 eb = fit_panel(ctx, obs_b, Lb, sec);
        const bool finite = std::isfinite(eb.value[0]) && std::isfinite(eb.value[1]);
        e::AnalyticSe seb{{kNaN, kNaN}, false};
        if (finite && !(eb.flags & (e::kFlagFlatSurface | e::kFlagNumeric))) {
            Timer t{&sec[kSe]};
            const double vb[2] = {eb.value[0], eb.value[1]};
            seb = e::analytic_se_scaled(rc::Objective{}, primary, obs_b.data(), T, g, vb);
        }
        for (int a = 0; a < 2; ++a) {
            if (!finite) continue;
            const double u = logit(eb.value[a]);
            ub[a].push_back(u);
            if (seb.ok) tb[a].push_back((u - u_hat[a]) / seb.se_u[a]);
        }
        if (b < kBartlettPanels) {
            s4b_panels.emplace_back(b, panel_hash(obs_b));  // the observations S-4b actually fits
            Timer t{&sec[kLrBoot]};
            const auto lr = e::profile_lr_statistic(rc::Objective{}, primary, obs_b.data(), T, g, Lb, eb, v_hat);
            for (int a = 0; a < 2; ++a) wstar[a].push_back(lr.w[a]);
        }
    }
    // The registration's requirement: S-4b's factor uses exactly S-10's own panels 0..198 (same seeds and
    // indices, same content), not a separate draw.
    bool same = s4b_panels.size() == kBartlettPanels && s4b_panels == s10_panels && wstar[0].size() == kBartlettPanels;
    for (std::uint32_t b = 0; same && b < kBartlettPanels; ++b) same = s4b_panels[b].first == b;
    if (!same) {
        std::fprintf(stderr, "S-4b did not use exactly S-10's bootstrap panels 0..%u (scenario %u, replicate %u)\n",
                     kBartlettPanels - 1, s.id, r);
        std::abort();
    }

    for (int a = 0; a < 2; ++a) {
        std::sort(ub[a].begin(), ub[a].end());
        std::sort(tb[a].begin(), tb[a].end());
        row.pct_left_out[a] = static_cast<std::int32_t>(kB - ub[a].size());
        row.pct_lo[a] = ub[a].empty() ? kNaN : inv_logit(resample::quantile_type7(ub[a], 0.025));
        row.pct_hi[a] = ub[a].empty() ? kNaN : inv_logit(resample::quantile_type7(ub[a], 0.975));
        row.stud_left_out[a] = static_cast<std::int32_t>(kB - tb[a].size());
        const bool stud_ok = se_hat.ok && static_cast<double>(row.stud_left_out[a]) <= kLeftOutMax * kB && !tb[a].empty();
        row.stud_lo[a] = stud_ok ? inv_logit(u_hat[a] - resample::quantile_type7(tb[a], 0.975) * se_hat.se_u[a]) : kNaN;
        row.stud_hi[a] = stud_ok ? inv_logit(u_hat[a] - resample::quantile_type7(tb[a], 0.025) * se_hat.se_u[a]) : kNaN;
        // S-4b's factor.
        double sum = 0.0;
        std::int32_t fin = 0;
        for (const double w : wstar[a]) {
            if (std::isfinite(w)) {
                sum += w;
                ++fin;
            }
        }
        row.wstar_finite[a] = fin;
        const bool on_bound = !(est.value[a] > g.axis[a].lo && est.value[a] < g.axis[a].hi);
        const bool too_few = static_cast<double>(kBartlettPanels - fin) > kLeftOutMax * kBartlettPanels;
        row.kb_fallback[a] = on_bound || too_few || fin == 0;
        row.kb[a] = row.kb_fallback[a] ? 1.0 : sum / fin;
        row.s4a_lo[a] = row.s4a_hi[a] = row.s4b_lo[a] = row.s4b_hi[a] = kNaN;
    }
    return row;
}

// Second pass (S-4 widths): the corrected intervals for replicates 0-199, S-4a with the scenario's
// in-sample factor k_a (the mean W0), S-4b with the replicate's own k_b.
void widths(Context& ctx, const rc::Scenario& s, std::uint32_t r, const double (&ka)[2], Row& row) {
    static const auto primary = quadrature::parity_rule();
    const Grid<2> g = rc::grid();
    Timer t{&row.seconds[kWidths]};
    const auto obs = as_obs(rc::panel(s, r), s.obligors);
    const auto T = static_cast<std::int64_t>(obs.size());
    std::vector<double> L;
    double scratch[kPhases] = {};
    const e::Estimate2 est = fit_panel(ctx, obs, L, scratch);
    for (int a = 0; a < 2; ++a) {
        if (std::isfinite(ka[a])) {
            const auto p = e::profile_intervals(rc::Objective{}, primary, obs.data(), T, g, L, est, 0.5 * kC * ka[a]);
            row.s4a_lo[a] = p.lo[a];
            row.s4a_hi[a] = p.hi[a];
            row.s4a_flags[a] = p.flags[a];
        }
        const auto q = e::profile_intervals(rc::Objective{}, primary, obs.data(), T, g, L, est, 0.5 * kC * row.kb[a]);
        row.s4b_lo[a] = q.lo[a];
        row.s4b_hi[a] = q.hi[a];
        row.s4b_flags[a] = q.flags[a];
    }
}

template <class F>
void run_pool(std::size_t jobs, unsigned threads, F&& f) {
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back([&] {
            for (;;) {
                const std::size_t j = next.fetch_add(1);
                if (j >= jobs) return;
                f(j);
            }
        });
    }
    for (auto& t : pool) t.join();
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
                if (end == c || id >= rc::kScenarios) {
                    std::fprintf(stderr, "--scenarios: ids 0..%u, comma-separated\n", rc::kScenarios - 1);
                    return 2;
                }
                ids.push_back(static_cast<std::uint32_t>(id));
                c = *end == ',' ? end + 1 : end;
            }
        } else {
            std::fprintf(stderr, "usage: study_parametric_bootstrap [--replicates R] [--scenarios ID,...] [--out FILE] "
                                 "[--threads N]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    if (threads == 0) threads = 1;
    const auto start = Clock::now();
    const auto elapsed = [&] { return std::chrono::duration<double>(Clock::now() - start).count(); };

    Context ctx;
    for (const auto id : ids) {
        const auto n = rc::scenario(id).obligors;
        if (ctx.jeffreys.count(n) == 0) {
            ctx.jeffreys.emplace(n, e::jeffreys_table(backends::CpuBackend{static_cast<int>(threads)}, rc::Objective{},
                                                      quadrature::parity_rule(), n, rc::grid()));
        }
    }
    std::fprintf(stderr, "Jeffreys tables, %.0f s\n", elapsed());

    // Pass 1: every replicate, heaviest scenarios first (T, then n).
    std::vector<std::pair<std::size_t, std::uint32_t>> jobs;
    for (std::size_t k = 0; k < ids.size(); ++k) {
        for (std::uint32_t r = 0; r < R; ++r) jobs.emplace_back(k, r);
    }
    std::stable_sort(jobs.begin(), jobs.end(), [&](const auto& a, const auto& b) {
        const rc::Scenario sa = rc::scenario(ids[a.first]), sb = rc::scenario(ids[b.first]);
        return sa.periods * sa.obligors > sb.periods * sb.obligors;
    });
    std::vector<Row> rows(ids.size() * R);
    std::mutex progress;
    std::size_t done = 0;
    run_pool(jobs.size(), threads, [&](std::size_t j) {
        const auto [k, r] = jobs[j];
        rows[k * R + r] = fit_one(ctx, rc::scenario(ids[k]), r);
        const std::lock_guard<std::mutex> lock(progress);
        if (++done % 100 == 0 || done == jobs.size()) std::fprintf(stderr, "fits %zu/%zu, %.0f s\n", done, jobs.size(), elapsed());
    });

    // S-4a's in-sample factor per scenario and parameter, then pass 2 (widths).
    std::vector<std::array<double, 2>> ka(ids.size());
    for (std::size_t k = 0; k < ids.size(); ++k) {
        for (int a = 0; a < 2; ++a) {
            double sum = 0.0;
            std::int64_t m = 0;
            for (std::uint32_t r = 0; r < R; ++r) {
                const double w = rows[k * R + r].w0[a];
                if (std::isfinite(w)) {
                    sum += w;
                    ++m;
                }
            }
            ka[k][a] = m > 0 ? sum / static_cast<double>(m) : kNaN;
        }
    }
    const std::uint32_t RW = std::min(R, kWidthReplicates);
    run_pool(ids.size() * RW, threads, [&](std::size_t j) {
        const std::size_t k = j / RW;
        const auto r = static_cast<std::uint32_t>(j % RW);
        const double kk[2] = {ka[k][0], ka[k][1]};
        widths(ctx, rc::scenario(ids[k]), r, kk, rows[k * R + r]);
    });
    std::printf("%u replicates x %zu scenarios in %.0f s (%u threads; %zu cached surface rows)\n", R, ids.size(),
                elapsed(), threads, ctx.cache.size());

    // Phase timings, CPU seconds summed over fits, per scenario and in total.
    double total[kPhases] = {};
    std::printf("scenario  T    n     |");
    for (const char* p : kPhaseNames) std::printf(" %s |", p);
    std::printf("\n");
    for (std::size_t k = 0; k < ids.size(); ++k) {
        const rc::Scenario s = rc::scenario(ids[k]);
        double sum[kPhases] = {};
        for (std::uint32_t r = 0; r < R; ++r) {
            for (int p = 0; p < kPhases; ++p) sum[p] += rows[k * R + r].seconds[p];
        }
        std::printf("%-9u %-4lld %-5lld |", s.id, static_cast<long long>(s.periods), static_cast<long long>(s.obligors));
        for (int p = 0; p < kPhases; ++p) {
            std::printf(" %.1f |", sum[p]);
            total[p] += sum[p];
        }
        std::printf("\n");
    }
    double all = 0.0;
    for (const double t : total) all += t;
    std::printf("total CPU s %.0f:", all);
    for (int p = 0; p < kPhases; ++p) std::printf(" %s %.0f (%.0f%%);", kPhaseNames[p], total[p], 100.0 * total[p] / all);
    std::printf("\n");

    if (!out_path.empty()) {
        std::ofstream out(out_path, std::ios::binary);
        out << "scenario,replicate,pd_true,rho_true,periods,obligors,pd,rho,flags";
        for (const char* p : {"pd", "rho"}) {
            out << ",se_calib_u_" << p << ",se_analytic_u_" << p << ",prof_lo_" << p << ",prof_hi_" << p << ",prof_flags_" << p
                << ",w0_" << p << ",jef_lo_" << p << ",jef_hi_" << p << ",pct_lo_" << p << ",pct_hi_" << p
                << ",pct_left_out_" << p << ",stud_lo_" << p << ",stud_hi_" << p << ",stud_left_out_" << p << ",kb_" << p
                << ",kb_fallback_" << p << ",wstar_finite_" << p << ",s4a_lo_" << p << ",s4a_hi_" << p << ",s4a_flags_" << p
                << ",s4b_lo_" << p << ",s4b_hi_" << p << ",s4b_flags_" << p;
        }
        out << ",jef_flags\n";
        char line[512];
        for (std::size_t k = 0; k < ids.size(); ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const Row& w = rows[k * R + r];
                std::snprintf(line, sizeof line, "%u,%u,%.17g,%.17g,%lld,%lld,%.17g,%.17g,%u", s.id, r, s.pd, s.rho,
                              static_cast<long long>(s.periods), static_cast<long long>(s.obligors), w.pd, w.rho, w.flags);
                out << line;
                for (int a = 0; a < 2; ++a) {
                    std::snprintf(line, sizeof line,
                                  ",%.17g,%.17g,%.17g,%.17g,%u,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%d,%.17g,%d,%d,"
                                  "%.17g,%.17g,%u,%.17g,%.17g,%u",
                                  w.se_calib_u[a], w.se_analytic_u[a], w.prof_lo[a], w.prof_hi[a], w.prof_flags[a], w.w0[a],
                                  w.jef_lo[a], w.jef_hi[a], w.pct_lo[a], w.pct_hi[a], w.pct_left_out[a], w.stud_lo[a],
                                  w.stud_hi[a], w.stud_left_out[a], w.kb[a], w.kb_fallback[a], w.wstar_finite[a],
                                  w.s4a_lo[a], w.s4a_hi[a], w.s4a_flags[a], w.s4b_lo[a], w.s4b_hi[a], w.s4b_flags[a]);
                    out << line;
                }
                out << ',' << w.jef_flags << '\n';
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
