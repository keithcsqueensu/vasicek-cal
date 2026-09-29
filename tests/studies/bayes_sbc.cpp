// SPDX-License-Identifier: Apache-2.0
//
// S-15: simulation-based calibration of the grid-Bayesian estimator (studies/bayes-sbc/PREDICTION.md,
// G1-G5; D-161, D-166). For each setting (n, T) and prior, N draws of (PD, rho) from the prior on
// the parity grid's box, each with a simulated panel and its grid posterior, with the resolution rule
// (the estimator) and without it (the diagnostic arm). The rank statistic is the marginal posterior
// CDF at the true value (engine::marginal_cdf), for PD and for rho. One row per (setting, prior,
// draw); studies/bayes-sbc/compare.py scores the predictions from these rows.
//
// The draw, as registered: a parity-grid cell with the prior's mass (the prior's cell masses,
// engine::posterior_detail::masses with a zero log-likelihood), then a point uniformly within the cell
// in the logit coordinate (half cells at the box's ends). Its uniforms come from Philox4x32-10 keyed
// by kSbcSeed, counter (setting, draw, prior, block), words to doubles as the DGP does
// (dgp::uniform_from_words): block 0 words 0-1 pick the cell by inverting the cumulative masses in
// grid order; block 0 words 2-3 and block 1 words 0-1 place the point on the PD and rho axes. The
// panel is the recovery DGP's (D-109) with scenario id 1000 + setting and replicate = draw.
//
//   study_bayes_sbc [--draws N] [--out FILE] [--threads N]
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

#include "dgp/philox.hpp"
#include "engine/posterior.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace e = vcal::engine;

inline constexpr std::uint64_t kSbcSeed = 0x5331355342434452ull;  // "S15SBCDR"; arbitrary, fixed

struct Setting {
    std::int64_t n, T;
    bool jeffreys;  // Jeffreys is run at (1,000, 20) only
};
inline constexpr Setting kSettings[3] = {{1000, 20, true}, {10000, 40, false}, {1000, 100, false}};

struct Draw {
    double pd, rho;
};

// The prior's cell masses on the parity grid, cumulated in grid order.
std::vector<double> cumulative_prior(e::Prior prior, const e::JeffreysTable* jt, const Grid<2>& g) {
    std::vector<double> m;
    if (!e::posterior_detail::masses(prior, jt, g, std::vector<double>(static_cast<std::size_t>(g.size()), 0.0), m)) {
        std::abort();
    }
    for (std::size_t k = 1; k < m.size(); ++k) m[k] += m[k - 1];
    return m;
}

Draw draw_from_prior(const std::vector<double>& cum, const Grid<2>& g, std::uint32_t setting, std::uint32_t draw,
                     std::uint32_t prior) {
    const dgp::PhiloxKey key{static_cast<std::uint32_t>(kSbcSeed), static_cast<std::uint32_t>(kSbcSeed >> 32)};
    const auto w0 = dgp::philox4x32({setting, draw, prior, 0}, key);
    const auto w1 = dgp::philox4x32({setting, draw, prior, 1}, key);
    const double uc = dgp::uniform_from_words(w0[0], w0[1]) * cum.back();
    const auto it = std::upper_bound(cum.begin(), cum.end(), uc);
    const std::int64_t k = std::min<std::int64_t>(it - cum.begin(), g.size() - 1);
    std::int32_t idx[2];
    g.unflatten(k, idx);
    const double ua[2] = {dgp::uniform_from_words(w0[2], w0[3]), dgp::uniform_from_words(w1[0], w1[1])};
    double v[2];
    for (int a = 0; a < 2; ++a) {
        const Axis& ax = g.axis[a];
        const double u = ax.scaled_at(idx[a]), h = ax.step();
        const double lo = idx[a] == 0 ? u : u - 0.5 * h, hi = idx[a] == ax.n - 1 ? u : u + 0.5 * h;
        v[a] = grid::from_scaled(ax.scale, lo + ua[a] * (hi - lo));
    }
    return {v[0], v[1]};
}

struct Row {
    std::uint32_t setting, prior, draw;
    double pd, rho;
    std::uint32_t flags;  // with the rule
    std::int32_t refinements;
    double rank[2];         // with the rule (NaN if refused)
    std::uint32_t flags_off;
    double rank_off[2];     // without the rule
};

Row sbc_one(std::uint32_t setting, std::uint32_t prior_id, std::uint32_t draw, const std::vector<double>& cum,
            const e::JeffreysTable* jt) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Setting& st = kSettings[setting];
    const Grid<2> g = rc::grid();
    const e::Prior prior = prior_id == 0 ? e::Prior::Flat : e::Prior::Jeffreys;
    const Draw th = draw_from_prior(cum, g, setting, draw, prior_id);
    std::vector<std::int64_t> n(static_cast<std::size_t>(st.T), st.n), d(n.size());
    const dgp::PanelSpec spec{rc::kSeed, 1000 + setting, draw, th.pd, th.rho, n.data(), st.T};
    if (dgp::simulate_panel(spec, d.data(), nullptr) != dgp::Status::Ok) std::abort();
    std::vector<rc::Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {st.n, d[t]};
    std::vector<double> L;
    e::Estimate2 est{};
    if (e::calibrate(backends::CpuBackend{1}, rc::Objective{}, primary, check, obs.data(), st.T, g, L, est) !=
        e::Status::Ok) {
        std::abort();
    }
    double se_u[2];
    for (int a = 0; a < 2; ++a) {
        se_u[a] = est.se[a] / grid::dvalue_dscaled(g.axis[a].scale, grid::to_scaled(g.axis[a].scale, est.value[a]));
    }
    const double truth_u[2] = {grid::to_scaled(g.axis[0].scale, th.pd), grid::to_scaled(g.axis[1].scale, th.rho)};
    Row row{setting, prior_id, draw, th.pd, th.rho, 0, 0, {NAN, NAN}, 0, {NAN, NAN}};
    const auto on = e::grid_posterior(backends::CpuBackend{1}, rc::Objective{}, primary, obs.data(), st.T, g, L, prior,
                                      jt, se_u, true);
    row.flags = on.flags;
    row.refinements = on.refinements;
    if (!(on.flags & (e::kPosteriorRefused | e::kPosteriorNumeric))) {
        for (int a = 0; a < 2; ++a) row.rank[a] = e::marginal_cdf(on.marginal[a], truth_u[a]);
    }
    const auto off = e::grid_posterior(backends::CpuBackend{1}, rc::Objective{}, primary, obs.data(), st.T, g, L,
                                       prior, jt, se_u, false);
    row.flags_off = off.flags;
    if (!(off.flags & e::kPosteriorNumeric)) {
        for (int a = 0; a < 2; ++a) row.rank_off[a] = e::marginal_cdf(off.marginal[a], truth_u[a]);
    }
    return row;
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t N = 1000;
    std::string out_path;
    unsigned threads = std::thread::hardware_concurrency();
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--draws") == 0 && i + 1 < argc) {
            N = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else {
            std::fprintf(stderr, "usage: study_bayes_sbc [--draws N] [--out FILE] [--threads N]\n");
            return 2;
        }
    }
    if (threads == 0) threads = 1;
    const auto start = std::chrono::steady_clock::now();
    const Grid<2> g = rc::grid();
    const e::JeffreysTable jt = e::jeffreys_table(backends::CpuBackend{static_cast<int>(threads)}, rc::Objective{},
                                                  quadrature::parity_rule(), 1000, g);
    const std::vector<double> cum_flat = cumulative_prior(e::Prior::Flat, nullptr, g);
    const std::vector<double> cum_jeffreys = cumulative_prior(e::Prior::Jeffreys, &jt, g);

    struct Job {
        std::uint32_t setting, prior, draw;
    };
    std::vector<Job> jobs;
    for (std::uint32_t s = 0; s < 3; ++s) {
        for (std::uint32_t p = 0; p < (kSettings[s].jeffreys ? 2u : 1u); ++p) {
            for (std::uint32_t d = 0; d < N; ++d) jobs.push_back({s, p, d});
        }
    }
    std::vector<Row> rows(jobs.size());
    std::atomic<std::size_t> next{0};
    std::mutex progress;
    std::size_t done = 0;
    auto worker = [&] {
        for (;;) {
            const std::size_t j = next.fetch_add(1);
            if (j >= jobs.size()) return;
            const Job& jb = jobs[j];
            rows[j] = sbc_one(jb.setting, jb.prior, jb.draw, jb.prior == 0 ? cum_flat : cum_jeffreys,
                              jb.prior == 0 ? nullptr : &jt);
            const std::lock_guard<std::mutex> lock(progress);
            if (++done % 500 == 0 || done == jobs.size()) {
                std::fprintf(stderr, "draws %zu/%zu, %.0f s\n", done, jobs.size(),
                             std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    std::printf("%zu draws in %.0f s (%u threads)\n", jobs.size(),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(), threads);
    if (!out_path.empty()) {
        std::ofstream out(out_path, std::ios::binary);
        out << "setting,n,periods,prior,draw,pd,rho,flags,refinements,rank_pd,rank_rho,flags_norule,rank_pd_norule,"
               "rank_rho_norule\n";
        char line[512];
        for (const Row& r : rows) {
            const Setting& st = kSettings[r.setting];
            std::snprintf(line, sizeof line, "%u,%lld,%lld,%s,%u,%.17g,%.17g,%u,%d,%.17g,%.17g,%u,%.17g,%.17g\n",
                          r.setting, static_cast<long long>(st.n), static_cast<long long>(st.T),
                          r.prior == 0 ? "flat" : "jeffreys", r.draw, r.pd, r.rho, r.flags, r.refinements, r.rank[0],
                          r.rank[1], r.flags_off, r.rank_off[0], r.rank_off[1]);
            out << line;
        }
        if (!out) {
            std::fprintf(stderr, "cannot write %s\n", out_path.c_str());
            return 1;
        }
        std::printf("wrote %s\n", out_path.c_str());
    }
    return 0;
}
