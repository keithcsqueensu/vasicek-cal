// SPDX-License-Identifier: Apache-2.0
//
// S-9: frequentist coverage of the grid-Bayesian credible intervals (studies/bayes-coverage/PREDICTION.md,
// F1-F5; D-161, D-166). Every recovery panel is fitted as the recovery harness fits it (calibrate on
// the parity grid, the surface rows shared through the row cache, D-167; the pinned profile
// intervals), and its grid posterior (engine::grid_posterior) is computed in four arms:
//   - prior: flat on the natural scale, or Jeffreys (the table for the scenario's n on the parity
//     grid, computed once per n);
//   - the resolution rule on (the estimator), or off (the diagnostic arm: the parity grid as it stands).
// One row per (scenario, replicate) with each arm's flags, refinements and 95% equal-tailed and HPD
// interval ends (natural scale), and the pinned profile interval's ends. Nothing is scored here:
// studies/bayes-coverage/compare.py scores the predictions from these rows.
//
//   study_bayes_coverage [--replicates R] [--scenarios ID,...] [--out FILE] [--threads N]
//
// Doubles are written to 17 significant digits (an exact round trip). Results do not depend on the
// thread count.
#include <algorithm>
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
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace e = vcal::engine;

struct Arm {
    std::uint32_t flags;
    std::int32_t refinements;
    double et_lo[2], et_hi[2], hpd_lo[2], hpd_hi[2];
    double seconds;
};

struct Row {
    double pd, rho;  // the MLE (refined grid estimate), for reference
    std::uint32_t fit_flags;
    double prof_lo[2], prof_hi[2];
    std::uint32_t prof_flags[2];
    Arm arm[4];  // flat+rule, jeffreys+rule, flat no rule, jeffreys no rule
};

const char* const kArmNames[4] = {"flat", "jeffreys", "flat_norule", "jeffreys_norule"};

Row fit_one(const rc::Scenario& s, std::uint32_t replicate, e::SurfaceRowCache& cache, const e::JeffreysTable& jt) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const std::vector<std::int64_t> d = rc::panel(s, replicate);
    std::vector<rc::Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {s.obligors, d[t]};
    const Grid<2> g = rc::grid();
    std::vector<double> L;
    e::Estimate2 est{};
    if (e::calibrate_cached(backends::CpuBackend{1}, cache, rc::Objective{}, primary, check, obs.data(), s.periods, g, L,
                            est) != e::Status::Ok) {
        std::abort();
    }
    const auto prof = e::profile_intervals(rc::Objective{}, primary, obs.data(), s.periods, g, L, est);
    Row row{};
    row.pd = est.value[0];
    row.rho = est.value[1];
    row.fit_flags = est.flags;
    for (int a = 0; a < 2; ++a) {
        row.prof_lo[a] = prof.lo[a];
        row.prof_hi[a] = prof.hi[a];
        row.prof_flags[a] = prof.flags[a];
    }
    double se_u[2];
    for (int a = 0; a < 2; ++a) {
        se_u[a] = est.se[a] / grid::dvalue_dscaled(g.axis[a].scale, grid::to_scaled(g.axis[a].scale, est.value[a]));
    }
    for (int k = 0; k < 4; ++k) {
        const e::Prior prior = (k % 2 == 0) ? e::Prior::Flat : e::Prior::Jeffreys;
        const bool rule = k < 2;
        const auto t0 = std::chrono::steady_clock::now();
        const auto p = e::grid_posterior(backends::CpuBackend{1}, rc::Objective{}, primary, obs.data(), s.periods, g, L,
                                         prior, &jt, se_u, rule);
        Arm& a = row.arm[k];
        a.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        a.flags = p.flags;
        a.refinements = p.refinements;
        for (int j = 0; j < 2; ++j) {
            a.et_lo[j] = p.et_lo[j];
            a.et_hi[j] = p.et_hi[j];
            a.hpd_lo[j] = p.hpd_lo[j];
            a.hpd_hi[j] = p.hpd_hi[j];
        }
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
                if (end == c || id >= rc::kScenarios) {
                    std::fprintf(stderr, "--scenarios: ids 0..%u, comma-separated\n", rc::kScenarios - 1);
                    return 2;
                }
                ids.push_back(static_cast<std::uint32_t>(id));
                c = *end == ',' ? end + 1 : end;
            }
        } else {
            std::fprintf(stderr,
                         "usage: study_bayes_coverage [--replicates R] [--scenarios ID,...] [--out FILE] [--threads N]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    if (threads == 0) threads = 1;
    const auto start = std::chrono::steady_clock::now();

    // The Jeffreys tables, one per n on the parity grid (all threads).
    std::map<std::int64_t, e::JeffreysTable> tables;
    for (const auto id : ids) {
        const auto n = rc::scenario(id).obligors;
        if (tables.count(n) == 0) {
            tables.emplace(n, e::jeffreys_table(backends::CpuBackend{static_cast<int>(threads)}, rc::Objective{},
                                                quadrature::parity_rule(), n, rc::grid()));
            std::fprintf(stderr, "Jeffreys table n = %lld, %.0f s\n", static_cast<long long>(n),
                         std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
        }
    }

    std::vector<std::pair<std::size_t, std::uint32_t>> jobs;
    for (std::size_t k = 0; k < ids.size(); ++k) {
        for (std::uint32_t r = 0; r < R; ++r) jobs.emplace_back(k, r);
    }
    std::stable_sort(jobs.begin(), jobs.end(), [&](const auto& a, const auto& b) {
        return rc::scenario(ids[a.first]).periods > rc::scenario(ids[b.first]).periods;
    });
    e::SurfaceRowCache cache;
    std::vector<Row> rows(ids.size() * R);
    std::atomic<std::size_t> next{0};
    std::mutex progress;
    std::size_t done = 0;
    auto worker = [&] {
        for (;;) {
            const std::size_t j = next.fetch_add(1);
            if (j >= jobs.size()) return;
            const auto [k, r] = jobs[j];
            const rc::Scenario s = rc::scenario(ids[k]);
            rows[k * R + r] = fit_one(s, r, cache, tables.at(s.obligors));
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
    std::printf("%u replicates x %zu scenarios in %.0f s (%u threads)\n", R, ids.size(),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(), threads);
    for (std::size_t k = 0; k < ids.size(); ++k) {
        const rc::Scenario s = rc::scenario(ids[k]);
        double secs[4] = {};
        std::int64_t refined = 0, refused = 0;
        for (std::uint32_t r = 0; r < R; ++r) {
            const Row& w = rows[k * R + r];
            for (int a = 0; a < 4; ++a) secs[a] += w.arm[a].seconds;
            refined += (w.arm[0].flags & e::kPosteriorRefined) != 0;
            refused += (w.arm[0].flags & e::kPosteriorRefused) != 0;
        }
        std::printf("%2u T %-3lld n %-5lld | CPU s: flat %.1f jeffreys %.1f, no rule %.2f %.2f | flat+rule refined %lld "
                    "refused %lld\n",
                    s.id, static_cast<long long>(s.periods), static_cast<long long>(s.obligors), secs[0], secs[1], secs[2],
                    secs[3], static_cast<long long>(refined), static_cast<long long>(refused));
    }
    if (!out_path.empty()) {
        std::ofstream out(out_path, std::ios::binary);
        out << "scenario,replicate,pd_true,rho_true,periods,obligors,mle_pd,mle_rho,fit_flags,prof_lo_pd,prof_hi_pd,"
               "prof_lo_rho,prof_hi_rho,prof_flags_pd,prof_flags_rho";
        for (const char* a : kArmNames) {
            out << ',' << a << "_flags," << a << "_refinements";
            for (const char* p : {"pd", "rho"}) {
                out << ',' << a << "_et_lo_" << p << ',' << a << "_et_hi_" << p << ',' << a << "_hpd_lo_" << p << ','
                    << a << "_hpd_hi_" << p;
            }
        }
        out << '\n';
        char line[512];
        for (std::size_t k = 0; k < ids.size(); ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            for (std::uint32_t r = 0; r < R; ++r) {
                const Row& w = rows[k * R + r];
                std::snprintf(line, sizeof line, "%u,%u,%.17g,%.17g,%lld,%lld,%.17g,%.17g,%u,%.17g,%.17g,%.17g,%.17g,%u,%u",
                              s.id, r, s.pd, s.rho, static_cast<long long>(s.periods),
                              static_cast<long long>(s.obligors), w.pd, w.rho, w.fit_flags, w.prof_lo[0], w.prof_hi[0],
                              w.prof_lo[1], w.prof_hi[1], w.prof_flags[0], w.prof_flags[1]);
                out << line;
                for (const Arm& a : w.arm) {
                    std::snprintf(line, sizeof line, ",%u,%d", a.flags, a.refinements);
                    out << line;
                    for (int j = 0; j < 2; ++j) {
                        std::snprintf(line, sizeof line, ",%.17g,%.17g,%.17g,%.17g", a.et_lo[j], a.et_hi[j], a.hpd_lo[j],
                                      a.hpd_hi[j]);
                        out << line;
                    }
                }
                out << '\n';
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
