// SPDX-License-Identifier: Apache-2.0
//
// P-9's measurement (perf/p9_repeat_rate.md): how often a run fits a panel it has already fitted. P-9
// would memoise fits of identical panels; it is built only if the repeat rate projects a real saving.
// Nothing is fitted beyond what is needed to regenerate the panels, and nothing here changes a result.
//
// Two populations of panels, regenerated exactly as the runs draw them:
//   - outer: the recovery panels (kSeed, scenario, replicate r = 0..999; D-109), fitted by the recovery
//     harness and by S-10/S-4;
//   - boot: S-10's parametric-bootstrap panels (S10PBOOT, scenario, replicate 1,000 r + b, b = 0..998),
//     simulated at the outer panel's estimate theta-hat_r, which this tool recomputes with the recovery
//     fit (estimate only; bit for bit the S-10 fit's estimate).
//
// A panel is identified by a 64-bit FNV-1a hash of its default counts (n is the same in every period
// and every panel of a (T, n) cell): ordered, the counts in period order; sorted, the counts sorted, the
// key the P-9 entry in STATE.md names. A fit of the ordered panel is bit for bit a fit of the same panel;
// a fit keyed by the sorted panel is not, since the log-likelihood sums periods in order. The expected
// number of 64-bit collisions among the at most 9 x 10^6 panels of a cell is below 3e-6.
//
// The repeats counted (a repeat is a panel identical to one counted earlier in the same scope):
//   - outer, per scenario and per (T, n) cell: an outer fit's estimate, SE, profile and q interval depend
//     only on the panel; its iid-bootstrap arm draws from a stream seeded by (scenario, replicate) and
//     so cannot be shared;
//   - boot, per replicate: the whole bootstrap fit (estimate, SE, and S-4b's W* at theta-hat_r for
//     b < 199), all of which a memo could share within the replicate;
//   - boot, per scenario and per cell: the estimate and SE only (W* depends on theta-hat_r).
//
//   p9_repeat_rate [--scenarios ID,...] [--boot-replicates R] [--boot-max-obligors N] [--threads N]
//                  [--out FILE] [--theta-out FILE]
//
// --boot-replicates R (default 1000) simulates the bootstrap panels of replicates 0..R-1;
// --boot-max-obligors N (default 10000) skips the bootstrap panels of scenarios with n > N (the
// simulation's cost grows with T x n; the n = 10,000 cells hold most of it). Results do not depend on the
// thread count.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <thread>
#include <vector>

#include "dgp/dgp.hpp"
#include "tests/harness/dgp_provenance.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;

inline constexpr std::uint64_t kBootSeed = 0x53313050424F4F54ull;  // "S10PBOOT", as S-10's tool
inline constexpr std::uint32_t kB = 999;
inline constexpr std::uint32_t kBartlettPanels = 199;

std::uint64_t fnv(std::vector<std::int64_t> d, bool sorted) {
    if (sorted) std::sort(d.begin(), d.end());
    std::uint64_t h = 1469598103934665603ull;
    for (const auto x : d) {
        for (int k = 0; k < 8; ++k) {
            h ^= static_cast<std::uint64_t>(x >> (8 * k)) & 0xFFu;
            h *= 1099511628211ull;
        }
    }
    return h;
}

// Panels minus distinct panels.
std::uint64_t repeats(std::vector<std::uint64_t> h) {
    std::sort(h.begin(), h.end());
    return h.size() - static_cast<std::uint64_t>(std::unique(h.begin(), h.end()) - h.begin());
}

template <class F>
void parallel(unsigned threads, std::uint32_t count, F f) {
    std::atomic<std::uint32_t> next{0};
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t) {
        pool.emplace_back([&] {
            for (std::uint32_t i; (i = next.fetch_add(1)) < count;) f(i);
        });
    }
    for (auto& t : pool) t.join();
}

struct Hashes {
    std::vector<std::uint64_t> h[2];  // ordered, sorted
};

struct Result {
    rc::Scenario s;
    std::uint64_t outer_rep[2] = {};             // per scenario
    bool boot = false;
    std::uint32_t boot_reps = 0;
    std::uint64_t boot_panels = 0;
    std::uint64_t boot_rep_within[2] = {};       // per replicate, all 999
    std::uint64_t boot_rep_within199[2] = {};    // per replicate, the first 199
    std::uint64_t boot_rep_scenario[2] = {};     // per scenario
    std::uint64_t boot_all_zero = 0;             // bootstrap panels with no defaults
    double seconds = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
    if (!vcal::test::announce_dgp_kernel("p9_repeat_rate")) return 2;
    std::vector<std::uint32_t> ids;
    std::uint32_t boot_r = rc::kReplicates;
    std::int64_t boot_max_n = 10000;
    unsigned threads = std::thread::hardware_concurrency();
    const char* out = "p9_repeat_rate.csv";
    const char* theta_out = nullptr;  // the estimates the bootstrap panels were simulated at, for a check
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            for (char* p = argv[++i]; *p;) {
                ids.push_back(static_cast<std::uint32_t>(std::strtoul(p, &p, 10)));
                if (*p == ',') ++p;
            }
        } else if (std::strcmp(argv[i], "--boot-replicates") == 0 && i + 1 < argc) {
            boot_r = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--boot-max-obligors") == 0 && i + 1 < argc) {
            boot_max_n = std::strtoll(argv[++i], nullptr, 10);
        } else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            threads = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out = argv[++i];
        } else if (std::strcmp(argv[i], "--theta-out") == 0 && i + 1 < argc) {
            theta_out = argv[++i];
        } else {
            std::fprintf(stderr, "usage: p9_repeat_rate [--scenarios ID,...] [--boot-replicates R] "
                                 "[--boot-max-obligors N] [--threads N] [--out FILE] [--theta-out FILE]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    if (threads == 0) threads = 1;
    if (boot_r > rc::kReplicates) boot_r = rc::kReplicates;

    engine::SurfaceRowCache cache;
    std::vector<Result> results;
    std::string thetas = "scenario,replicate,pd,rho\n";
    // Per (T, n) cell: every scenario's outer hashes and bootstrap hashes, for the cell scope.
    std::map<std::pair<std::int64_t, std::int64_t>, Hashes> cell_outer, cell_boot;
    for (const std::uint32_t id : ids) {
        const auto t0 = std::chrono::steady_clock::now();
        Result res;
        res.s = rc::scenario(id);
        const rc::Scenario& s = res.s;
        const auto T = static_cast<std::size_t>(s.periods);
        const std::uint32_t R = rc::kReplicates;
        std::vector<std::uint64_t> outer[2] = {std::vector<std::uint64_t>(R), std::vector<std::uint64_t>(R)};
        std::vector<double> theta(2 * static_cast<std::size_t>(R));
        res.boot = s.obligors <= boot_max_n && boot_r > 0;
        res.boot_reps = res.boot ? boot_r : 0;
        std::vector<std::uint64_t> boot[2];
        for (auto& v : boot) v.resize(static_cast<std::size_t>(res.boot_reps) * kB);
        std::vector<std::uint64_t> within[2] = {std::vector<std::uint64_t>(res.boot_reps), std::vector<std::uint64_t>(res.boot_reps)};
        std::vector<std::uint64_t> within199[2] = {within[0], within[1]};
        std::atomic<std::uint64_t> all_zero{0};
        parallel(threads, R, [&](std::uint32_t r) {
            const auto d = rc::panel(s, r);
            for (int k = 0; k < 2; ++k) outer[k][r] = fnv(d, k == 1);
            if (r >= res.boot_reps) return;
            const rc::Fit f = rc::fit(s, r, &cache, 0);  // the estimate S-10 simulates at
            theta[2 * r] = f.value[0];
            theta[2 * r + 1] = f.value[1];
            std::vector<std::int64_t> n(T, s.obligors), db(T);
            std::vector<std::uint64_t> mine[2];
            std::uint64_t zero = 0;
            for (std::uint32_t b = 0; b < kB; ++b) {
                const dgp::PanelSpec spec{kBootSeed, s.id, 1000u * r + b, f.value[0], f.value[1], n.data(), s.periods};
                if (dgp::simulate_panel(spec, db.data(), nullptr) != dgp::Status::Ok) std::abort();
                bool z = true;
                for (const auto x : db) z = z && x == 0;
                zero += z ? 1 : 0;
                for (int k = 0; k < 2; ++k) {
                    const std::uint64_t h = fnv(db, k == 1);
                    boot[k][static_cast<std::size_t>(r) * kB + b] = h;
                    mine[k].push_back(h);
                }
            }
            all_zero += zero;
            for (int k = 0; k < 2; ++k) {
                within[k][r] = repeats(mine[k]);
                within199[k][r] = repeats(std::vector<std::uint64_t>(mine[k].begin(), mine[k].begin() + kBartlettPanels));
            }
        });
        const auto cell = std::make_pair(s.periods, s.obligors);
        for (int k = 0; k < 2; ++k) {
            res.outer_rep[k] = repeats(outer[k]);
            auto& co = cell_outer[cell].h[k];
            co.insert(co.end(), outer[k].begin(), outer[k].end());
            if (!res.boot) continue;
            for (std::uint32_t r = 0; r < res.boot_reps; ++r) {
                res.boot_rep_within[k] += within[k][r];
                res.boot_rep_within199[k] += within199[k][r];
            }
            res.boot_rep_scenario[k] = repeats(boot[k]);
            auto& cb = cell_boot[cell].h[k];
            cb.insert(cb.end(), boot[k].begin(), boot[k].end());
        }
        for (std::uint32_t r = 0; r < res.boot_reps; ++r) {
            char buf[96];
            std::snprintf(buf, sizeof buf, "%u,%u,%.17g,%.17g\n", id, r, theta[2 * r], theta[2 * r + 1]);
            thetas += buf;
        }
        res.boot_panels = static_cast<std::uint64_t>(res.boot_reps) * kB;
        res.boot_all_zero = all_zero;
        res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::fprintf(stderr, "scenario %2u (T %3lld, n %5lld, PD %.3f, rho %.2f): outer repeats %llu/%llu, boot %llu "
                             "panels, within-replicate %llu/%llu, within-scenario %llu/%llu (ordered/sorted); %.0f s\n",
                     id, static_cast<long long>(s.periods), static_cast<long long>(s.obligors), s.pd, s.rho,
                     static_cast<unsigned long long>(res.outer_rep[0]), static_cast<unsigned long long>(res.outer_rep[1]),
                     static_cast<unsigned long long>(res.boot_panels),
                     static_cast<unsigned long long>(res.boot_rep_within[0]),
                     static_cast<unsigned long long>(res.boot_rep_within[1]),
                     static_cast<unsigned long long>(res.boot_rep_scenario[0]),
                     static_cast<unsigned long long>(res.boot_rep_scenario[1]), res.seconds);
        results.push_back(res);
    }

    std::string text = "# P-9 repeat rate (perf/p9_repeat_rate.md): repeated panels per scope, ordered and sorted keys\n"
                       "scope,scenario,periods,obligors,pd,rho,panels,repeats_ordered,repeats_sorted,all_zero\n";
    char line[256];
    for (const auto& r : results) {
        const auto& s = r.s;
        const auto row = [&](const char* scope, std::uint64_t panels, const std::uint64_t* rep, std::uint64_t zero) {
            std::snprintf(line, sizeof line, "%s,%u,%lld,%lld,%.17g,%.17g,%llu,%llu,%llu,%llu\n", scope, s.id,
                          static_cast<long long>(s.periods), static_cast<long long>(s.obligors), s.pd, s.rho,
                          static_cast<unsigned long long>(panels), static_cast<unsigned long long>(rep[0]),
                          static_cast<unsigned long long>(rep[1]), static_cast<unsigned long long>(zero));
            text += line;
        };
        row("outer_scenario", rc::kReplicates, r.outer_rep, 0);
        if (!r.boot) continue;
        row("boot_replicate", r.boot_panels, r.boot_rep_within, r.boot_all_zero);
        const std::uint64_t p199 = static_cast<std::uint64_t>(r.boot_reps) * kBartlettPanels;
        row("boot_replicate_199", p199, r.boot_rep_within199, 0);
        row("boot_scenario", r.boot_panels, r.boot_rep_scenario, r.boot_all_zero);
    }
    for (const auto* m : {&cell_outer, &cell_boot}) {
        for (const auto& [cell, hs] : *m) {
            if (hs.h[0].empty()) continue;
            std::snprintf(line, sizeof line, "%s,,%lld,%lld,,,%zu,%llu,%llu,\n", m == &cell_outer ? "outer_cell" : "boot_cell",
                          static_cast<long long>(cell.first), static_cast<long long>(cell.second), hs.h[0].size(),
                          static_cast<unsigned long long>(repeats(hs.h[0])),
                          static_cast<unsigned long long>(repeats(hs.h[1])));
            text += line;
        }
    }
    std::ofstream f(out, std::ios::binary);
    f << text;
    f.close();
    if (!f) {
        std::fprintf(stderr, "cannot write %s\n", out);
        return 1;
    }
    if (theta_out != nullptr) {
        std::ofstream g(theta_out, std::ios::binary);
        g << thetas;
        g.close();
        if (!g) return 1;
    }
    return 0;
}
