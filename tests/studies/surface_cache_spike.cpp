// SPDX-License-Identifier: Apache-2.0
//
// Spike: caching per-period surface rows by observation (n, d) across the replicates of a scenario,
// extending D-122 (equal observations share one row within a panel) across panels. Not a study and not
// merged without review.
//
// For each scenario of the chosen set:
//   1. simulate every replicate's panel (the recovery DGP; recovery::panel);
//   2. UNCACHED: fit every replicate as recovery::fit does, instrumented by phase: the surface
//      (evaluate_surface, D-122 within the panel), the rest of calibrate (argmax, refinement,
//      quadrature check, Hessian), the profile intervals, the bootstrap (weights, W x L reduction,
//      percentile ends), and S-23's q interval and q bootstrap ends;
//   3. CACHED: compute each distinct (n, d) of the scenario's panels once, as a row over the grid
//      (log_contrib at every grid point, exactly as evaluate_surface computes it), then fit every
//      replicate from its rows copied out of the cache (engine::calibrate_from_surface), the other
//      phases unchanged;
//   4. compare every field of every replicate's Fit, cached and uncached, with recovery::fit: bit for
//      bit (the padding bytes excluded).
// Reports surface rows computed, and per phase the CPU time summed over replicates (each fit is
// serial) and the wall-clock time of the whole pass.
//
//   surface_cache_spike [--replicates R] [--scenarios ID,...]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
using Clock = std::chrono::steady_clock;
using Objective = rc::Objective;

struct Phases {
    double surface = 0, calibrate = 0, profile = 0, bootstrap = 0, q = 0;
    Phases& operator+=(const Phases& o) {
        surface += o.surface;
        calibrate += o.calibrate;
        profile += o.profile;
        bootstrap += o.bootstrap;
        q += o.q;
        return *this;
    }
    double total() const { return surface + calibrate + profile + bootstrap + q; }
};

double since(Clock::time_point& t0) {
    const auto t1 = Clock::now();
    const double s = std::chrono::duration<double>(t1 - t0).count();
    t0 = t1;
    return s;
}

using RowCache = std::map<std::int64_t, std::vector<double>>;  // d -> K values (n is the scenario's)

// recovery::fit, phase by phase; with cache != nullptr the surface is assembled from cached rows.
rc::Fit fit_timed(const rc::Scenario& s, std::uint32_t replicate, const std::vector<std::int64_t>& d,
                  const RowCache* cache, Phases& ph, std::int64_t& rows) {
    std::vector<Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {s.obligors, d[t]};
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    const std::int64_t K = g.size(), T = s.periods;
    const backends::CpuBackend serial{1};
    auto t0 = Clock::now();
    std::vector<double> L(static_cast<std::size_t>(T * K));
    if (cache == nullptr) {
        engine::evaluate_surface(serial, Objective{}, primary, obs.data(), T, g, L.data());
        rows += static_cast<std::int64_t>(std::set<std::int64_t>(d.begin(), d.end()).size());
    } else {
        for (std::int64_t t = 0; t < T; ++t) {
            const auto& row = cache->at(d[static_cast<std::size_t>(t)]);
            std::copy(row.begin(), row.end(), L.begin() + t * K);
        }
    }
    ph.surface += since(t0);
    engine::Estimate2 est{};
    if (engine::calibrate_from_surface(serial, Objective{}, primary, check, obs.data(), T, g, L, est) != engine::Status::Ok) {
        std::abort();
    }
    ph.calibrate += since(t0);
    const auto prof = engine::profile_intervals(Objective{}, primary, obs.data(), T, g, L, est);
    ph.profile += since(t0);
    const auto idx = resample::bootstrap_indices(rc::bootstrap_seed(s.id, replicate), resample::Scheme::IidBootstrap,
                                                 rc::kBootstrapReplicates, T);
    const auto W = resample::weights_from_indices(idx.data(), rc::kBootstrapReplicates, T, T);
    std::vector<resample::Replicate2> reps(rc::kBootstrapReplicates);
    resample::replicate_estimates(serial, g, L.data(), T, W.data(), rc::kBootstrapReplicates, reps.data());
    std::uint32_t edge = 0;
    for (const auto& r : reps) edge += (r.flags & engine::kFlagGridEdge) ? 1u : 0u;
    const auto b0 = resample::percentile_interval(reps.data(), rc::kBootstrapReplicates, 0);
    const auto b1 = resample::percentile_interval(reps.data(), rc::kBootstrapReplicates, 1);
    ph.bootstrap += since(t0);
    const auto qi = engine::conditional_pd_interval(Objective{}, primary, obs.data(), T, g, L, est, prof);
    const auto qg = engine::conditional_pd_logit_gradient(est.value[0], est.value[1]);
    double se_u[2];
    for (int a = 0; a < 2; ++a) {
        se_u[a] = est.se[a] / grid::dvalue_dscaled(g.axis[a].scale, grid::to_scaled(g.axis[a].scale, est.value[a]));
    }
    const double q_var = qg.ds_du[0] * qg.ds_du[0] * se_u[0] * se_u[0] + qg.ds_du[1] * qg.ds_du[1] * se_u[1] * se_u[1] +
                         2.0 * qg.ds_du[0] * qg.ds_du[1] * est.corr * se_u[0] * se_u[1];
    std::vector<double> q_reps;
    q_reps.reserve(reps.size());
    for (const auto& r : reps) {
        const double q = std::isfinite(r.value[0]) && std::isfinite(r.value[1]) ? engine::conditional_pd(r.value[0], r.value[1])
                                                                                : std::nan("");
        if (std::isfinite(q)) q_reps.push_back(q);
    }
    std::sort(q_reps.begin(), q_reps.end());
    const double nan = std::nan("");
    const double q_boot_lo = q_reps.empty() ? nan : resample::quantile_type7(q_reps, 0.025);
    const double q_boot_hi = q_reps.empty() ? nan : resample::quantile_type7(q_reps, 0.975);
    ph.q += since(t0);
    return {{est.value[0], est.value[1]},
            {est.se[0], est.se[1]},
            est.loglik,
            est.flags,
            {prof.lo[0], prof.lo[1]},
            {prof.hi[0], prof.hi[1]},
            {prof.flags[0], prof.flags[1]},
            prof.residual_max,
            {b0.lo, b1.lo},
            {b0.hi, b1.hi},
            edge,
            static_cast<std::uint32_t>(b0.excluded > b1.excluded ? b0.excluded : b1.excluded),
            qi.estimate,
            q_var >= 0.0 ? std::sqrt(q_var) : nan,
            qi.lo,
            qi.hi,
            qi.flags,
            qi.residual_max,
            q_boot_lo,
            q_boot_hi};
}

template <class T>
bool same(const T& a, const T& b) {
    return std::memcmp(&a, &b, sizeof(T)) == 0;
}

// Field by field (Fit has padding, whose bytes are unspecified).
bool same_fit(const rc::Fit& a, const rc::Fit& b) {
    return same(a.value, b.value) && same(a.se, b.se) && same(a.loglik, b.loglik) && a.flags == b.flags &&
           same(a.prof_lo, b.prof_lo) && same(a.prof_hi, b.prof_hi) && same(a.prof_flags, b.prof_flags) &&
           same(a.prof_residual, b.prof_residual) && same(a.boot_lo, b.boot_lo) && same(a.boot_hi, b.boot_hi) &&
           a.boot_edge == b.boot_edge && a.boot_excluded == b.boot_excluded && same(a.q_hat, b.q_hat) &&
           same(a.q_se_s, b.q_se_s) && same(a.q_prof_lo, b.q_prof_lo) && same(a.q_prof_hi, b.q_prof_hi) &&
           a.q_prof_flags == b.q_prof_flags && same(a.q_prof_residual, b.q_prof_residual) &&
           same(a.q_boot_lo, b.q_boot_lo) && same(a.q_boot_hi, b.q_boot_hi);
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids = {29, 37, 72, 4, 68, 49, 7, 43, 51};  // the study subset (D-150)
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            ids.clear();
            for (const char* c = argv[++i]; *c != '\0';) {
                char* end = nullptr;
                ids.push_back(static_cast<std::uint32_t>(std::strtoul(c, &end, 10)));
                c = *end == ',' ? end + 1 : end;
            }
        } else {
            std::fprintf(stderr, "usage: surface_cache_spike [--replicates R] [--scenarios ID,...]\n");
            return 2;
        }
    }
    const backends::CpuBackend all{};
    static const auto primary = quadrature::parity_rule();
    const Grid<2> g = rc::grid();
    const std::int64_t K = g.size();
    Phases un_total, ca_total;
    double un_wall = 0, ca_wall = 0, build_wall = 0, build_cpu = 0;
    std::int64_t un_rows = 0, ca_rows = 0, mismatches = 0, compared = 0;
    std::printf("scenario  T    n     | rows uncached  cached | CPU s uncached: surface  rest | cached: build  rest | "
                "identical\n");
    for (const std::uint32_t id : ids) {
        const rc::Scenario s = rc::scenario(id);
        std::vector<std::vector<std::int64_t>> panels(R);
        for (std::uint32_t r = 0; r < R; ++r) panels[r] = rc::panel(s, r);
        // UNCACHED.
        std::vector<Phases> ph(R);
        std::vector<std::int64_t> rows(R, 0);
        std::vector<rc::Fit> un(R), ca(R), ref(R);
        auto t0 = Clock::now();
        all.parallel_for(R, [&](std::int64_t r) {
            un[static_cast<std::size_t>(r)] = fit_timed(s, static_cast<std::uint32_t>(r), panels[static_cast<std::size_t>(r)],
                                                        nullptr, ph[static_cast<std::size_t>(r)], rows[static_cast<std::size_t>(r)]);
        });
        const double w_un = since(t0);
        Phases p_un;
        std::int64_t r_un = 0;
        for (std::uint32_t r = 0; r < R; ++r) {
            p_un += ph[r];
            r_un += rows[r];
        }
        // CACHED: the scenario's distinct counts, each row once (parallel over rows x grid points).
        std::set<std::int64_t> distinct;
        for (const auto& p : panels) distinct.insert(p.begin(), p.end());
        const std::vector<std::int64_t> ds(distinct.begin(), distinct.end());
        RowCache cache;
        for (const auto d : ds) cache[d].assign(static_cast<std::size_t>(K), 0.0);
        t0 = Clock::now();
        std::vector<double> cpu(static_cast<std::size_t>(ds.size() * static_cast<std::size_t>(K)));
        all.parallel_for(static_cast<std::int64_t>(ds.size()) * K, [&](std::int64_t j) {
            const auto t1 = Clock::now();
            const std::int64_t d = ds[static_cast<std::size_t>(j / K)], k = j % K;
            double v[2];
            g.values(k, v);
            cache[d][static_cast<std::size_t>(k)] = Objective{}.log_contrib({s.obligors, d}, Objective::theta(v), primary);
            cpu[static_cast<std::size_t>(j)] = std::chrono::duration<double>(Clock::now() - t1).count();
        });
        const double w_build = since(t0);
        double c_build = 0;
        for (const double c : cpu) c_build += c;
        std::vector<Phases> ph2(R);
        std::vector<std::int64_t> rows2(R, 0);
        all.parallel_for(R, [&](std::int64_t r) {
            ca[static_cast<std::size_t>(r)] = fit_timed(s, static_cast<std::uint32_t>(r), panels[static_cast<std::size_t>(r)],
                                                        &cache, ph2[static_cast<std::size_t>(r)], rows2[static_cast<std::size_t>(r)]);
        });
        const double w_ca = since(t0);
        Phases p_ca;
        for (std::uint32_t r = 0; r < R; ++r) p_ca += ph2[r];
        // The reference: the library's own recovery::fit.
        all.parallel_for(R, [&](std::int64_t r) { ref[static_cast<std::size_t>(r)] = rc::fit(s, static_cast<std::uint32_t>(r)); });
        std::int64_t bad = 0;
        for (std::uint32_t r = 0; r < R; ++r) {
            bad += !same_fit(ca[r], ref[r]);
            bad += !same_fit(un[r], ref[r]);
        }
        compared += 2 * R;
        mismatches += bad;
        un_total += p_un;
        ca_total += p_ca;
        un_wall += w_un;
        ca_wall += w_ca;
        build_wall += w_build;
        build_cpu += c_build;
        un_rows += r_un;
        ca_rows += static_cast<std::int64_t>(ds.size());
        std::printf("%-8u  %-4lld %-5lld | %8lld  %8zu | %22.1f %5.1f | %13.2f %5.1f | %s\n", id,
                    static_cast<long long>(s.periods), static_cast<long long>(s.obligors), static_cast<long long>(r_un),
                    ds.size(), p_un.surface, p_un.total() - p_un.surface, c_build, p_ca.total() - p_ca.surface,
                    bad == 0 ? "yes" : "NO");
    }
    const auto line = [](const char* name, const Phases& p) {
        std::printf("  %-9s surface %8.1f  calibrate %7.1f  profile %7.1f  bootstrap %7.1f  q %7.1f  total %8.1f\n", name,
                    p.surface, p.calibrate, p.profile, p.bootstrap, p.q, p.total());
    };
    std::printf("\nsurface rows computed: uncached %lld, cached %lld (%.1fx fewer)\n", static_cast<long long>(un_rows),
                static_cast<long long>(ca_rows), static_cast<double>(un_rows) / static_cast<double>(ca_rows));
    std::printf("CPU seconds summed over replicates (each fit serial):\n");
    line("uncached", un_total);
    line("cached", ca_total);
    std::printf("  cache build (CPU) %.2f; cached total incl. build %.1f (%.2fx the uncached)\n", build_cpu,
                ca_total.total() + build_cpu, (ca_total.total() + build_cpu) / un_total.total());
    std::printf("wall clock, all threads: uncached %.1f s; cached %.1f s (build %.1f s + fits %.1f s)\n", un_wall,
                build_wall + ca_wall, build_wall, ca_wall);
    std::printf("fits compared with recovery::fit: %lld, differing: %lld\n", static_cast<long long>(compared),
                static_cast<long long>(mismatches));
    return mismatches == 0 ? 0 : 1;
}
