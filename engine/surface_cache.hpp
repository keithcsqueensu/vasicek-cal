// SPDX-License-Identifier: Apache-2.0
//
// Surface rows shared across panels (D-167). A row of the per-period surface, l(obs, theta_k) at every
// grid point k, depends on nothing but the observation and the configuration that evaluates it: the
// objective (its type, which carries its precision policy and, for the rate objective, its zero-rate
// treatment; and its state, if it has any), the integrator (its type and every table and count, N
// included), and the grid (each axis's bounds, point count and scale). D-122 shares rows between equal
// observations within a panel; this cache shares them across panels: replicates, scenarios with the
// same n, bootstrap panels.
//
// The key is (configuration, observation), with the configuration spelled out field by field, so no
// configuration change can hit a row computed under another (tested: unit_surface_cache). Rows are
// filled on demand, the first time an observation appears, so the range of counts need not be known
// in advance (large n, S-27). A row is computed exactly as evaluate_surface computes it, so a surface
// assembled from the cache is identical to evaluate_surface's bit for bit. Thread-safe: concurrent
// requests for one row compute it once.
#pragma once

#include <atomic>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include "core/grid.hpp"
#include "engine/argmax_tiles.hpp"
#include "core/quadrature/composite_legendre.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "core/quadrature/integrator.hpp"

namespace vcal::engine {

namespace cache_detail {

template <class T>
void put(std::string& k, const T& v) {
    static_assert(std::is_arithmetic_v<T> || std::is_pointer_v<T> || std::is_enum_v<T>, "put scalars only");
    k.append(reinterpret_cast<const char*>(&v), sizeof v);
}

// Integrator signatures: the type's tag and every field that changes its values (the tables'
// addresses identify them within the process; the counts are explicit).
inline void signature(std::string& k, const quadrature::GaussHermiteRule& r) {
    put(k, r.node);
    put(k, r.log_weight);
    put(k, r.n);
}
inline void signature(std::string& k, const quadrature::GaussLegendreRule& r) {
    put(k, r.node);
    put(k, r.weight);
    put(k, r.n);
}
template <class P>
void signature(std::string& k, const quadrature::GaussHermiteFixed<P>& i) {
    k += "GHF";
    signature(k, i.rule);
}
template <class P>
void signature(std::string& k, const quadrature::GaussHermiteAdaptive<P>& i) {
    k += "GHA";
    signature(k, i.rule);
}
template <class P>
void signature(std::string& k, const quadrature::CompositeLegendre<P>& i) {
    k += "CGL";
    signature(k, i.rule);
    put(k, i.panels);
}
template <class P>
void signature(std::string& k, const quadrature::SplitRule<P>& i) {
    k += "SPL";
    signature(k, i.interior);
    signature(k, i.one_sided);
}

template <class O, class = void>
struct has_cache_signature : std::false_type {};
template <class O>
struct has_cache_signature<O, std::void_t<decltype(std::declval<const O&>().cache_signature(std::declval<std::string&>()))>>
    : std::true_type {};

}  // namespace cache_detail

// The configuration part of a row's key.
template <class Objective, class Integrator, int D>
std::string surface_configuration(const Objective& objective, const Integrator& integrator, const Grid<D>& grid) {
    namespace cd = cache_detail;
    std::string k;
    k += typeid(Objective).name();  // the type: its precision policy and any template treatment
    k += '|';
    if constexpr (cd::has_cache_signature<Objective>::value) {
        objective.cache_signature(k);  // an objective with state spells it out
    } else {
        static_assert(std::is_empty_v<Objective>,
                      "vcal: an objective with state must declare cache_signature(std::string&) to be cached");
    }
    k += '|';
    k += typeid(Integrator).name();
    cd::signature(k, integrator);
    k += '|';
    for (int a = 0; a < D; ++a) {
        cd::put(k, grid.axis[a].lo);
        cd::put(k, grid.axis[a].hi);
        cd::put(k, grid.axis[a].n);
        cd::put(k, grid.axis[a].scale);
    }
    return k;
}

class SurfaceRowCache {
public:
    // Fills by the calling thread, across every cache (D-193): the seconds spent computing rows on first
    // request, and how many. Timing instrumentation for runtime projections, which count a fill once
    // rather than per replicate; no result depends on it.
    struct FillTally {
        double seconds = 0.0;
        std::int64_t rows = 0;
    };
    static FillTally& thread_fills() {
        thread_local FillTally tally;
        return tally;
    }

    // A cached row with what is computed from it once, at fill time (P-27, D-410): on a two-dimensional grid its
    // tile maxima for the bounded argmax (engine/argmax_tiles.hpp), and whether every value is finite. Never
    // changed after the fill, and never moved (entries are held by pointer), so callers may keep pointers into it.
    struct CachedRow {
        std::vector<double> row;
        std::vector<double> tile_max;  // empty unless D == 2
        bool finite = false;
    };

    // The row for observation y under this configuration (config = surface_configuration(...)),
    // computed on first request with the backend's parallel_for over grid points.
    template <class Backend, class Objective, class Integrator, int D>
    const std::vector<double>& row(const Backend& backend, const std::string& config, const Objective& objective,
                                   const Integrator& integrator, const Grid<D>& grid, const typename Objective::Obs& y) {
        return entry(backend, config, objective, integrator, grid, y).row;
    }

    // The same row with its tile maxima and finiteness (P-27).
    template <class Backend, class Objective, class Integrator, int D>
    const CachedRow& entry(const Backend& backend, const std::string& config, const Objective& objective,
                           const Integrator& integrator, const Grid<D>& grid, const typename Objective::Obs& y) {
        using Obs = typename Objective::Obs;
        static_assert(std::is_trivially_copyable_v<Obs> && sizeof(Obs) % 8 == 0,
                      "vcal: an observation must be trivially copyable, with 8-byte fields only, to key the cache");
        std::string key = config;
        key += '#';
        key.append(reinterpret_cast<const char*>(&y), sizeof(Obs));
        Entry* e = nullptr;
        {
            const std::lock_guard<std::mutex> lock(mutex_);
            auto& slot = rows_[key];
            if (!slot) slot = std::make_unique<Entry>();
            e = slot.get();
        }
        bool computed = false;
        std::call_once(e->once, [&] {
            const auto t0 = std::chrono::steady_clock::now();
            const std::int64_t K = grid.size();
            e->row.assign(static_cast<std::size_t>(K), 0.0);
            backend.parallel_for(K, [&](std::int64_t k) {
                double v[D];
                grid.values(k, v);
                e->row[static_cast<std::size_t>(k)] = objective.log_contrib(y, Objective::theta(v), integrator);
            });
            e->finite = true;
            for (const double x : e->row) e->finite = e->finite && std::isfinite(x);
            if constexpr (D == 2) {
                e->tile_max.resize(static_cast<std::size_t>(argmax_tiles(grid)));
                tile_maxima(grid, e->row.data(), e->tile_max.data());
            }
            computed = true;
            FillTally& tally = thread_fills();
            tally.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            ++tally.rows;
        });
        (computed ? misses_ : hits_).fetch_add(1, std::memory_order_relaxed);
        return *e;
    }

    std::int64_t hits() const { return hits_.load(); }
    std::int64_t misses() const { return misses_.load(); }
    std::size_t size() const {
        const std::lock_guard<std::mutex> lock(mutex_);
        return rows_.size();
    }

private:
    struct Entry : CachedRow {
        std::once_flag once;
    };
    mutable std::mutex mutex_;
    std::unordered_map<std::string, std::unique_ptr<Entry>> rows_;
    std::atomic<std::int64_t> hits_{0}, misses_{0};
};

// evaluate_surface through the cache: L[t][k] for every period, each distinct observation's row
// taken from the cache (computed on first request). Identical to evaluate_surface bit for bit. Writes every
// element of L (periods x K), whatever it held before, so callers need not clear it (P-26, D-409).
template <class Backend, class Objective, class Integrator, int D>
void evaluate_surface_cached(const Backend& backend, SurfaceRowCache& cache, const Objective& objective,
                             const Integrator& integrator, const typename Objective::Obs* obs, std::int64_t periods,
                             const Grid<D>& grid, double* L) {
    const std::string config = surface_configuration(objective, integrator, grid);
    const std::int64_t K = grid.size();
    for (std::int64_t t = 0; t < periods; ++t) {
        std::int64_t u = 0;
        while (!(obs[u] == obs[t])) ++u;
        if (u < t) {
            std::memcpy(L + t * K, L + u * K, static_cast<std::size_t>(K) * sizeof(double));
            continue;
        }
        const auto& r = cache.row(backend, config, objective, integrator, grid, obs[t]);
        std::memcpy(L + t * K, r.data(), static_cast<std::size_t>(K) * sizeof(double));
    }
}

// A panel's surface as pointers into the row cache instead of a T x K copy (P-27, D-410): row[t] is
// L[t][.] as evaluate_surface_cached would write it, tile_max[t] that row's tile maxima, and finite whether every
// value of every row is finite. Valid while the cache lives.
struct SurfaceRows {
    std::vector<const double*> row;
    std::vector<const double*> tile_max;
    bool finite = true;
    std::int64_t periods = 0;
    std::int64_t K = 0;
};

// The rows of evaluate_surface_cached, as pointers: the same lookups (one per distinct observation, in period
// order), so the cache's fills, hits and misses are as before.
template <class Backend, class Objective, class Integrator>
void surface_rows_cached(const Backend& backend, SurfaceRowCache& cache, const Objective& objective,
                         const Integrator& integrator, const typename Objective::Obs* obs, std::int64_t periods,
                         const Grid<2>& grid, SurfaceRows& out) {
    const std::string config = surface_configuration(objective, integrator, grid);
    out.periods = periods;
    out.K = grid.size();
    out.finite = true;
    out.row.resize(static_cast<std::size_t>(periods));
    out.tile_max.resize(static_cast<std::size_t>(periods));
    for (std::int64_t t = 0; t < periods; ++t) {
        std::int64_t u = 0;
        while (!(obs[u] == obs[t])) ++u;
        if (u < t) {
            out.row[static_cast<std::size_t>(t)] = out.row[static_cast<std::size_t>(u)];
            out.tile_max[static_cast<std::size_t>(t)] = out.tile_max[static_cast<std::size_t>(u)];
            continue;
        }
        const auto& e = cache.entry(backend, config, objective, integrator, grid, obs[t]);
        out.row[static_cast<std::size_t>(t)] = e.row.data();
        out.tile_max[static_cast<std::size_t>(t)] = e.tile_max.data();
        out.finite = out.finite && e.finite;
    }
}

// The contiguous surface the rows stand for, as evaluate_surface_cached writes it (L holds periods x K values).
inline void copy_rows(const SurfaceRows& s, double* L) {
    for (std::int64_t t = 0; t < s.periods; ++t) {
        std::memcpy(L + t * s.K, s.row[static_cast<std::size_t>(t)], static_cast<std::size_t>(s.K) * sizeof(double));
    }
}

}  // namespace vcal::engine
