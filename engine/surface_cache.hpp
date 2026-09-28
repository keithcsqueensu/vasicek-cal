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
    // The row for observation y under this configuration (config = surface_configuration(...)),
    // computed on first request with the backend's parallel_for over grid points.
    template <class Backend, class Objective, class Integrator, int D>
    const std::vector<double>& row(const Backend& backend, const std::string& config, const Objective& objective,
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
            const std::int64_t K = grid.size();
            e->row.assign(static_cast<std::size_t>(K), 0.0);
            backend.parallel_for(K, [&](std::int64_t k) {
                double v[D];
                grid.values(k, v);
                e->row[static_cast<std::size_t>(k)] = objective.log_contrib(y, Objective::theta(v), integrator);
            });
            computed = true;
        });
        (computed ? misses_ : hits_).fetch_add(1, std::memory_order_relaxed);
        return e->row;
    }

    std::int64_t hits() const { return hits_.load(); }
    std::int64_t misses() const { return misses_.load(); }
    std::size_t size() const {
        const std::lock_guard<std::mutex> lock(mutex_);
        return rows_.size();
    }

private:
    struct Entry {
        std::once_flag once;
        std::vector<double> row;
    };
    mutable std::mutex mutex_;
    std::unordered_map<std::string, std::unique_ptr<Entry>> rows_;
    std::atomic<std::int64_t> hits_{0}, misses_{0};
};

// evaluate_surface through the cache: L[t][k] for every period, each distinct observation's row
// taken from the cache (computed on first request). Identical to evaluate_surface bit for bit.
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

}  // namespace vcal::engine
