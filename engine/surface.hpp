// SPDX-License-Identifier: Apache-2.0
//
// Surface engine (ARCHITECTURE.md §3.7, §5; D-008, D-009).
//
//   Stage 1  evaluate_surface: L[t][k] = log_contrib(obs_t, theta(grid point k)), T x K row-major.
//            Computed once per (data, grid, model). l_t depends only on obs_t and theta, so
//            periods with equal observations are evaluated once and the row copied (D-122):
//            the copy is bitwise what evaluating it again would give. Low-default panels
//            repeat heavily (a T = 100 panel with n = 100 and PD 0.1% has 3 distinct d).
//   Stage 2  reduce_weighted:  for each weight row b, fold s_b[k] = sum_t W[b][t] L[t][k] into a
//            Reducer state. Fused and tiled over k, so the B x K matrix is never stored.
//
// Determinism (§6): parallelism is only across grid points (stage 1) and across fixed-size
// k-tiles (stage 2). Each sum over t runs serially, in ascending t, inside one call of
// weighted_sum. Tiles have a fixed size independent of the thread count, and their states
// are merged in tile order. Results are therefore identical for any number of threads.
//
// A Backend provides `void parallel_for(std::int64_t count, const F& f) const`, calling
// f(i) exactly once for each i in [0, count). The engine includes no backend.
#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "core/grid.hpp"
#include "core/precision.hpp"

namespace vcal::engine {

inline constexpr std::int64_t kReduceTileSize = 1024;

// sum_t w[t] * L[t * K + k] in ascending t, Neumaier-compensated; periods with w[t] == 0 are
// skipped, so a zero weight never meets a non-finite L. The single definition every reducer
// and refinement uses, so they all see the same values bit for bit.
VCAL_HD double weighted_sum(const double* L, std::int64_t T, std::int64_t K, const double* w, std::int64_t k) {
    double sum = 0.0;
    double comp = 0.0;
    double nonfinite = 0.0;  // +-inf / NaN terms bypass compensation, which they would poison
    bool any_nonfinite = false;
    for (std::int64_t t = 0; t < T; ++t) {
        if (w[t] == 0.0) continue;
        const double x = w[t] * L[t * K + k];
        if (!std::isfinite(x)) {
            nonfinite += x;
            any_nonfinite = true;
            continue;
        }
        const double s = sum + x;
        comp += std::fabs(sum) >= std::fabs(x) ? (sum - s) + x : (x - s) + sum;
        sum = s;
    }
    return any_nonfinite ? nonfinite + (sum + comp) : sum + comp;
}

template <class Backend, class Objective, class Integrator, int D>
void evaluate_surface(const Backend& backend, const Objective& objective, const Integrator& integrator,
                      const typename Objective::Obs* obs, std::int64_t periods, const Grid<D>& grid, double* L) {
    static_assert(D == Objective::n_params, "grid dimension must match Objective::n_params");
    const std::int64_t K = grid.size();
    // first[t]: the earliest period with the same observation as t. O(T^2) comparisons, which is
    // negligible next to T x K integrals for any realistic number of periods.
    std::vector<std::int64_t> first(static_cast<std::size_t>(periods));
    for (std::int64_t t = 0; t < periods; ++t) {
        std::int64_t u = 0;
        while (!(obs[u] == obs[t])) ++u;
        first[static_cast<std::size_t>(t)] = u;
    }
    backend.parallel_for(K, [&](std::int64_t k) {
        double v[D];
        grid.values(k, v);
        const auto theta = Objective::theta(v);
        for (std::int64_t t = 0; t < periods; ++t) {
            const std::int64_t u = first[static_cast<std::size_t>(t)];
            L[t * K + k] = u == t ? objective.log_contrib(obs[t], theta, integrator) : L[u * K + k];
        }
    });
}

// out[b] = reduction of s_b over all k, for b in [0, B). W is B x T row-major.
template <class Backend, class Reducer>
void reduce_weighted(const Backend& backend, const Reducer& reducer, const double* L, std::int64_t T,
                     std::int64_t K, const double* W, std::int64_t B, typename Reducer::State* out) {
    const std::int64_t tiles = (K + kReduceTileSize - 1) / kReduceTileSize;
    std::vector<typename Reducer::State> partial(static_cast<std::size_t>(B * tiles));
    backend.parallel_for(B * tiles, [&](std::int64_t job) {
        const std::int64_t b = job / tiles;
        const std::int64_t tile = job % tiles;
        auto state = reducer.init();
        const std::int64_t end = tile * kReduceTileSize + kReduceTileSize < K ? tile * kReduceTileSize + kReduceTileSize : K;
        for (std::int64_t k = tile * kReduceTileSize; k < end; ++k) {
            reducer.push(state, k, weighted_sum(L, T, K, W + b * T, k));
        }
        partial[static_cast<std::size_t>(job)] = state;
    });
    for (std::int64_t b = 0; b < B; ++b) {
        auto state = reducer.init();
        for (std::int64_t tile = 0; tile < tiles; ++tile) reducer.merge(state, partial[static_cast<std::size_t>(b * tiles + tile)]);
        out[b] = state;
    }
}

}  // namespace vcal::engine
