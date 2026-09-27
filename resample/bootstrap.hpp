// SPDX-License-Identifier: Apache-2.0
//
// Replicate estimates from W x L, and bootstrap percentile intervals (M2b, D-134, D-135).
//
// replicate_estimates: for each row b of W, the estimate that calibrate would give on the
// reweighted panel, without ever storing the B x K surfaces (ARCHITECTURE.md §5.2):
//   1. engine::reduce_weighted with ArgMax: fused and tiled over (b, k-tile), each surface value
//      s_b[k] = sum_t W[b][t] L[t][k] computed on the fly (ascending t, compensated), tile states
//      merged in tile order. Memory O(B x tiles).
//   2. per replicate, engine::refine_2d on the 3x3 stencil around the replicate's argmax, the
//      nine stencil values recomputed from L and W[b] by the same weighted_sum.
// Results do not depend on the thread count. Standard errors are not computed per replicate
// (the bootstrap distribution is the uncertainty); each replicate carries its flags (grid edge,
// flat surface, refinement rejected, numeric).
//
// percentile_interval: the (1 - level)/2 and (1 + level)/2 sample quantiles of the replicate
// estimates, by linear interpolation between order statistics (Hyndman and Fan's type 7, the
// default of numpy.percentile and R's quantile). Percentile intervals are invariant under
// monotone reparametrisation, so taking them on the logit scale, as recommended for bounded
// parameters, gives the same interval. Replicates on the grid edge are kept at their grid value
// (they are estimates, not failures) and counted; non-finite ones are excluded and counted.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "core/grid.hpp"
#include "core/reducers/argmax.hpp"
#include "engine/refine.hpp"
#include "engine/surface.hpp"

namespace vcal::resample {

struct Replicate2 {
    double value[2];      // natural scale, grid axis order
    double surface_max;   // s_b at the grid argmax (the replicate's log-likelihood on the grid)
    std::int64_t grid_index;
    std::uint32_t flags;  // engine::kFlag* (edge, flat, rejected, numeric)
};

template <class Backend>
void replicate_estimates(const Backend& backend, const Grid<2>& grid, const double* L, std::int64_t periods,
                         const double* W, std::int64_t replicates, Replicate2* out) {
    const std::int64_t K = grid.size();
    std::vector<reducers::ArgMax::State> best(static_cast<std::size_t>(replicates));
    engine::reduce_weighted(backend, reducers::ArgMax{}, L, periods, K, W, replicates, best.data());
    backend.parallel_for(replicates, [&](std::int64_t b) {
        const auto& s = best[static_cast<std::size_t>(b)];
        Replicate2& r = out[b];
        const double nan = std::numeric_limits<double>::quiet_NaN();
        if (s.k < 0) {
            r = {{nan, nan}, nan, -1, engine::kFlagNumeric};
            return;
        }
        const double* w = W + b * periods;
        const auto value_at = [&](std::int32_t i0, std::int32_t i1) {
            const std::int32_t idx[2] = {i0, i1};
            return engine::weighted_sum(L, periods, K, w, grid.flatten(idx));
        };
        const engine::Refined2 f = engine::refine_2d(grid, s.k, value_at);
        r = {{f.value[0], f.value[1]}, s.best, s.k,
             f.flags | (s.nan_count > 0 ? static_cast<std::uint32_t>(engine::kFlagNumeric) : 0u)};
    });
}

struct PercentileInterval {
    double lo;
    double hi;
    std::int64_t used;      // finite replicate values
    std::int64_t excluded;  // non-finite replicate values
};

// Type-7 sample quantile of sorted finite values x (n >= 1), p in [0, 1].
inline double quantile_type7(const std::vector<double>& x, double p) {
    const double h = static_cast<double>(x.size() - 1) * p;
    const auto i = static_cast<std::size_t>(std::floor(h));
    if (i + 1 >= x.size()) return x.back();
    return x[i] + (h - static_cast<double>(i)) * (x[i + 1] - x[i]);
}

inline PercentileInterval percentile_interval(const Replicate2* reps, std::int64_t replicates, int param,
                                              double level = 0.95) {
    std::vector<double> x;
    x.reserve(static_cast<std::size_t>(replicates));
    for (std::int64_t b = 0; b < replicates; ++b) {
        if (std::isfinite(reps[b].value[param])) x.push_back(reps[b].value[param]);
    }
    const auto excluded = replicates - static_cast<std::int64_t>(x.size());
    if (x.empty()) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan, nan, 0, excluded};
    }
    std::sort(x.begin(), x.end());
    return {quantile_type7(x, 0.5 * (1.0 - level)), quantile_type7(x, 0.5 * (1.0 + level)),
            static_cast<std::int64_t>(x.size()), excluded};
}

}  // namespace vcal::resample
