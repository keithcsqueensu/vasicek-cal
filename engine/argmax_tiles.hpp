// SPDX-License-Identifier: Apache-2.0
//
// The bounded argmax's tiling (D-172): tiles of kArgmaxTile x kArgmaxTile grid indices, each row's maximum over a
// tile, and the rounding allowance on a tile's bound. Shared by the bootstrap's search (resample/bootstrap.hpp) and
// by calibrate on cached rows (P-27, D-410), which keeps each row's tile maxima in the row cache. Moved here from
// resample/bootstrap.hpp unchanged.
#pragma once

#include <cmath>
#include <cstdint>

#include "core/grid.hpp"

namespace vcal::engine {

inline constexpr std::int32_t kArgmaxTile = 4;
// Rounding allowance on a tile's bound, relative to sum_j |M[b][j] U[j][T]|. A computed s_b[k] differs from its exact
// value by a few eps times sum_j |M[b][j] R[j][k]|; a point with |R[j][k]| much larger than |U[j][T]| lies far below
// the bound, so the allowance only needs to cover points near it. 1e-9 is about 1e7 eps: generous, and it costs only
// an occasional extra tile.
inline constexpr double kArgmaxBoundMargin = 1e-9;

inline std::int64_t argmax_tiles(const Grid<2>& grid) {
    const std::int32_t t0 = (grid.axis[0].n + kArgmaxTile - 1) / kArgmaxTile;
    const std::int32_t t1 = (grid.axis[1].n + kArgmaxTile - 1) / kArgmaxTile;
    return static_cast<std::int64_t>(t0) * t1;
}

// U[tile] = max over the tile of row[k], for one row of K = grid.size() values, as the bootstrap's search computes
// it (fmax from -inf, indices in ascending order).
inline void tile_maxima(const Grid<2>& grid, const double* row, double* U) {
    const std::int32_t n0 = grid.axis[0].n, n1 = grid.axis[1].n;
    const std::int32_t t1 = (n1 + kArgmaxTile - 1) / kArgmaxTile;
    const std::int64_t tiles = argmax_tiles(grid);
    for (std::int64_t i = 0; i < tiles; ++i) U[i] = -HUGE_VAL;
    for (std::int32_t i0 = 0; i0 < n0; ++i0) {
        for (std::int32_t i1 = 0; i1 < n1; ++i1) {
            const std::int32_t idx[2] = {i0, i1};
            const std::int64_t tile = static_cast<std::int64_t>(i0 / kArgmaxTile) * t1 + i1 / kArgmaxTile;
            U[tile] = std::fmax(U[tile], row[grid.flatten(idx)]);
        }
    }
}

}  // namespace vcal::engine
