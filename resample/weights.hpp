// SPDX-License-Identifier: Apache-2.0
//
// Resampling weight matrices W (B x T, row-major) for the W x L engine (D-009, M2b, D-132..D-134).
// A replicate's surface is s_b[k] = sum_t W[b][t] L[t][k], so every scheme below is a way of
// writing W; the quadrature behind L runs once.
//
// Random draws (D-132). Philox4x32-10, as in dgp/, but in a key domain of its own:
//     key     = seed XOR kResampleKeyTag, split into (low 32 bits, high 32 bits),
//     counter = (scheme, replicate, j / 2, 0),
// and draw j uses words (0, 1) of that block if j is even, (2, 3) if odd, turned into a uniform
// u = (k + 1/2) 2^-52 exactly as the DGP does (dgp.md). Because the key differs from the DGP's
// for the same seed, resampling a simulated panel never reuses or shifts the stream that
// generated it. An index in [0, m) is floor(u m), clamped to m - 1 (u m can round up to m when u
// is within half an ulp of 1); its bias is below m 2^-52.
//
// Schemes:
//   iid bootstrap         draw j = 0..T-1 is the index of the period in position j.
//   moving-block bootstrap (Kunsch 1989; D-133) with block length l, 1 <= l <= T: blocks
//                         i = 0..ceil(T/l)-1 start at floor(u_i (T - l + 1)) (draw i), are laid
//                         end to end and cut at T positions. Non-circular.
//   jackknife             B = T, W[b][t] = 1 - [t == b] (delete-one).
//   walk-forward          window b covers periods [end_b - w, end_b), end_b = first_end + b step,
//                         for every end_b <= T; w = 0 means expanding windows [0, end_b).
//   external indices      a B x m matrix of period indices from any tool (D-134); row b becomes
//                         counts W[b][t] = #{j : index[b][j] = t}. Supplying the indices that
//                         another tool drew reproduces its replicates one for one.
// Indices are returned alongside weights so that any replicate can be reproduced elsewhere.
#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "dgp/philox.hpp"

namespace vcal::resample {

enum class Scheme : std::uint32_t { IidBootstrap = 1, BlockBootstrap = 2 };

// "RSMPBOOT" in ASCII: the resampling key domain (D-132).
inline constexpr std::uint64_t kResampleKeyTag = 0x52534D50424F4F54ull;

inline dgp::PhiloxKey resample_key(std::uint64_t seed) {
    const std::uint64_t k = seed ^ kResampleKeyTag;
    return {static_cast<std::uint32_t>(k), static_cast<std::uint32_t>(k >> 32)};
}

// Draw j of (seed, scheme, replicate): a uniform in (0, 1), exactly representable.
inline double draw_uniform(std::uint64_t seed, Scheme scheme, std::uint32_t replicate, std::uint32_t j) {
    const dgp::PhiloxCounter w =
        dgp::philox4x32({static_cast<std::uint32_t>(scheme), replicate, j / 2u, 0u}, resample_key(seed));
    const std::uint32_t a = (j % 2u == 0u) ? w[0] : w[2];
    const std::uint32_t b = (j % 2u == 0u) ? w[1] : w[3];
    const std::uint64_t k = (static_cast<std::uint64_t>(a >> 6) << 26) | static_cast<std::uint64_t>(b >> 6);
    return (static_cast<double>(k) + 0.5) * 0x1p-52;
}

// An index in [0, m) from draw j.
inline std::int64_t draw_index(std::uint64_t seed, Scheme scheme, std::uint32_t replicate, std::uint32_t j,
                               std::int64_t m) {
    const auto i = static_cast<std::int64_t>(draw_uniform(seed, scheme, replicate, j) * static_cast<double>(m));
    return i < m ? i : m - 1;
}

// B x T period indices of the iid (block_length == 1 with scheme IidBootstrap) or moving-block
// bootstrap, row-major.
inline std::vector<std::int64_t> bootstrap_indices(std::uint64_t seed, Scheme scheme, std::uint32_t replicates,
                                                   std::int64_t periods, std::int64_t block_length = 1) {
    if (periods < 1 || periods > (std::int64_t{1} << 31)) throw std::invalid_argument("periods out of range");
    if (scheme == Scheme::IidBootstrap) block_length = 1;
    if (block_length < 1 || block_length > periods) throw std::invalid_argument("need 1 <= block_length <= periods");
    std::vector<std::int64_t> idx(static_cast<std::size_t>(replicates) * static_cast<std::size_t>(periods));
    for (std::uint32_t b = 0; b < replicates; ++b) {
        std::int64_t* row = idx.data() + static_cast<std::size_t>(b) * static_cast<std::size_t>(periods);
        if (scheme == Scheme::IidBootstrap) {
            for (std::int64_t j = 0; j < periods; ++j) row[j] = draw_index(seed, scheme, b, static_cast<std::uint32_t>(j), periods);
            continue;
        }
        std::int64_t pos = 0;
        for (std::uint32_t i = 0; pos < periods; ++i) {
            const std::int64_t start = draw_index(seed, scheme, b, i, periods - block_length + 1);
            for (std::int64_t s = 0; s < block_length && pos < periods; ++s) row[pos++] = start + s;
        }
    }
    return idx;
}

// Counts from a B x m index matrix (each index in [0, T)).
inline std::vector<double> weights_from_indices(const std::int64_t* indices, std::int64_t replicates,
                                                std::int64_t draws, std::int64_t periods) {
    if (replicates < 1 || draws < 1 || periods < 1) throw std::invalid_argument("empty index matrix");
    std::vector<double> w(static_cast<std::size_t>(replicates * periods), 0.0);
    for (std::int64_t b = 0; b < replicates; ++b) {
        for (std::int64_t j = 0; j < draws; ++j) {
            const std::int64_t t = indices[b * draws + j];
            if (t < 0 || t >= periods) {
                throw std::invalid_argument("index " + std::to_string(t) + " outside [0, " + std::to_string(periods) + ")");
            }
            w[static_cast<std::size_t>(b * periods + t)] += 1.0;
        }
    }
    return w;
}

inline std::vector<double> jackknife_weights(std::int64_t periods) {
    if (periods < 2) throw std::invalid_argument("the jackknife needs at least 2 periods");
    std::vector<double> w(static_cast<std::size_t>(periods * periods), 1.0);
    for (std::int64_t b = 0; b < periods; ++b) w[static_cast<std::size_t>(b * periods + b)] = 0.0;
    return w;
}

// Windows ending at first_end, first_end + step, ... <= T; `window` periods long, or expanding
// from period 0 when window == 0. Returns B x T 0/1 weights (B = number of windows).
inline std::vector<double> walk_forward_weights(std::int64_t periods, std::int64_t window, std::int64_t first_end,
                                                std::int64_t step) {
    if (step < 1 || first_end < 1 || first_end > periods || window < 0 || window > first_end) {
        throw std::invalid_argument("need step >= 1, 1 <= first_end <= periods, 0 <= window <= first_end");
    }
    std::vector<double> w;
    for (std::int64_t end = first_end; end <= periods; end += step) {
        const std::int64_t begin = window == 0 ? 0 : end - window;
        for (std::int64_t t = 0; t < periods; ++t) w.push_back(t >= begin && t < end ? 1.0 : 0.0);
    }
    return w;
}

}  // namespace vcal::resample
