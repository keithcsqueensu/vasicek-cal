// SPDX-License-Identifier: Apache-2.0
//
// Philox4x32 counter-based generator (Salmon, Moraes, Dror & Shaw, SC'11), as in Random123:
// each round multiplies counter words 0 and 2 by fixed 32-bit constants, swaps and mixes the
// 64-bit products with the other words and the key, then bumps the key by the Weyl
// constants. Integer arithmetic only, so the output is identical on every platform.
// Verified against Random123's published known answers
// (tests/golden/dgp/random123_philox4x32_kat.txt).
#pragma once

#include <array>
#include <cstdint>

namespace vcal::dgp {

using PhiloxCounter = std::array<std::uint32_t, 4>;
using PhiloxKey = std::array<std::uint32_t, 2>;

inline PhiloxCounter philox4x32(PhiloxCounter c, PhiloxKey k, int rounds = 10) {
    constexpr std::uint64_t kM0 = 0xD2511F53u;
    constexpr std::uint64_t kM1 = 0xCD9E8D57u;
    constexpr std::uint32_t kW0 = 0x9E3779B9u;
    constexpr std::uint32_t kW1 = 0xBB67AE85u;
    for (int r = 0; r < rounds; ++r) {
        if (r > 0) {
            k[0] += kW0;
            k[1] += kW1;
        }
        const std::uint64_t p0 = kM0 * c[0];
        const std::uint64_t p1 = kM1 * c[2];
        c = {static_cast<std::uint32_t>(p1 >> 32) ^ c[1] ^ k[0], static_cast<std::uint32_t>(p1),
             static_cast<std::uint32_t>(p0 >> 32) ^ c[3] ^ k[1], static_cast<std::uint32_t>(p0)};
    }
    return c;
}

}  // namespace vcal::dgp
