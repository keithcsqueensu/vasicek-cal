// SPDX-License-Identifier: Apache-2.0
// core/reducers/argmax (M1.4; ARCHITECTURE.md §3.6, §6).
#include "core/reducers/argmax.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "tests/harness/vcal_test.hpp"

namespace {

using vcal::reducers::ArgMax;

bool same_state(const ArgMax::State& a, const ArgMax::State& b) {
    return std::memcmp(&a.best, &b.best, sizeof a.best) == 0 && a.k == b.k && a.nan_count == b.nan_count;
}

// Deterministic values with repeated maxima, NaNs and -inf (a 64-bit LCG; no <random> needed).
std::vector<double> sample_values(std::size_t n) {
    std::vector<double> v(n);
    std::uint64_t x = 0x9E3779B97F4A7C15ull;
    for (auto& e : v) {
        x = x * 6364136223846793005ull + 1442695040888963407ull;
        const auto r = static_cast<int>(x >> 58);  // 0..63
        e = r == 0 ? std::numeric_limits<double>::quiet_NaN()
            : r == 1 ? -std::numeric_limits<double>::infinity()
                     : static_cast<double>(r % 7);  // max value 6 repeats many times
    }
    return v;
}

}  // namespace

VCAL_TEST(argmax_prefers_larger_value_then_smaller_index) {
    const ArgMax r;
    auto s = r.init();
    r.push(s, 5, 1.0);
    r.push(s, 2, 1.0);  // tie: smaller k wins
    r.push(s, 9, 0.5);
    VCAL_CHECK_EQ(s.k, std::int64_t{2});
    r.push(s, 7, 2.0);
    VCAL_CHECK_EQ(s.k, std::int64_t{7});
    VCAL_CHECK_EQ(s.best, 2.0);
}

VCAL_TEST(argmax_counts_nan_and_never_selects_it) {
    const ArgMax r;
    auto s = r.init();
    r.push(s, 0, std::numeric_limits<double>::quiet_NaN());
    VCAL_CHECK_EQ(s.k, std::int64_t{-1});
    r.push(s, 1, -std::numeric_limits<double>::infinity());
    VCAL_CHECK_EQ(s.k, std::int64_t{1});  // -inf is a value; it beats "nothing yet"
    r.push(s, 2, std::numeric_limits<double>::quiet_NaN());
    VCAL_CHECK_EQ(s.nan_count, std::int64_t{2});
    VCAL_CHECK_EQ(s.k, std::int64_t{1});
}

// Any split into chunks, merged in any order, gives the identical state (§6).
VCAL_TEST(argmax_merge_is_order_and_split_independent) {
    const ArgMax r;
    const auto v = sample_values(5000);
    auto serial = r.init();
    for (std::size_t k = 0; k < v.size(); ++k) r.push(serial, static_cast<std::int64_t>(k), v[k]);
    VCAL_CHECK(serial.k >= 0);

    for (const std::size_t chunk : {1u, 7u, 64u, 1000u, 4999u}) {
        std::vector<ArgMax::State> parts;
        for (std::size_t lo = 0; lo < v.size(); lo += chunk) {
            auto s = r.init();
            for (std::size_t k = lo; k < std::min(v.size(), lo + chunk); ++k) r.push(s, static_cast<std::int64_t>(k), v[k]);
            parts.push_back(s);
        }
        auto forward = r.init();
        for (const auto& p : parts) r.merge(forward, p);
        auto backward = r.init();
        for (auto it = parts.rbegin(); it != parts.rend(); ++it) r.merge(backward, *it);
        VCAL_CHECK(same_state(forward, serial));
        VCAL_CHECK(same_state(backward, serial));
    }
}
