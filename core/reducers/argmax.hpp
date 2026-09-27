// SPDX-License-Identifier: Apache-2.0
//
// ArgMax reducer (ARCHITECTURE.md §3.6): folds weighted-surface values s[k] into the grid
// index of the maximum. Order: larger value wins; equal values -> smaller k. That is a total
// order, so merge() is exactly associative and commutative and the result cannot depend on
// how the grid was split across threads or blocks (§6). NaN values never win; they are
// counted so the caller can flag a numerically broken surface.
#pragma once

#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "core/precision.hpp"

namespace vcal::reducers {

struct ArgMax {
    struct State {
        double best;
        std::int64_t k;  // -1 until a non-NaN value has been pushed
        std::int64_t nan_count;
    };

    VCAL_HD State init() const { return {-HUGE_VAL, -1, 0}; }

    VCAL_HD void push(State& s, std::int64_t k, double v) const {
        if (std::isnan(v)) {
            ++s.nan_count;
        } else if (beats(v, k, s)) {
            s.best = v;
            s.k = k;
        }
    }

    VCAL_HD void merge(State& s, const State& other) const {
        s.nan_count += other.nan_count;
        if (other.k >= 0 && beats(other.best, other.k, s)) {
            s.best = other.best;
            s.k = other.k;
        }
    }

private:
    VCAL_HD static bool beats(double v, std::int64_t k, const State& s) {
        return s.k < 0 || v > s.best || (v == s.best && k < s.k);
    }
};

// --- contract (C++17 trait + per-clause static_asserts, D-063) -------------------------------

namespace detail {
template <class R, class = void>
struct has_reducer_members : std::false_type {};
template <class R>
struct has_reducer_members<
    R, std::void_t<typename R::State, decltype(std::declval<const R&>().init()),
                   decltype(std::declval<const R&>().push(std::declval<typename R::State&>(), std::int64_t{}, 0.0)),
                   decltype(std::declval<const R&>().merge(std::declval<typename R::State&>(),
                                                           std::declval<const typename R::State&>()))>>
    : std::true_type {};
}  // namespace detail

template <class R>
constexpr bool check_reducer() {
    static_assert(detail::has_reducer_members<R>::value,
                  "vcal: a Reducer must declare State, init(), push(State&, int64_t k, double v) and "
                  "merge(State&, const State&)");
    return true;
}

static_assert(check_reducer<ArgMax>());

}  // namespace vcal::reducers
