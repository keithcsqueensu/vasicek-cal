// SPDX-License-Identifier: Apache-2.0
//
// Precision policies and the host/device function qualifier (ARCHITECTURE.md §3.1).
//
//   eval_t  : arithmetic inside the integrand
//   accum_t : log-sum-exp and the Σ_t accumulation; always FP64 (D-039)
//
// PrecisionMixed is native-profile only (D-035, D-039); parity is FP64 everywhere.
//
// Device-visible, so C++17 (D-058): the contract is a trait, not a concept.
//   is_precision_v<P>        boolean, for SFINAE and tests
//   check_precision<P>()     static_asserts with a message per clause (D-063)
#pragma once

#include <limits>
#include <type_traits>

#if defined(__CUDACC__)
#define VCAL_HD __host__ __device__ inline
#else
#define VCAL_HD inline
#endif

static_assert(std::numeric_limits<double>::is_iec559, "vcal requires IEEE 754 binary64 double");
static_assert(std::numeric_limits<float>::is_iec559, "vcal requires IEEE 754 binary32 float");

namespace vcal {

struct PrecisionF64 {
    using eval_t = double;
    using accum_t = double;
};

struct PrecisionMixed {
    using eval_t = float;
    using accum_t = double;
};

namespace detail {
template <class P, class = void>
struct has_precision_members : std::false_type {};
template <class P>
struct has_precision_members<P, std::void_t<typename P::eval_t, typename P::accum_t>>
    : std::true_type {};
}  // namespace detail

template <class P, class = void>
struct is_precision : std::false_type {};
template <class P>
struct is_precision<P, std::enable_if_t<detail::has_precision_members<P>::value>>
    : std::bool_constant<std::is_floating_point_v<typename P::eval_t> &&
                         std::is_same_v<typename P::accum_t, double>> {};
template <class P>
inline constexpr bool is_precision_v = is_precision<P>::value;

template <class P>
constexpr bool check_precision() {
    static_assert(detail::has_precision_members<P>::value,
                  "vcal: a Precision policy must declare member types eval_t and accum_t");
    if constexpr (detail::has_precision_members<P>::value) {
        static_assert(std::is_floating_point_v<typename P::eval_t>,
                      "vcal: Precision::eval_t must be a floating-point type");
        static_assert(std::is_same_v<typename P::accum_t, double>,
                      "vcal: Precision::accum_t must be double (D-039: accumulation is always FP64)");
    }
    return true;
}

static_assert(check_precision<PrecisionF64>());
static_assert(check_precision<PrecisionMixed>());

}  // namespace vcal
