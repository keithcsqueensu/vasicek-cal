// SPDX-License-Identifier: Apache-2.0
//
// Minimal in-house test runner (D-029).
//
// - VCAL_TEST(name) defines a self-registering test.
// - VCAL_CHECK* are non-fatal: the failure is reported and the test continues.
// - VCAL_REQUIRE aborts the current test on failure.
// - Floating-point values are printed with max_digits10 significant digits, so a
//   failure message carries enough to reproduce the comparison exactly.
// - VCAL_CHECK_NEAR / VCAL_CHECK_REL fail on NaN in either argument.
//
// The runner (runner.cpp) provides main(): no arguments runs every test in name
// order; `--list` prints test names; any other arguments are exact test names.
#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>

namespace vcal::test {

using TestFn = void (*)();

struct Registrar {
    Registrar(const char* name, TestFn fn, const char* file, int line);
};

// Thrown by VCAL_REQUIRE to abort the current test; caught by the runner.
struct RequireFailed {};

void report_failure(const char* file, int line, const std::string& message);

// |actual - expected| <= abs_tol. False if either value is NaN.
bool near_abs(double actual, double expected, double abs_tol);

// |actual - expected| <= rel_tol * |expected|. False if either value is NaN.
// With expected == 0 this requires actual == 0; use near_abs there.
bool near_rel(double actual, double expected, double rel_tol);

// Number of representable doubles between a and b. +0 and -0 are 0 apart, equal
// infinities are 0 apart, and NaN in either argument gives UINT64_MAX.
std::uint64_t ulp_distance(double a, double b);

// Informational line in the test output (observed error levels, etc.).
void note(const std::string& message);

// Value of an environment variable, or "" if unset.
std::string env_or_empty(const char* name);

// Worst ULP error over a sweep, with the input that produced it.
struct UlpStats {
    std::uint64_t max_ulps = 0;
    double worst_input = 0.0;
    double worst_actual = 0.0;
    double worst_expected = 0.0;
    std::size_t count = 0;

    void add(double input, double actual, double expected);
    std::string summary(const std::string& label) const;
};

template <class T>
std::string describe(const T& value) {
    std::ostringstream os;
    if constexpr (std::is_floating_point_v<T>) {
        os.precision(std::numeric_limits<T>::max_digits10);
    }
    os << value;
    return os.str();
}

}  // namespace vcal::test

#define VCAL_TEST(name)                                                                  \
    static void vcal_test_fn_##name();                                                   \
    static const ::vcal::test::Registrar vcal_test_reg_##name(#name, &vcal_test_fn_##name, \
                                                              __FILE__, __LINE__);       \
    static void vcal_test_fn_##name()

#define VCAL_CHECK(cond)                                                                 \
    do {                                                                                 \
        if (!(cond)) ::vcal::test::report_failure(__FILE__, __LINE__, "VCAL_CHECK(" #cond ")"); \
    } while (0)

#define VCAL_REQUIRE(cond)                                                               \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            ::vcal::test::report_failure(__FILE__, __LINE__, "VCAL_REQUIRE(" #cond ")"); \
            throw ::vcal::test::RequireFailed{};                                         \
        }                                                                                \
    } while (0)

#define VCAL_CHECK_EQ(actual, expected)                                                  \
    do {                                                                                 \
        const auto& vcal_a_ = (actual);                                                  \
        const auto& vcal_e_ = (expected);                                                \
        if (!(vcal_a_ == vcal_e_))                                                       \
            ::vcal::test::report_failure(                                                \
                __FILE__, __LINE__,                                                      \
                "VCAL_CHECK_EQ(" #actual ", " #expected "): " +                          \
                    ::vcal::test::describe(vcal_a_) + " != " + ::vcal::test::describe(vcal_e_)); \
    } while (0)

#define VCAL_CHECK_NEAR(actual, expected, abs_tol)                                       \
    do {                                                                                 \
        const double vcal_a_ = (actual);                                                 \
        const double vcal_e_ = (expected);                                               \
        const double vcal_t_ = (abs_tol);                                                \
        if (!::vcal::test::near_abs(vcal_a_, vcal_e_, vcal_t_))                          \
            ::vcal::test::report_failure(                                                \
                __FILE__, __LINE__,                                                      \
                "VCAL_CHECK_NEAR(" #actual ", " #expected "): " +                        \
                    ::vcal::test::describe(vcal_a_) + " vs " + ::vcal::test::describe(vcal_e_) + \
                    ", abs_tol " + ::vcal::test::describe(vcal_t_));                     \
    } while (0)

// Notes the sweep's worst case, and fails at the caller's line if it exceeds limit or the
// sweep saw no points.
#define VCAL_CHECK_ULP_STATS(stats, limit, label)                                     \
    do {                                                                                 \
        const ::vcal::test::UlpStats& vcal_s_ = (stats);                                 \
        const std::string vcal_l_ = (label);                                             \
        const std::uint64_t vcal_m_ = (limit);                                           \
        ::vcal::test::note(vcal_s_.summary(vcal_l_));                                    \
        if (vcal_s_.count == 0) {                                                        \
            ::vcal::test::report_failure(__FILE__, __LINE__, vcal_l_ + ": no points");   \
        } else if (vcal_s_.max_ulps > vcal_m_) {                                         \
            ::vcal::test::report_failure(__FILE__, __LINE__,                             \
                                         vcal_s_.summary(vcal_l_) + " exceeds " +        \
                                             std::to_string(vcal_m_) + " ulp");          \
        }                                                                                \
    } while (0)

#define VCAL_CHECK_REL(actual, expected, rel_tol)                                      \
    do {                                                                                 \
        const double vcal_a_ = (actual);                                                 \
        const double vcal_e_ = (expected);                                               \
        const double vcal_t_ = (rel_tol);                                                \
        if (!::vcal::test::near_rel(vcal_a_, vcal_e_, vcal_t_))                          \
            ::vcal::test::report_failure(                                                \
                __FILE__, __LINE__,                                                      \
                "VCAL_CHECK_REL(" #actual ", " #expected "): " +                         \
                    ::vcal::test::describe(vcal_a_) + " vs " + ::vcal::test::describe(vcal_e_) + \
                    ", rel_tol " + ::vcal::test::describe(vcal_t_));                     \
    } while (0)
