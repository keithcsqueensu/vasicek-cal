// SPDX-License-Identifier: Apache-2.0
// Passing-path checks for the runner's comparison semantics.
#include "tests/harness/vcal_test.hpp"

#include <limits>

namespace {
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
}  // namespace

VCAL_TEST(near_abs_semantics) {
    VCAL_CHECK(vcal::test::near_abs(1.0, 1.0, 0.0));
    VCAL_CHECK(vcal::test::near_abs(1.0 + 1e-12, 1.0, 1e-11));
    VCAL_CHECK(!vcal::test::near_abs(1.0 + 1e-10, 1.0, 1e-11));
    VCAL_CHECK(vcal::test::near_abs(kInf, kInf, 0.0));
    VCAL_CHECK(!vcal::test::near_abs(kInf, -kInf, 1e300));
    VCAL_CHECK(!vcal::test::near_abs(kInf, 1.0, 1e300));
    VCAL_CHECK(!vcal::test::near_abs(kNaN, kNaN, kInf));
    VCAL_CHECK(!vcal::test::near_abs(0.0, kNaN, kInf));
}

VCAL_TEST(near_rel_semantics) {
    VCAL_CHECK(vcal::test::near_rel(1e-300 * (1 + 1e-13), 1e-300, 1e-12));
    VCAL_CHECK(!vcal::test::near_rel(1e-300 * (1 + 1e-11), 1e-300, 1e-12));
    VCAL_CHECK(vcal::test::near_rel(0.0, 0.0, 0.0));
    VCAL_CHECK(!vcal::test::near_rel(1e-300, 0.0, 1.0));  // expected 0 requires exact 0
    VCAL_CHECK(!vcal::test::near_rel(kNaN, 1.0, kInf));
}

VCAL_TEST(describe_round_trips_doubles) {
    VCAL_CHECK_EQ(vcal::test::describe(0.1), std::string("0.10000000000000001"));
    VCAL_CHECK_EQ(vcal::test::describe(0.1f), std::string("0.100000001"));
    VCAL_CHECK_EQ(vcal::test::describe(42), std::string("42"));
}

VCAL_TEST(check_macros_pass) {
    VCAL_CHECK(true);
    VCAL_REQUIRE(true);
    VCAL_CHECK_EQ(3, 3);
    VCAL_CHECK_NEAR(1.0, 1.0 + 1e-15, 1e-14);
    VCAL_CHECK_REL(2.0, 2.0 * (1 + 1e-15), 1e-14);
}
