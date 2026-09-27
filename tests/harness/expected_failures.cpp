// SPDX-License-Identifier: Apache-2.0
// Every test here must fail with exactly one failed check. CTest asserts the
// runner's summary line and non-zero exit (see tests/CMakeLists.txt).
#include "tests/harness/vcal_test.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

VCAL_TEST(check_fails) {
    int one = 1;  // non-const on purpose: MSVC 14.44 warns (C4127) on a constant condition
    VCAL_CHECK(one + one == 3);
}

VCAL_TEST(check_eq_fails) { VCAL_CHECK_EQ(1, 2); }

VCAL_TEST(check_near_fails_on_nan) {
    VCAL_CHECK_NEAR(std::numeric_limits<double>::quiet_NaN(), 0.0,
                    std::numeric_limits<double>::infinity());
}

VCAL_TEST(check_rel_fails) { VCAL_CHECK_REL(1.0 + 1e-10, 1.0, 1e-12); }

VCAL_TEST(require_aborts_test) {
    VCAL_REQUIRE(false);
    VCAL_CHECK(false);  // must not be reached
}

VCAL_TEST(exception_is_a_failure) { throw std::runtime_error("boom"); }

VCAL_TEST(ulp_stats_exceed_limit) {
    vcal::test::UlpStats s;
    s.add(0.0, 1.0, std::nextafter(std::nextafter(1.0, 2.0), 2.0));  // 2 ulp apart
    VCAL_CHECK_ULP_STATS(s, 1, "two-ulp sweep");
}
