// SPDX-License-Identifier: Apache-2.0
// core/grid (M1.4; ARCHITECTURE.md §3.3, D-089).
#include "core/grid.hpp"

#include <cmath>
#include <cstdint>
#include <string>

#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

const vcal::AxisScale kScales[] = {vcal::AxisScale::Linear, vcal::AxisScale::Log, vcal::AxisScale::Logit,
                                   vcal::AxisScale::Probit};

std::string name(vcal::AxisScale s) {
    switch (s) {
        case vcal::AxisScale::Linear: return "linear";
        case vcal::AxisScale::Log: return "log";
        case vcal::AxisScale::Logit: return "logit";
        case vcal::AxisScale::Probit: return "probit";
    }
    return "?";
}

}  // namespace

VCAL_TEST(axis_endpoints_are_exact_and_values_increase) {
    for (const auto s : kScales) {
        const vcal::Axis a{1e-4, 0.5, 37, s};
        VCAL_CHECK(vcal::axis_error(a) == nullptr);
        VCAL_CHECK_EQ(a.value_at(0), 1e-4);
        VCAL_CHECK_EQ(a.value_at(36), 0.5);
        for (std::int32_t i = 1; i < a.n; ++i) {
            if (!(a.value_at(i - 1) < a.value_at(i))) {
                vcal::test::report_failure(__FILE__, __LINE__, name(s) + " axis not increasing at " + std::to_string(i));
            }
        }
    }
}

VCAL_TEST(scale_round_trip) {
    vcal::test::UlpStats stats;
    for (const auto s : kScales) {
        for (const double v : {1e-12, 1e-6, 1e-3, 0.02, 0.12, 0.3, 0.5, 0.75, 0.999}) {
            stats.add(v, vcal::grid::from_scaled(s, vcal::grid::to_scaled(s, v)), v);
        }
    }
    VCAL_CHECK_ULP_STATS(stats, vcal::tol::TOL_GRID_ROUND_TRIP_ULP, "from_scaled(to_scaled(v))");
}

// dv/du against a central difference in u (the delta method's Jacobian, D-095).
VCAL_TEST(dvalue_dscaled_matches_central_difference) {
    double worst = 0.0;
    for (const auto s : kScales) {
        for (const double v : {1e-6, 0.01, 0.2, 0.5, 0.8}) {
            const double u = vcal::grid::to_scaled(s, v);
            const double h = 1e-5 * (1.0 + std::fabs(u));
            const double fd =
                (vcal::grid::from_scaled(s, u + h) - vcal::grid::from_scaled(s, u - h)) / (2.0 * h);
            worst = std::fmax(worst, std::fabs(vcal::grid::dvalue_dscaled(s, u) / fd - 1.0));
            VCAL_CHECK_REL(vcal::grid::dvalue_dscaled(s, u), fd, vcal::tol::TOL_GRID_JACOBIAN_FD_REL);
        }
    }
    vcal::test::note("dv/du vs central difference: worst relative deviation " + vcal::test::describe(worst));
}

VCAL_TEST(logit_inverse_does_not_overflow) {
    VCAL_CHECK_EQ(vcal::grid::from_scaled(vcal::AxisScale::Logit, 800.0), 1.0);
    VCAL_CHECK_EQ(vcal::grid::from_scaled(vcal::AxisScale::Logit, -800.0), 0.0);
    VCAL_CHECK(std::isfinite(vcal::grid::from_scaled(vcal::AxisScale::Logit, -700.0)));
}

VCAL_TEST(flatten_is_row_major_and_round_trips) {
    const vcal::Grid<2> g{{{0.0, 1.0, 5, vcal::AxisScale::Linear}, {0.0, 1.0, 7, vcal::AxisScale::Linear}}};
    VCAL_CHECK_EQ(g.size(), std::int64_t{35});
    const std::int32_t i[2] = {2, 3};
    VCAL_CHECK_EQ(g.flatten(i), std::int64_t{2 * 7 + 3});  // axis 0 slowest
    for (std::int64_t k = 0; k < g.size(); ++k) {
        std::int32_t j[2];
        g.unflatten(k, j);
        VCAL_CHECK_EQ(g.flatten(j), k);
    }
}

VCAL_TEST(axis_validation) {
    VCAL_CHECK(vcal::axis_error({0.0, 1.0, 2, vcal::AxisScale::Linear}) != nullptr);  // < 3 points
    VCAL_CHECK(vcal::axis_error({1.0, 1.0, 5, vcal::AxisScale::Linear}) != nullptr);  // lo == hi
    VCAL_CHECK(vcal::axis_error({0.0, 1.0, 5, vcal::AxisScale::Log}) != nullptr);     // log needs lo > 0
    VCAL_CHECK(vcal::axis_error({0.1, 1.0, 5, vcal::AxisScale::Logit}) != nullptr);   // needs hi < 1
    VCAL_CHECK(vcal::axis_error({0.0, 0.5, 5, vcal::AxisScale::Probit}) != nullptr);  // needs lo > 0
    VCAL_CHECK(vcal::axis_error({1e-3, vcal::kDefaultRhoUpper, 5, vcal::AxisScale::Logit}) == nullptr);
    const vcal::Grid<2> bad{{{1e-3, 0.5, 5, vcal::AxisScale::Logit}, {0.0, 0.5, 5, vcal::AxisScale::Log}}};
    VCAL_CHECK(vcal::grid_error(bad) != nullptr);
}
