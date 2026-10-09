// SPDX-License-Identifier: Apache-2.0
// resample/bartlett.hpp (S-4b as a native option, D-182, D-184): the Bartlett-corrected profile interval
// reproduces study S-4's run bit for bit; it is deterministic and independent of the thread count; its
// factor is the mean of the per-panel profile statistics; a factor of 1 gives the plain profile interval.
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "resample/bartlett.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

namespace e = vcal::engine;
namespace rs = vcal::resample;
namespace rc = vcal::recovery;
using Backend = vcal::backends::CpuBackend;

constexpr std::uint64_t kS4Seed = 0x53313050424F4F54ull;  // "S10PBOOT", the S-4 run's bootstrap seed

bool same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

// The S-4 golden was written by MSVC 19.51 (D-181). As the recovery replay (D-121): bit for bit on that
// compiler; elsewhere a platform's libm can move the fits' polished maxima inside the solvers'
// tolerances, so the factor and the ends are compared at TOL_PROFILE_CROSS_PLATFORM_REL, and every
// flag and count exactly.
#if defined(_MSC_VER) && _MSC_VER == 1951
constexpr bool kGoldenCompiler = true;
#else
constexpr bool kGoldenCompiler = false;
#endif

bool matches_golden(double got, double want) {
    if (kGoldenCompiler) return same(got, want);
    return got == want || std::fabs(got - want) <= vcal::tol::TOL_PROFILE_CROSS_PLATFORM_REL * std::fabs(want);
}

struct Fitted {
    std::vector<rc::Objective::Obs> obs;
    std::vector<double> L;
    e::Estimate2 est{};
};

Fitted fit_replicate(std::uint32_t scenario, std::uint32_t replicate) {
    static const auto primary = vcal::quadrature::parity_rule();
    static const auto check = vcal::quadrature::parity_rule(true);
    const rc::Scenario s = rc::scenario(scenario);
    const auto d = rc::panel(s, replicate);
    Fitted f;
    f.obs.resize(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) f.obs[t] = {s.obligors, d[t]};
    VCAL_REQUIRE(e::calibrate(Backend{1}, rc::Objective{}, primary, check, f.obs.data(),
                              static_cast<std::int64_t>(f.obs.size()), rc::grid(), f.L, f.est) == e::Status::Ok);
    return f;
}

rs::BartlettIntervals2 bartlett(const Fitted& f, const rs::BartlettSpec& spec, int threads = 1) {
    static const auto primary = vcal::quadrature::parity_rule();
    static const auto check = vcal::quadrature::parity_rule(true);
    rs::BartlettIntervals2 out{};
    VCAL_REQUIRE(rs::bartlett_profile_intervals(Backend{threads}, rc::Objective{}, primary, check, f.obs.data(),
                                                static_cast<std::int64_t>(f.obs.size()), rc::grid(), f.L, f.est, spec,
                                                out) == rs::BartlettStatus::Ok);
    return out;
}

bool same_result(const rs::BartlettIntervals2& a, const rs::BartlettIntervals2& b) {
    bool ok = a.panels == b.panels && a.panels_quad_flagged == b.panels_quad_flagged;
    for (int k = 0; k < 2; ++k) {
        ok = ok && same(a.factor[k], b.factor[k]) && a.factor_flags[k] == b.factor_flags[k] &&
             a.finite[k] == b.finite[k] && same(a.intervals.lo[k], b.intervals.lo[k]) &&
             same(a.intervals.hi[k], b.intervals.hi[k]) && a.intervals.flags[k] == b.intervals.flags[k];
    }
    return ok;
}

}  // namespace

// The golden is five replicates of the S-4 run (D-181): T = 20 at n = 10^4 and 100, T = 100 at n = 100,
// and two group-A replicates whose factor falls back on the grid bound. 995 bootstrap fits: slow.
VCAL_TEST(bartlett_reproduces_the_s4_run) {
    const auto g = vcal::test::read_golden_csv("bartlett/s4_run.csv");
    VCAL_REQUIRE(g.rows.size() == 5u);
    const char* p[2] = {"pd", "rho"};
    for (const auto& row : g.rows) {
        const auto sc = static_cast<std::uint32_t>(vcal::test::parse_int(row[g.column("scenario")]));
        const auto r = static_cast<std::uint32_t>(vcal::test::parse_int(row[g.column("replicate")]));
        const auto b = bartlett(fit_replicate(sc, r), {kS4Seed, sc, 1000u * r, rs::kBartlettPanels});
        VCAL_CHECK_EQ(b.panels, std::int64_t{rs::kBartlettPanels});
        for (int a = 0; a < 2; ++a) {
            const std::string s = std::string("_") + p[a];
            const double kb = vcal::test::parse_double(row[g.column("kb" + s + "_hex")]);
            const bool fallback = vcal::test::parse_int(row[g.column("kb_fallback" + s)]) != 0;
            VCAL_CHECK(matches_golden(b.factor[a], kb));
            VCAL_CHECK_EQ(b.factor_flags[a] != 0u, fallback);
            VCAL_CHECK_EQ(b.finite[a], vcal::test::parse_int(row[g.column("wstar_finite" + s)]));
            VCAL_CHECK(matches_golden(b.intervals.lo[a], vcal::test::parse_double(row[g.column("s4b_lo" + s + "_hex")])));
            VCAL_CHECK(matches_golden(b.intervals.hi[a], vcal::test::parse_double(row[g.column("s4b_hi" + s + "_hex")])));
            VCAL_CHECK_EQ(std::int64_t{b.intervals.flags[a]}, vcal::test::parse_int(row[g.column("s4b_flags" + s)]));
        }
    }
}

// Scenario 63 (T = 20, n = 100), 9 panels: the factor is the mean of the finite W*, in panel order, of the
// panels the DGP gives for (seed, stream, first_panel + b) fitted as the original; the intervals are the
// profile intervals at threshold k c / 2.
VCAL_TEST(bartlett_factor_is_the_mean_profile_statistic) {
    static const auto primary = vcal::quadrature::parity_rule();
    static const auto check = vcal::quadrature::parity_rule(true);
    const Fitted f = fit_replicate(63, 0);
    const rs::BartlettSpec spec{0x1234ull, 7u, 50u, 9u};
    const auto b = bartlett(f, spec);
    const auto T = static_cast<std::int64_t>(f.obs.size());
    const double v_hat[2] = {f.est.value[0], f.est.value[1]};
    double sum[2] = {0.0, 0.0};
    std::int64_t fin[2] = {0, 0};
    std::vector<std::int64_t> n(f.obs.size(), f.obs[0].n), d(n.size());
    for (std::uint32_t k = 0; k < spec.panels; ++k) {
        const vcal::dgp::PanelSpec ps{spec.seed, spec.stream, spec.first_panel + k, v_hat[0], v_hat[1], n.data(), T};
        VCAL_REQUIRE(vcal::dgp::simulate_panel(ps, d.data(), nullptr) == vcal::dgp::Status::Ok);
        std::vector<rc::Objective::Obs> ob(n.size());
        for (std::size_t t = 0; t < n.size(); ++t) ob[t] = {n[t], d[t]};
        std::vector<double> Lb;
        e::Estimate2 eb{};
        VCAL_REQUIRE(e::calibrate(Backend{1}, rc::Objective{}, primary, check, ob.data(), T, rc::grid(), Lb, eb) ==
                     e::Status::Ok);
        const auto lr = e::profile_lr_statistic(rc::Objective{}, primary, ob.data(), T, rc::grid(), Lb, eb, v_hat);
        for (int a = 0; a < 2; ++a) {
            if (std::isfinite(lr.w[a])) {
                sum[a] += lr.w[a];
                ++fin[a];
            }
        }
    }
    for (int a = 0; a < 2; ++a) {
        VCAL_CHECK_EQ(b.finite[a], fin[a]);
        VCAL_CHECK_EQ(b.factor_flags[a], 0u);
        VCAL_CHECK(same(b.factor[a], sum[a] / static_cast<double>(fin[a])));
        const auto p = e::profile_intervals(rc::Objective{}, primary, f.obs.data(), T, rc::grid(), f.L, f.est,
                                            e::kProfileThreshold95 * b.factor[a]);
        VCAL_CHECK(same(b.intervals.lo[a], p.lo[a]) && same(b.intervals.hi[a], p.hi[a]));
        VCAL_CHECK_EQ(b.intervals.flags[a], p.flags[a]);
    }
    VCAL_CHECK_EQ(b.panels, std::int64_t{9});
}

VCAL_TEST(bartlett_is_deterministic_and_independent_of_the_thread_count) {
    const Fitted f = fit_replicate(37, 3);
    const rs::BartlettSpec spec{0xABCDull, 1u, 0u, 7u};
    const auto a = bartlett(f, spec, 1);
    VCAL_CHECK(same_result(a, bartlett(f, spec, 1)));
    VCAL_CHECK(same_result(a, bartlett(f, spec, 4)));
    const auto other = bartlett(f, {0xABCEull, 1u, 0u, 7u});  // another seed: other panels, another factor
    VCAL_CHECK(!same(a.factor[1], other.factor[1]));
}

// Scenario 0, replicate 1 (group A): both estimates on the grid bound, so k = 1 for both and the intervals
// are the plain profile intervals bit for bit.
VCAL_TEST(bartlett_factor_on_the_bound_gives_the_profile_interval) {
    static const auto primary = vcal::quadrature::parity_rule();
    const Fitted f = fit_replicate(0, 1);
    const auto b = bartlett(f, {kS4Seed, 0u, 1000u, 5u});
    const auto p = e::profile_intervals(rc::Objective{}, primary, f.obs.data(), static_cast<std::int64_t>(f.obs.size()),
                                        rc::grid(), f.L, f.est);
    for (int a = 0; a < 2; ++a) {
        VCAL_CHECK_EQ(b.factor_flags[a], std::uint32_t{rs::kBartlettOnBound});
        VCAL_CHECK(same(b.factor[a], 1.0));
        VCAL_CHECK(same(b.intervals.lo[a], p.lo[a]) && same(b.intervals.hi[a], p.hi[a]));
        VCAL_CHECK_EQ(b.intervals.flags[a], p.flags[a]);
    }
}

VCAL_TEST(bartlett_rejects_an_invalid_spec) {
    static const auto primary = vcal::quadrature::parity_rule();
    static const auto check = vcal::quadrature::parity_rule(true);
    const Fitted f = fit_replicate(63, 0);
    const auto T = static_cast<std::int64_t>(f.obs.size());
    rs::BartlettIntervals2 out{};
    VCAL_CHECK(rs::bartlett_profile_intervals(Backend{1}, rc::Objective{}, primary, check, f.obs.data(), T, rc::grid(),
                                              f.L, f.est, {1u, 0u, 0u, 0u}, out) == rs::BartlettStatus::InvalidSpec);
    VCAL_CHECK(rs::bartlett_profile_intervals(Backend{1}, rc::Objective{}, primary, check, f.obs.data(), T, rc::grid(),
                                              f.L, f.est, {1u, 0u, 0xFFFFFFFFu, 2u}, out) ==
               rs::BartlettStatus::InvalidSpec);
    VCAL_CHECK(rs::bartlett_profile_intervals(Backend{1}, rc::Objective{}, primary, check, f.obs.data(), T, rc::grid(),
                                              f.L, f.est, {1u, 0u, 0xFFFFFFFFu, 1u}, out) == rs::BartlettStatus::Ok);
}
