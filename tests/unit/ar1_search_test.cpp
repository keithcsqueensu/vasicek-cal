// SPDX-License-Identifier: Apache-2.0
// The AR(1)-aware MLE's search and profile intervals (M6; D-195, D-196; benchmark_models/ar1/search.hpp), on panels
// from the AR(1) DGP (dgp::simulate_panel_ar1) with known phi and on hand-made ones:
//   - the estimate is a maximum of ref's independent dense likelihood, and no fixed-phi maximum on a fine phi grid
//     beats it (the scan does not miss the global maximum);
//   - at phi = 0 the inner search is the binomial-mixture MLE (calibrate and its polished maximum);
//   - on long informative panels it recovers PD, rho and phi within 4 SEs, unflagged;
//   - at the box: phi held at +-kPhiBound with its score pointing out, flagged; an all-zero panel;
//   - determinism; and the profile ends against cold-started maximisations at the end points.
// Seed 0x4D36524543 ("M6REC"): fresh, disjoint from S-35's (recovery::kSeed), the validation run's 0x7A11DA7E, the
// expected information's 0x1F020261010, the bootstrap comparator's 0xB0071A61 and the scratch 0xA11CE5EED.
#include "benchmark_models/ar1/search.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
#include "engine/profile.hpp"
#include "ref/ref.hpp"
#include "tests/harness/vcal_test.hpp"

namespace {

namespace ar1 = vcal::benchmark_models::ar1;
using vcal::test::describe;

constexpr std::uint64_t kSearchSeed = 0x4D36524543ull;
const ar1::Box kBox{1e-4, 0.2, 1e-3, 0.5};  // the recovery box (D-089, D-115)

struct Panel {
    std::vector<ar1::Obs> obs;
    std::string name;
    std::int64_t T() const { return static_cast<std::int64_t>(obs.size()); }
};

Panel dgp_panel(std::uint32_t scenario, double pd, double rho, double phi, std::int64_t T, std::int64_t n) {
    std::vector<std::int64_t> ns(static_cast<std::size_t>(T), n), d(ns.size());
    const vcal::dgp::PanelSpec spec{kSearchSeed, scenario, 0u, pd, rho, ns.data(), T};
    VCAL_REQUIRE(vcal::dgp::simulate_panel_ar1(spec, phi, d.data(), nullptr) == vcal::dgp::Status::Ok);
    Panel p;
    for (std::size_t t = 0; t < d.size(); ++t) p.obs.push_back({n, d[t]});
    p.name = "scenario " + std::to_string(scenario) + " (PD " + describe(pd) + ", rho " + describe(rho) + ", phi " +
             describe(phi) + ", T " + std::to_string(T) + ", n " + std::to_string(n) + ")";
    return p;
}

ar1::Estimate3 fit(const Panel& p) {
    ar1::Estimate3 e{};
    VCAL_REQUIRE(ar1::fit(p.obs.data(), p.T(), kBox, 0.01, 0.12, e));
    return e;
}

double ref_ll(const Panel& p, double pd, double rho, double phi) {
    std::vector<vcalref::Period> q;
    for (const auto& o : p.obs) q.push_back({o.n, o.d});
    return vcalref::ar1_loglik(pd, rho, phi, q);
}

// The 2-D maximum over (PD, rho) at fixed phi, from the box's centre (a cold start).
double fixed_phi_maximum(const Panel& p, double phi, double (&u)[3]) {
    ar1::detail::Search s(p.obs.data(), p.T(), kBox);
    u[0] = 0.5 * (s.lo_[0] + s.hi_[0]);
    u[1] = 0.5 * (s.lo_[1] + s.hi_[1]);
    u[2] = phi;
    const bool pd_rho[3] = {true, true, false};
    bool conv = false;
    const double l = s.newton(u, pd_rho, ar1::kStepTolerance, conv);
    VCAL_CHECK(conv);
    return l;
}

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

}  // namespace

// The estimate is a maximum of an independent likelihood: ref's dense filter (no node placement, no quadrature in
// common) at the estimate is at least its value a step of 1e-3 away in every interior coordinate direction and
// along the diagonals. The polish locates the maximum to 1e-10 in u, so the first-order change over the step is
// negligible against the second-order drop (about 1e-6 times the curvature); the two filters agree to 1e-13.
VCAL_TEST(ar1_fit_is_a_maximum_of_the_dense_reference) {
    const std::vector<Panel> panels = {dgp_panel(0, 0.01, 0.12, 0.6, 12, 1000), dgp_panel(1, 0.02, 0.2, -0.3, 10, 500),
                                       dgp_panel(2, 0.05, 0.1, 0.0, 12, 200)};
    for (const auto& p : panels) {
        const auto e = fit(p);
        const double u_hat[3] = {ar1::detail::logit(e.value[0]), ar1::detail::logit(e.value[1]), e.value[2]};
        const double at = ref_ll(p, e.value[0], e.value[1], e.value[2]);
        double least_drop = HUGE_VAL;
        for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
                for (const double si : {-1.0, 1.0}) {
                    for (const double sj : {-1.0, 1.0}) {
                        if (i == j && sj < 0.0) continue;
                        double u[3] = {u_hat[0], u_hat[1], u_hat[2]};
                        u[i] += si * 1e-3;
                        if (j != i) u[j] += sj * 1e-3;
                        if (std::fabs(u[2]) > ar1::kPhiBound) continue;
                        const double drop =
                            at - ref_ll(p, ar1::detail::inv_logit(u[0]), ar1::detail::inv_logit(u[1]), u[2]);
                        least_drop = std::fmin(least_drop, drop);
                        if (!(drop > 0.0)) {
                            vcal::test::report_failure(__FILE__, __LINE__,
                                                       p.name + ": ref is higher a step away in (" + std::to_string(i) +
                                                           ", " + std::to_string(j) + "): drop " + describe(drop));
                        }
                    }
                }
            }
        }
        VCAL_CHECK_EQ(e.flags, 0u);
        vcal::test::note(p.name + ": phi^ " + describe(e.value[2]) + ", least drop of ref " + describe(least_drop));
    }
}

// The scan does not miss the global maximum: the 2-D maxima at phi = -0.95, -0.90, ..., 0.95, each cold-started
// from the box's centre, never exceed the estimate's log-likelihood (beyond the polish's tolerance), on panels
// where phi is weakly identified (short, few defaults: the profile in phi is flat) and where it is not. Slow.
VCAL_TEST(ar1_fit_is_not_beaten_on_a_fine_phi_grid) {
    const std::vector<Panel> panels = {dgp_panel(3, 0.01, 0.12, 0.6, 20, 100), dgp_panel(4, 0.01, 0.12, -0.6, 20, 1000),
                                       dgp_panel(5, 0.001, 0.05, 0.3, 20, 1000), dgp_panel(6, 0.05, 0.24, 0.9, 40, 1000)};
    for (const auto& p : panels) {
        const auto e = fit(p);
        double best = -HUGE_VAL, best_phi = 0.0;
        for (int k = -19; k <= 19; ++k) {
            const double phi = k == -19 ? -ar1::kPhiBound : k == 19 ? ar1::kPhiBound : 0.05 * k;
            double u[3];
            const double l = fixed_phi_maximum(p, phi, u);
            if (l > best) {
                best = l;
                best_phi = phi;
            }
        }
        VCAL_CHECK(best <= e.loglik + 1e-9);
        vcal::test::note(p.name + ": estimate " + describe(e.loglik) + " at phi " + describe(e.value[2]) +
                         "; best on the grid " + describe(best) + " at phi " + describe(best_phi) + "; local maxima " +
                         std::to_string(e.scan_local_maxima));
    }
}

// At phi = 0 the filter is the per-period likelihood bit for bit (unit_ar1_filter), so the search's 2-D maximum
// there is the binomial-mixture MLE: calibrate on the recovery grid, polished by engine::profile_intervals
// (max_at). The two searches differ (grid then polish against this quasi-Newton), so they agree to the polish's
// location tolerance, and their maxima's values to far below it.
VCAL_TEST(ar1_search_at_phi_zero_is_the_binomial_mle) {
    using Objective = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
    const auto primary = vcal::quadrature::parity_rule();
    const auto check = vcal::quadrature::parity_rule(true);
    const vcal::Grid<2> g{{{1e-4, 0.2, 61, vcal::AxisScale::Logit}, {1e-3, 0.5, 41, vcal::AxisScale::Logit}}};
    for (const auto& p : {dgp_panel(7, 0.01, 0.12, 0.0, 20, 1000), dgp_panel(8, 0.05, 0.24, 0.6, 40, 100),
                          dgp_panel(9, 0.001, 0.12, -0.3, 20, 10000)}) {
        std::vector<Objective::Obs> obs;
        for (const auto& o : p.obs) obs.push_back({o.n, o.d});
        std::vector<double> L;
        vcal::engine::Estimate2 est{};
        VCAL_REQUIRE(vcal::engine::calibrate(vcal::backends::CpuBackend{1}, Objective{}, primary, check, obs.data(), p.T(),
                                             g, L, est) == vcal::engine::Status::Ok);
        const auto prof = vcal::engine::profile_intervals(Objective{}, primary, obs.data(), p.T(), g, L, est);
        double u[3];
        const double l = fixed_phi_maximum(p, 0.0, u);
        ar1::FilterWork w;
        const double l_binomial =
            ar1::log_likelihood(p.obs.data(), p.T(), {prof.max_at[0], prof.max_at[1], 0.0}, primary, w).value;
        const double gap_pd = std::fabs(u[0] - ar1::detail::logit(prof.max_at[0]));
        const double gap_rho = std::fabs(u[1] - ar1::detail::logit(prof.max_at[1]));
        VCAL_CHECK(gap_pd <= 1e-5);
        VCAL_CHECK(gap_rho <= 1e-5);
        VCAL_CHECK(std::fabs(l - l_binomial) <= 1e-9);
        vcal::test::note(p.name + ": logit gaps " + describe(gap_pd) + ", " + describe(gap_rho) + "; log-likelihoods " +
                         describe(l - l_binomial) + " apart");
    }
}

// Recovery where phi is identified (the expected information, information_m6.csv: long panels with many defaults):
// T = 200, n = 10,000, one panel at each of phi = -0.5, 0.3 and 0.8. Every coordinate within 4 SEs of the truth
// (observed information), no flag, one local maximum along the scan. Slow: T = 200.
VCAL_TEST(ar1_fit_recovers_parameters_on_long_informative_panels) {
    std::uint32_t scenario = 10;
    for (const double phi : {-0.5, 0.3, 0.8}) {
        const auto p = dgp_panel(scenario++, 0.02, 0.1, phi, 200, 10000);
        const auto e = fit(p);
        const double truth[3] = {0.02, 0.1, phi};
        std::string line = p.name + ":";
        for (int i = 0; i < 3; ++i) {
            const double z = (e.value[i] - truth[i]) / e.se[i];
            VCAL_CHECK(std::fabs(z) <= 4.0);
            line += " " + describe(e.value[i]) + " (" + describe(z) + " SE)";
        }
        VCAL_CHECK_EQ(e.flags, 0u);
        VCAL_CHECK_EQ(e.scan_local_maxima, 1);
        vcal::test::note(line);
    }
}

// At the box: a panel from phi = +-0.999 (outside |phi| <= 0.95) with many defaults puts the maximum beyond the
// bound. The estimate holds phi exactly at the bound, the score's phi component there points out of the box, the
// fit is flagged at a bound (kFlagGridEdge) and nothing else is wrong with it. An all-zero panel: the likelihood
// rises towards the lower PD bound; the fit stays finite and is flagged at a bound.
VCAL_TEST(ar1_fit_holds_a_bound_and_flags_it) {
    std::uint32_t scenario = 20;
    for (const double phi : {0.999, -0.999}) {
        const auto p = dgp_panel(scenario++, 0.02, 0.1, phi, 100, 10000);
        const auto e = fit(p);
        VCAL_CHECK(bits_equal(e.value[2], phi > 0.0 ? ar1::kPhiBound : -ar1::kPhiBound));
        VCAL_CHECK((e.flags & vcal::engine::kFlagGridEdge) != 0u);
        VCAL_CHECK_EQ(e.flags & ~vcal::engine::kFlagGridEdge, 0u);
        ar1::detail::Search s(p.obs.data(), p.T(), kBox);
        const double u[3] = {ar1::detail::logit(e.value[0]), ar1::detail::logit(e.value[1]), e.value[2]};
        double grad[3];
        s.value_and_gradient(u, grad);
        VCAL_CHECK(grad[2] * phi > 0.0);
        vcal::test::note(p.name + ": phi^ " + describe(e.value[2]) + ", score in phi " + describe(grad[2]));
    }
    Panel zeros;
    for (int t = 0; t < 20; ++t) zeros.obs.push_back({1000, 0});
    const auto e = fit(zeros);
    VCAL_CHECK(std::isfinite(e.loglik));
    VCAL_CHECK((e.flags & vcal::engine::kFlagGridEdge) != 0u);
    VCAL_CHECK(bits_equal(e.value[0], ar1::detail::inv_logit(ar1::detail::logit(kBox.pd_lo))));
    vcal::test::note("all-zero panel: PD^ " + describe(e.value[0]) + ", rho^ " + describe(e.value[1]) + ", phi^ " +
                     describe(e.value[2]) + ", flags " + std::to_string(e.flags));
}

// The search is deterministic: a second fit of the same panel returns the same estimate, checks, diagnostics and
// work counts bit for bit (one thread; the search has no parallel reduction).
VCAL_TEST(ar1_fit_is_deterministic) {
    const auto p = dgp_panel(30, 0.01, 0.12, 0.6, 20, 1000);
    const auto a = fit(p), b = fit(p);
    for (int i = 0; i < 3; ++i) {
        VCAL_CHECK(bits_equal(a.value[i], b.value[i]));
        VCAL_CHECK(bits_equal(a.se[i], b.se[i]));
    }
    VCAL_CHECK(bits_equal(a.loglik, b.loglik));
    VCAL_CHECK(bits_equal(a.doubled_check, b.doubled_check));
    VCAL_CHECK(bits_equal(a.placement_check, b.placement_check));
    VCAL_CHECK(bits_equal(a.second_maximum_gap, b.second_maximum_gap));
    VCAL_CHECK_EQ(a.scan_local_maxima, b.scan_local_maxima);
    VCAL_CHECK_EQ(a.scan_unconverged, b.scan_unconverged);
    VCAL_CHECK_EQ(a.evaluations, b.evaluations);
    VCAL_CHECK_EQ(a.scores, b.scores);
    VCAL_CHECK_EQ(a.flags, b.flags);
}

// Each solved profile end is where the profile crosses l_max - 1.9207: the 2-D maximum over the other two
// parameters at the end, the best of four cold starts spread over the box (not the profile's warm start), is within
// 1e-8 of the threshold; a truncated end is the bound itself with the profile still above the threshold there. One
// start is not enough: at an end the inner problem can have two maxima (scenario 31's upper PD end: phi at 0.95
// with rho 0.34 on the threshold, and rho at 0.5 with phi -0.73 some 16 below it). The best of several also checks
// the other way, that the profile's warm start has not followed a lower branch.
VCAL_TEST(ar1_profile_ends_are_where_the_profile_crosses_the_threshold) {
    constexpr double kStarts[4][3] = {{0.005, 0.03, -0.6}, {0.005, 0.3, 0.0}, {0.02, 0.03, 0.6}, {0.02, 0.3, 0.9}};
    for (const auto& p : {dgp_panel(31, 0.01, 0.12, 0.6, 20, 1000), dgp_panel(32, 0.05, 0.24, 0.3, 40, 100)}) {
        const auto e = fit(p);
        ar1::ProfileIntervals3 pi{};
        const bool which[3] = {true, true, true};
        ar1::profile_intervals(p.obs.data(), p.T(), kBox, e, which, pi);
        const double target = e.loglik - vcal::engine::kProfileThreshold95;
        std::string line = p.name + ":";
        for (int a = 0; a < 3; ++a) {
            VCAL_CHECK((pi.flags[a] & vcal::engine::kIntervalNotComputed) == 0u);
            for (const int side : {-1, 1}) {
                const double end = side < 0 ? pi.lo[a] : pi.hi[a];
                const bool truncated = (pi.flags[a] & (side < 0 ? vcal::engine::kIntervalLowerTruncated
                                                                : vcal::engine::kIntervalUpperTruncated)) != 0u;
                VCAL_CHECK(side < 0 ? end < e.value[a] : end > e.value[a]);
                // The best of cold-started maxima over the other two parameters with parameter a at the end.
                double l = -HUGE_VAL;
                for (const auto& st : kStarts) {
                    ar1::detail::Search s(p.obs.data(), p.T(), kBox);
                    double u[3] = {ar1::detail::logit(st[0]), ar1::detail::logit(st[1]), st[2]};
                    u[a] = a < 2 ? ar1::detail::logit(end) : end;
                    bool others[3] = {true, true, true};
                    others[a] = false;
                    bool conv = false;
                    l = std::fmax(l, s.newton(u, others, ar1::kStepTolerance, conv));
                    VCAL_CHECK(conv);
                }
                if (truncated) {
                    VCAL_CHECK(l >= target);
                } else {
                    VCAL_CHECK(std::fabs(l - target) <= 1e-8);
                }
                line += " " + describe(end) + (truncated ? " (truncated)" : "") + " [" + describe(l - target) + "]";
            }
        }
        vcal::test::note(line);
    }
}
