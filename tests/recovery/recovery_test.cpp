// SPDX-License-Identifier: Apache-2.0
//
// M1.8 recovery harness tests (D-121). The full run (R = 1000 replicates of 81 scenarios) is too
// slow for CI, so it is a tool (recovery_harness --write) whose output is committed. These tests:
//   - check the summary statistics on a constructed sample (flagged replicates stay in bias and
//     RMSE but out of coverage; the verdict rule; the interval);
//   - enforce the exit criteria on the committed summary (D-124). Engine correctness must hold:
//     RMSE falling with T, no quadrature or numeric flags. Coverage verdicts must be what the
//     recorded counts give under the current tolerances, with no UNREVIEWED verdict: the set of
//     out-of-band verdicts must equal the reviewed list kKnownFindings exactly;
//   - replay replicates 0..kReplayReplicates-1 of every scenario (slow; Release only) against the
//     committed replay values, at the replay tolerances (a D-107 comparison of its own).
#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/recovery/recovery.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace rc = vcal::recovery;
namespace e = vcal::engine;
namespace tol = vcal::tol;
using vcal::test::describe;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// A fit whose estimate sits `k` interval half-widths from the truth in the logit coordinate, on
// both axes, with scaled-coordinate SE `se_u`.
rc::Fit fit_at(const rc::Scenario& s, double k, double se_u, std::uint32_t flags) {
    const auto g = rc::grid();
    const double truth[2] = {s.pd, s.rho};
    rc::Fit f{};
    for (int a = 0; a < 2; ++a) {
        const auto sc = g.axis[a].scale;
        const double u = vcal::grid::to_scaled(sc, truth[a]) + k * rc::kZ975 * se_u;
        f.value[a] = vcal::grid::from_scaled(sc, u);
        f.se[a] = (flags & e::kFlagGridEdge) ? kNaN : se_u * vcal::grid::dvalue_dscaled(sc, u);
        // A profile interval equal to the Wald one, except on the grid edge, where it runs from the
        // lower bound of the box (truncated) to just past the truth.
        const double half = rc::kZ975 * se_u;
        if (flags & e::kFlagGridEdge) {
            f.prof_lo[a] = g.axis[a].lo;
            f.prof_hi[a] = vcal::grid::from_scaled(sc, vcal::grid::to_scaled(sc, truth[a]) + 0.1);
            f.prof_flags[a] = e::kIntervalLowerTruncated;
        } else {
            f.prof_lo[a] = vcal::grid::from_scaled(sc, u - half);
            f.prof_hi[a] = vcal::grid::from_scaled(sc, u + half);
            f.prof_flags[a] = 0u;
        }
        f.boot_lo[a] = f.prof_lo[a];  // the bootstrap interval, here equal to the profile one
        f.boot_hi[a] = f.prof_hi[a];
    }
    f.loglik = -1.0;
    f.flags = flags;
    return f;
}

}  // namespace

// --- the statistics on a constructed sample --------------------------------------------------------

VCAL_TEST(recovery_summary_statistics_on_a_constructed_sample) {
    const rc::Scenario s = rc::scenario(40);  // PD 1%, rho 0.12, T 40, n 1000
    VCAL_CHECK_EQ(s.pd, 0.01);
    VCAL_CHECK_EQ(s.rho, 0.12);
    VCAL_CHECK_EQ(s.periods, 40);
    VCAL_CHECK_EQ(s.obligors, 1000);

    // Inside and outside the interval, just either side of its edge.
    VCAL_CHECK(rc::covers(rc::grid(), 0, fit_at(s, 1.0 - 1e-9, 0.2, 0), s.pd));
    VCAL_CHECK(!rc::covers(rc::grid(), 0, fit_at(s, 1.0 + 1e-9, 0.2, 0), s.pd));
    VCAL_CHECK(rc::covers(rc::grid(), 1, fit_at(s, -(1.0 - 1e-9), 0.3, 0), s.rho));
    VCAL_CHECK(!rc::covers(rc::grid(), 1, fit_at(s, -(1.0 + 1e-9), 0.3, 0), s.rho));

    // Four replicates: covered, not covered, on the grid edge (no SE), near a bound (would cover).
    const rc::Fit fits[4] = {fit_at(s, 0.0, 0.2, 0), fit_at(s, 2.0, 0.2, 0), fit_at(s, 5.0, 0.2, e::kFlagGridEdge),
                             fit_at(s, 0.0, 0.2, e::kFlagNearBound)};
    const rc::Summary sum = rc::summarise(s, fits, 4);
    VCAL_CHECK_EQ(sum.edge, 1);
    VCAL_CHECK_EQ(sum.near_bound, 1);
    VCAL_CHECK_EQ(sum.flagged, 2);
    VCAL_CHECK_EQ(sum.unflagged, 2);
    VCAL_CHECK_EQ(sum.flagged_fraction, 0.5);
    const double truth[2] = {s.pd, s.rho};
    for (int a = 0; a < 2; ++a) {
        const auto& p = sum.param[a];
        VCAL_CHECK_EQ(p.covered, 1);  // the near-bound replicate would cover, but is not counted
        VCAL_CHECK_EQ(p.coverage, 0.5);
        VCAL_CHECK(p.verdict == rc::Verdict::Deferred);
        // Bias and RMSE use all four, the edge replicate included.
        double mean = 0.0, sq = 0.0;
        for (const auto& f : fits) {
            mean += f.value[a] / 4.0;
            sq += (f.value[a] - truth[a]) * (f.value[a] - truth[a]) / 4.0;
        }
        VCAL_CHECK_REL(p.bias, mean - truth[a], 1e-12);
        VCAL_CHECK_REL(p.rmse, std::sqrt(sq), 1e-12);
        // The SE diagnostic uses the two unflagged replicates only: u offsets 0 and 2 z 0.2.
        VCAL_CHECK_REL(p.sd_u, 2.0 * rc::kZ975 * 0.2 / std::sqrt(2.0), 1e-9);
        VCAL_CHECK_REL(p.rms_se_u, 0.2, 1e-12);
        // t(T-1), T = 40: t = 2.02 < 2 z, so the replicate 2 half-widths out is still not covered.
        VCAL_CHECK_EQ(p.t_covered, 1);
        // Profile: over ALL replicates. Covered: the one at the truth, the edge one (truncated at the
        // bound, reaching past the truth) and the near-bound one; not the one 2 half-widths out.
        VCAL_CHECK_EQ(p.profile_covered, 3);
        VCAL_CHECK_EQ(p.profile_truncated, 1);
        VCAL_CHECK_EQ(p.profile_not_computed, 0);
        VCAL_CHECK_EQ(p.profile_coverage, 0.75);
        VCAL_CHECK(p.profile_verdict == rc::Verdict::Pass);  // 4 replicates: the band is 0.95 +- 0.36
        VCAL_CHECK_EQ(p.boot_covered, 3);
        VCAL_CHECK_EQ(p.boot_degenerate, 0);
    }

    // The verdict rule: DEFERRED from the flagged-fraction threshold up, else the Monte Carlo band.
    double lo = 0.0, hi = 0.0;
    rc::coverage_band(1000, lo, hi);
    VCAL_CHECK_REL(hi - 0.95, tol::TOL_RECOVERY_COVERAGE_BAND_Z * std::sqrt(0.95 * 0.05 / 1000.0), 1e-15);
    const double f_max = tol::TOL_RECOVERY_MAX_FLAGGED_FRACTION;
    for (const bool reviewed : {false, true}) {
        const rc::Verdict outside = reviewed ? rc::Verdict::KnownFinding : rc::Verdict::Unreviewed;
        VCAL_CHECK(rc::verdict(0.0, 1000, 0.95, lo, hi, reviewed) == rc::Verdict::Pass);
        VCAL_CHECK(rc::verdict(0.0, 1000, lo - 1e-9, lo, hi, reviewed) == outside);
        VCAL_CHECK(rc::verdict(0.0, 1000, hi + 1e-9, lo, hi, reviewed) == outside);
        VCAL_CHECK(rc::verdict(f_max, 1000, 0.95, lo, hi, reviewed) == rc::Verdict::Deferred);
        VCAL_CHECK(rc::verdict(f_max, 1000, 0.5, lo, hi, reviewed) == rc::Verdict::Deferred);
        VCAL_CHECK(rc::verdict(f_max * 0.99, 1000, 0.95, lo, hi, reviewed) == rc::Verdict::Pass);
    }
    // The profile verdict is never DEFERRED, and each reviewed kind holds only on its own side.
    const rc::KnownProfileFinding conservative{0, 0, rc::ProfileDiagnosis::TruncationConservative};
    const rc::KnownProfileFinding under{0, 0, rc::ProfileDiagnosis::SmallTUndercoverage};
    VCAL_CHECK(rc::profile_verdict(0.95, lo, hi, nullptr) == rc::Verdict::Pass);
    VCAL_CHECK(rc::profile_verdict(0.95, lo, hi, &conservative) == rc::Verdict::Pass);
    VCAL_CHECK(rc::profile_verdict(hi + 1e-9, lo, hi, nullptr) == rc::Verdict::Unreviewed);
    VCAL_CHECK(rc::profile_verdict(hi + 1e-9, lo, hi, &conservative) == rc::Verdict::Conservative);
    VCAL_CHECK(rc::profile_verdict(hi + 1e-9, lo, hi, &under) == rc::Verdict::Unreviewed);
    VCAL_CHECK(rc::profile_verdict(lo - 1e-9, lo, hi, &under) == rc::Verdict::KnownFinding);
    VCAL_CHECK(rc::profile_verdict(lo - 1e-9, lo, hi, &conservative) == rc::Verdict::Unreviewed);
    // An interval that was not computed does not cover.
    rc::Fit missing = fit_at(s, 0.0, 0.2, 0);
    missing.prof_flags[0] = e::kIntervalNotComputed;
    VCAL_CHECK(!rc::profile_covers(0, missing, s.pd));
    VCAL_CHECK(rc::profile_covers(0, fit_at(s, 0.0, 0.2, 0), s.pd));

    // RMSE must fall strictly with T for each (PD, rho, n) and parameter.
    std::vector<rc::Summary> all(rc::kScenarios);
    for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
        all[id].s = rc::scenario(id);
        for (auto& p : all[id].param) p.rmse = 1.0 / static_cast<double>(all[id].s.periods);
    }
    VCAL_CHECK(rc::rmse_violations(all).empty());
    all[4].param[1].rmse = all[1].param[1].rmse;  // scenario 4 is scenario 1 with T = 40 instead of 20
    const auto bad = rc::rmse_violations(all);
    VCAL_REQUIRE(bad.size() == 1u);
    vcal::test::note("constructed violation detected as expected: " + bad[0]);
}

// S-23: q's summary on a constructed sample. Five replicates around the true q, with profile intervals
// covering, entirely below, entirely above, not computed and box-limited at the lower end.
VCAL_TEST(recovery_q_summary_on_a_constructed_sample) {
    const rc::Scenario s = rc::scenario(40);
    const double q = rc::q_truth(s);
    VCAL_CHECK_REL(q, e::conditional_pd(0.01, 0.12), 0.0);
    rc::Fit fits[5] = {fit_at(s, 0.0, 0.2, 0), fit_at(s, 0.0, 0.2, 0), fit_at(s, 0.0, 0.2, 0),
                       fit_at(s, 0.0, 0.2, 0), fit_at(s, 0.0, 0.2, e::kFlagNearBound)};
    const double lo[5] = {0.9 * q, 0.5 * q, 1.1 * q, kNaN, 0.8 * q};
    const double hi[5] = {1.2 * q, 0.9 * q, 1.5 * q, kNaN, 1.3 * q};
    const double hat[5] = {q, 0.7 * q, 1.3 * q, q, 1.4 * q};  // the last lies outside its own interval
    const std::uint32_t flags[5] = {0u, 0u, 0u, e::kIntervalNotComputed, e::kIntervalLowerBoxLimited};
    for (int r = 0; r < 5; ++r) {
        fits[r].q_hat = hat[r];
        fits[r].q_se_s = 0.1;
        fits[r].q_prof_lo = lo[r];
        fits[r].q_prof_hi = hi[r];
        fits[r].q_prof_flags = flags[r];
        fits[r].q_boot_lo = lo[r];
        fits[r].q_boot_hi = hi[r];
    }
    const rc::Summary sum = rc::summarise(s, fits, 5);
    const auto& qs = sum.q;
    VCAL_CHECK_EQ(qs.profile_covered, 2);
    VCAL_CHECK_EQ(qs.profile_miss_below, 1);
    VCAL_CHECK_EQ(qs.profile_miss_above, 1);
    VCAL_CHECK_EQ(qs.profile_not_computed, 1);
    VCAL_CHECK_EQ(qs.profile_box_limited, 1);
    VCAL_CHECK_EQ(qs.profile_box_limited_lo, 1);
    VCAL_CHECK_EQ(qs.profile_truncated, 0);
    VCAL_CHECK_EQ(qs.profile_excludes_estimate, 1);
    VCAL_CHECK_EQ(qs.profile_coverage, 0.4);
    VCAL_CHECK_EQ(qs.boot_covered, 2);  // NaN ends do not cover
    VCAL_CHECK_EQ(qs.boot_miss_below, 1);
    VCAL_CHECK_EQ(qs.boot_miss_above, 1);
    // Wald on the four unflagged: logit(q_hat) +- 1.96 x 0.1 contains logit(q) only for the first
    // and the fourth (q_hat = q); 0.7 q lies below, 1.3 q above.
    VCAL_CHECK_EQ(qs.wald_covered, 2);
    VCAL_CHECK_EQ(qs.wald_miss_below, 1);
    VCAL_CHECK_EQ(qs.wald_miss_above, 1);
    VCAL_CHECK_EQ(qs.wald_coverage, 0.5);
    // The median of q_hat / q - 1 over all five: {0, -0.3, 0.3, 0, 0.4} -> 0.
    VCAL_CHECK(std::fabs(qs.median_rel_err) <= 1e-15);
    // The median log width over the four computed intervals: log(1.2/0.9), log(0.9/0.5), log(1.5/1.1)
    // and log(1.3/0.8), whose middle two are the third and the fourth.
    VCAL_CHECK_REL(qs.profile_log_width_median, 0.5 * (std::log(1.5 / 1.1) + std::log(1.3 / 0.8)), 1e-12);
}

// S-23 acceptance: q's interval contains q_hat. Two replicates whose estimate is the box's corner
// (PD 1e-4, rho 1e-3), where q_hat is the lower limit of q and the lower end is truncated there;
// in the first full-matrix run the end was 1 ulp above q_hat (a round trip through logit).
VCAL_TEST(recovery_q_interval_contains_the_estimate_at_the_corner) {
    const std::uint32_t cases[2][2] = {{6, 534}, {22, 793}};
    for (const auto& c : cases) {
        const rc::Fit f = rc::fit(rc::scenario(c[0]), c[1]);
        const auto g = rc::grid();
        vcal::test::note("scenario " + std::to_string(c[0]) + " replicate " + std::to_string(c[1]) + ": estimate (" +
                         describe(f.value[0]) + ", " + describe(f.value[1]) + "), q_hat " + describe(f.q_hat) +
                         ", interval [" + describe(f.q_prof_lo) + ", " + describe(f.q_prof_hi) + "], flags " +
                         std::to_string(f.q_prof_flags));
        VCAL_CHECK_EQ(f.value[0], g.axis[0].lo);
        VCAL_CHECK_EQ(f.value[1], g.axis[1].lo);
        VCAL_CHECK(f.q_prof_flags & e::kIntervalLowerTruncated);
        VCAL_CHECK_EQ(f.q_prof_lo, f.q_hat);
        VCAL_CHECK(f.q_prof_lo <= f.q_hat && f.q_hat <= f.q_prof_hi);
    }
}

// --- the committed goldens meet the exit criteria ------------------------------------------------

VCAL_TEST(recovery_goldens_meet_the_exit_criteria) {
    // Every output of recovery_harness --write exists and is non-empty.
    for (const std::string& f : {vcal::test::golden_path("recovery/MANIFEST.json"),
                                 vcal::test::golden_path("recovery/replay.csv"),
                                 vcal::test::golden_path("recovery/replay_panels.csv"),
                                 vcal::test::golden_path("recovery/summary.csv"),
                                 std::string(VCAL_SOURCE_DIR) + "/docs/methodology/recovery_results.md"}) {
        std::error_code ec;
        const auto size = std::filesystem::file_size(f, ec);
        VCAL_CHECK(!ec && size > 0);
    }
    const auto t = vcal::test::read_golden_csv("recovery/summary.csv");
    VCAL_REQUIRE(t.rows.size() == rc::kScenarios);
    std::vector<rc::Summary> all(rc::kScenarios);
    int count[rc::kVerdicts] = {};
    int pcount[rc::kVerdicts] = {};
    int bcount[rc::kVerdicts] = {};
    std::string unreviewed, stale;
    const char* prefix[2] = {"pd", "rho"};
    for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
        const auto& r = t.rows[id];
        const auto col = [&](const std::string& c) { return r[t.column(c)]; };
        const auto num = [&](const std::string& c) { return vcal::test::parse_int(col(c)); };
        const rc::Scenario s = rc::scenario(id);
        VCAL_REQUIRE(num("scenario") == id);
        VCAL_CHECK_EQ(num("periods"), s.periods);
        VCAL_CHECK_EQ(num("obligors"), s.obligors);
        VCAL_CHECK_EQ(num("replicates"), static_cast<std::int64_t>(rc::kReplicates));
        VCAL_CHECK_EQ(num("quad_unconverged"), 0);  // the parity rule converges on every replicate
        VCAL_CHECK_EQ(num("numeric"), 0);
        const std::int64_t flagged = num("flagged"), unflagged = num("unflagged");
        VCAL_CHECK_EQ(flagged + unflagged, static_cast<std::int64_t>(rc::kReplicates));
        const double fraction = static_cast<double>(flagged) / static_cast<double>(rc::kReplicates);
        double lo = 0.0, hi = 0.0;
        rc::coverage_band(unflagged, lo, hi);
        rc::Summary& sum = all[id];
        sum.s = s;
        for (int a = 0; a < 2; ++a) {
            const std::string p = prefix[a];
            sum.param[a].rmse = vcal::test::parse_double(col(p + "_rmse_hex"));
            const double coverage = unflagged > 0 ? static_cast<double>(num(p + "_covered")) / static_cast<double>(unflagged)
                                                  : kNaN;
            // The recorded verdict must be what the current tolerances and the reviewed list give for
            // the recorded counts.
            const bool reviewed = rc::known_finding(id, a) != nullptr;
            const rc::Verdict v = rc::verdict(fraction, unflagged, coverage, lo, hi, reviewed);
            VCAL_CHECK_EQ(col(p + "_verdict"), std::string(rc::verdict_text(v)));
            ++count[static_cast<int>(v)];
            const std::string at = " " + std::to_string(id) + "/" + p + " (coverage " + describe(coverage) +
                                   ", band [" + describe(lo) + ", " + describe(hi) + "])";
            if (v == rc::Verdict::Unreviewed) unreviewed += at;
            if (reviewed && v != rc::Verdict::KnownFinding) stale += at;  // listed, but no longer outside the band
            // Profile likelihood: coverage over all replicates, never deferred (D-131).
            double plo = 0.0, phi = 0.0;
            rc::coverage_band(static_cast<std::int64_t>(rc::kReplicates), plo, phi);
            const double pcov = static_cast<double>(num(p + "_profile_covered")) / static_cast<double>(rc::kReplicates);
            const auto* previewed = rc::known_profile_finding(id, a);
            const rc::Verdict pv = rc::profile_verdict(pcov, plo, phi, previewed);
            VCAL_CHECK_EQ(col(p + "_profile_verdict"), std::string(rc::verdict_text(pv)));
            ++pcount[static_cast<int>(pv)];
            const std::string pat = " " + std::to_string(id) + "/" + p + " (profile coverage " + describe(pcov) + ")";
            if (pv == rc::Verdict::Unreviewed) unreviewed += pat;
            if (previewed && pv == rc::Verdict::Pass) stale += pat;  // listed, but now inside the band
            // iid bootstrap percentile interval: same rules as the profile interval (D-137).
            const double bcov = static_cast<double>(num(p + "_boot_covered")) / static_cast<double>(rc::kReplicates);
            const auto* breviewed = rc::known_bootstrap_finding(id, a);
            const rc::Verdict bv = rc::bootstrap_verdict(bcov, plo, phi, breviewed);
            VCAL_CHECK_EQ(col(p + "_boot_verdict"), std::string(rc::verdict_text(bv)));
            ++bcount[static_cast<int>(bv)];
            const std::string bat = " " + std::to_string(id) + "/" + p + " (bootstrap coverage " + describe(bcov) + ")";
            if (bv == rc::Verdict::Unreviewed) unreviewed += bat;
            if (breviewed && bv == rc::Verdict::Pass) stale += bat;
            // Engine correctness: every fit's endpoint residual within its tolerance (D-131).
            VCAL_CHECK(vcal::test::parse_double(col("profile_residual_max_hex")) <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
        }
    }
    vcal::test::note("Wald: " + std::to_string(count[0]) + " PASS, " + std::to_string(count[1]) + " DEFERRED, " +
                     std::to_string(count[2]) + " KNOWN FINDING, " + std::to_string(count[3]) +
                     " UNREVIEWED (scenario x parameter)");
    vcal::test::note("profile: " + std::to_string(pcount[0]) + " PASS, " +
                     std::to_string(pcount[static_cast<int>(rc::Verdict::Conservative)]) + " CONSERVATIVE, " +
                     std::to_string(pcount[2]) + " KNOWN FINDING, " + std::to_string(pcount[3]) + " UNREVIEWED");
    VCAL_CHECK_EQ(pcount[static_cast<int>(rc::Verdict::Conservative)],
                  static_cast<int>(rc::known_profile_finding_count(rc::ProfileDiagnosis::TruncationConservative)));
    VCAL_CHECK_EQ(pcount[2], static_cast<int>(rc::known_profile_finding_count(rc::ProfileDiagnosis::SmallTUndercoverage)));
    vcal::test::note("bootstrap: " + std::to_string(bcount[0]) + " PASS, " +
                     std::to_string(bcount[static_cast<int>(rc::Verdict::Conservative)]) + " CONSERVATIVE, " +
                     std::to_string(bcount[2]) + " KNOWN FINDING, " + std::to_string(bcount[3]) + " UNREVIEWED");
    VCAL_CHECK_EQ(bcount[static_cast<int>(rc::Verdict::Conservative)],
                  static_cast<int>(rc::known_bootstrap_finding_count(rc::BootstrapDiagnosis::Conservative)));
    VCAL_CHECK_EQ(bcount[2], static_cast<int>(rc::known_bootstrap_finding_count(rc::BootstrapDiagnosis::BoundaryBreakdown) +
                                              rc::known_bootstrap_finding_count(rc::BootstrapDiagnosis::NoBiasSkewCorrection)));
    // S-23: q's three verdict families, by the same rules; its acceptance checks on every scenario.
    int qp[rc::kVerdicts] = {}, qw[rc::kVerdicts] = {}, qb[rc::kVerdicts] = {};
    for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
        const auto& r = t.rows[id];
        const auto col = [&](const std::string& c) { return r[t.column(c)]; };
        const auto num = [&](const std::string& c) { return vcal::test::parse_int(col(c)); };
        const double R = static_cast<double>(rc::kReplicates);
        double plo = 0.0, phi = 0.0, lo = 0.0, hi = 0.0;
        rc::coverage_band(static_cast<std::int64_t>(rc::kReplicates), plo, phi);
        const std::int64_t unflagged = num("unflagged");
        rc::coverage_band(unflagged, lo, hi);
        VCAL_CHECK_EQ(num("q_profile_excludes_estimate"), 0);
        VCAL_CHECK(vcal::test::parse_double(col("q_profile_residual_max_hex")) <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
        VCAL_CHECK_EQ(vcal::test::parse_double(col("q_true_hex")), rc::q_truth(rc::scenario(id)));
        const double pcov = static_cast<double>(num("q_profile_covered")) / R;
        const auto* pr = rc::known_q_finding(rc::kKnownQProfileFindings, id);
        const rc::Verdict pv = rc::profile_verdict(pcov, plo, phi, pr);
        VCAL_CHECK_EQ(col("q_profile_verdict"), std::string(rc::verdict_text(pv)));
        ++qp[static_cast<int>(pv)];
        const double wcov = unflagged > 0 ? static_cast<double>(num("q_wald_covered")) / static_cast<double>(unflagged) : kNaN;
        const bool wr = rc::known_q_finding(rc::kKnownQWaldFindings, id) != nullptr;
        const rc::Verdict wv = rc::verdict(static_cast<double>(num("flagged")) / R, unflagged, wcov, lo, hi, wr);
        VCAL_CHECK_EQ(col("q_wald_verdict"), std::string(rc::verdict_text(wv)));
        ++qw[static_cast<int>(wv)];
        const double bcov = static_cast<double>(num("q_boot_covered")) / R;
        const auto* br = rc::known_q_finding(rc::kKnownQBootstrapFindings, id);
        const rc::Verdict bv = rc::bootstrap_verdict(bcov, plo, phi, br);
        VCAL_CHECK_EQ(col("q_boot_verdict"), std::string(rc::verdict_text(bv)));
        ++qb[static_cast<int>(bv)];
        const std::string at = " " + std::to_string(id) + "/q";
        if (pv == rc::Verdict::Unreviewed) unreviewed += at + " (profile coverage " + describe(pcov) + ")";
        if (wv == rc::Verdict::Unreviewed) unreviewed += at + " (Wald coverage " + describe(wcov) + ")";
        if (bv == rc::Verdict::Unreviewed) unreviewed += at + " (bootstrap coverage " + describe(bcov) + ")";
        if (pr && pv == rc::Verdict::Pass) stale += at + " (profile)";
        if (wr && wv != rc::Verdict::KnownFinding) stale += at + " (Wald)";
        if (br && bv == rc::Verdict::Pass) stale += at + " (bootstrap)";
    }
    const auto line = [](const int* c) {
        return std::to_string(c[0]) + " PASS, " + std::to_string(c[static_cast<int>(rc::Verdict::Conservative)]) +
               " CONSERVATIVE, " + std::to_string(c[2]) + " KNOWN FINDING, " + std::to_string(c[1]) + " DEFERRED, " +
               std::to_string(c[3]) + " UNREVIEWED";
    };
    vcal::test::note("q profile: " + line(qp) + "; q Wald: " + line(qw) + "; q bootstrap: " + line(qb));
    VCAL_CHECK_EQ(qp[static_cast<int>(rc::Verdict::Conservative)] + qp[2], static_cast<int>(rc::kKnownQProfileFindings.size()));
    VCAL_CHECK_EQ(qw[2], static_cast<int>(rc::kKnownQWaldFindings.size()));
    VCAL_CHECK_EQ(qb[static_cast<int>(rc::Verdict::Conservative)] + qb[2], static_cast<int>(rc::kKnownQBootstrapFindings.size()));
    if (!unreviewed.empty()) vcal::test::note("unreviewed, outside the band:" + unreviewed);
    if (!stale.empty()) vcal::test::note("in a reviewed list but not outside the band:" + stale);
    VCAL_CHECK(unreviewed.empty());
    VCAL_CHECK(stale.empty());
    VCAL_CHECK_EQ(count[2], static_cast<int>(sizeof rc::kKnownFindings / sizeof rc::kKnownFindings[0]));
    const auto bad = rc::rmse_violations(all);
    for (const auto& b : bad) vcal::test::note(b);
    VCAL_CHECK(bad.empty());
}

// --- replay: this platform re-fits the recorded replicates ----------------------------------------

// The panels the scipy cross-check reads (S-23) are the DGP's replay panels, count for count.
VCAL_TEST(recovery_replay_panels_match_the_dgp) {
    const auto t = vcal::test::read_golden_csv("recovery/replay_panels.csv");
    VCAL_REQUIRE(t.rows.size() == static_cast<std::size_t>(rc::kScenarios) * rc::kReplayReplicates);
    int mismatches = 0;
    for (const auto& r : t.rows) {
        const rc::Scenario s = rc::scenario(static_cast<std::uint32_t>(vcal::test::parse_int(r[t.column("scenario")])));
        const auto d = rc::panel(s, static_cast<std::uint32_t>(vcal::test::parse_int(r[t.column("replicate")])));
        std::string want;
        for (std::size_t k = 0; k < d.size(); ++k) want += (k ? " " : "") + std::to_string(d[k]);
        if (vcal::test::parse_int(r[t.column("n")]) != s.obligors || r[t.column("d")] != want) ++mismatches;
    }
    VCAL_CHECK_EQ(mismatches, 0);
}

VCAL_TEST(recovery_replay) {
    const auto t = vcal::test::read_golden_csv("recovery/replay.csv");
    VCAL_REQUIRE(t.rows.size() == static_cast<std::size_t>(rc::kScenarios) * rc::kReplayReplicates);
    std::vector<rc::Fit> got(t.rows.size());
    vcal::backends::CpuBackend{}.parallel_for(static_cast<std::int64_t>(t.rows.size()), [&](std::int64_t i) {
        const auto& r = t.rows[static_cast<std::size_t>(i)];
        const auto id = static_cast<std::uint32_t>(vcal::test::parse_int(r[t.column("scenario")]));
        const auto rep = static_cast<std::uint32_t>(vcal::test::parse_int(r[t.column("replicate")]));
        got[static_cast<std::size_t>(i)] = rc::fit(rc::scenario(id), rep);
    });
    double worst = 0.0, worst_se = 0.0, worst_profile = 0.0, worst_boot = 0.0, worst_q = 0.0, worst_q_se = 0.0;
    std::string worst_at, worst_se_at;
    int flag_mismatches = 0;
    for (std::size_t i = 0; i < t.rows.size(); ++i) {
        const auto& r = t.rows[i];
        const rc::Fit& f = got[i];
        const std::string at = "scenario " + r[t.column("scenario")] + " replicate " + r[t.column("replicate")];
        if (static_cast<std::int64_t>(f.flags) != vcal::test::parse_int(r[t.column("flags")])) ++flag_mismatches;
        if (static_cast<std::int64_t>(f.prof_flags[0] | (f.prof_flags[1] << 4)) !=
            vcal::test::parse_int(r[t.column("profile_flags")])) {
            ++flag_mismatches;
        }
        const double bends[4] = {f.boot_lo[0], f.boot_hi[0], f.boot_lo[1], f.boot_hi[1]};
        const char* bcols[4] = {"boot_pd_lo_hex", "boot_pd_hi_hex", "boot_rho_lo_hex", "boot_rho_hi_hex"};
        for (int k = 0; k < 4; ++k) {
            const double want = vcal::test::parse_double(r[t.column(bcols[k])]);
            if (std::isnan(want) && std::isnan(bends[k])) continue;
            worst_boot = std::fmax(worst_boot, std::fabs(bends[k] / want - 1.0));
        }
        if (static_cast<std::int64_t>(f.boot_edge) != vcal::test::parse_int(r[t.column("boot_edge")])) ++flag_mismatches;
        // S-23: q_hat and its SE, q's profile ends and flags, q's bootstrap ends.
        if (static_cast<std::int64_t>(f.q_prof_flags) != vcal::test::parse_int(r[t.column("q_profile_flags")])) {
            ++flag_mismatches;
        }
        const auto rel_to = [&](double value, const char* c) {
            const double want = vcal::test::parse_double(r[t.column(c)]);
            return std::isnan(want) && std::isnan(value) ? 0.0 : std::fabs(value / want - 1.0);
        };
        worst_q = std::fmax(worst_q, rel_to(f.q_hat, "q_hat_hex"));
        worst_q_se = std::fmax(worst_q_se, rel_to(f.q_se_s, "q_se_s_hex"));
        worst_profile = std::fmax(worst_profile, std::fmax(rel_to(f.q_prof_lo, "q_lo_hex"), rel_to(f.q_prof_hi, "q_hi_hex")));
        worst_boot = std::fmax(worst_boot, std::fmax(rel_to(f.q_boot_lo, "q_boot_lo_hex"), rel_to(f.q_boot_hi, "q_boot_hi_hex")));
        const double ends[4] = {f.prof_lo[0], f.prof_hi[0], f.prof_lo[1], f.prof_hi[1]};
        const char* end_cols[4] = {"pd_lo_hex", "pd_hi_hex", "rho_lo_hex", "rho_hi_hex"};
        for (int k = 0; k < 4; ++k) {
            const double want = vcal::test::parse_double(r[t.column(end_cols[k])]);
            if (std::isnan(want) && std::isnan(ends[k])) continue;
            worst_profile = std::fmax(worst_profile, std::fabs(ends[k] / want - 1.0));
        }
        const double values[5] = {f.value[0], f.value[1], f.se[0], f.se[1], f.loglik};
        const char* cols[5] = {"pd_hex", "rho_hex", "se_pd_hex", "se_rho_hex", "loglik_hex"};
        for (int k = 0; k < 5; ++k) {
            const double want = vcal::test::parse_double(r[t.column(cols[k])]);
            if (std::isnan(want) && std::isnan(values[k])) continue;
            const double rel = std::fabs(values[k] / want - 1.0);
            const bool is_se = k == 2 || k == 3;
            double& w = is_se ? worst_se : worst;
            std::string& w_at = is_se ? worst_se_at : worst_at;
            if (!(rel <= w)) {
                w = rel;
                w_at = std::string(cols[k]) + " of " + at;
            }
        }
    }
    vcal::test::note(std::to_string(t.rows.size()) + " replicates re-fitted: estimates and loglik worst relative "
                     "difference " + describe(worst) + (worst_at.empty() ? "" : " (" + worst_at + ")") + "; SEs " +
                     describe(worst_se) + (worst_se_at.empty() ? "" : " (" + worst_se_at + ")") +
                     "; q_hat " + describe(worst_q) + ", its SE " + describe(worst_q_se) +
                     "; profile endpoints (PD, rho, q) " + describe(worst_profile) + "; bootstrap ends " + describe(worst_boot) +
                     "; " + std::to_string(flag_mismatches) +
                     " flag mismatches");
    VCAL_CHECK_EQ(flag_mismatches, 0);
    VCAL_CHECK(worst <= tol::TOL_RECOVERY_REPLAY_REL);
    VCAL_CHECK(worst_se <= tol::TOL_RECOVERY_REPLAY_SE_REL);
    VCAL_CHECK(worst_q <= tol::TOL_RECOVERY_REPLAY_REL);
    VCAL_CHECK(worst_q_se <= tol::TOL_RECOVERY_REPLAY_SE_REL);
    VCAL_CHECK(worst_profile <= tol::TOL_PROFILE_CROSS_PLATFORM_REL);
    VCAL_CHECK(worst_boot <= tol::TOL_BOOTSTRAP_CROSS_PLATFORM_REL);
}

// A fit with only some arms (a study asking for what it scores) computes those arms exactly as a
// full fit does, and leaves the others NaN and flagged not computed.
VCAL_TEST(recovery_fit_arms_compute_only_what_is_requested) {
    const auto s = rc::scenario(37);
    const rc::Fit all = rc::fit(s, 3);
    const auto same = [](double x, double y) { return std::memcmp(&x, &y, sizeof x) == 0; };
    const rc::Fit prof = rc::fit(s, 3, nullptr, rc::kArmProfile);
    VCAL_CHECK(same(prof.prof_lo[0], all.prof_lo[0]) && same(prof.prof_hi[1], all.prof_hi[1]));
    VCAL_CHECK(prof.prof_flags[0] == all.prof_flags[0] && prof.prof_flags[1] == all.prof_flags[1]);
    VCAL_CHECK(std::isnan(prof.boot_lo[0]) && std::isnan(prof.q_prof_lo) && std::isnan(prof.q_se_s));
    VCAL_CHECK(prof.q_prof_flags == vcal::engine::kIntervalNotComputed);
    const rc::Fit boot = rc::fit(s, 3, nullptr, rc::kArmBootstrap);
    VCAL_CHECK(same(boot.boot_lo[0], all.boot_lo[0]) && same(boot.boot_hi[1], all.boot_hi[1]));
    VCAL_CHECK(boot.boot_edge == all.boot_edge);
    VCAL_CHECK(std::isnan(boot.prof_lo[0]) && boot.prof_flags[0] == vcal::engine::kIntervalNotComputed);
    VCAL_CHECK(std::isnan(boot.q_boot_lo));
    const rc::Fit pq = rc::fit(s, 3, nullptr, rc::kArmProfile | rc::kArmQ);
    VCAL_CHECK(same(pq.q_prof_lo, all.q_prof_lo) && same(pq.q_prof_hi, all.q_prof_hi) && same(pq.q_se_s, all.q_se_s));
    VCAL_CHECK(std::isnan(pq.q_boot_lo) && std::isnan(pq.boot_lo[1]));
    VCAL_CHECK(same(pq.value[0], all.value[0]) && same(pq.se[1], all.se[1]));  // the estimate is always there
}

// D-170 regression: scenario 20, replicate 388. q's lower end is held where the curve q = c meets
// the PD axis's edge, a bound that moves with q, so the envelope slope there is wrong. A root
// finder that stopped on a short Newton step returned a residual of 6.2e-4; convergence is now by
// bracketing.
VCAL_TEST(q_interval_end_on_a_moving_bound_meets_the_residual) {
    const rc::Fit f = rc::fit(rc::scenario(20), 388, nullptr, rc::kArmProfile | rc::kArmQ);
    VCAL_CHECK(f.q_prof_residual <= vcal::tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL);
    VCAL_CHECK(f.q_prof_lo <= f.q_hat && f.q_hat <= f.q_prof_hi);
    VCAL_CHECK(f.q_prof_flags & vcal::engine::kIntervalLowerBoxLimited);
}
