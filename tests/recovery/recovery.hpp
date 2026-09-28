// SPDX-License-Identifier: Apache-2.0
//
// M1.8 recovery harness (D-121): does the parity binomial-mixture MLE recover the parameters of
// the model that generated the data? Shared by the harness executable (which writes the goldens)
// and the tests (which replay a subset and enforce the exit criteria).
//
// Scenario matrix: PD in {0.1%, 1%, 5%} x rho in {0.02, 0.12, 0.24} x T in {20, 40, 100} x
// n in {100, 1000, 10000}, n constant over the T periods. Scenario id
//     id = ((i_pd * 3 + i_rho) * 3 + i_T) * 3 + i_n,
// is also the DGP scenario number, so replicate r of scenario id is the panel
// (kSeed, id, r, PD, rho, n, T) of dgp/ (D-109): reproducible from those numbers alone.
//
// Each replicate is fitted exactly as a user would: engine::calibrate with the parity rule and
// its doubled check (D-118), on the grid below, SEs from the Hessian at the estimate (D-119).
//
// Per scenario and parameter, three numbers, none of which quietly drops replicates (D-121):
//   1. bias and RMSE over ALL replicates (grid-edge estimates included, at the grid value);
//   2. the fraction of replicates without a reliable Wald interval ("flagged"): on the grid edge,
//      within 2 SEs of a bound (kFlagNearBound), or a flat surface;
//   3. coverage of the 95% Wald interval among the unflagged replicates, labelled conditional.
//      The interval is symmetric in the axis's scaled (logit) coordinate, where the Hessian is
//      taken: u_hat +- z_0.975 * se_u, se_u = se / (dv/du). It covers when it contains u_true.
// Two kinds of criterion, kept apart (owner, D-121 as amended):
//   - engine correctness, which must hold: RMSE falls with T; no quadrature or numeric flags;
//     replay agrees across platforms (and, in unit_xref, estimates and SEs agree with ref);
//   - interval coverage, a property of the Wald interval method rather than of the engine, which
//     is reported. Verdict per parameter:
//       DEFERRED       flagged fraction >= TOL_RECOVERY_MAX_FLAGGED_FRACTION: wait for M2's
//                      profile-likelihood intervals;
//       PASS           conditional coverage inside 0.95 +- TOL_RECOVERY_COVERAGE_BAND_Z *
//                      sqrt(0.95 * 0.05 / m), m the number of unflagged replicates;
//       KNOWN FINDING  outside the band, and reviewed: listed in kKnownFindings with its diagnosis;
//       UNREVIEWED     outside the band and not listed. The tests allow none, so any change in
//                      the set of out-of-band verdicts needs a review.
// M2 adds, side by side with Wald (D-131): the profile-likelihood interval (engine/profile.hpp),
// with coverage over ALL replicates (nothing deferred: a truncated interval is still an interval),
// verdict PASS / CONSERVATIVE / KNOWN FINDING / UNREVIEWED against kKnownProfileFindings; and the t(T-1) Wald
// interval, a `native` comparison reported without a verdict.
// Evidence for the diagnoses, per scenario and parameter among unflagged replicates, in the logit
// coordinate: sd(u_hat) against the root-mean-square Hessian SE. A ratio near 1 means the SEs
// match the actual spread of the estimates.
// S-23 adds a third quantity, the 99.9% conditional PD q (engine/conditional_pd.hpp), with three
// intervals and a verdict family each, under the same band and policy: the profile interval over
// ALL replicates; the delta-method Wald interval in logit(q), conditional on the same unflagged
// replicates as PD and rho, DEFERRED at the same flagged fraction; and the bootstrap percentile
// interval of q at the same B = 999 replicate estimates, over ALL replicates. Pre-registered in
// studies/derived-quantity-intervals/PREDICTION.md.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
#include "engine/conditional_pd.hpp"
#include "engine/profile.hpp"
#include "resample/bootstrap.hpp"
#include "resample/weights.hpp"
#include "tests/tolerances.hpp"

namespace vcal::recovery {

using Objective = objectives::BinomialMixture<PrecisionF64>;

inline constexpr std::uint64_t kSeed = 0x4D31385245434F56ull;  // "M18RECOV"; arbitrary, fixed
inline constexpr std::uint32_t kReplicates = 1000;            // R per scenario for the goldens
inline constexpr std::uint32_t kReplayReplicates = 2;         // replicates 0..1 replayed in CI

inline constexpr double kPds[3] = {0.001, 0.01, 0.05};
inline constexpr double kRhos[3] = {0.02, 0.12, 0.24};
inline constexpr std::int64_t kPeriods[3] = {20, 40, 100};
inline constexpr std::int64_t kObligors[3] = {100, 1000, 10000};
inline constexpr std::uint32_t kScenarios = 81;
// iid bootstrap of periods per panel (M2b, D-136): B replicates, percentile intervals.
inline constexpr std::uint32_t kBootstrapReplicates = 999;

// Each panel's bootstrap stream (in the resampling key domain, D-132).
inline std::uint64_t bootstrap_seed(std::uint32_t scenario, std::uint32_t replicate) {
    return kSeed ^ ((static_cast<std::uint64_t>(scenario) << 32) | replicate);
}

// Phi^-1(0.975), the two-sided 95% normal quantile (mpmath, 17 significant digits).
inline constexpr double kZ975 = 1.959963984540054;
// Student t 0.975 quantiles on T - 1 = 19, 39, 99 degrees of freedom (mpmath, 17 significant
// digits), for the t(T-1) Wald interval, a `native` comparison only (D-131).
inline constexpr double kT975[3] = {2.0930240544083098, 2.0226909200367611, 1.9842169515864175};
inline constexpr double kNominalCoverage = 0.95;

struct Scenario {
    std::uint32_t id;
    double pd;
    double rho;
    std::int64_t periods;
    std::int64_t obligors;
};

inline Scenario scenario(std::uint32_t id) {
    const std::uint32_t i_n = id % 3, i_t = (id / 3) % 3, i_rho = (id / 9) % 3, i_pd = id / 27;
    return {id, kPds[i_pd], kRhos[i_rho], kPeriods[i_t], kObligors[i_n]};
}

// The estimation grid: the M1.7 box (D-115). PD [1e-4, 0.2] and rho [1e-3, 0.5] (D-089), logit axes.
inline Grid<2> grid() {
    return {{{1e-4, 0.2, 61, AxisScale::Logit}, {1e-3, kDefaultRhoUpper, 41, AxisScale::Logit}}};
}

struct Fit {
    double value[2];  // PD, rho
    double se[2];
    double loglik;
    std::uint32_t flags;
    // Profile-likelihood 95% interval (M2, D-128..D-130), natural scale.
    double prof_lo[2];
    double prof_hi[2];
    std::uint32_t prof_flags[2];
    double prof_residual;  // max |P(endpoint) - (l_max - c)| over solved endpoints
    // iid bootstrap percentile interval (M2b, D-135, D-136), natural scale.
    double boot_lo[2];
    double boot_hi[2];
    std::uint32_t boot_edge;      // bootstrap replicates on the grid edge
    std::uint32_t boot_excluded;  // bootstrap replicates with non-finite estimates
    // S-23: the 99.9% conditional PD q.
    double q_hat;            // q at the estimate
    double q_se_s;           // delta-method SE of logit(q) (NaN where the Hessian SEs are)
    double q_prof_lo;        // profile-likelihood 95% interval, natural scale
    double q_prof_hi;
    std::uint32_t q_prof_flags;  // engine::kInterval*: truncated, box-limited, not computed
    double q_prof_residual;  // max |P_q(endpoint) - (l_max - c)| over solved endpoints
    double q_boot_lo;        // bootstrap percentile interval of q, natural scale
    double q_boot_hi;
};

// Replicate r of scenario s: the default counts d_t, n = s.obligors in every period (D-109).
inline std::vector<std::int64_t> panel(const Scenario& s, std::uint32_t replicate) {
    std::vector<std::int64_t> n(static_cast<std::size_t>(s.periods), s.obligors);
    std::vector<std::int64_t> d(n.size());
    const dgp::PanelSpec spec{kSeed, s.id, replicate, s.pd, s.rho, n.data(), s.periods};
    if (dgp::simulate_panel(spec, d.data(), nullptr) != dgp::Status::Ok) std::abort();
    return d;
}

// The interval arms a fit computes beyond the estimate (whose Wald and t intervals are free): a
// study requests only what it scores. Arms not requested are NaN, with their interval flags
// engine::kIntervalNotComputed, never a stale value. The q interval needs the profile. The recovery
// harness requests every arm.
enum Arms : std::uint32_t {
    kArmProfile = 1u << 0,    // PD and rho profile-likelihood intervals
    kArmBootstrap = 1u << 1,  // iid bootstrap percentile intervals (and q's bootstrap ends with kArmQ)
    kArmQ = 1u << 2,          // S-23: q's profile interval and its delta-method SE (needs kArmProfile)
    kArmAll = kArmProfile | kArmBootstrap | kArmQ,
};

// Simulates and fits one replicate. The fit is serial (the caller parallelises over replicates);
// results do not depend on the thread count (§6). With a row cache (D-167) the surface's rows come from
// it, shared across replicates and scenarios; the result is identical bit for bit. Each requested arm
// is computed exactly as with kArmAll.
inline Fit fit(const Scenario& s, std::uint32_t replicate, engine::SurfaceRowCache* cache = nullptr,
               std::uint32_t arms = kArmAll) {
    if ((arms & kArmQ) && !(arms & kArmProfile)) std::abort();  // q's interval is built on the profile
    const std::vector<std::int64_t> d = panel(s, replicate);
    std::vector<Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {s.obligors, d[t]};
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    std::vector<double> L;
    engine::Estimate2 est{};
    const Grid<2> g = grid();
    const auto status = cache != nullptr ? engine::calibrate_cached(backends::CpuBackend{1}, *cache, Objective{}, primary,
                                                                    check, obs.data(), s.periods, g, L, est)
                                         : engine::calibrate(backends::CpuBackend{1}, Objective{}, primary, check,
                                                             obs.data(), s.periods, g, L, est);
    if (status != engine::Status::Ok) std::abort();
    const double nan = std::nan("");
    engine::ProfileIntervals2 prof{{nan, nan}, {nan, nan}, {engine::kIntervalNotComputed, engine::kIntervalNotComputed},
                                   nan, {nan, nan}, nan, 0};
    if (arms & kArmProfile) prof = engine::profile_intervals(Objective{}, primary, obs.data(), s.periods, g, L, est);
    std::vector<resample::Replicate2> reps;
    std::uint32_t edge = 0;
    resample::PercentileInterval b0{nan, nan, 0, 0}, b1{nan, nan, 0, 0};
    if (arms & kArmBootstrap) {
        const auto idx = resample::bootstrap_indices(bootstrap_seed(s.id, replicate), resample::Scheme::IidBootstrap,
                                                     kBootstrapReplicates, s.periods);
        const auto W = resample::weights_from_indices(idx.data(), kBootstrapReplicates, s.periods, s.periods);
        reps.resize(kBootstrapReplicates);
        resample::replicate_estimates(backends::CpuBackend{1}, g, L.data(), s.periods, W.data(), kBootstrapReplicates,
                                      reps.data());
        for (const auto& r : reps) edge += (r.flags & engine::kFlagGridEdge) ? 1u : 0u;
        b0 = resample::percentile_interval(reps.data(), kBootstrapReplicates, 0);
        b1 = resample::percentile_interval(reps.data(), kBootstrapReplicates, 1);
    }
    // S-23: q's profile interval, its delta-method SE in logit(q), and q at each bootstrap estimate.
    engine::ConditionalPdInterval qi{nan, nan, nan, engine::kIntervalNotComputed, nan, 0};
    if (arms & kArmQ) qi = engine::conditional_pd_interval(Objective{}, primary, obs.data(), s.periods, g, L, est, prof);
    const auto qg = engine::conditional_pd_logit_gradient(est.value[0], est.value[1]);
    double se_u[2];
    for (int a = 0; a < 2; ++a) {
        se_u[a] = est.se[a] / grid::dvalue_dscaled(g.axis[a].scale, grid::to_scaled(g.axis[a].scale, est.value[a]));
    }
    const double q_var = qg.ds_du[0] * qg.ds_du[0] * se_u[0] * se_u[0] + qg.ds_du[1] * qg.ds_du[1] * se_u[1] * se_u[1] +
                         2.0 * qg.ds_du[0] * qg.ds_du[1] * est.corr * se_u[0] * se_u[1];
    std::vector<double> q_reps;
    q_reps.reserve(reps.size());
    for (const auto& r : reps) {
        const double q = std::isfinite(r.value[0]) && std::isfinite(r.value[1])
                             ? engine::conditional_pd(r.value[0], r.value[1])
                             : std::nan("");
        if (std::isfinite(q)) q_reps.push_back(q);
    }
    std::sort(q_reps.begin(), q_reps.end());
    if (!(arms & kArmQ)) q_reps.clear();  // q's bootstrap ends belong to the q arm
    const double q_boot_lo = q_reps.empty() ? nan : resample::quantile_type7(q_reps, 0.025);
    const double q_boot_hi = q_reps.empty() ? nan : resample::quantile_type7(q_reps, 0.975);
    return {{est.value[0], est.value[1]},
            {est.se[0], est.se[1]},
            est.loglik,
            est.flags,
            {prof.lo[0], prof.lo[1]},
            {prof.hi[0], prof.hi[1]},
            {prof.flags[0], prof.flags[1]},
            prof.residual_max,
            {b0.lo, b1.lo},
            {b0.hi, b1.hi},
            edge,
            static_cast<std::uint32_t>(b0.excluded > b1.excluded ? b0.excluded : b1.excluded),
            qi.estimate,
            (arms & kArmQ) && q_var >= 0.0 ? std::sqrt(q_var) : nan,
            qi.lo,
            qi.hi,
            qi.flags,
            qi.residual_max,
            q_boot_lo,
            q_boot_hi};
}

inline constexpr std::uint32_t kNoReliableInterval =
    engine::kFlagGridEdge | engine::kFlagFlatSurface | engine::kFlagNearBound;

// Does the 95% Wald interval, symmetric in axis a's scaled coordinate, contain the true value?
// `quantile` is z_0.975 (parity) or t_0.975 on T - 1 degrees of freedom (the native comparison).
inline bool covers(const Grid<2>& g, int a, const Fit& f, double truth, double quantile = kZ975) {
    const AxisScale sc = g.axis[a].scale;
    const double u_hat = grid::to_scaled(sc, f.value[a]);
    const double se_u = f.se[a] / grid::dvalue_dscaled(sc, u_hat);
    return std::fabs(u_hat - grid::to_scaled(sc, truth)) <= quantile * se_u;
}

// Does the profile-likelihood interval contain the true value? An interval that was not
// computed does not cover (it is counted, never dropped).
inline bool profile_covers(int a, const Fit& f, double truth) {
    if (f.prof_flags[a] & engine::kIntervalNotComputed) return false;
    return f.prof_lo[a] <= truth && truth <= f.prof_hi[a];
}

// Does the bootstrap percentile interval contain the true value? Non-finite ends do not cover.
inline bool bootstrap_covers(int a, const Fit& f, double truth) {
    return f.boot_lo[a] <= truth && truth <= f.boot_hi[a];
}

inline double t_quantile(std::int64_t periods) {
    for (int i = 0; i < 3; ++i) {
        if (kPeriods[i] == periods) return kT975[i];
    }
    return std::nan("");
}

// Conservative: reviewed and ABOVE the band (a safe-side property). KnownFinding: reviewed and
// outside it (for the profile interval, below it: a genuine limitation). D-124, D-131.
enum class Verdict { Pass, Deferred, KnownFinding, Unreviewed, Conservative };
inline constexpr int kVerdicts = 5;

inline const char* verdict_text(Verdict v) {
    switch (v) {
        case Verdict::Pass: return "PASS";
        case Verdict::Deferred: return "DEFERRED";
        case Verdict::KnownFinding: return "KNOWN FINDING";
        case Verdict::Unreviewed: return "UNREVIEWED";
        case Verdict::Conservative: return "CONSERVATIVE";
    }
    return "?";
}

// Reviewed out-of-band coverage verdicts (D-124), from the R = 1000 run of 2026-09-26. Each is a
// property of Wald intervals, not of the engine; the diagnostic ratio sd(u_hat) / rms(se_u) is in
// recovery_results.md. All are re-assessed with M2's profile-likelihood intervals.
// The last two are S-23's, for the delta-method interval of q.
enum class Diagnosis { SmallTUndercoverage, SkewedOvercoverage, BiasedLowUndercoverage, SeOverstatedOvercoverage };

inline const char* diagnosis_text(Diagnosis d) {
    switch (d) {
        case Diagnosis::SmallTUndercoverage:
            return "small-T Wald undercoverage: the Hessian SE understates the spread of the estimates";
        case Diagnosis::SkewedOvercoverage:
            return "overcoverage of a skewed estimate: the Hessian SE does not understate the spread (SE ratio <= 1), "
                   "but u_hat is skewed in the logit coordinate, so the symmetric interval is miscalibrated";
        case Diagnosis::BiasedLowUndercoverage:
            return "low estimates with narrow intervals: q-hat is biased low (rho-hat's downward bias carried into "
                   "q), and its delta-method SE moves with it (error and SE correlate at +0.64 to +0.88), so the "
                   "lowest estimates get the narrowest intervals and the misses fall almost all below the truth; "
                   "removing the mean bias alone would leave coverage at 0.916-0.947";
        case Diagnosis::SeOverstatedOvercoverage:
            return "overcoverage: the delta-method SE overstates the spread of logit(q-hat) (SE ratio 0.94)";
    }
    return "?";
}

struct KnownFinding {
    std::uint32_t scenario;
    int param;  // 0 = PD, 1 = rho
    Diagnosis diagnosis;
};

inline constexpr KnownFinding kKnownFindings[] = {
    {11, 0, Diagnosis::SmallTUndercoverage}, {29, 0, Diagnosis::SmallTUndercoverage},
    {37, 0, Diagnosis::SmallTUndercoverage}, {38, 0, Diagnosis::SmallTUndercoverage},
    {47, 0, Diagnosis::SmallTUndercoverage}, {49, 0, Diagnosis::SmallTUndercoverage},
    {55, 0, Diagnosis::SmallTUndercoverage}, {58, 0, Diagnosis::SmallTUndercoverage},
    {65, 0, Diagnosis::SmallTUndercoverage}, {68, 0, Diagnosis::SmallTUndercoverage},
    {73, 0, Diagnosis::SmallTUndercoverage}, {74, 0, Diagnosis::SmallTUndercoverage},
    {29, 1, Diagnosis::SmallTUndercoverage}, {74, 1, Diagnosis::SmallTUndercoverage},
    {13, 1, Diagnosis::SkewedOvercoverage}, {31, 1, Diagnosis::SkewedOvercoverage},
    {42, 1, Diagnosis::SkewedOvercoverage}, {51, 1, Diagnosis::SkewedOvercoverage},
    {63, 1, Diagnosis::SkewedOvercoverage},
};

inline const KnownFinding* known_finding(std::uint32_t scenario, int param) {
    for (const auto& k : kKnownFindings) {
        if (k.scenario == scenario && k.param == param) return &k;
    }
    return nullptr;
}

// Reviewed out-of-band verdicts for the profile-likelihood interval (D-131), from the R = 1000
// run of 2026-09-26. Coverage over ALL replicates, nothing deferred. Two kinds, each valid only on
// its own side of the band: a reviewed verdict that moves to the other side is UNREVIEWED again.
enum class ProfileDiagnosis { TruncationConservative, SmallTUndercoverage };

inline const char* profile_diagnosis_text(ProfileDiagnosis d) {
    return d == ProfileDiagnosis::TruncationConservative
               ? "conservative: near-uninformative data, many intervals truncated at a bound of the box (wider), "
                 "so they contain the truth more often than advertised"
               : "small-T undercoverage of the likelihood-ratio interval; for rho, misses are mostly below the "
                 "truth (the downward small-T bias of rho-hat)";
}

struct KnownProfileFinding {
    std::uint32_t scenario;
    int param;  // 0 = PD, 1 = rho
    ProfileDiagnosis diagnosis;
};

inline constexpr KnownProfileFinding kKnownProfileFindings[] = {
    // 20 conservative: n = 100 at PD 0.1%, or rho = 0.02 with little data; 25-100% truncated.
    {0, 0, ProfileDiagnosis::TruncationConservative}, {0, 1, ProfileDiagnosis::TruncationConservative},
    {1, 1, ProfileDiagnosis::TruncationConservative}, {3, 0, ProfileDiagnosis::TruncationConservative},
    {3, 1, ProfileDiagnosis::TruncationConservative}, {4, 1, ProfileDiagnosis::TruncationConservative},
    {6, 0, ProfileDiagnosis::TruncationConservative}, {6, 1, ProfileDiagnosis::TruncationConservative},
    {9, 0, ProfileDiagnosis::TruncationConservative}, {9, 1, ProfileDiagnosis::TruncationConservative},
    {12, 0, ProfileDiagnosis::TruncationConservative}, {12, 1, ProfileDiagnosis::TruncationConservative},
    {15, 1, ProfileDiagnosis::TruncationConservative}, {18, 0, ProfileDiagnosis::TruncationConservative},
    {18, 1, ProfileDiagnosis::TruncationConservative}, {21, 0, ProfileDiagnosis::TruncationConservative},
    {21, 1, ProfileDiagnosis::TruncationConservative}, {27, 1, ProfileDiagnosis::TruncationConservative},
    {30, 1, ProfileDiagnosis::TruncationConservative}, {33, 1, ProfileDiagnosis::TruncationConservative},
    // 6 undercovering, T = 20-40: coverage 0.920-0.927 against a band from 0.927.
    {29, 0, ProfileDiagnosis::SmallTUndercoverage}, {29, 1, ProfileDiagnosis::SmallTUndercoverage},
    {55, 1, ProfileDiagnosis::SmallTUndercoverage}, {68, 0, ProfileDiagnosis::SmallTUndercoverage},
    {72, 1, ProfileDiagnosis::SmallTUndercoverage}, {74, 1, ProfileDiagnosis::SmallTUndercoverage},
};

inline const KnownProfileFinding* known_profile_finding(std::uint32_t scenario, int param) {
    for (const auto& k : kKnownProfileFindings) {
        if (k.scenario == scenario && k.param == param) return &k;
    }
    return nullptr;
}

inline std::size_t known_profile_finding_count(ProfileDiagnosis d) {
    std::size_t n = 0;
    for (const auto& k : kKnownProfileFindings) n += k.diagnosis == d ? 1 : 0;
    return n;
}

// Reviewed out-of-band verdicts for the iid bootstrap percentile interval (D-137), from the
// R = 1000 run of 2026-09-26: 125 of 162, all below the band, against predictions registered
// beforehand (docs/methodology/bootstrap_predictions.md). Same rules as the profile list:
// coverage over all replicates, each kind valid only on its side of the band. Pinned per verdict
// for CI; the documents describe them at the method level.
enum class BootstrapDiagnosis { Conservative, BoundaryBreakdown, NoBiasSkewCorrection };

inline const char* bootstrap_diagnosis_text(BootstrapDiagnosis d) {
    switch (d) {
        case BootstrapDiagnosis::Conservative: return "conservative";
        case BootstrapDiagnosis::BoundaryBreakdown:
            return "boundary breakdown: near a bound of the box the bootstrap is inconsistent (Andrews 2000); "
                   "resampled estimates pile on the bound and intervals collapse";
        case BootstrapDiagnosis::NoBiasSkewCorrection:
            return "the percentile interval corrects neither the downward bias nor the skew of rho-hat (and at "
                   "small T the resampling variance is (T-1)/T of the true one)";
    }
    return "?";
}
struct KnownBootstrapFinding {
    std::uint32_t scenario;
    int param;
    BootstrapDiagnosis diagnosis;
};
inline constexpr KnownBootstrapFinding kKnownBootstrapFindings[] = {
    {0, 0, BootstrapDiagnosis::BoundaryBreakdown}, {0, 1, BootstrapDiagnosis::BoundaryBreakdown}, {1, 0, BootstrapDiagnosis::BoundaryBreakdown},
    {1, 1, BootstrapDiagnosis::BoundaryBreakdown}, {2, 0, BootstrapDiagnosis::BoundaryBreakdown}, {2, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {3, 0, BootstrapDiagnosis::BoundaryBreakdown}, {3, 1, BootstrapDiagnosis::BoundaryBreakdown}, {4, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {5, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {6, 1, BootstrapDiagnosis::BoundaryBreakdown}, {7, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {8, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {9, 0, BootstrapDiagnosis::BoundaryBreakdown}, {9, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {10, 0, BootstrapDiagnosis::BoundaryBreakdown}, {10, 1, BootstrapDiagnosis::BoundaryBreakdown}, {11, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {11, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {12, 0, BootstrapDiagnosis::BoundaryBreakdown}, {12, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {13, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {13, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {14, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {14, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {15, 0, BootstrapDiagnosis::BoundaryBreakdown}, {15, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {16, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {17, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {18, 0, BootstrapDiagnosis::BoundaryBreakdown},
    {18, 1, BootstrapDiagnosis::BoundaryBreakdown}, {19, 0, BootstrapDiagnosis::BoundaryBreakdown}, {19, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {20, 0, BootstrapDiagnosis::BoundaryBreakdown}, {20, 1, BootstrapDiagnosis::BoundaryBreakdown}, {21, 0, BootstrapDiagnosis::BoundaryBreakdown},
    {21, 1, BootstrapDiagnosis::BoundaryBreakdown}, {22, 0, BootstrapDiagnosis::BoundaryBreakdown}, {22, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {23, 0, BootstrapDiagnosis::BoundaryBreakdown}, {23, 1, BootstrapDiagnosis::BoundaryBreakdown}, {24, 0, BootstrapDiagnosis::BoundaryBreakdown},
    {24, 1, BootstrapDiagnosis::BoundaryBreakdown}, {25, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {25, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {27, 0, BootstrapDiagnosis::BoundaryBreakdown}, {27, 1, BootstrapDiagnosis::BoundaryBreakdown}, {28, 0, BootstrapDiagnosis::BoundaryBreakdown},
    {28, 1, BootstrapDiagnosis::BoundaryBreakdown}, {29, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {29, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {30, 0, BootstrapDiagnosis::BoundaryBreakdown}, {30, 1, BootstrapDiagnosis::BoundaryBreakdown}, {31, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {32, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {33, 1, BootstrapDiagnosis::BoundaryBreakdown}, {34, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {35, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {36, 0, BootstrapDiagnosis::BoundaryBreakdown}, {36, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {37, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {37, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {38, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {38, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {39, 0, BootstrapDiagnosis::BoundaryBreakdown}, {39, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {40, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {40, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {41, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {41, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {42, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {43, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {43, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {45, 0, BootstrapDiagnosis::BoundaryBreakdown}, {45, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {46, 0, BootstrapDiagnosis::BoundaryBreakdown}, {46, 1, BootstrapDiagnosis::BoundaryBreakdown}, {47, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {47, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {48, 0, BootstrapDiagnosis::BoundaryBreakdown}, {48, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {49, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {49, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {50, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {50, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {51, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {52, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {53, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {54, 0, BootstrapDiagnosis::BoundaryBreakdown}, {54, 1, BootstrapDiagnosis::BoundaryBreakdown},
    {55, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {55, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {56, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {57, 0, BootstrapDiagnosis::BoundaryBreakdown}, {57, 1, BootstrapDiagnosis::BoundaryBreakdown}, {58, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {59, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {60, 1, BootstrapDiagnosis::BoundaryBreakdown}, {62, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {63, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {63, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {64, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {64, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {65, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {65, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {66, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {67, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {68, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {68, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {70, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {72, 0, BootstrapDiagnosis::BoundaryBreakdown},
    {72, 1, BootstrapDiagnosis::BoundaryBreakdown}, {73, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {73, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {74, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {74, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {75, 0, BootstrapDiagnosis::NoBiasSkewCorrection},
    {75, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {76, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {76, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {77, 0, BootstrapDiagnosis::NoBiasSkewCorrection}, {77, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {78, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
    {79, 1, BootstrapDiagnosis::NoBiasSkewCorrection}, {80, 1, BootstrapDiagnosis::NoBiasSkewCorrection},
};

inline const KnownBootstrapFinding* known_bootstrap_finding(std::uint32_t scenario, int param) {
    for (const auto& k : kKnownBootstrapFindings) {
        if (k.scenario == scenario && k.param == param) return &k;
    }
    return nullptr;
}

inline std::size_t known_bootstrap_finding_count(BootstrapDiagnosis d) {
    std::size_t n = 0;
    for (const auto& k : kKnownBootstrapFindings) n += k.diagnosis == d ? 1 : 0;
    return n;
}

inline Verdict bootstrap_verdict(double coverage, double lo, double hi, const KnownBootstrapFinding* reviewed) {
    if (coverage >= lo && coverage <= hi) return Verdict::Pass;
    if (reviewed && reviewed->diagnosis == BootstrapDiagnosis::Conservative && coverage > hi) return Verdict::Conservative;
    if (reviewed && reviewed->diagnosis != BootstrapDiagnosis::Conservative && coverage < lo) return Verdict::KnownFinding;
    return Verdict::Unreviewed;
}

// --- S-23: the 99.9% conditional PD q -------------------------------------------------------------

inline double q_truth(const Scenario& s) { return engine::conditional_pd(s.pd, s.rho); }

// Does the delta-method Wald interval, symmetric in logit(q), contain the true q?
inline bool q_wald_covers(const Fit& f, double truth) {
    const double s_hat = grid::to_scaled(AxisScale::Logit, f.q_hat);
    return std::fabs(s_hat - grid::to_scaled(AxisScale::Logit, truth)) <= kZ975 * f.q_se_s;
}

// Does q's profile interval contain the true q? One that was not computed does not cover.
inline bool q_profile_covers(const Fit& f, double truth) {
    if (f.q_prof_flags & engine::kIntervalNotComputed) return false;
    return f.q_prof_lo <= truth && truth <= f.q_prof_hi;
}

inline bool q_bootstrap_covers(const Fit& f, double truth) { return f.q_boot_lo <= truth && truth <= f.q_boot_hi; }

// Reviewed out-of-band verdicts for q (S-23), with param = 2, one list per interval; same rules as
// the lists for PD and rho, each kind valid only on its own side of the band.
// Reviewed from the S-23 full-matrix run of 2026-09-27 (studies/README.md, S-23; D-156).
inline const std::vector<KnownProfileFinding> kKnownQProfileFindings = {
    // 15 conservative: near-uninformative data, the lower end box-limited in 52-100% of replicates.
    {0, 2, ProfileDiagnosis::TruncationConservative},  {1, 2, ProfileDiagnosis::TruncationConservative},
    {3, 2, ProfileDiagnosis::TruncationConservative},  {4, 2, ProfileDiagnosis::TruncationConservative},
    {6, 2, ProfileDiagnosis::TruncationConservative},  {9, 2, ProfileDiagnosis::TruncationConservative},
    {12, 2, ProfileDiagnosis::TruncationConservative}, {15, 2, ProfileDiagnosis::TruncationConservative},
    {18, 2, ProfileDiagnosis::TruncationConservative}, {21, 2, ProfileDiagnosis::TruncationConservative},
    {24, 2, ProfileDiagnosis::TruncationConservative}, {27, 2, ProfileDiagnosis::TruncationConservative},
    {30, 2, ProfileDiagnosis::TruncationConservative}, {33, 2, ProfileDiagnosis::TruncationConservative},
    {54, 2, ProfileDiagnosis::TruncationConservative},
    // 2 undercovering at T = 20, misses mostly below the truth, as for rho in the same scenarios.
    {72, 2, ProfileDiagnosis::SmallTUndercoverage}, {74, 2, ProfileDiagnosis::SmallTUndercoverage},
};
inline const std::vector<KnownFinding> kKnownQWaldFindings = {
    {11, 2, Diagnosis::BiasedLowUndercoverage}, {29, 2, Diagnosis::BiasedLowUndercoverage},
    {32, 2, Diagnosis::BiasedLowUndercoverage}, {35, 2, Diagnosis::BiasedLowUndercoverage},
    {37, 2, Diagnosis::BiasedLowUndercoverage}, {38, 2, Diagnosis::BiasedLowUndercoverage},
    {47, 2, Diagnosis::BiasedLowUndercoverage}, {49, 2, Diagnosis::BiasedLowUndercoverage},
    {55, 2, Diagnosis::BiasedLowUndercoverage}, {56, 2, Diagnosis::BiasedLowUndercoverage},
    {59, 2, Diagnosis::BiasedLowUndercoverage}, {64, 2, Diagnosis::BiasedLowUndercoverage},
    {65, 2, Diagnosis::BiasedLowUndercoverage}, {67, 2, Diagnosis::BiasedLowUndercoverage},
    {68, 2, Diagnosis::BiasedLowUndercoverage}, {73, 2, Diagnosis::BiasedLowUndercoverage},
    {74, 2, Diagnosis::BiasedLowUndercoverage}, {13, 2, Diagnosis::SeOverstatedOvercoverage},
};
// 73 below the band, none above. Each takes the diagnosis of rho's bootstrap finding in the same
// scenario (else PD's): q inherits them.
inline const std::vector<KnownBootstrapFinding> kKnownQBootstrapFindings = {
    {0, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {1, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {2, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {3, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {4, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {5, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {6, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {7, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {8, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {9, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {10, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {11, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {12, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {13, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {14, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {15, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {16, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {17, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {18, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {19, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {20, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {21, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {22, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {23, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {24, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {25, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {27, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {28, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {29, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {30, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {31, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {32, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {33, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {35, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {36, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {37, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {38, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {39, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {40, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {41, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {42, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {43, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {45, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {46, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {47, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {48, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {49, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {50, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {51, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {52, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {53, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {54, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {55, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {56, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {57, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {58, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {59, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {60, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {63, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {64, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {65, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {66, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {67, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {68, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {70, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {72, 2, BootstrapDiagnosis::BoundaryBreakdown},
    {73, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {74, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {75, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {76, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {77, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {78, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
    {80, 2, BootstrapDiagnosis::NoBiasSkewCorrection},
};

template <class K>
const K* known_q_finding(const std::vector<K>& list, std::uint32_t scenario) {
    for (const auto& k : list) {
        if (k.scenario == scenario) return &k;
    }
    return nullptr;
}

struct QSummary {
    double truth;
    double median_rel_err;  // median(q_hat / q - 1) over all replicates
    double rmse;            // sqrt(mean((q_hat - q)^2)) over all replicates
    // Profile-likelihood interval, over ALL replicates.
    std::int64_t profile_covered;
    std::int64_t profile_truncated;       // either end at the limit of q in the box
    std::int64_t profile_box_limited;     // either end box-limited (truncation included)
    std::int64_t profile_box_limited_lo;  // the lower end box-limited
    std::int64_t profile_not_computed;    // counted as not covering
    std::int64_t profile_miss_below;      // intervals entirely below the truth
    std::int64_t profile_miss_above;
    std::int64_t profile_excludes_estimate;  // intervals not containing q_hat: an acceptance check, none allowed
    double profile_log_width_median;         // median log(hi / lo) over computed intervals
    double profile_residual_max;
    double profile_coverage;
    Verdict profile_verdict;
    // Delta-method Wald in logit(q), among the unflagged replicates (as for PD and rho).
    std::int64_t wald_covered;
    std::int64_t wald_miss_below;
    std::int64_t wald_miss_above;
    double wald_coverage;
    Verdict wald_verdict;
    // Bootstrap percentile, over ALL replicates.
    std::int64_t boot_covered;
    std::int64_t boot_miss_below;
    std::int64_t boot_miss_above;
    double boot_coverage;
    Verdict boot_verdict;
};

inline double median_of(std::vector<double> x) {
    if (x.empty()) return std::nan("");
    std::sort(x.begin(), x.end());
    const std::size_t m = x.size() / 2;
    return x.size() % 2 == 1 ? x[m] : 0.5 * (x[m - 1] + x[m]);
}

struct ParamSummary {
    double mean;
    double bias;
    double bias_mcse;  // Monte Carlo SE of the bias: sd / sqrt(R)
    double rmse;
    std::int64_t covered;  // among the unflagged
    double coverage;       // covered / unflagged; NaN if none unflagged
    double sd_u;           // sd of u_hat among the unflagged (logit coordinate)
    double rms_se_u;       // root mean square of se_u among the unflagged
    Verdict verdict;
    // t(T-1) Wald interval, native comparison (D-131): coverage among the same unflagged replicates.
    std::int64_t t_covered;
    double t_coverage;
    // Profile-likelihood interval (D-128..D-131): coverage among ALL replicates.
    std::int64_t profile_covered;
    std::int64_t profile_truncated;     // replicates with either end at a bound of the box
    std::int64_t profile_not_computed;  // counted as not covering
    double profile_coverage;
    Verdict profile_verdict;  // PASS, KNOWN FINDING or UNREVIEWED; never DEFERRED
    // iid bootstrap percentile interval (D-136, D-137): coverage among ALL replicates.
    std::int64_t boot_covered;
    std::int64_t boot_degenerate;  // intervals of zero width (every bootstrap estimate equal)
    double boot_coverage;
    Verdict boot_verdict;
};

struct Summary {
    Scenario s;
    std::int64_t replicates;
    std::int64_t edge, near_bound, flat, rejected, quad_unconverged, numeric;
    std::int64_t flagged;    // edge | flat | near-bound: no reliable Wald interval
    std::int64_t unflagged;  // m
    double flagged_fraction;
    double band_lo, band_hi;  // Monte Carlo band for the conditional coverage
    double profile_band_lo, profile_band_hi;  // Monte Carlo band for all R replicates
    double profile_residual_max;              // worst self-reported endpoint residual, log-likelihood units
    double boot_edge_fraction;                // mean fraction of bootstrap replicates on the grid edge
    ParamSummary param[2];
    QSummary q;  // S-23
};

// Neumaier-compensated sum in the order given, so summaries do not depend on anything but the fits.
struct CompensatedSum {
    double sum = 0.0, comp = 0.0;
    void add(double x) {
        const double t = sum + x;
        comp += std::fabs(sum) >= std::fabs(x) ? (sum - t) + x : (x - t) + sum;
        sum = t;
    }
    double value() const { return sum + comp; }
};

inline Verdict verdict(double flagged_fraction, std::int64_t unflagged, double coverage, double lo, double hi,
                       bool reviewed) {
    if (!(flagged_fraction < tol::TOL_RECOVERY_MAX_FLAGGED_FRACTION) || unflagged == 0) return Verdict::Deferred;
    if (coverage >= lo && coverage <= hi) return Verdict::Pass;
    return reviewed ? Verdict::KnownFinding : Verdict::Unreviewed;
}

inline Verdict profile_verdict(double coverage, double lo, double hi, const KnownProfileFinding* reviewed) {
    if (coverage >= lo && coverage <= hi) return Verdict::Pass;
    if (reviewed && reviewed->diagnosis == ProfileDiagnosis::TruncationConservative && coverage > hi) {
        return Verdict::Conservative;
    }
    if (reviewed && reviewed->diagnosis == ProfileDiagnosis::SmallTUndercoverage && coverage < lo) {
        return Verdict::KnownFinding;
    }
    return Verdict::Unreviewed;
}

inline void coverage_band(std::int64_t m, double& lo, double& hi) {
    const double half = m > 0 ? tol::TOL_RECOVERY_COVERAGE_BAND_Z *
                                    std::sqrt(kNominalCoverage * (1.0 - kNominalCoverage) / static_cast<double>(m))
                              : std::nan("");
    lo = kNominalCoverage - half;
    hi = kNominalCoverage + half;
}

inline Summary summarise(const Scenario& s, const Fit* fits, std::int64_t R) {
    const Grid<2> g = grid();
    const double truth[2] = {s.pd, s.rho};
    Summary out{};
    out.s = s;
    out.replicates = R;
    for (std::int64_t r = 0; r < R; ++r) {
        const std::uint32_t f = fits[r].flags;
        out.edge += (f & engine::kFlagGridEdge) ? 1 : 0;
        out.near_bound += (f & engine::kFlagNearBound) ? 1 : 0;
        out.flat += (f & engine::kFlagFlatSurface) ? 1 : 0;
        out.rejected += (f & engine::kFlagRefinementRejected) ? 1 : 0;
        out.quad_unconverged += (f & engine::kFlagQuadratureUnconverged) ? 1 : 0;
        out.numeric += (f & engine::kFlagNumeric) ? 1 : 0;
        out.flagged += (f & kNoReliableInterval) ? 1 : 0;
    }
    out.unflagged = R - out.flagged;
    out.flagged_fraction = static_cast<double>(out.flagged) / static_cast<double>(R);
    coverage_band(out.unflagged, out.band_lo, out.band_hi);
    coverage_band(R, out.profile_band_lo, out.profile_band_hi);
    for (std::int64_t r = 0; r < R; ++r) out.profile_residual_max = std::fmax(out.profile_residual_max, fits[r].prof_residual);
    {
        CompensatedSum edge;
        for (std::int64_t r = 0; r < R; ++r) edge.add(static_cast<double>(fits[r].boot_edge) / static_cast<double>(kBootstrapReplicates));
        out.boot_edge_fraction = edge.value() / static_cast<double>(R);
    }
    const double tq = t_quantile(s.periods);
    for (int a = 0; a < 2; ++a) {
        CompensatedSum sum;
        for (std::int64_t r = 0; r < R; ++r) sum.add(fits[r].value[a]);
        const double mean = sum.value() / static_cast<double>(R);
        CompensatedSum dev2, err2, u_sum, se2_sum;
        std::int64_t covered = 0, t_covered = 0, prof_covered = 0, prof_truncated = 0, prof_missing = 0;
        std::int64_t boot_covered = 0, boot_degenerate = 0;
        const AxisScale sc = g.axis[a].scale;
        for (std::int64_t r = 0; r < R; ++r) {
            const double v = fits[r].value[a];
            dev2.add((v - mean) * (v - mean));
            err2.add((v - truth[a]) * (v - truth[a]));
            if (profile_covers(a, fits[r], truth[a])) ++prof_covered;
            if (fits[r].prof_flags[a] & (engine::kIntervalLowerTruncated | engine::kIntervalUpperTruncated)) ++prof_truncated;
            if (fits[r].prof_flags[a] & engine::kIntervalNotComputed) ++prof_missing;
            if (bootstrap_covers(a, fits[r], truth[a])) ++boot_covered;
            if (fits[r].boot_lo[a] == fits[r].boot_hi[a]) ++boot_degenerate;
            if (fits[r].flags & kNoReliableInterval) continue;
            if (covers(g, a, fits[r], truth[a])) ++covered;
            if (covers(g, a, fits[r], truth[a], tq)) ++t_covered;
            const double u = grid::to_scaled(sc, v);
            const double se_u = fits[r].se[a] / grid::dvalue_dscaled(sc, u);
            u_sum.add(u);
            se2_sum.add(se_u * se_u);
        }
        const double m = static_cast<double>(out.unflagged);
        const double u_mean = u_sum.value() / m;
        CompensatedSum u_dev2;
        for (std::int64_t r = 0; r < R; ++r) {
            if (fits[r].flags & kNoReliableInterval) continue;
            const double u = grid::to_scaled(sc, fits[r].value[a]);
            u_dev2.add((u - u_mean) * (u - u_mean));
        }
        ParamSummary& p = out.param[a];
        p.mean = mean;
        p.bias = mean - truth[a];
        p.bias_mcse = R > 1 ? std::sqrt(dev2.value() / static_cast<double>(R - 1) / static_cast<double>(R)) : std::nan("");
        p.rmse = std::sqrt(err2.value() / static_cast<double>(R));
        p.covered = covered;
        p.coverage = out.unflagged > 0 ? static_cast<double>(covered) / m : std::nan("");
        p.sd_u = out.unflagged > 1 ? std::sqrt(u_dev2.value() / (m - 1.0)) : std::nan("");
        p.rms_se_u = out.unflagged > 0 ? std::sqrt(se2_sum.value() / m) : std::nan("");
        p.verdict = verdict(out.flagged_fraction, out.unflagged, p.coverage, out.band_lo, out.band_hi,
                            known_finding(s.id, a) != nullptr);
        p.t_covered = t_covered;
        p.t_coverage = out.unflagged > 0 ? static_cast<double>(t_covered) / m : std::nan("");
        p.profile_covered = prof_covered;
        p.profile_truncated = prof_truncated;
        p.profile_not_computed = prof_missing;
        p.profile_coverage = static_cast<double>(prof_covered) / static_cast<double>(R);
        p.profile_verdict = profile_verdict(p.profile_coverage, out.profile_band_lo, out.profile_band_hi,
                                            known_profile_finding(s.id, a));
        p.boot_covered = boot_covered;
        p.boot_degenerate = boot_degenerate;
        p.boot_coverage = static_cast<double>(boot_covered) / static_cast<double>(R);
        p.boot_verdict = bootstrap_verdict(p.boot_coverage, out.profile_band_lo, out.profile_band_hi,
                                           known_bootstrap_finding(s.id, a));
    }

    // S-23: q.
    QSummary& q = out.q;
    q.truth = q_truth(s);
    std::vector<double> rel_err, log_width;
    CompensatedSum q_err2;
    for (std::int64_t r = 0; r < R; ++r) {
        const Fit& f = fits[r];
        rel_err.push_back(f.q_hat / q.truth - 1.0);
        q_err2.add((f.q_hat - q.truth) * (f.q_hat - q.truth));
        const std::uint32_t pf = f.q_prof_flags;
        if (q_profile_covers(f, q.truth)) ++q.profile_covered;
        if (pf & (engine::kIntervalLowerTruncated | engine::kIntervalUpperTruncated)) ++q.profile_truncated;
        if (pf & (engine::kIntervalLowerBoxLimited | engine::kIntervalUpperBoxLimited)) ++q.profile_box_limited;
        if (pf & engine::kIntervalLowerBoxLimited) ++q.profile_box_limited_lo;
        if (pf & engine::kIntervalNotComputed) {
            ++q.profile_not_computed;
        } else {
            if (f.q_prof_hi < q.truth) ++q.profile_miss_below;
            if (f.q_prof_lo > q.truth) ++q.profile_miss_above;
            if (!(f.q_prof_lo <= f.q_hat && f.q_hat <= f.q_prof_hi)) ++q.profile_excludes_estimate;
            log_width.push_back(std::log(f.q_prof_hi / f.q_prof_lo));
        }
        q.profile_residual_max = std::fmax(q.profile_residual_max, f.q_prof_residual);
        if (q_bootstrap_covers(f, q.truth)) ++q.boot_covered;
        if (f.q_boot_hi < q.truth) ++q.boot_miss_below;
        if (f.q_boot_lo > q.truth) ++q.boot_miss_above;
        if (f.flags & kNoReliableInterval) continue;
        if (q_wald_covers(f, q.truth)) {
            ++q.wald_covered;
        } else {
            const double s_hat = grid::to_scaled(AxisScale::Logit, f.q_hat);
            (s_hat < grid::to_scaled(AxisScale::Logit, q.truth) ? q.wald_miss_below : q.wald_miss_above) += 1;
        }
    }
    q.median_rel_err = median_of(rel_err);
    q.rmse = std::sqrt(q_err2.value() / static_cast<double>(R));
    q.profile_log_width_median = median_of(log_width);
    q.profile_coverage = static_cast<double>(q.profile_covered) / static_cast<double>(R);
    q.profile_verdict = profile_verdict(q.profile_coverage, out.profile_band_lo, out.profile_band_hi,
                                        known_q_finding(kKnownQProfileFindings, s.id));
    q.wald_coverage = out.unflagged > 0 ? static_cast<double>(q.wald_covered) / static_cast<double>(out.unflagged)
                                        : std::nan("");
    q.wald_verdict = verdict(out.flagged_fraction, out.unflagged, q.wald_coverage, out.band_lo, out.band_hi,
                             known_q_finding(kKnownQWaldFindings, s.id) != nullptr);
    q.boot_coverage = static_cast<double>(q.boot_covered) / static_cast<double>(R);
    q.boot_verdict = bootstrap_verdict(q.boot_coverage, out.profile_band_lo, out.profile_band_hi,
                                       known_q_finding(kKnownQBootstrapFindings, s.id));
    return out;
}

// Exit criterion: RMSE falls with T, strictly, for each (PD, rho, n) and parameter. Returns one
// line per violation (empty when the criterion holds). `all` is indexed by scenario id.
inline std::vector<std::string> rmse_violations(const std::vector<Summary>& all) {
    static const char* names[2] = {"PD", "rho"};
    std::vector<std::string> out;
    for (std::uint32_t i_pd = 0; i_pd < 3; ++i_pd) {
        for (std::uint32_t i_rho = 0; i_rho < 3; ++i_rho) {
            for (std::uint32_t i_n = 0; i_n < 3; ++i_n) {
                for (int a = 0; a < 2; ++a) {
                    for (std::uint32_t i_t = 0; i_t + 1 < 3; ++i_t) {
                        const std::uint32_t id = ((i_pd * 3 + i_rho) * 3 + i_t) * 3 + i_n;
                        const Summary& shorter = all[id];
                        const Summary& longer = all[id + 3];
                        if (!(longer.param[a].rmse < shorter.param[a].rmse)) {
                            out.push_back(std::string(names[a]) + " RMSE does not fall from T = " +
                                          std::to_string(shorter.s.periods) + " to T = " +
                                          std::to_string(longer.s.periods) + " (scenarios " + std::to_string(id) +
                                          ", " + std::to_string(id + 3) + ")");
                        }
                    }
                }
            }
        }
    }
    return out;
}

}  // namespace vcal::recovery
