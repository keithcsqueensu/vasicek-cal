// SPDX-License-Identifier: Apache-2.0
//
// M1.8 recovery harness executable (D-121, D-124). Fits every replicate of every scenario in
// tests/recovery/recovery.hpp and summarises them.
//
//   recovery_harness [--replicates R]   print the summary for R replicates (default kReplicates)
//   recovery_harness --write            run R = kReplicates and write the goldens:
//       tests/golden/recovery/summary.csv    per-scenario summary (exact hex + decimal for readers)
//       tests/golden/recovery/replay.csv     replicates 0..kReplayReplicates-1 of every scenario,
//                                            re-fitted in CI (unit_recovery: recovery_replay)
//       tests/golden/recovery/replay_panels.csv  those replicates' default counts, for the scipy
//                                            cross-check of q's profile interval (S-23)
//       tests/golden/recovery/MANIFEST.json  provenance
//       docs/methodology/recovery_results.md the tables, generated; do not edit by hand
//     then checks that every file exists and is non-empty.
//   --replay-dir DIR                    also write replay.csv and replay_panels.csv to DIR (every
//                                        scenario; any R >= kReplayReplicates), e.g. for the scipy
//                                        cross-check before a full run.
//   --summary-out FILE                  also write the summary (summary.csv's format) to FILE, one row
//                                        per scenario run (with --scenarios, those only), --load-fits
//                                        included: for reviewing verdicts before a --write run, and for
//                                        studies at other R (S-13).
//   --results-md FILE                   also write the results page (recovery_results.md's format) to
//                                        FILE, from any run over every scenario, --load-fits included:
//                                        the page is a function of the fits and the reviewed labels,
//                                        so a change of labels or diagnoses needs no refit.
//   --scenarios ID,ID,...               fit and summarise only these scenarios (exploration, e.g. the
//                                        study subset of D-150); not with --write or --check.
//   --save-fits FILE / --load-fits FILE  developer convenience: write every fit to a binary file,
//                                        or read them back instead of refitting, to re-summarise
//                                        after a review. Goldens are only committed from a run that
//                                        fitted (the harness refuses --write with --load-fits).
//   recovery_harness --check            run R = kReplicates and compare with the committed summary:
//       on the platform that wrote it (same compiler string as MANIFEST.json), every count and
//       verdict exactly and every value bitwise; elsewhere, counts within
//       TOL_RECOVERY_CROSS_PLATFORM_COUNT replicates, values at the replay tolerances
//       (TOL_RECOVERY_REPLAY_*), and verdicts wherever the counts they rest on agree exactly.
//
// Replicates run in parallel (CpuBackend over (replicate, scenario) jobs); each fit is serial,
// and results do not depend on the thread count, so the goldens do not either.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tests/harness/golden.hpp"
#include "tests/recovery/recovery.hpp"

#define VCAL_STR2(x) #x
#define VCAL_STR(x) VCAL_STR2(x)
#if defined(_MSC_VER)
#define VCAL_COMPILER "MSVC " VCAL_STR(_MSC_VER)
#elif defined(__clang__)
#define VCAL_COMPILER "Clang " __clang_version__
#elif defined(__GNUC__)
#define VCAL_COMPILER "GCC " __VERSION__
#else
#define VCAL_COMPILER "unknown compiler"
#endif

namespace {

namespace rc = vcal::recovery;
namespace tol = vcal::tol;
using vcal::test::to_hex;

std::string fmt(const char* format, double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, format, v);
    return buf;
}

std::string pct(double v) { return fmt("%.1f%%", 100.0 * v); }

std::string seed_hex() {
    char buf[32];
    std::snprintf(buf, sizeof buf, "0x%016llX", static_cast<unsigned long long>(rc::kSeed));
    return buf;
}

const char* kParamName[2] = {"pd", "rho"};
const char* kParamLabel[2] = {"PD", "ρ"};

double ratio(const rc::ParamSummary& p) { return p.sd_u / p.rms_se_u; }

// --- writing ---------------------------------------------------------------------------------------

// Opens `path` for writing, creating its directory; a file that cannot be written is an error,
// never a silent no-op.
std::ofstream open_out(const std::string& path) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot write " + path);
    return stream;
}

void close_out(std::ofstream& out, const std::string& path) {
    out.close();
    if (!out) throw std::runtime_error("error writing " + path);
}

void write_summary_csv(const std::string& path, const std::vector<rc::Summary>& all) {
    std::ofstream out = open_out(path);
    out << "# Recovery summary (M1.8 D-121, D-124; M2 D-131): one row per scenario, R replicates each.\n"
           "# Written by recovery_harness --write (" VCAL_COMPILER "); provenance in MANIFEST.json.\n"
           "# Bias and RMSE use all replicates. Coverage, sd_u and rms_se_u are conditional on the replicate\n"
           "# being unflagged (not on the grid edge, not within 2 SEs of a bound, surface not flat); sd_u and\n"
           "# rms_se_u are in the logit coordinate. t_*: the t(T-1) Wald interval on the same unflagged\n"
           "# replicates (native comparison). profile_*: the profile-likelihood interval over ALL replicates.\n"
           "# boot_*: the iid bootstrap percentile interval (B = 999) over ALL replicates.\n"
           "# q_*: the 99.9% conditional PD (S-23): profile interval and bootstrap over ALL replicates, the\n"
           "# delta-method Wald interval in logit(q) on the unflagged replicates; miss_below / miss_above count\n"
           "# intervals entirely below / above the true q; log_width is log(hi / lo) of the profile interval.\n"
           "# Doubles in *_hex columns are exact; the decimal columns are for readers.\n"
           "scenario,pd,rho,periods,obligors,replicates,edge,near_bound,flat,rejected,quad_unconverged,numeric,"
           "flagged,unflagged,band_lo_hex,band_hi_hex,profile_band_lo_hex,profile_band_hi_hex,"
           "profile_residual_max_hex,boot_edge_fraction_hex";
    for (const char* p : kParamName) {
        out << ',' << p << "_mean_hex," << p << "_bias_hex," << p << "_bias_mcse_hex," << p << "_rmse_hex," << p
            << "_sd_u_hex," << p << "_rms_se_u_hex," << p << "_covered," << p << "_verdict," << p << "_bias," << p
            << "_rmse," << p << "_coverage," << p << "_se_ratio," << p << "_t_covered," << p << "_profile_covered,"
            << p << "_profile_truncated," << p << "_profile_not_computed," << p << "_profile_verdict," << p
            << "_t_coverage," << p << "_profile_coverage," << p << "_boot_covered," << p << "_boot_degenerate," << p
            << "_boot_verdict," << p << "_boot_coverage";
    }
    out << ",q_true_hex,q_median_rel_err_hex,q_rmse_hex,q_profile_covered,q_profile_truncated,q_profile_box_limited,"
           "q_profile_box_limited_lo,q_profile_not_computed,q_profile_miss_below,q_profile_miss_above,"
           "q_profile_excludes_estimate,q_profile_log_width_median_hex,q_profile_residual_max_hex,q_profile_verdict,"
           "q_wald_covered,q_wald_miss_below,q_wald_miss_above,q_wald_verdict,q_boot_covered,q_boot_miss_below,"
           "q_boot_miss_above,q_boot_verdict,q_true,q_median_rel_err,q_profile_coverage,q_profile_width_ratio_median,"
           "q_wald_coverage,q_boot_coverage";
    out << '\n';
    for (const auto& s : all) {
        out << s.s.id << ',' << fmt("%.17g", s.s.pd) << ',' << fmt("%.17g", s.s.rho) << ',' << s.s.periods << ','
            << s.s.obligors << ',' << s.replicates << ',' << s.edge << ',' << s.near_bound << ',' << s.flat << ','
            << s.rejected << ',' << s.quad_unconverged << ',' << s.numeric << ',' << s.flagged << ',' << s.unflagged
            << ',' << to_hex(s.band_lo) << ',' << to_hex(s.band_hi) << ',' << to_hex(s.profile_band_lo) << ','
            << to_hex(s.profile_band_hi) << ',' << to_hex(s.profile_residual_max) << ',' << to_hex(s.boot_edge_fraction);
        for (const auto& p : s.param) {
            out << ',' << to_hex(p.mean) << ',' << to_hex(p.bias) << ',' << to_hex(p.bias_mcse) << ','
                << to_hex(p.rmse) << ',' << to_hex(p.sd_u) << ',' << to_hex(p.rms_se_u) << ',' << p.covered << ','
                << rc::verdict_text(p.verdict) << ',' << fmt("%.6g", p.bias) << ',' << fmt("%.6g", p.rmse) << ','
                << fmt("%.4f", p.coverage) << ',' << fmt("%.4f", ratio(p)) << ',' << p.t_covered << ','
                << p.profile_covered << ',' << p.profile_truncated << ',' << p.profile_not_computed << ','
                << rc::verdict_text(p.profile_verdict) << ',' << fmt("%.4f", p.t_coverage) << ','
                << fmt("%.4f", p.profile_coverage) << ',' << p.boot_covered << ',' << p.boot_degenerate << ','
                << rc::verdict_text(p.boot_verdict) << ',' << fmt("%.4f", p.boot_coverage);
        }
        const auto& q = s.q;
        out << ',' << to_hex(q.truth) << ',' << to_hex(q.median_rel_err) << ',' << to_hex(q.rmse) << ','
            << q.profile_covered << ',' << q.profile_truncated << ',' << q.profile_box_limited << ','
            << q.profile_box_limited_lo << ',' << q.profile_not_computed << ',' << q.profile_miss_below << ','
            << q.profile_miss_above << ',' << q.profile_excludes_estimate << ',' << to_hex(q.profile_log_width_median)
            << ',' << to_hex(q.profile_residual_max) << ',' << rc::verdict_text(q.profile_verdict) << ','
            << q.wald_covered << ',' << q.wald_miss_below << ',' << q.wald_miss_above << ','
            << rc::verdict_text(q.wald_verdict) << ',' << q.boot_covered << ',' << q.boot_miss_below << ','
            << q.boot_miss_above << ',' << rc::verdict_text(q.boot_verdict) << ',' << fmt("%.6g", q.truth) << ','
            << fmt("%.4f", q.median_rel_err) << ',' << fmt("%.4f", q.profile_coverage) << ','
            << fmt("%.4f", std::exp(q.profile_log_width_median)) << ',' << fmt("%.4f", q.wald_coverage) << ','
            << fmt("%.4f", q.boot_coverage);
        out << '\n';
    }
    close_out(out, path);
}

void write_replay_csv(const std::string& path, const std::vector<rc::Fit>& fits) {
    std::ofstream out = open_out(path);
    out << "# M1.8 replay values (D-121): replicates 0.." << rc::kReplayReplicates - 1
        << " of every scenario, as fitted by recovery_harness --write (" VCAL_COMPILER ").\n"
           "# Re-fitted in CI by unit_recovery (recovery_replay) and compared at TOL_RECOVERY_REPLAY_REL\n"
           "# (estimates, loglik) and TOL_RECOVERY_REPLAY_SE_REL (SEs); flags exactly.\n"
           "# Profile endpoints at TOL_PROFILE_CROSS_PLATFORM_REL; profile flags exactly. Bootstrap interval\n"
           "# ends at TOL_BOOTSTRAP_CROSS_PLATFORM_REL. S-23's q: q_hat at TOL_RECOVERY_REPLAY_REL, its delta-method\n"
           "# SE at TOL_RECOVERY_REPLAY_SE_REL, profile ends and flags as above, bootstrap ends as above.\n"
           "scenario,replicate,pd_hex,rho_hex,se_pd_hex,se_rho_hex,loglik_hex,flags,pd_lo_hex,pd_hi_hex,rho_lo_hex,"
           "rho_hi_hex,profile_flags,boot_pd_lo_hex,boot_pd_hi_hex,boot_rho_lo_hex,boot_rho_hi_hex,boot_edge,"
           "q_hat_hex,q_se_s_hex,q_lo_hex,q_hi_hex,q_profile_flags,q_boot_lo_hex,q_boot_hi_hex\n";
    for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
        for (std::uint32_t r = 0; r < rc::kReplayReplicates; ++r) {
            const rc::Fit& f = fits[static_cast<std::size_t>(r) * rc::kScenarios + id];
            out << id << ',' << r << ',' << to_hex(f.value[0]) << ',' << to_hex(f.value[1]) << ',' << to_hex(f.se[0])
                << ',' << to_hex(f.se[1]) << ',' << to_hex(f.loglik) << ',' << f.flags << ',' << to_hex(f.prof_lo[0])
                << ',' << to_hex(f.prof_hi[0]) << ',' << to_hex(f.prof_lo[1]) << ',' << to_hex(f.prof_hi[1]) << ','
                << (f.prof_flags[0] | (f.prof_flags[1] << 4)) << ',' << to_hex(f.boot_lo[0]) << ','
                << to_hex(f.boot_hi[0]) << ',' << to_hex(f.boot_lo[1]) << ',' << to_hex(f.boot_hi[1]) << ','
                << f.boot_edge << ',' << to_hex(f.q_hat) << ',' << to_hex(f.q_se_s) << ',' << to_hex(f.q_prof_lo) << ','
                << to_hex(f.q_prof_hi) << ',' << f.q_prof_flags << ',' << to_hex(f.q_boot_lo) << ','
                << to_hex(f.q_boot_hi) << '\n';
        }
    }
    close_out(out, path);
}

void write_replay_panels_csv(const std::string& path) {
    std::ofstream out = open_out(path);
    out << "# The recovery replay panels (S-23): replicates 0.." << rc::kReplayReplicates - 1
        << " of every scenario, from dgp/ (seed, scenario, replicate).\n"
           "# n obligors in every period; d is the default count per period, in period order, space-separated.\n"
           "# Read by validation/scipy/conditional_pd_profile.py; unit_recovery (recovery_replay) checks it against dgp/.\n"
           "scenario,replicate,n,d\n";
    for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
        for (std::uint32_t r = 0; r < rc::kReplayReplicates; ++r) {
            const rc::Scenario s = rc::scenario(id);
            out << id << ',' << r << ',' << s.obligors << ',';
            const auto d = rc::panel(s, r);
            for (std::size_t k = 0; k < d.size(); ++k) out << (k ? " " : "") << d[k];
            out << '\n';
        }
    }
    close_out(out, path);
}

void write_manifest(const std::string& path, double seconds, std::uint32_t R) {
    const auto g = rc::grid();
    std::ofstream out = open_out(path);
    out << "{\n"
           "  \"description\": \"M1.8 recovery harness goldens (D-121, D-124)\",\n"
           "  \"generated_by\": \"recovery_harness --write\",\n"
           "  \"compiler\": \"" VCAL_COMPILER "\",\n"
           "  \"dgp_seed\": \""
        << seed_hex() << "\",\n"
        << "  \"replicates_per_scenario\": " << R << ",\n"
        << "  \"replay_replicates\": " << rc::kReplayReplicates << ",\n"
        << "  \"scenarios\": \"id = ((i_pd*3 + i_rho)*3 + i_T)*3 + i_n; PD {0.001, 0.01, 0.05}, rho {0.02, 0.12, "
           "0.24}, T {20, 40, 100}, n {100, 1000, 10000}\",\n"
        << "  \"grid\": \"PD logit [" << fmt("%.17g", g.axis[0].lo) << ", " << fmt("%.17g", g.axis[0].hi) << "] x "
        << g.axis[0].n << "; rho logit [" << fmt("%.17g", g.axis[1].lo) << ", " << fmt("%.17g", g.axis[1].hi)
        << "] x " << g.axis[1].n << "\",\n"
        << "  \"integrator\": \"parity rule (D-118): adaptive GH 128 for 0 < d < n, composite Gauss-Legendre 16 x 16 "
           "for d in {0, n}; check: the rule doubled\",\n"
        << "  \"standard_errors\": \"Hessian of the objective at the estimate, step "
        << fmt("%.17g", vcal::engine::kHessianStepFraction) << " x stencil SE (D-119)\",\n"
        << "  \"interval\": \"95% Wald, symmetric in the logit coordinate, z = " << fmt("%.16g", rc::kZ975) << "\",\n"
        << "  \"coverage_band_z\": " << fmt("%.17g", tol::TOL_RECOVERY_COVERAGE_BAND_Z) << ",\n"
        << "  \"max_flagged_fraction\": " << fmt("%.17g", tol::TOL_RECOVERY_MAX_FLAGGED_FRACTION) << ",\n"
        << "  \"known_findings\": " << (sizeof rc::kKnownFindings / sizeof rc::kKnownFindings[0]) << ",\n"
        << "  \"profile_interval\": \"profile likelihood (D-128), threshold chi2_1(0.95)/2 = "
        << fmt("%.17g", vcal::engine::kProfileThreshold95) << ", truncated at the box, never extrapolated\",\n"
        << "  \"t_interval\": \"Wald with t(T-1) quantiles 19/39/99 df, native comparison (D-131)\",\n"
        << "  \"bootstrap\": \"iid bootstrap of periods, B = " << rc::kBootstrapReplicates
        << ", percentile interval (type 7), seed = recovery seed XOR (scenario << 32 | replicate), resampling key "
           "domain (D-132, D-135, D-136)\",\n"
        << "  \"profile_findings\": \"" << rc::known_profile_finding_count(rc::ProfileDiagnosis::TruncationConservative)
        << " conservative, " << rc::known_profile_finding_count(rc::ProfileDiagnosis::SmallTUndercoverage)
        << " undercoverage (reviewed)\",\n"
        << "  \"conditional_pd\": \"S-23: q = Phi((Phi^-1(PD) + sqrt(rho) z) / sqrt(1 - rho)), z = Phi^-1(0.999) = "
        << fmt("%.17g", vcal::engine::kAdverseZ999)
        << "; profile interval along q = c (engine/conditional_pd.hpp), delta-method Wald in logit(q), bootstrap "
           "percentile of q at the same replicates; reviewed findings: "
        << rc::kKnownQProfileFindings.size() << " profile, " << rc::kKnownQWaldFindings.size() << " Wald, "
        << rc::kKnownQBootstrapFindings.size() << " bootstrap\",\n"
        << "  \"wall_seconds\": " << fmt("%.0f", seconds) << "\n"
        << "}\n";
    close_out(out, path);
}

// S-23: the 99.9% conditional PD. Verdict counts per interval, then one row per scenario.
void write_conditional_pd_md(std::ofstream& out, const std::vector<rc::Summary>& all, std::uint32_t R) {
    int pc[rc::kVerdicts] = {}, wc[rc::kVerdicts] = {}, bc[rc::kVerdicts] = {};
    std::int64_t excludes = 0;
    double residual = 0.0;
    for (const auto& s : all) {
        ++pc[static_cast<int>(s.q.profile_verdict)];
        ++wc[static_cast<int>(s.q.wald_verdict)];
        ++bc[static_cast<int>(s.q.boot_verdict)];
        excludes += s.q.profile_excludes_estimate;
        residual = std::fmax(residual, s.q.profile_residual_max);
    }
    const auto counts = [&](const int* c) {
        return std::to_string(c[0]) + " PASS, " + std::to_string(c[static_cast<int>(rc::Verdict::Conservative)]) +
               " CONSERVATIVE, " + std::to_string(c[2]) + " KNOWN FINDING, " + std::to_string(c[1]) + " DEFERRED, " +
               std::to_string(c[3]) + " UNREVIEWED";
    };
    const auto frac = [&](std::int64_t n) { return pct(static_cast<double>(n) / static_cast<double>(R)); };
    out << "## The 99.9% conditional PD (S-23)\n\n"
        << "q = Φ((Φ⁻¹(PD) + √ρ·Φ⁻¹(0.999))/√(1 − ρ)), the PD at the 0.1% adverse factor level. Study, predictions "
           "and the comparison with them: S-23 in the [study index](../../studies/) and its "
           "[pre-registration](../../studies/derived-quantity-intervals/). Same band and policy as PD "
           "and ρ.\n\n"
        << "- **Profile likelihood, over all replicates:** " << counts(pc) << ".\n"
        << "- **Delta-method Wald in logit(q), among the unflagged replicates:** " << counts(wc) << ".\n"
        << "- **Bootstrap percentile of q, over all replicates:** " << counts(bc) << ".\n"
        << "- **Acceptance checks:** profile intervals not containing q̂: " << excludes
        << "; worst end-point residual " << fmt("%.2g", residual) << " (tolerance "
        << fmt("%.0e", tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL) << ", asserted on every fit).\n\n"
        << "\"Box-limited\": an end whose inner maximiser lies on a bound of the box, or the limit of q in the box "
           "(\"truncated\"). \"Below\": intervals entirely below the true q. Width: median of the profile interval's "
           "hi / lo.\n\n"
        << "| scenario | PD | ρ | T | n | q | median q̂/q − 1 | profile cov. | box-limited (lower) | truncated | "
           "below / above | width | profile verdict | Wald cov. | Wald verdict | bootstrap cov. | bootstrap "
           "verdict |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n";
    for (const auto& s : all) {
        const auto& q = s.q;
        out << "| " << s.s.id << " | " << fmt("%g", s.s.pd) << " | " << fmt("%g", s.s.rho) << " | " << s.s.periods
            << " | " << s.s.obligors << " | " << fmt("%.4g", q.truth) << " | " << fmt("%+.3f", q.median_rel_err)
            << " | " << fmt("%.3f", q.profile_coverage) << " | " << frac(q.profile_box_limited) << " ("
            << frac(q.profile_box_limited_lo) << ") | " << frac(q.profile_truncated) << " | " << q.profile_miss_below
            << " / " << q.profile_miss_above << " | " << fmt("%.2f", std::exp(q.profile_log_width_median)) << " | "
            << rc::verdict_text(q.profile_verdict) << " | "
            << (s.unflagged > 0 ? fmt("%.3f", q.wald_coverage) : std::string("—")) << " | "
            << rc::verdict_text(q.wald_verdict) << " | " << fmt("%.3f", q.boot_coverage) << " | "
            << rc::verdict_text(q.boot_verdict) << " |\n";
    }
    out << "\n### The 99.9% conditional PD: reviewed findings\n\n"
        << "| interval | scenario | coverage | verdict | diagnosis |\n|---|---|---|---|---|\n";
    for (const auto& k : rc::kKnownQProfileFindings) {
        const auto& q = all[k.scenario].q;
        out << "| profile | " << k.scenario << " | " << fmt("%.3f", q.profile_coverage) << " | "
            << rc::verdict_text(q.profile_verdict) << " | " << rc::profile_diagnosis_text(k.diagnosis) << " |\n";
    }
    for (const auto& k : rc::kKnownQWaldFindings) {
        const auto& q = all[k.scenario].q;
        out << "| delta-method Wald | " << k.scenario << " | " << fmt("%.3f", q.wald_coverage) << " | "
            << rc::verdict_text(q.wald_verdict) << " | " << rc::diagnosis_text(k.diagnosis) << " |\n";
    }
    int boundary = 0, no_correction = 0;
    for (const auto& k : rc::kKnownQBootstrapFindings) {
        (k.diagnosis == rc::BootstrapDiagnosis::BoundaryBreakdown ? boundary : no_correction) += 1;
    }
    out << "\nBootstrap: " << rc::kKnownQBootstrapFindings.size() << " reviewed findings, all below the band: "
        << boundary << " " << rc::bootstrap_diagnosis_text(rc::BootstrapDiagnosis::BoundaryBreakdown) << "; "
        << no_correction << " " << rc::bootstrap_diagnosis_text(rc::BootstrapDiagnosis::NoBiasSkewCorrection)
        << ". Each takes the diagnosis of ρ's bootstrap finding in the same scenario (else PD's).\n\n";
}

void write_results_md(const std::string& path, const std::vector<rc::Summary>& all,
                      const std::vector<std::string>& rmse_bad, std::uint32_t R) {
    std::ofstream out = open_out(path);
    int count[rc::kVerdicts] = {};
    int pcount[rc::kVerdicts] = {};
    int bcount[rc::kVerdicts] = {};
    for (const auto& s : all) {
        for (const auto& p : s.param) {
            ++count[static_cast<int>(p.verdict)];
            ++pcount[static_cast<int>(p.profile_verdict)];
            ++bcount[static_cast<int>(p.boot_verdict)];
        }
    }
    out << "<!-- Generated by recovery_harness --write from tests/golden/recovery/summary.csv. Do not edit. -->\n"
           "# Recovery results (M1.8, M2)\n\n"
           "Definitions, design and interpretation: [recovery.md](recovery.md) (D-121, D-124). "
        << R << " replicates per scenario, " << rc::kScenarios
        << " scenarios. Bias ± its Monte Carlo SE; bias and RMSE over **all** replicates. \"Flagged\": no "
           "reliable Wald interval (grid edge, within 2 SEs of a bound, or flat); \"edge\" is the subset on the "
           "grid edge. Coverage is of the 95% Wald interval **among unflagged replicates only (conditional)**, "
           "with the Monte Carlo band for that count. \"SE ratio\" is sd(û) / rms(se_u) among the unflagged, "
           "in the logit coordinate: above 1, the Hessian SEs understate the actual spread of the estimates. "
           "The Monte Carlo band for profile coverage over all " + std::to_string(R) + " replicates is "
        << fmt("%.3f", all[0].profile_band_lo) << "–" << fmt("%.3f", all[0].profile_band_hi) << ".\n\n"
        << "## Engine correctness\n\n"
        << "- **RMSE falls with T:** "
        << (rmse_bad.empty() ? std::string("holds for every (PD, ρ, n) and both parameters.")
                             : std::to_string(rmse_bad.size()) + " violation(s):")
        << "\n";
    for (const auto& v : rmse_bad) out << "  - " << v << "\n";
    std::int64_t quad = 0, numeric = 0;
    for (const auto& s : all) {
        quad += s.quad_unconverged;
        numeric += s.numeric;
    }
    out << "- **Quadrature check:** " << quad << " of " << static_cast<std::int64_t>(R) * rc::kScenarios
        << " fits flagged. **Non-finite likelihoods:** " << numeric << ".\n\n"
        << "## Interval coverage (a property of the interval method)\n\n"
        << "**Profile likelihood (M2, D-128..D-131), coverage over all replicates:** " << pcount[0] << " PASS, "
        << pcount[static_cast<int>(rc::Verdict::Conservative)] << " CONSERVATIVE (reviewed; above the band, a "
        << "safe-side property), " << pcount[2] << " KNOWN FINDING (reviewed; below the band), " << pcount[3]
        << " UNREVIEWED.\n\n"
        << "**iid bootstrap percentile interval (M2b, B = " << rc::kBootstrapReplicates
        << ", D-135–D-137), coverage over all replicates:** " << bcount[0] << " PASS, "
        << bcount[static_cast<int>(rc::Verdict::Conservative)] << " CONSERVATIVE, " << bcount[2]
        << " KNOWN FINDING, " << bcount[3] << " UNREVIEWED. Predictions made before the run: "
        << "[bootstrap_predictions.md](bootstrap_predictions.md).\n\n"
        << "**The band.** 0.95 ± " << fmt("%.2f", tol::TOL_RECOVERY_COVERAGE_BAND_Z)
        << " × √(0.95 × 0.05 / m): a two-sided 0.1% band, so of 162 verdicts about 0.16 would fall outside "
        << "by chance alone. The out-of-band verdicts are real effects, not noise.\n\n"
        << "**Wald, coverage among unflagged replicates (conditional):** " << count[0] << " PASS, " << count[1]
        << " DEFERRED (at least " << pct(tol::TOL_RECOVERY_MAX_FLAGGED_FRACTION) << " flagged), " << count[2]
        << " KNOWN FINDING (reviewed; listed below), " << count[3] << " UNREVIEWED. The t(T − 1) Wald interval "
        << "is a `native` comparison on the same replicates, reported without a verdict.\n\n"
        << "### Wald known findings, and what the profile interval does there\n\n"
        << "| scenario | PD | ρ | T | n | parameter | flagged | Wald coverage | band | SE ratio | t(T−1) coverage | "
           "profile coverage (all) | profile verdict | diagnosis |\n"
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n";
    for (const auto& k : rc::kKnownFindings) {
        const rc::Summary& s = all[k.scenario];
        const auto& p = s.param[k.param];
        out << "| " << s.s.id << " | " << fmt("%g", s.s.pd) << " | " << fmt("%g", s.s.rho) << " | " << s.s.periods
            << " | " << s.s.obligors << " | " << kParamLabel[k.param] << " | " << pct(s.flagged_fraction) << " | "
            << fmt("%.3f", p.coverage) << " | " << fmt("%.3f", s.band_lo) << "–" << fmt("%.3f", s.band_hi) << " | "
            << fmt("%.3f", ratio(p)) << " | " << fmt("%.3f", p.t_coverage) << " | " << fmt("%.3f", p.profile_coverage)
            << " | " << rc::verdict_text(p.profile_verdict) << " | " << rc::diagnosis_text(k.diagnosis) << " |\n";
    }
    out << "\n";

    // Wald vs profile, cross-tabulated (the headline of M2a).
    const rc::Verdict order[] = {rc::Verdict::Pass, rc::Verdict::Conservative, rc::Verdict::KnownFinding,
                                 rc::Verdict::Unreviewed};
    int cross[rc::kVerdicts][rc::kVerdicts] = {};
    for (const auto& s : all) {
        for (const auto& p : s.param) ++cross[static_cast<int>(p.verdict)][static_cast<int>(p.profile_verdict)];
    }
    out << "### Wald verdict → profile verdict (scenario × parameter)\n\n| Wald \\ profile |";
    for (const auto v : order) out << " " << rc::verdict_text(v) << " |";
    out << " total |\n|---|---|---|---|---|---|\n";
    for (const auto w : {rc::Verdict::Pass, rc::Verdict::Deferred, rc::Verdict::KnownFinding, rc::Verdict::Unreviewed}) {
        int total = 0;
        out << "| " << rc::verdict_text(w) << " |";
        for (const auto v : order) {
            out << " " << cross[static_cast<int>(w)][static_cast<int>(v)] << " |";
            total += cross[static_cast<int>(w)][static_cast<int>(v)];
        }
        out << " " << total << " |\n";
    }
    int pb[rc::kVerdicts][rc::kVerdicts] = {};
    for (const auto& s : all) {
        for (const auto& p : s.param) ++pb[static_cast<int>(p.profile_verdict)][static_cast<int>(p.boot_verdict)];
    }
    out << "\n### Profile verdict → bootstrap verdict (scenario × parameter)\n\n| profile \\ bootstrap |";
    for (const auto v : order) out << " " << rc::verdict_text(v) << " |";
    out << " total |\n|---|---|---|---|---|---|\n";
    for (const auto w : order) {
        int total = 0;
        out << "| " << rc::verdict_text(w) << " |";
        for (const auto v : order) {
            out << " " << pb[static_cast<int>(w)][static_cast<int>(v)] << " |";
            total += pb[static_cast<int>(w)][static_cast<int>(v)];
        }
        out << " " << total << " |\n";
    }
    // Bootstrap coverage by group, the groups of the pre-registered predictions: A = Wald DEFERRED
    // (near a bound or near-uninformative); otherwise B, C, D by T = 20, 40, 100.
    out << "\n### Bootstrap coverage by group (groups as in [bootstrap_predictions.md](bootstrap_predictions.md))\n\n"
        << "| group | parameter | scenarios | PASS | below the band | coverage range | zero-width intervals | "
           "profile coverage range |\n|---|---|---|---|---|---|---|---|\n";
    const char* group_name[4] = {"A: near a bound / near-uninformative", "B: T = 20", "C: T = 40", "D: T = 100"};
    for (int g = 0; g < 4; ++g) {
        for (int a = 0; a < 2; ++a) {
            int n = 0, pass = 0, below = 0;
            std::int64_t degenerate = 0;
            double lo = 1.0, hi = 0.0, plo = 1.0, phi = 0.0;
            for (const auto& s : all) {
                const bool in_a = s.param[0].verdict == rc::Verdict::Deferred;
                const int group = in_a ? 0 : s.s.periods == 20 ? 1 : s.s.periods == 40 ? 2 : 3;
                if (group != g) continue;
                const auto& p = s.param[a];
                ++n;
                pass += p.boot_verdict == rc::Verdict::Pass;
                below += p.boot_coverage < s.profile_band_lo;
                degenerate += p.boot_degenerate;
                lo = std::fmin(lo, p.boot_coverage);
                hi = std::fmax(hi, p.boot_coverage);
                plo = std::fmin(plo, p.profile_coverage);
                phi = std::fmax(phi, p.profile_coverage);
            }
            out << "| " << group_name[g] << " | " << kParamLabel[a] << " | " << n << " | " << pass << " | " << below
                << " | " << fmt("%.3f", lo) << "–" << fmt("%.3f", hi) << " | " << degenerate << " | " << fmt("%.3f", plo)
                << "–" << fmt("%.3f", phi) << " |\n";
        }
    }
    out << "\nEvery out-of-band bootstrap verdict is pinned in `kKnownBootstrapFindings` "
           "(tests/recovery/recovery.hpp) with its diagnosis: boundary breakdown (group A) or no bias/skew "
           "correction (the rest).\n";
    out << "\n### Profile findings (reviewed)\n\n"
        << "| scenario | PD | ρ | T | n | parameter | profile coverage | band | truncated | verdict | Wald verdict | "
           "diagnosis |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n";
    for (const auto& k : rc::kKnownProfileFindings) {
        const rc::Summary& s = all[k.scenario];
        const auto& p = s.param[k.param];
        out << "| " << s.s.id << " | " << fmt("%g", s.s.pd) << " | " << fmt("%g", s.s.rho) << " | " << s.s.periods
            << " | " << s.s.obligors << " | " << kParamLabel[k.param] << " | " << fmt("%.3f", p.profile_coverage)
            << " | " << fmt("%.3f", s.profile_band_lo) << "–" << fmt("%.3f", s.profile_band_hi) << " | "
            << pct(static_cast<double>(p.profile_truncated) / static_cast<double>(s.replicates)) << " | "
            << rc::verdict_text(p.profile_verdict) << " | " << rc::verdict_text(p.verdict) << " | "
            << rc::profile_diagnosis_text(k.diagnosis) << " |\n";
    }
    double residual = 0.0;
    for (const auto& s : all) residual = std::fmax(residual, s.profile_residual_max);
    out << "\nWorst profile endpoint residual over all " << static_cast<std::int64_t>(R) * rc::kScenarios
        << " fits: " << fmt("%.2g", residual) << " in log-likelihood (tolerance "
        << fmt("%.0e", tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL) << ", asserted on every fit).\n\n";
    write_conditional_pd_md(out, all, R);
    for (std::uint32_t i_pd = 0; i_pd < 3; ++i_pd) {
        out << "## PD = " << fmt("%g", rc::kPds[i_pd] * 100) << "%\n\n"
            << "| ρ | T | n | PD bias | PD RMSE | ρ bias | ρ RMSE | flagged (edge) | Wald cov. PD / ρ | "
               "t(T−1) cov. PD / ρ | Wald verdict PD / ρ | profile cov. PD / ρ (all) | profile truncated PD / ρ | "
               "profile verdict PD / ρ | bootstrap cov. PD / ρ (all) | bootstrap verdict PD / ρ |\n"
               "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n";
        for (std::uint32_t id = i_pd * 27; id < (i_pd + 1) * 27; ++id) {
            const rc::Summary& s = all[id];
            const auto& p = s.param[0];
            const auto& q = s.param[1];
            const bool any = s.unflagged > 0;
            out << "| " << fmt("%g", s.s.rho) << " | " << s.s.periods << " | " << s.s.obligors << " | "
                << fmt("%.2e", p.bias) << " ± " << fmt("%.1e", p.bias_mcse) << " | " << fmt("%.2e", p.rmse) << " | "
                << fmt("%+.4f", q.bias) << " ± " << fmt("%.4f", q.bias_mcse) << " | " << fmt("%.4f", q.rmse) << " | "
                << pct(s.flagged_fraction) << " (" << pct(static_cast<double>(s.edge) / static_cast<double>(s.replicates))
                << ") | " << (any ? fmt("%.3f", p.coverage) + " / " + fmt("%.3f", q.coverage) : "—") << " | "
                << (any ? fmt("%.3f", p.t_coverage) + " / " + fmt("%.3f", q.t_coverage) : "—") << " | "
                << rc::verdict_text(p.verdict) << " / " << rc::verdict_text(q.verdict) << " | "
                << fmt("%.3f", p.profile_coverage) << " / " << fmt("%.3f", q.profile_coverage) << " | "
                << pct(static_cast<double>(p.profile_truncated) / static_cast<double>(s.replicates)) << " / "
                << pct(static_cast<double>(q.profile_truncated) / static_cast<double>(s.replicates)) << " | "
                << rc::verdict_text(p.profile_verdict) << " / " << rc::verdict_text(q.profile_verdict) << " | "
                << fmt("%.3f", p.boot_coverage) << " / " << fmt("%.3f", q.boot_coverage) << " | "
                << rc::verdict_text(p.boot_verdict) << " / " << rc::verdict_text(q.boot_verdict) << " |\n";
        }
        out << "\n";
    }
    close_out(out, path);
}

// After writing: every output must exist and be non-empty, so a silent failure cannot recur.
void verify_written(const std::vector<std::string>& paths) {
    for (const auto& p : paths) {
        std::error_code ec;
        const auto size = std::filesystem::file_size(p, ec);
        if (ec || size == 0) throw std::runtime_error("output missing or empty after writing: " + p);
    }
}

// --- checking --------------------------------------------------------------------------------------

std::string manifest_compiler() {
    std::ifstream in(vcal::test::golden_path("recovery/MANIFEST.json"), std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string text = ss.str();
    const std::string key = "\"compiler\": \"";
    const auto a = text.find(key);
    if (a == std::string::npos) throw std::runtime_error("MANIFEST.json has no compiler");
    const auto b = text.find('"', a + key.size());
    return text.substr(a + key.size(), b - a - key.size());
}

// Returns the number of mismatches, printing each.
int check_against_goldens(const std::vector<rc::Summary>& all) {
    const auto t = vcal::test::read_golden_csv("recovery/summary.csv");
    const bool same_platform = manifest_compiler() == VCAL_COMPILER;
    std::printf("checking against tests/golden/recovery/summary.csv: %s\n",
                same_platform ? "same platform, exact" : "other platform, within tolerances");
    int bad = 0;
    const auto report = [&](std::uint32_t id, const std::string& what) {
        std::printf("  scenario %u: %s\n", id, what.c_str());
        ++bad;
    };
    if (t.rows.size() != all.size()) throw std::runtime_error("summary.csv has the wrong number of rows");
    for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
        const auto& r = t.rows[id];
        const auto num = [&](const std::string& c) { return vcal::test::parse_int(r[t.column(c)]); };
        const auto dbl = [&](const std::string& c) { return vcal::test::parse_double(r[t.column(c)]); };
        const rc::Summary& s = all[id];
        const auto count = [&](const std::string& c, std::int64_t got) {
            const std::int64_t want = num(c);
            const std::int64_t allowed = same_platform ? 0 : tol::TOL_RECOVERY_CROSS_PLATFORM_COUNT;
            if (std::llabs(got - want) > allowed) report(id, c + " " + std::to_string(got) + " vs " + std::to_string(want));
            return got == want;
        };
        const auto value = [&](const std::string& c, double got, double rel) {
            const double want = dbl(c);
            if (std::isnan(want) && std::isnan(got)) return;
            const bool ok = same_platform ? to_hex(got) == to_hex(want) : std::fabs(got / want - 1.0) <= rel;
            if (!ok) report(id, c + " " + to_hex(got) + " vs " + to_hex(want));
        };
        // Every count is compared (and reported), so no short-circuiting here.
        int unequal = 0;
        unequal += count("edge", s.edge) ? 0 : 1;
        unequal += count("near_bound", s.near_bound) ? 0 : 1;
        unequal += count("flat", s.flat) ? 0 : 1;
        unequal += count("flagged", s.flagged) ? 0 : 1;
        unequal += count("quad_unconverged", s.quad_unconverged) ? 0 : 1;
        unequal += count("numeric", s.numeric) ? 0 : 1;
        count("rejected", s.rejected);
        const bool counts_equal = unequal == 0;
        for (int a = 0; a < 2; ++a) {
            const std::string p = kParamName[a];
            const auto& q = s.param[a];
            const bool covered_equal = count(p + "_covered", q.covered);
            count(p + "_t_covered", q.t_covered);
            const bool profile_covered_equal = count(p + "_profile_covered", q.profile_covered);
            const bool profile_truncated_equal = count(p + "_profile_truncated", q.profile_truncated);
            const bool profile_equal = profile_covered_equal && profile_truncated_equal;
            const bool boot_equal = count(p + "_boot_covered", q.boot_covered);
            if ((same_platform || boot_equal) && r[t.column(p + "_boot_verdict")] != rc::verdict_text(q.boot_verdict)) {
                report(id, p + "_boot_verdict " + rc::verdict_text(q.boot_verdict) + " vs " + r[t.column(p + "_boot_verdict")]);
            }
            if ((same_platform || profile_equal) &&
                r[t.column(p + "_profile_verdict")] != rc::verdict_text(q.profile_verdict)) {
                report(id, p + "_profile_verdict " + rc::verdict_text(q.profile_verdict) + " vs " +
                               r[t.column(p + "_profile_verdict")]);
            }
            value(p + "_mean_hex", q.mean, tol::TOL_RECOVERY_REPLAY_REL);
            value(p + "_rmse_hex", q.rmse, tol::TOL_RECOVERY_REPLAY_REL);
            value(p + "_sd_u_hex", q.sd_u, tol::TOL_RECOVERY_REPLAY_REL);
            value(p + "_rms_se_u_hex", q.rms_se_u, tol::TOL_RECOVERY_REPLAY_SE_REL);
            // A verdict rests on the flagged and covered counts; compare it wherever they agree.
            if ((same_platform || (counts_equal && covered_equal)) &&
                r[t.column(p + "_verdict")] != rc::verdict_text(q.verdict)) {
                report(id, p + "_verdict " + rc::verdict_text(q.verdict) + " vs " + r[t.column(p + "_verdict")]);
            }
        }
        // S-23: q.
        const auto& q = s.q;
        const bool qp_equal = count("q_profile_covered", q.profile_covered);
        count("q_profile_truncated", q.profile_truncated);
        count("q_profile_box_limited", q.profile_box_limited);
        count("q_profile_not_computed", q.profile_not_computed);
        count("q_profile_excludes_estimate", q.profile_excludes_estimate);
        const bool qw_equal = count("q_wald_covered", q.wald_covered);
        const bool qb_equal = count("q_boot_covered", q.boot_covered);
        value("q_median_rel_err_hex", q.median_rel_err, tol::TOL_RECOVERY_REPLAY_REL);
        value("q_rmse_hex", q.rmse, tol::TOL_RECOVERY_REPLAY_REL);
        const auto q_verdict = [&](const std::string& c, rc::Verdict v, bool equal) {
            if ((same_platform || equal) && r[t.column(c)] != rc::verdict_text(v)) {
                report(id, c + " " + rc::verdict_text(v) + " vs " + r[t.column(c)]);
            }
        };
        q_verdict("q_profile_verdict", q.profile_verdict, qp_equal);
        q_verdict("q_wald_verdict", q.wald_verdict, counts_equal && qw_equal);
        q_verdict("q_boot_verdict", q.boot_verdict, qb_equal);
    }
    std::printf("%d mismatch(es)\n", bad);
    return bad;
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t R = rc::kReplicates;
    bool write = false, check = false;
    std::string save_fits, load_fits;
    std::vector<std::uint32_t> ids;
    std::string replay_dir, summary_out, results_md;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--write") == 0) {
            write = true;
        } else if (std::strcmp(argv[i], "--check") == 0) {
            check = true;
        } else if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--save-fits") == 0 && i + 1 < argc) {
            save_fits = argv[++i];
        } else if (std::strcmp(argv[i], "--load-fits") == 0 && i + 1 < argc) {
            load_fits = argv[++i];
        } else if (std::strcmp(argv[i], "--results-md") == 0 && i + 1 < argc) {
            results_md = argv[++i];
        } else if (std::strcmp(argv[i], "--summary-out") == 0 && i + 1 < argc) {
            summary_out = argv[++i];
        } else if (std::strcmp(argv[i], "--replay-dir") == 0 && i + 1 < argc) {
            replay_dir = argv[++i];
        } else if (std::strcmp(argv[i], "--scenarios") == 0 && i + 1 < argc) {
            for (const char* c = argv[++i]; *c != '\0';) {
                char* end = nullptr;
                const unsigned long id = std::strtoul(c, &end, 10);
                if (end == c || id >= rc::kScenarios) {
                    std::fprintf(stderr, "--scenarios: ids 0..%u, comma-separated\n", rc::kScenarios - 1);
                    return 2;
                }
                ids.push_back(static_cast<std::uint32_t>(id));
                c = *end == ',' ? end + 1 : end;
            }
        } else {
            std::fprintf(stderr, "usage: recovery_harness [--replicates R] [--write | --check] [--scenarios ID,...] "
                                 "[--save-fits FILE | --load-fits FILE] [--summary-out FILE] [--results-md FILE] "
                                 "[--replay-dir DIR]\n");
            return 2;
        }
    }
    if ((write || check) && (R != rc::kReplicates || (write && check))) {
        std::fprintf(stderr, "--write and --check are exclusive and require R = %u\n", rc::kReplicates);
        return 2;
    }
    if (!load_fits.empty() && (write || check)) {
        std::fprintf(stderr, "--load-fits cannot be combined with --write or --check: goldens come from fitting\n");
        return 2;
    }
    const bool subset = !ids.empty();
    if (subset && (write || check || !replay_dir.empty() || !results_md.empty())) {
        std::fprintf(stderr, "--scenarios cannot be combined with --write, --check, --replay-dir or --results-md: "
                             "those cover every scenario\n");
        return 2;
    }
    if (!subset) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    if (R < rc::kReplayReplicates) {
        std::fprintf(stderr, "need R >= %u\n", rc::kReplayReplicates);
        return 2;
    }

    // fits[r * kScenarios + id]; jobs run in batches of replicates so progress can be reported.
    std::vector<rc::Fit> fits(static_cast<std::size_t>(R) * rc::kScenarios);
    const auto start = std::chrono::steady_clock::now();
    const std::uint32_t batch = 25;
    const std::size_t bytes = fits.size() * sizeof(rc::Fit);
    if (!load_fits.empty()) {
        std::ifstream in(load_fits, std::ios::binary);
        in.read(reinterpret_cast<char*>(fits.data()), static_cast<std::streamsize>(bytes));
        if (!in || in.peek() != std::char_traits<char>::eof()) {
            std::fprintf(stderr, "cannot read %zu bytes of fits from %s\n", bytes, load_fits.c_str());
            return 1;
        }
        R = 0;  // skip the fitting loop
    }
    for (std::uint32_t r0 = 0; r0 < R; r0 += batch) {
        const std::uint32_t r1 = r0 + batch < R ? r0 + batch : R;
        const auto per = static_cast<std::int64_t>(ids.size());
        const std::int64_t jobs = static_cast<std::int64_t>(r1 - r0) * per;
        vcal::backends::CpuBackend{}.parallel_for(jobs, [&](std::int64_t j) {
            const std::uint32_t r = r0 + static_cast<std::uint32_t>(j / per);
            const std::uint32_t id = ids[static_cast<std::size_t>(j % per)];
            fits[static_cast<std::size_t>(r) * rc::kScenarios + id] = rc::fit(rc::scenario(id), r);
        });
        const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::fprintf(stderr, "replicates %u/%u, %.0f s\n", r1, R, s);
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const auto selected = [&](std::size_t i) {
        return std::find(ids.begin(), ids.end(), static_cast<std::uint32_t>(i % rc::kScenarios)) != ids.end();
    };
    // Engine correctness (D-131): every fit's profile endpoints solve the threshold to tolerance.
    for (std::size_t i = 0; i < fits.size(); ++i) {
        if (!selected(i)) continue;
        if (!(fits[i].prof_residual <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL)) {
            std::fprintf(stderr, "profile endpoint residual %.3g exceeds %.3g (scenario %zu, replicate %zu)\n",
                         fits[i].prof_residual, tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL, i % rc::kScenarios,
                         i / rc::kScenarios);
            return 1;
        }
    }
    R = static_cast<std::uint32_t>(fits.size() / rc::kScenarios);
    if (!save_fits.empty()) {
        std::ofstream out(save_fits, std::ios::binary);
        out.write(reinterpret_cast<const char*>(fits.data()), static_cast<std::streamsize>(bytes));
        if (!out) {
            std::fprintf(stderr, "cannot write %s\n", save_fits.c_str());
            return 1;
        }
    }

    // S-23 acceptance conditions, checked after saving so the fits survive for diagnosis: every
    // solved end point of q's interval meets the residual tolerance, and every interval contains q_hat.
    int q_bad = 0;
    for (std::size_t i = 0; i < fits.size(); ++i) {
        if (!selected(i)) continue;
        const rc::Fit& f = fits[i];
        const bool computed = !(f.q_prof_flags & vcal::engine::kIntervalNotComputed);
        if (!(f.q_prof_residual <= tol::TOL_PROFILE_ENDPOINT_RESIDUAL_LL) ||
            (computed && !(f.q_prof_lo <= f.q_hat && f.q_hat <= f.q_prof_hi))) {
            if (++q_bad <= 20) {
                std::fprintf(stderr, "q interval [%.17g, %.17g] flags %u, q_hat %.17g, residual %.3g (scenario %zu, "
                                     "replicate %zu)\n", f.q_prof_lo, f.q_prof_hi, f.q_prof_flags, f.q_hat,
                             f.q_prof_residual, i % rc::kScenarios, i / rc::kScenarios);
            }
        }
    }

    std::vector<rc::Summary> all(rc::kScenarios);
    std::vector<rc::Fit> scenario_fits(R);
    for (const std::uint32_t id : ids) {
        for (std::uint32_t r = 0; r < R; ++r) scenario_fits[r] = fits[static_cast<std::size_t>(r) * rc::kScenarios + id];
        all[id] = rc::summarise(rc::scenario(id), scenario_fits.data(), R);
    }
    const auto rmse_bad = subset ? std::vector<std::string>{} : rc::rmse_violations(all);

    for (const std::uint32_t id : ids) {
        const rc::Summary& s = all[id];
        std::printf("%2u PD %-5g rho %-4g T %-3lld n %-5lld | PD bias %+.2e rmse %.2e | rho bias %+.4f rmse %.4f | "
                    "flagged %5.1f%% edge %5.1f%% | cov %.3f %.3f band [%.3f, %.3f] | ratio %.3f %.3f | %s / %s | "
                    "profile %.3f %.3f %s / %s | t %.3f %.3f | boot %.3f %.3f %s / %s | quad %lld\n",
                    s.s.id, s.s.pd, s.s.rho, static_cast<long long>(s.s.periods), static_cast<long long>(s.s.obligors),
                    s.param[0].bias, s.param[0].rmse, s.param[1].bias, s.param[1].rmse, 100.0 * s.flagged_fraction,
                    100.0 * static_cast<double>(s.edge) / static_cast<double>(R), s.param[0].coverage,
                    s.param[1].coverage, s.band_lo, s.band_hi, ratio(s.param[0]), ratio(s.param[1]),
                    rc::verdict_text(s.param[0].verdict), rc::verdict_text(s.param[1].verdict),
                    s.param[0].profile_coverage, s.param[1].profile_coverage,
                    rc::verdict_text(s.param[0].profile_verdict), rc::verdict_text(s.param[1].profile_verdict),
                    s.param[0].t_coverage, s.param[1].t_coverage, s.param[0].boot_coverage, s.param[1].boot_coverage,
                    rc::verdict_text(s.param[0].boot_verdict), rc::verdict_text(s.param[1].boot_verdict),
                    static_cast<long long>(s.quad_unconverged));
        const auto& q = s.q;
        std::printf("   q %.5f | median rel err %+.4f | profile %.3f %s box-limited %.1f%% (lo %.1f%%) truncated %.1f%% "
                    "miss below/above %lld/%lld width ratio %.3f | wald %.3f %s | boot %.3f %s\n",
                    q.truth, q.median_rel_err, q.profile_coverage, rc::verdict_text(q.profile_verdict),
                    100.0 * static_cast<double>(q.profile_box_limited) / static_cast<double>(R),
                    100.0 * static_cast<double>(q.profile_box_limited_lo) / static_cast<double>(R),
                    100.0 * static_cast<double>(q.profile_truncated) / static_cast<double>(R),
                    static_cast<long long>(q.profile_miss_below), static_cast<long long>(q.profile_miss_above),
                    std::exp(q.profile_log_width_median), q.wald_coverage, rc::verdict_text(q.wald_verdict),
                    q.boot_coverage, rc::verdict_text(q.boot_verdict));
    }
    for (const auto& v : rmse_bad) std::printf("RMSE: %s\n", v.c_str());
    std::printf("%u replicates x %zu scenarios in %.0f s\n", R, ids.size(), seconds);
    if (q_bad > 0) {
        std::fprintf(stderr, "%d fits fail the S-23 acceptance conditions for q's profile interval\n", q_bad);
        return 1;
    }

    try {
        if (write) {
            const std::vector<std::string> paths = {
                vcal::test::golden_path("recovery/summary.csv"), vcal::test::golden_path("recovery/replay.csv"),
                vcal::test::golden_path("recovery/MANIFEST.json"),
                std::string(VCAL_SOURCE_DIR) + "/docs/methodology/recovery_results.md",
                vcal::test::golden_path("recovery/replay_panels.csv")};
            write_summary_csv(paths[0], all);
            write_replay_csv(paths[1], fits);
            write_manifest(paths[2], seconds, R);
            write_results_md(paths[3], all, rmse_bad, R);
            write_replay_panels_csv(paths[4]);
            verify_written(paths);
            for (const auto& p : paths) std::printf("wrote %s\n", p.c_str());
        }
        if (!results_md.empty()) {
            write_results_md(results_md, all, rmse_bad, R);
            verify_written({results_md});
            std::printf("wrote %s\n", results_md.c_str());
        }
        if (!summary_out.empty()) {
            std::vector<rc::Summary> chosen;
            for (const std::uint32_t id : ids) chosen.push_back(all[id]);
            write_summary_csv(summary_out, chosen);
            verify_written({summary_out});
            std::printf("wrote %s\n", summary_out.c_str());
        }
        if (!replay_dir.empty()) {
            const std::vector<std::string> paths = {replay_dir + "/replay.csv", replay_dir + "/replay_panels.csv"};
            write_replay_csv(paths[0], fits);
            write_replay_panels_csv(paths[1]);
            verify_written(paths);
            for (const auto& p : paths) std::printf("wrote %s\n", p.c_str());
        }
        if (check && check_against_goldens(all) != 0) return 1;
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "recovery_harness: %s\n", ex.what());
        return 1;
    }
    return 0;
}
