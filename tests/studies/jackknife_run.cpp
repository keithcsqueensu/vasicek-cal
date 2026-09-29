// SPDX-License-Identifier: Apache-2.0
//
// The shared jackknife run: studies S-3, S-5 and S-21, with the bias-corrected Wald arm for q
// (studies/jackknife-bias-rho/PREDICTION.md, J1-J18; D-150, D-155).
//
// Every recovery panel is fitted exactly as recovery_harness fits it (tests/recovery/recovery.hpp:
// the same grid, parity rule, bootstrap seed and B), and then:
//   - the T delete-one estimates from W x L (jackknife_weights, replicate_estimates);
//   - S-3: rho-tilde = T rho-hat - (T - 1) mean rho-hat_(-t), natural scale, set to a bound of the
//     box when outside it; the parity profile interval for rho shifted in logit rho by
//     logit rho-tilde - logit rho-hat, an end truncated at the box staying there;
//   - S-5: BCa intervals (resample/jackknife.hpp) for PD and rho in the logit coordinate;
//   - S-21: leave-one- and leave-two-period-out influence on rho-hat and PD-hat in Hessian SE units,
//     and whether the most influential period is the one with the most extreme realised factor;
//   - J15-J18: the Wald interval for logit(q) centred at logit q(PD-hat, rho-tilde), with the
//     delta-method SE (a) or the jackknife SE (b);
//   - J14: for replicates 0-4, each delete-one panel is refitted and polished off the grid (the
//     profile code's maximum); the refined estimate's distance from it is reported in SE units.
//
// Output (--summary-out): one row per scenario of counts and statistics, with each verdict family's
// raw position against the band (PASS, BELOW, ABOVE or DEFERRED). Reviewed labels are not written
// here: they live in studies/jackknife-bias-rho/reviewed.csv, so a review never needs a refit
// (D-155). Results do not depend on the thread count.
//
//   study_jackknife_run [--replicates R] [--scenarios ID,...] [--summary-out FILE] [--replicates-out FILE]
//                       [--polished]
//
// --polished adds the arm of the addendum of 2026-09-28: S-3's corrected rho and shifted interval,
// and the q Wald sub-arms, from exact off-grid maxima (one profile per delete-one panel; costly).
//
// --replicates-out writes one CSV row per (scenario, replicate) with every JkFit field, doubles
// to 17 significant digits (exact round trip), for conversion to a committed Parquet file.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "resample/jackknife.hpp"
#include "tests/harness/golden.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
using vcal::test::to_hex;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr std::uint32_t kExactRefitReplicates = 5;

double logit(double v) { return std::log(v) - std::log1p(-v); }
double inv_logit(double u) { return grid::from_scaled(AxisScale::Logit, u); }

// One replicate's results. Coverage flags are 1 (covers), 0 (misses, interval entirely above or
// below: see the `below` flags), or -1 (not computed / not assessed).
struct JkFit {
    double pd, rho, rho_tilde, q_true_err_a;  // q_true_err_a: logit q(PD-hat, rho-tilde) - logit q
    std::uint32_t flags;                      // the parity fit's flags
    bool clamped;
    std::int8_t parity_rho_cover, shifted_rho_cover;
    std::int8_t bca_cover[2], bca_below[2], pct_cover[2];
    std::int8_t qa_cover, qa_below, qb_cover, qb_below, qb_cover_all;
    double se_delta, se_jack;       // SEs of logit q (a, b)
    double infl_rho, infl_pd;       // largest |delete-one change| in SE units (NaN when flagged)
    double infl_pair_rho;           // largest |delete-two change| of rho in SE units
    std::int8_t extreme_z_is_top;   // most influential period for rho == most extreme |Z_t|
    double refit_gap;               // J14: max_t |refined - polished| of rho in SE units (NaN if unchecked)
    double refit_gap_clean;         // the same over delete-one fits without flags (refinement accepted)
    // The polished arm (addendum, 2026-09-28; --polished only): the same S-3 and q quantities from
    // the exact off-grid maxima of the full and delete-one panels instead of the grid refinement.
    bool polished;
    double rho_tilde_p;
    std::int8_t shifted_p_cover, qa_p_cover, qa_p_below, qb_p_cover, qb_p_below;
};

double clamp(double v, double lo, double hi) { return v < lo ? lo : v > hi ? hi : v; }

engine::SurfaceRowCache g_cache;  // shared across the run's replicates and scenarios (D-167)

JkFit jk_fit(const rc::Scenario& s, std::uint32_t replicate, bool polished) {
    using Objective = rc::Objective;
    const std::int64_t T = s.periods;
    std::vector<std::int64_t> n(static_cast<std::size_t>(T), s.obligors), d(static_cast<std::size_t>(T));
    std::vector<double> z(static_cast<std::size_t>(T));
    const dgp::PanelSpec spec{rc::kSeed, s.id, replicate, s.pd, s.rho, n.data(), T};
    if (dgp::simulate_panel(spec, d.data(), z.data()) != dgp::Status::Ok) std::abort();
    std::vector<Objective::Obs> obs(static_cast<std::size_t>(T));
    for (std::size_t t = 0; t < obs.size(); ++t) obs[t] = {s.obligors, d[t]};
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    std::vector<double> L;
    engine::Estimate2 est{};
    if (engine::calibrate_cached(backends::CpuBackend{1}, g_cache, Objective{}, primary, check, obs.data(), T, g, L,
                                 est) != engine::Status::Ok) {
        std::abort();
    }
    const auto prof = engine::profile_intervals(Objective{}, primary, obs.data(), T, g, L, est);
    const backends::CpuBackend serial{1};

    // Bootstrap replicates, exactly as recovery.hpp's fit.
    const auto idx = resample::bootstrap_indices(rc::bootstrap_seed(s.id, replicate), resample::Scheme::IidBootstrap,
                                                 rc::kBootstrapReplicates, T);
    const auto W = resample::weights_from_indices(idx.data(), rc::kBootstrapReplicates, T, T);
    std::vector<resample::Replicate2> boot(rc::kBootstrapReplicates);
    resample::replicate_estimates_compact(serial, g, obs.data(), L.data(), T, W.data(), rc::kBootstrapReplicates, boot.data());
    // Delete-one and delete-two estimates.
    const auto J1 = resample::jackknife_weights(T);
    std::vector<resample::Replicate2> jk(static_cast<std::size_t>(T));
    resample::replicate_estimates_compact(serial, g, obs.data(), L.data(), T, J1.data(), T, jk.data());
    const auto J2 = resample::delete_two_weights(T);
    const std::int64_t pairs = T * (T - 1) / 2;
    std::vector<resample::Replicate2> jk2(static_cast<std::size_t>(pairs));
    resample::replicate_estimates_compact(serial, g, obs.data(), L.data(), T, J2.data(), pairs, jk2.data());

    JkFit out{};
    out.pd = est.value[0];
    out.rho = est.value[1];
    out.flags = est.flags;
    const bool unflagged = !(est.flags & rc::kNoReliableInterval);
    const double truth[2] = {s.pd, s.rho};

    // S-3.
    std::vector<double> rho_minus(static_cast<std::size_t>(T)), pd_minus(static_cast<std::size_t>(T));
    for (std::int64_t t = 0; t < T; ++t) {
        rho_minus[static_cast<std::size_t>(t)] = jk[static_cast<std::size_t>(t)].value[1];
        pd_minus[static_cast<std::size_t>(t)] = jk[static_cast<std::size_t>(t)].value[0];
    }
    const double raw = resample::jackknife_bias_corrected(est.value[1], rho_minus.data(), T);
    out.rho_tilde = clamp(raw, g.axis[1].lo, g.axis[1].hi);
    out.clamped = out.rho_tilde != raw;
    const bool prof_ok = !(prof.flags[1] & engine::kIntervalNotComputed);
    out.parity_rho_cover = static_cast<std::int8_t>(prof_ok ? (prof.lo[1] <= s.rho && s.rho <= prof.hi[1]) : -1);
    if (prof_ok) {
        const double shift = logit(out.rho_tilde) - logit(est.value[1]);
        const double lo = (prof.flags[1] & engine::kIntervalLowerTruncated) ? prof.lo[1] : inv_logit(logit(prof.lo[1]) + shift);
        const double hi = (prof.flags[1] & engine::kIntervalUpperTruncated) ? prof.hi[1] : inv_logit(logit(prof.hi[1]) + shift);
        out.shifted_rho_cover = lo <= s.rho && s.rho <= hi;
    } else {
        out.shifted_rho_cover = -1;
    }

    // S-5: BCa in logit coordinates, and the percentile interval from the same replicates.
    for (int a = 0; a < 2; ++a) {
        std::vector<double> u_boot(boot.size()), u_jk(static_cast<std::size_t>(T));
        for (std::size_t b = 0; b < boot.size(); ++b) u_boot[b] = logit(boot[b].value[a]);
        for (std::int64_t t = 0; t < T; ++t) u_jk[static_cast<std::size_t>(t)] = logit(jk[static_cast<std::size_t>(t)].value[a]);
        const auto bca = resample::bca_interval(u_boot.data(), static_cast<std::int64_t>(u_boot.size()), logit(est.value[a]),
                                                u_jk.data(), T);
        if (bca.computed) {
            const double lo = inv_logit(bca.lo), hi = inv_logit(bca.hi);
            out.bca_cover[a] = lo <= truth[a] && truth[a] <= hi;
            out.bca_below[a] = hi < truth[a];
        } else {
            out.bca_cover[a] = -1;
        }
        const auto pct = resample::percentile_interval(boot.data(), static_cast<std::int64_t>(boot.size()), a);
        out.pct_cover[a] = pct.lo <= truth[a] && truth[a] <= pct.hi;
    }

    // J15-J18: Wald for logit(q) centred at logit q(PD-hat, rho-tilde).
    const double s_true = logit(engine::conditional_pd(s.pd, s.rho));
    const double s_tilde = logit(engine::conditional_pd(est.value[0], out.rho_tilde));
    out.q_true_err_a = s_tilde - s_true;
    {
        const auto qg = engine::conditional_pd_logit_gradient(est.value[0], est.value[1]);
        double se_u[2];
        for (int a = 0; a < 2; ++a) se_u[a] = est.se[a] / grid::dvalue_dscaled(AxisScale::Logit, logit(est.value[a]));
        const double var = qg.ds_du[0] * qg.ds_du[0] * se_u[0] * se_u[0] + qg.ds_du[1] * qg.ds_du[1] * se_u[1] * se_u[1] +
                           2.0 * qg.ds_du[0] * qg.ds_du[1] * est.corr * se_u[0] * se_u[1];
        out.se_delta = var >= 0.0 ? std::sqrt(var) : kNaN;
        double mean = 0.0;
        std::vector<double> sm(static_cast<std::size_t>(T));
        for (std::int64_t t = 0; t < T; ++t) {
            sm[static_cast<std::size_t>(t)] = logit(engine::conditional_pd(pd_minus[static_cast<std::size_t>(t)],
                                                                           rho_minus[static_cast<std::size_t>(t)]));
            mean += sm[static_cast<std::size_t>(t)];
        }
        mean /= static_cast<double>(T);
        double ss = 0.0;
        for (const double v : sm) ss += (v - mean) * (v - mean);
        out.se_jack = std::sqrt(static_cast<double>(T - 1) / static_cast<double>(T) * ss);
        const auto covers = [&](double se) { return std::fabs(out.q_true_err_a) <= rc::kZ975 * se; };
        out.qb_cover_all = static_cast<std::int8_t>(std::isfinite(out.se_jack) ? covers(out.se_jack) : -1);
        if (unflagged) {
            out.qa_cover = covers(out.se_delta);
            out.qa_below = !out.qa_cover && out.q_true_err_a < 0.0;
            out.qb_cover = covers(out.se_jack);
            out.qb_below = !out.qb_cover && out.q_true_err_a < 0.0;
        } else {
            out.qa_cover = out.qb_cover = -1;
        }
    }

    // S-21: influence in SE units, on replicates without a Wald flag.
    out.infl_rho = out.infl_pd = out.infl_pair_rho = kNaN;
    out.extreme_z_is_top = -1;
    if (unflagged) {
        std::int64_t top = 0, extreme = 0;
        double best = -1.0, best_pd = 0.0;
        for (std::int64_t t = 0; t < T; ++t) {
            const double dr = std::fabs(rho_minus[static_cast<std::size_t>(t)] - est.value[1]) / est.se[1];
            const double dp = std::fabs(pd_minus[static_cast<std::size_t>(t)] - est.value[0]) / est.se[0];
            if (dr > best) {
                best = dr;
                top = t;
            }
            best_pd = std::fmax(best_pd, dp);
            if (std::fabs(z[static_cast<std::size_t>(t)]) > std::fabs(z[static_cast<std::size_t>(extreme)])) extreme = t;
        }
        out.infl_rho = best;
        out.infl_pd = best_pd;
        out.extreme_z_is_top = top == extreme;
        double best2 = 0.0;
        for (const auto& r : jk2) best2 = std::fmax(best2, std::fabs(r.value[1] - est.value[1]) / est.se[1]);
        out.infl_pair_rho = best2;
    }

    // J14: exact refits of the delete-one panels for the first replicates.
    out.refit_gap = out.refit_gap_clean = kNaN;
    if (replicate < kExactRefitReplicates && unflagged) {
        double gap = 0.0, gap_clean = 0.0;
        std::vector<Objective::Obs> sub;
        std::vector<double> Ls;
        for (std::int64_t t = 0; t < T; ++t) {
            sub.clear();
            for (std::int64_t u = 0; u < T; ++u) {
                if (u != t) sub.push_back(obs[static_cast<std::size_t>(u)]);
            }
            engine::Estimate2 e{};
            if (engine::calibrate(serial, Objective{}, primary, check, sub.data(), T - 1, g, Ls, e) != engine::Status::Ok) {
                std::abort();
            }
            const auto p = engine::profile_intervals(Objective{}, primary, sub.data(), T - 1, g, Ls, e);
            if (!std::isfinite(p.max_at[1])) continue;
            const double gt = std::fabs(e.value[1] - p.max_at[1]) / est.se[1];
            gap = std::fmax(gap, gt);
            if (e.flags == 0) gap_clean = std::fmax(gap_clean, gt);
        }
        out.refit_gap = gap;
        out.refit_gap_clean = gap_clean;
    }
    // The polished arm: each delete-one panel's exact maximum, from the profile code on that panel
    // (its per-period surfaces are L without row t), and the full panel's from prof.max_at.
    out.polished = polished;
    if (polished) {
        out.shifted_p_cover = out.qa_p_cover = out.qb_p_cover = -1;
        const std::int64_t K = g.size();
        std::vector<double> rho_p(static_cast<std::size_t>(T)), pd_p(static_cast<std::size_t>(T)), Ls;
        std::vector<Objective::Obs> sub;
        bool ok = std::isfinite(prof.max_at[0]) && std::isfinite(prof.max_at[1]);
        for (std::int64_t t = 0; t < T && ok; ++t) {
            sub.clear();
            Ls.clear();
            for (std::int64_t u = 0; u < T; ++u) {
                if (u == t) continue;
                sub.push_back(obs[static_cast<std::size_t>(u)]);
                Ls.insert(Ls.end(), L.begin() + u * K, L.begin() + (u + 1) * K);
            }
            const auto& e = jk[static_cast<std::size_t>(t)];
            engine::Estimate2 et{};
            et.value[0] = e.value[0];
            et.value[1] = e.value[1];
            et.flags = e.flags & (engine::kFlagFlatSurface | engine::kFlagNumeric);
            et.loglik = -std::numeric_limits<double>::infinity();  // the polish finds the maximum itself
            const auto pt = engine::profile_intervals(Objective{}, primary, sub.data(), T - 1, g, Ls, et);
            ok = std::isfinite(pt.max_at[0]) && std::isfinite(pt.max_at[1]);
            pd_p[static_cast<std::size_t>(t)] = pt.max_at[0];
            rho_p[static_cast<std::size_t>(t)] = pt.max_at[1];
        }
        if (ok) {
            const double raw_p = resample::jackknife_bias_corrected(prof.max_at[1], rho_p.data(), T);
            out.rho_tilde_p = clamp(raw_p, g.axis[1].lo, g.axis[1].hi);
            if (prof_ok) {
                const double shift = logit(out.rho_tilde_p) - logit(prof.max_at[1]);
                const double lo = (prof.flags[1] & engine::kIntervalLowerTruncated) ? prof.lo[1] : inv_logit(logit(prof.lo[1]) + shift);
                const double hi = (prof.flags[1] & engine::kIntervalUpperTruncated) ? prof.hi[1] : inv_logit(logit(prof.hi[1]) + shift);
                out.shifted_p_cover = lo <= s.rho && s.rho <= hi;
            }
            const double err = logit(engine::conditional_pd(prof.max_at[0], out.rho_tilde_p)) - s_true;
            double mean = 0.0;
            std::vector<double> sm(static_cast<std::size_t>(T));
            for (std::int64_t t = 0; t < T; ++t) {
                sm[static_cast<std::size_t>(t)] = logit(engine::conditional_pd(pd_p[static_cast<std::size_t>(t)],
                                                                               rho_p[static_cast<std::size_t>(t)]));
                mean += sm[static_cast<std::size_t>(t)];
            }
            mean /= static_cast<double>(T);
            double ss = 0.0;
            for (const double v : sm) ss += (v - mean) * (v - mean);
            const double se_jack_p = std::sqrt(static_cast<double>(T - 1) / static_cast<double>(T) * ss);
            if (unflagged) {
                out.qa_p_cover = std::fabs(err) <= rc::kZ975 * out.se_delta;
                out.qa_p_below = !out.qa_p_cover && err < 0.0;
                out.qb_p_cover = std::fabs(err) <= rc::kZ975 * se_jack_p;
                out.qb_p_below = !out.qb_p_cover && err < 0.0;
            }
        }
    }
    return out;
}

const char* band_class(double coverage, double lo, double hi) {
    if (!(coverage == coverage)) return "DEFERRED";
    return coverage < lo ? "BELOW" : coverage > hi ? "ABOVE" : "PASS";
}

double median(std::vector<double> x) {
    x.erase(std::remove_if(x.begin(), x.end(), [](double v) { return !std::isfinite(v); }), x.end());
    if (x.empty()) return kNaN;
    std::sort(x.begin(), x.end());
    const std::size_t m = x.size() / 2;
    return x.size() % 2 ? x[m] : 0.5 * (x[m - 1] + x[m]);
}

double quantile(std::vector<double> x, double p) {
    x.erase(std::remove_if(x.begin(), x.end(), [](double v) { return !std::isfinite(v); }), x.end());
    if (x.empty()) return kNaN;
    std::sort(x.begin(), x.end());
    return resample::quantile_type7(x, p);
}

double correlation(const std::vector<double>& a, const std::vector<double>& b) {
    double ma = 0.0, mb = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        ma += a[i];
        mb += b[i];
    }
    ma /= static_cast<double>(a.size());
    mb /= static_cast<double>(b.size());
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sab += (a[i] - ma) * (b[i] - mb);
        saa += (a[i] - ma) * (a[i] - ma);
        sbb += (b[i] - mb) * (b[i] - mb);
    }
    return sab / std::sqrt(saa * sbb);
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids;
    std::string summary_out, replicates_out;
    bool polished = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--polished") == 0) {
            polished = true;
        } else if (std::strcmp(argv[i], "--replicates-out") == 0 && i + 1 < argc) {
            replicates_out = argv[++i];
        } else if (std::strcmp(argv[i], "--summary-out") == 0 && i + 1 < argc) {
            summary_out = argv[++i];
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
            std::fprintf(stderr, "usage: study_jackknife_run [--replicates R] [--scenarios ID,...] [--summary-out FILE] "
                                 "[--replicates-out FILE] [--polished]\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    const auto per = static_cast<std::int64_t>(ids.size());
    std::vector<JkFit> fits(static_cast<std::size_t>(R) * ids.size());
    const auto start = std::chrono::steady_clock::now();
    const std::uint32_t batch = 25;
    for (std::uint32_t r0 = 0; r0 < R; r0 += batch) {
        const std::uint32_t r1 = r0 + batch < R ? r0 + batch : R;
        backends::CpuBackend{}.parallel_for(static_cast<std::int64_t>(r1 - r0) * per, [&](std::int64_t j) {
            const std::uint32_t r = r0 + static_cast<std::uint32_t>(j / per);
            const std::size_t k = static_cast<std::size_t>(j % per);
            fits[static_cast<std::size_t>(r) * ids.size() + k] = jk_fit(rc::scenario(ids[k]), r, polished);
        });
        std::fprintf(stderr, "replicates %u/%u, %.0f s\n", r1, R,
                     std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
    }

    double band_lo = 0.0, band_hi = 0.0;
    rc::coverage_band(static_cast<std::int64_t>(R), band_lo, band_hi);
    std::string csv =
        "# The shared jackknife run (S-3, S-5, S-21; J1-J18): one row per scenario, R replicates each.\n"
        "# Written by study_jackknife_run. *_class is the raw position against the band (PASS, BELOW,\n"
        "# ABOVE, DEFERRED); reviewed labels live in studies/jackknife-bias-rho/reviewed.csv (D-155).\n"
        "scenario,replicates,unflagged,rho_bias_hex,rho_rmse_hex,rho_tilde_bias_hex,rho_tilde_rmse_hex,clamped,"
        "parity_rho_covered,shifted_rho_covered,shifted_gained,shifted_lost,shifted_class,"
        "pd_bca_covered,pd_bca_below,pd_bca_not_computed,pd_bca_class,pd_pct_covered,"
        "rho_bca_covered,rho_bca_below,rho_bca_not_computed,rho_bca_class,rho_pct_covered,"
        "qa_covered,qa_below,qa_class,qb_covered,qb_below,qb_class,qb_covered_all,corr_err_se_delta_hex,"
        "corr_err_se_jack_hex,infl_rho_median_hex,infl_rho_p90_hex,infl_pd_median_hex,extreme_z_top,extreme_z_assessed,"
        "pair_ratio_median_hex,refit_gap_max_hex,refit_gap_clean_max_hex,shifted_coverage,rho_bca_coverage,pd_bca_coverage,qa_coverage,"
        "qb_coverage,rho_tilde_p_bias_hex,rho_tilde_p_rmse_hex,shifted_p_covered,shifted_p_class,qa_p_covered,qa_p_below,"
        "qa_p_class,qb_p_covered,qb_p_below,qb_p_class\n";
    for (std::size_t k = 0; k < ids.size(); ++k) {
        const rc::Scenario s = rc::scenario(ids[k]);
        std::int64_t unflagged = 0, clamped = 0, par = 0, sh = 0, gained = 0, lost = 0;
        std::int64_t bca[2] = {}, bca_below[2] = {}, bca_nc[2] = {}, pct[2] = {};
        std::int64_t qa = 0, qa_below = 0, qb = 0, qb_below = 0, qb_all = 0, top = 0, assessed = 0;
        rc::CompensatedSum b_hat, b_tilde, e_hat, e_tilde;
        std::vector<double> infl_rho, infl_pd, ratio, err, sd, sj;
        double refit = 0.0, refit_clean = 0.0;
        rc::CompensatedSum b_p, e_p;
        std::int64_t sh_p = 0, qa_p = 0, qa_p_below = 0, qb_p = 0, qb_p_below = 0;
        for (std::uint32_t r = 0; r < R; ++r) {
            const JkFit& f = fits[static_cast<std::size_t>(r) * ids.size() + k];
            b_hat.add(f.rho - s.rho);
            b_tilde.add(f.rho_tilde - s.rho);
            e_hat.add((f.rho - s.rho) * (f.rho - s.rho));
            e_tilde.add((f.rho_tilde - s.rho) * (f.rho_tilde - s.rho));
            clamped += f.clamped;
            par += f.parity_rho_cover == 1;
            sh += f.shifted_rho_cover == 1;
            gained += f.shifted_rho_cover == 1 && f.parity_rho_cover != 1;
            lost += f.shifted_rho_cover != 1 && f.parity_rho_cover == 1;
            for (int a = 0; a < 2; ++a) {
                bca[a] += f.bca_cover[a] == 1;
                bca_below[a] += f.bca_cover[a] == 0 && f.bca_below[a];
                bca_nc[a] += f.bca_cover[a] == -1;
                pct[a] += f.pct_cover[a] == 1;
            }
            qb_all += f.qb_cover_all == 1;
            if (f.polished && std::isfinite(f.rho_tilde_p)) {
                b_p.add(f.rho_tilde_p - s.rho);
                e_p.add((f.rho_tilde_p - s.rho) * (f.rho_tilde_p - s.rho));
            }
            sh_p += f.polished && f.shifted_p_cover == 1;
            if (f.polished && f.qa_p_cover >= 0) {
                qa_p += f.qa_p_cover;
                qa_p_below += f.qa_p_below;
                qb_p += f.qb_p_cover;
                qb_p_below += f.qb_p_below;
            }
            if (f.qa_cover < 0) continue;
            ++unflagged;
            qa += f.qa_cover;
            qa_below += f.qa_below;
            qb += f.qb_cover;
            qb_below += f.qb_below;
            err.push_back(f.q_true_err_a);
            sd.push_back(f.se_delta);
            sj.push_back(f.se_jack);
            infl_rho.push_back(f.infl_rho);
            infl_pd.push_back(f.infl_pd);
            ratio.push_back(f.infl_pair_rho / f.infl_rho);
            top += f.extreme_z_is_top == 1;
            assessed += f.extreme_z_is_top >= 0;
            if (std::isfinite(f.refit_gap)) refit = std::fmax(refit, f.refit_gap);
            if (std::isfinite(f.refit_gap_clean)) refit_clean = std::fmax(refit_clean, f.refit_gap_clean);
        }
        const double Rd = static_cast<double>(R);
        double m_lo = 0.0, m_hi = 0.0;
        rc::coverage_band(unflagged, m_lo, m_hi);
        const double fraction = static_cast<double>(static_cast<std::int64_t>(R) - unflagged) / Rd;
        const bool deferred = !(fraction < tol::TOL_RECOVERY_MAX_FLAGGED_FRACTION) || unflagged == 0;
        const double sh_cov = static_cast<double>(sh) / Rd;
        const double bca_cov[2] = {static_cast<double>(bca[0]) / Rd, static_cast<double>(bca[1]) / Rd};
        const double qa_cov = deferred ? kNaN : static_cast<double>(qa) / static_cast<double>(unflagged);
        const double qb_cov = deferred ? kNaN : static_cast<double>(qb) / static_cast<double>(unflagged);
        char line[4096];
        std::snprintf(line, sizeof line,
                      "%u,%u,%lld,%s,%s,%s,%s,%lld,%lld,%lld,%lld,%lld,%s,%lld,%lld,%lld,%s,%lld,%lld,%lld,%lld,%s,%lld,"
                      "%lld,%lld,%s,%lld,%lld,%s,%lld,%s,%s,%s,%s,%s,%lld,%lld,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                      s.id, R, static_cast<long long>(unflagged), to_hex(b_hat.value() / Rd).c_str(),
                      to_hex(std::sqrt(e_hat.value() / Rd)).c_str(), to_hex(b_tilde.value() / Rd).c_str(),
                      to_hex(std::sqrt(e_tilde.value() / Rd)).c_str(), static_cast<long long>(clamped),
                      static_cast<long long>(par), static_cast<long long>(sh), static_cast<long long>(gained),
                      static_cast<long long>(lost), band_class(sh_cov, band_lo, band_hi), static_cast<long long>(bca[0]),
                      static_cast<long long>(bca_below[0]), static_cast<long long>(bca_nc[0]),
                      band_class(bca_cov[0], band_lo, band_hi), static_cast<long long>(pct[0]),
                      static_cast<long long>(bca[1]), static_cast<long long>(bca_below[1]),
                      static_cast<long long>(bca_nc[1]), band_class(bca_cov[1], band_lo, band_hi),
                      static_cast<long long>(pct[1]), static_cast<long long>(qa), static_cast<long long>(qa_below),
                      band_class(qa_cov, m_lo, m_hi), static_cast<long long>(qb), static_cast<long long>(qb_below),
                      band_class(qb_cov, m_lo, m_hi), static_cast<long long>(qb_all),
                      to_hex(err.size() > 2 ? correlation(err, sd) : kNaN).c_str(),
                      to_hex(err.size() > 2 ? correlation(err, sj) : kNaN).c_str(), to_hex(median(infl_rho)).c_str(),
                      to_hex(quantile(infl_rho, 0.9)).c_str(), to_hex(median(infl_pd)).c_str(),
                      static_cast<long long>(top), static_cast<long long>(assessed), to_hex(median(ratio)).c_str(),
                      to_hex(refit).c_str(), to_hex(refit_clean).c_str(), sh_cov, bca_cov[1], bca_cov[0], qa_cov, qb_cov);
        csv += line;
        csv.pop_back();  // the row continues with the polished arm's columns
        const double shp_cov = polished ? static_cast<double>(sh_p) / Rd : kNaN;
        const double qap_cov = polished && !deferred ? static_cast<double>(qa_p) / static_cast<double>(unflagged) : kNaN;
        const double qbp_cov = polished && !deferred ? static_cast<double>(qb_p) / static_cast<double>(unflagged) : kNaN;
        std::snprintf(line, sizeof line, ",%s,%s,%lld,%s,%lld,%lld,%s,%lld,%lld,%s\n",
                      to_hex(polished ? b_p.value() / Rd : kNaN).c_str(),
                      to_hex(polished ? std::sqrt(e_p.value() / Rd) : kNaN).c_str(), static_cast<long long>(sh_p),
                      polished ? band_class(shp_cov, band_lo, band_hi) : "-", static_cast<long long>(qa_p),
                      static_cast<long long>(qa_p_below), polished ? band_class(qap_cov, m_lo, m_hi) : "-",
                      static_cast<long long>(qb_p), static_cast<long long>(qb_p_below),
                      polished ? band_class(qbp_cov, m_lo, m_hi) : "-");
        csv += line;
        if (polished) {
            std::printf("   polished: rho-tilde bias %+.4f rmse %.4f | shifted %.3f %s | q Wald a %.3f %s b %.3f %s\n",
                        b_p.value() / Rd, std::sqrt(e_p.value() / Rd), shp_cov, band_class(shp_cov, band_lo, band_hi),
                        qap_cov, band_class(qap_cov, m_lo, m_hi), qbp_cov, band_class(qbp_cov, m_lo, m_hi));
        }
        std::printf("%2u T %-3lld n %-5lld | rho bias %+.4f -> %+.4f rmse %.4f -> %.4f clamp %lld | profile %.3f -> "
                    "shifted %.3f %s | BCa PD %.3f %s rho %.3f %s (pct %.3f %.3f) | q Wald a %.3f %s b %.3f %s | "
                    "infl rho %.2f pair ratio %.2f top=extreme %.2f | refit gap %.3g (unflagged %.3g)\n",
                    s.id, static_cast<long long>(s.periods), static_cast<long long>(s.obligors), b_hat.value() / Rd,
                    b_tilde.value() / Rd, std::sqrt(e_hat.value() / Rd), std::sqrt(e_tilde.value() / Rd),
                    static_cast<long long>(clamped), static_cast<double>(par) / Rd, sh_cov,
                    band_class(sh_cov, band_lo, band_hi), bca_cov[0], band_class(bca_cov[0], band_lo, band_hi),
                    bca_cov[1], band_class(bca_cov[1], band_lo, band_hi), static_cast<double>(pct[0]) / Rd,
                    static_cast<double>(pct[1]) / Rd, qa_cov, band_class(qa_cov, m_lo, m_hi), qb_cov,
                    band_class(qb_cov, m_lo, m_hi), median(infl_rho), median(ratio),
                    assessed > 0 ? static_cast<double>(top) / static_cast<double>(assessed) : kNaN, refit, refit_clean);
    }
    std::printf("%u replicates x %zu scenarios in %.0f s\n", R, ids.size(),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
    if (!replicates_out.empty()) {
        std::ofstream out(replicates_out, std::ios::binary);
        out << "scenario,replicate,pd,rho,rho_tilde,q_err_a,flags,clamped,parity_rho_cover,shifted_rho_cover,"
               "pd_bca_cover,pd_bca_below,rho_bca_cover,rho_bca_below,pd_pct_cover,rho_pct_cover,qa_cover,qa_below,"
               "qb_cover,qb_below,qb_cover_all,se_delta,se_jack,infl_rho,infl_pd,infl_pair_rho,extreme_z_is_top,"
               "refit_gap,refit_gap_clean,rho_tilde_p,shifted_p_cover,qa_p_cover,qa_p_below,qb_p_cover,qb_p_below\n";
        char line[1024];
        for (std::size_t k = 0; k < ids.size(); ++k) {
            for (std::uint32_t r = 0; r < R; ++r) {
                const JkFit& f = fits[static_cast<std::size_t>(r) * ids.size() + k];
                std::snprintf(line, sizeof line,
                              "%u,%u,%.17g,%.17g,%.17g,%.17g,%u,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g,"
                              "%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%d,%d,%d,%d,%d\n",
                              ids[k], r, f.pd, f.rho, f.rho_tilde, f.q_true_err_a, f.flags, f.clamped ? 1 : 0,
                              f.parity_rho_cover, f.shifted_rho_cover, f.bca_cover[0], f.bca_below[0], f.bca_cover[1],
                              f.bca_below[1], f.pct_cover[0], f.pct_cover[1], f.qa_cover, f.qa_below, f.qb_cover,
                              f.qb_below, f.qb_cover_all, f.se_delta, f.se_jack, f.infl_rho, f.infl_pd, f.infl_pair_rho,
                              f.extreme_z_is_top, f.refit_gap, f.refit_gap_clean, f.polished ? f.rho_tilde_p : kNaN,
                              f.polished ? f.shifted_p_cover : -1, f.polished ? f.qa_p_cover : -1,
                              f.polished ? f.qa_p_below : -1, f.polished ? f.qb_p_cover : -1,
                              f.polished ? f.qb_p_below : -1);
                out << line;
            }
        }
        if (!out) {
            std::fprintf(stderr, "cannot write %s\n", replicates_out.c_str());
            return 1;
        }
        std::printf("wrote %s\n", replicates_out.c_str());
    }
    if (!summary_out.empty()) {
        std::ofstream out(summary_out, std::ios::binary);
        out << csv;
        if (!out) {
            std::fprintf(stderr, "cannot write %s\n", summary_out.c_str());
            return 1;
        }
        std::printf("wrote %s\n", summary_out.c_str());
    }
    return 0;
}
