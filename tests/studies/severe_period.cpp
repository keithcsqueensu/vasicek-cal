// SPDX-License-Identifier: Apache-2.0
//
// S-34, sensitivity to severe new periods (studies/severe-period-sensitivity/PREDICTION.md, K1-K10;
// D-152, D-160).
//
// Every recovery panel is fitted as recovery_harness fits it (tests/recovery/recovery.hpp: the same
// grid, parity rule and check), and then refitted four more times, extended by k = 1 or 2 periods
// at a 1-in-100 or 1-in-1,000 adverse factor level. Each added period has the scenario's n obligors
// and the median default count of Binomial(n, p_a), p_a the conditional PD at that level
// (engine::conditional_pd, the engine's convention), so the added periods are the same for every
// replicate of a scenario. For the original and each extended panel:
//   - PD-hat and rho-hat are the exact off-grid maximum (the profile code's polished maximum), not
//     the grid refinement (D-159), and q-hat = q(PD-hat, rho-hat), the 99.9% conditional PD;
//   - the parity profile intervals for PD and rho, and S-23's for q.
// Shifts are in the original panel's SE units: its Hessian SEs for PD and rho (D-119), S-23's
// delta-method SE of logit q for logit q, on replicates without a Wald flag.
//
// Output (--summary-out): one row per scenario and variant (variant 0 is the original panel).
// --replicates-out: one row per scenario, replicate and variant, doubles to 17 significant digits,
// for conversion to a committed Parquet file. --counts prints each scenario's added default counts
// and exits, reading no panel. Results do not depend on the thread count.
//
//   study_severe_period [--replicates R] [--scenarios ID,...] [--summary-out FILE] [--replicates-out FILE]
//   study_severe_period --counts
#include <algorithm>
#include <array>
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

#include "tests/harness/golden.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
using vcal::test::to_hex;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
// Phi^-1(0.99), the 1-in-100 adverse factor level (mpmath: 2.32634787404084110...); the 1-in-1,000
// level is engine::kAdverseZ999. Adverse in the engine's convention: a higher Z is better.
constexpr double kAdverseZ99 = 2.3263478740408408;

struct Variant {
    int added;     // k, the number of added periods
    double z_a;    // the adverse factor level, as a positive quantile
    const char* name;
};
constexpr int kVariants = 4;
constexpr Variant kVariant[kVariants] = {{1, kAdverseZ99, "1x1-in-100"},
                                         {2, kAdverseZ99, "2x1-in-100"},
                                         {1, engine::kAdverseZ999, "1x1-in-1000"},
                                         {2, engine::kAdverseZ999, "2x1-in-1000"}};

// The median of Binomial(n, p): the smallest d with P(D <= d) >= 1/2, summed in log space.
std::int64_t binomial_median(std::int64_t n, double p) {
    const double lg_n = std::lgamma(static_cast<double>(n) + 1.0), lp = std::log(p), lq = std::log1p(-p);
    double cdf = 0.0;
    for (std::int64_t d = 0; d < n; ++d) {
        const double dd = static_cast<double>(d);
        cdf += std::exp(lg_n - std::lgamma(dd + 1.0) - std::lgamma(static_cast<double>(n - d) + 1.0) + dd * lp +
                        static_cast<double>(n - d) * lq);
        if (cdf >= 0.5) return d;
    }
    return n;
}

std::int64_t added_count(const rc::Scenario& s, const Variant& v) {
    return binomial_median(s.obligors, engine::conditional_pd(s.pd, s.rho, v.z_a));
}

// One panel's fit: estimates at the polished maximum, and the three profile intervals.
// cover: 1 covers the truth, 0 misses, -1 not computed (counted as not covering, D-131).
struct Side {
    double value[3];  // PD-hat, rho-hat, q-hat
    double lo[3], hi[3];
    std::uint32_t flags;  // the refined fit's flags (kNoReliableInterval marks a Wald flag)
    std::int8_t cover[3];
};

struct SpFit {
    Side side[1 + kVariants];  // [0] the original panel, [1..4] the variants
    double se[3];              // the original panel's SEs: PD, rho (natural scale), logit q
};

// se, when given, receives the SEs at the refined estimate as recovery.hpp's fit computes them: the
// Hessian SEs of PD and rho (D-119) and S-23's delta-method SE of logit q.
Side fit_panel(const rc::Scenario& s, const std::vector<rc::Objective::Obs>& obs, double* se = nullptr) {
    using Objective = rc::Objective;
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    const auto T = static_cast<std::int64_t>(obs.size());
    std::vector<double> L;
    engine::Estimate2 est{};
    if (engine::calibrate(backends::CpuBackend{1}, Objective{}, primary, check, obs.data(), T, g, L, est) !=
        engine::Status::Ok) {
        std::abort();
    }
    const auto prof = engine::profile_intervals(Objective{}, primary, obs.data(), T, g, L, est);
    const auto qi = engine::conditional_pd_interval(Objective{}, primary, obs.data(), T, g, L, est, prof);
    if (se != nullptr) {
        const auto qg = engine::conditional_pd_logit_gradient(est.value[0], est.value[1]);
        double se_u[2];
        for (int a = 0; a < 2; ++a) {
            se_u[a] = est.se[a] / grid::dvalue_dscaled(g.axis[a].scale, grid::to_scaled(g.axis[a].scale, est.value[a]));
        }
        const double q_var = qg.ds_du[0] * qg.ds_du[0] * se_u[0] * se_u[0] +
                             qg.ds_du[1] * qg.ds_du[1] * se_u[1] * se_u[1] +
                             2.0 * qg.ds_du[0] * qg.ds_du[1] * est.corr * se_u[0] * se_u[1];
        se[0] = est.se[0];
        se[1] = est.se[1];
        se[2] = q_var >= 0.0 ? std::sqrt(q_var) : kNaN;
    }
    Side out{};
    out.value[0] = prof.max_at[0];
    out.value[1] = prof.max_at[1];
    out.value[2] = std::isfinite(prof.max_at[0]) && std::isfinite(prof.max_at[1])
                       ? engine::conditional_pd(prof.max_at[0], prof.max_at[1])
                       : kNaN;
    out.flags = est.flags;
    const double truth[3] = {s.pd, s.rho, engine::conditional_pd(s.pd, s.rho)};
    for (int a = 0; a < 3; ++a) {
        out.lo[a] = a < 2 ? prof.lo[a] : qi.lo;
        out.hi[a] = a < 2 ? prof.hi[a] : qi.hi;
        const bool computed = !((a < 2 ? prof.flags[a] : qi.flags) & engine::kIntervalNotComputed);
        out.cover[a] = static_cast<std::int8_t>(computed ? (out.lo[a] <= truth[a] && truth[a] <= out.hi[a]) : -1);
    }
    return out;
}

double logit(double v) { return std::log(v) - std::log1p(-v); }

SpFit sp_fit(const rc::Scenario& s, std::uint32_t replicate, const std::int64_t (&counts)[kVariants]) {
    const std::vector<std::int64_t> d = rc::panel(s, replicate);
    std::vector<rc::Objective::Obs> obs(d.size());
    for (std::size_t t = 0; t < d.size(); ++t) obs[t] = {s.obligors, d[t]};
    SpFit out{};
    out.side[0] = fit_panel(s, obs, out.se);
    for (int v = 0; v < kVariants; ++v) {
        std::vector<rc::Objective::Obs> ext = obs;
        for (int k = 0; k < kVariant[v].added; ++k) ext.push_back({s.obligors, counts[v]});
        out.side[1 + v] = fit_panel(s, ext);
    }
    return out;
}

double median(std::vector<double> x) {
    std::vector<double> y;
    for (const double v : x) {
        if (std::isfinite(v)) y.push_back(v);
    }
    if (y.empty()) return kNaN;
    std::sort(y.begin(), y.end());
    const std::size_t m = y.size() / 2;
    return y.size() % 2 ? y[m] : 0.5 * (y[m - 1] + y[m]);
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t R = rc::kReplicates;
    std::vector<std::uint32_t> ids;
    std::string summary_out, replicates_out;
    bool counts_only = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--counts") == 0) {
            counts_only = true;
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
            std::fprintf(stderr, "usage: study_severe_period [--replicates R] [--scenarios ID,...] [--summary-out FILE] "
                                 "[--replicates-out FILE]\n       study_severe_period --counts\n");
            return 2;
        }
    }
    if (ids.empty()) {
        for (std::uint32_t id = 0; id < rc::kScenarios; ++id) ids.push_back(id);
    }
    std::vector<std::array<std::int64_t, kVariants>> counts(ids.size());
    for (std::size_t k = 0; k < ids.size(); ++k) {
        for (int v = 0; v < kVariants; ++v) counts[k][static_cast<std::size_t>(v)] = added_count(rc::scenario(ids[k]), kVariant[v]);
    }
    if (counts_only) {
        std::printf("scenario,count_1in100,count_1in1000,median_count_z0\n");
        for (std::size_t k = 0; k < ids.size(); ++k) {
            const rc::Scenario s = rc::scenario(ids[k]);
            std::printf("%u,%lld,%lld,%lld\n", s.id, static_cast<long long>(counts[k][0]),
                        static_cast<long long>(counts[k][2]),
                        static_cast<long long>(binomial_median(s.obligors, engine::conditional_pd(s.pd, s.rho, 0.0))));
        }
        return 0;
    }

    const auto per = static_cast<std::int64_t>(ids.size());
    std::vector<SpFit> fits(static_cast<std::size_t>(R) * ids.size());
    const auto start = std::chrono::steady_clock::now();
    const std::uint32_t batch = 25;
    for (std::uint32_t r0 = 0; r0 < R; r0 += batch) {
        const std::uint32_t r1 = r0 + batch < R ? r0 + batch : R;
        backends::CpuBackend{}.parallel_for(static_cast<std::int64_t>(r1 - r0) * per, [&](std::int64_t j) {
            const std::uint32_t r = r0 + static_cast<std::uint32_t>(j / per);
            const std::size_t k = static_cast<std::size_t>(j % per);
            std::int64_t c[kVariants];
            for (int v = 0; v < kVariants; ++v) c[v] = counts[k][static_cast<std::size_t>(v)];
            fits[static_cast<std::size_t>(r) * ids.size() + k] = sp_fit(rc::scenario(ids[k]), r, c);
        });
        std::fprintf(stderr, "replicates %u/%u, %.0f s\n", r1, R,
                     std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
    }

    std::string csv =
        "# S-34, sensitivity to severe new periods (K1-K10): one row per scenario and variant, R replicates each.\n"
        "# Written by study_severe_period. Variant none is the original panel. Shifts are extended minus\n"
        "# original, in the original panel's SE units (medians over replicates without a Wald flag; logit\n"
        "# scale for q); rel_* are medians of extended/original - 1 and *_ratio of extended/original ends,\n"
        "# over all replicates. *_covered counts intervals covering the truth (not computed = not covering).\n"
        "scenario,variant,added,added_count,replicates,unflagged,flagged_ext,q_not_computed,"
        "shift_pd_se_hex,shift_rho_se_hex,shift_q_se_hex,rel_pd_hex,rel_rho_hex,rel_q_hex,"
        "pd_lo_ratio_hex,pd_hi_ratio_hex,rho_lo_ratio_hex,rho_hi_ratio_hex,q_lo_ratio_hex,q_hi_ratio_hex,"
        "q_above_old_upper,old_q_below_new_lower,pd_covered,rho_covered,q_covered\n";
    for (std::size_t k = 0; k < ids.size(); ++k) {
        const rc::Scenario s = rc::scenario(ids[k]);
        std::int64_t q_covered_before = 0;
        for (std::uint32_t r = 0; r < R; ++r) q_covered_before += fits[static_cast<std::size_t>(r) * ids.size() + k].side[0].cover[2] == 1;
        for (int v = 0; v <= kVariants; ++v) {
            std::int64_t unflagged = 0, flagged_ext = 0, q_nc = 0, above = 0, below = 0, covered[3] = {};
            std::vector<double> shift[3], rel[3], lo_ratio[3], hi_ratio[3];
            for (std::uint32_t r = 0; r < R; ++r) {
                const SpFit& f = fits[static_cast<std::size_t>(r) * ids.size() + k];
                const Side& o = f.side[0];
                const Side& e = f.side[v];
                flagged_ext += (e.flags & rc::kNoReliableInterval) ? 1 : 0;
                q_nc += e.cover[2] < 0;
                for (int a = 0; a < 3; ++a) {
                    covered[a] += e.cover[a] == 1;
                    rel[a].push_back(e.value[a] / o.value[a] - 1.0);
                    lo_ratio[a].push_back(e.lo[a] / o.lo[a]);
                    hi_ratio[a].push_back(e.hi[a] / o.hi[a]);
                }
                above += e.value[2] > o.hi[2];
                below += o.value[2] < e.lo[2];
                if (o.flags & rc::kNoReliableInterval) continue;
                ++unflagged;
                shift[0].push_back((e.value[0] - o.value[0]) / f.se[0]);
                shift[1].push_back((e.value[1] - o.value[1]) / f.se[1]);
                shift[2].push_back((logit(e.value[2]) - logit(o.value[2])) / f.se[2]);
            }
            char line[2048];
            std::snprintf(line, sizeof line,
                          "%u,%s,%d,%lld,%u,%lld,%lld,%lld,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%lld,%lld,%lld,%lld,%lld\n",
                          s.id, v == 0 ? "none" : kVariant[v - 1].name, v == 0 ? 0 : kVariant[v - 1].added,
                          static_cast<long long>(v == 0 ? 0 : counts[k][static_cast<std::size_t>(v - 1)]), R,
                          static_cast<long long>(unflagged), static_cast<long long>(flagged_ext),
                          static_cast<long long>(q_nc), to_hex(median(shift[0])).c_str(),
                          to_hex(median(shift[1])).c_str(), to_hex(median(shift[2])).c_str(),
                          to_hex(median(rel[0])).c_str(), to_hex(median(rel[1])).c_str(),
                          to_hex(median(rel[2])).c_str(), to_hex(median(lo_ratio[0])).c_str(),
                          to_hex(median(hi_ratio[0])).c_str(), to_hex(median(lo_ratio[1])).c_str(),
                          to_hex(median(hi_ratio[1])).c_str(), to_hex(median(lo_ratio[2])).c_str(),
                          to_hex(median(hi_ratio[2])).c_str(), static_cast<long long>(above),
                          static_cast<long long>(below), static_cast<long long>(covered[0]),
                          static_cast<long long>(covered[1]), static_cast<long long>(covered[2]));
            csv += line;
            if (v == 0) continue;
            std::printf("%2u T %-3lld n %-5lld %-11s d=%-5lld | shift SE: PD %+.2f rho %+.2f q %+.2f | rel q %+.3f | "
                        "q above old %.3f | cover q %.3f -> %.3f\n",
                        s.id, static_cast<long long>(s.periods), static_cast<long long>(s.obligors), kVariant[v - 1].name,
                        static_cast<long long>(counts[k][static_cast<std::size_t>(v - 1)]), median(shift[0]),
                        median(shift[1]), median(shift[2]), median(rel[2]), static_cast<double>(above) / R,
                        static_cast<double>(q_covered_before) / R, static_cast<double>(covered[2]) / R);
        }
    }
    std::printf("%u replicates x %zu scenarios x %d variants in %.0f s\n", R, ids.size(), kVariants,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
    if (!replicates_out.empty()) {
        std::ofstream out(replicates_out, std::ios::binary);
        out << "scenario,replicate,variant,pd,rho,q,pd_lo,pd_hi,rho_lo,rho_hi,q_lo,q_hi,flags,pd_cover,rho_cover,"
               "q_cover,se_pd,se_rho,se_q\n";
        char line[1024];
        for (std::size_t k = 0; k < ids.size(); ++k) {
            for (std::uint32_t r = 0; r < R; ++r) {
                const SpFit& f = fits[static_cast<std::size_t>(r) * ids.size() + k];
                for (int v = 0; v <= kVariants; ++v) {
                    const Side& e = f.side[v];
                    std::snprintf(line, sizeof line,
                                  "%u,%u,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%u,%d,%d,%d,%.17g,%.17g,%.17g\n",
                                  ids[k], r, v, e.value[0], e.value[1], e.value[2], e.lo[0], e.hi[0], e.lo[1], e.hi[1],
                                  e.lo[2], e.hi[2], e.flags, e.cover[0], e.cover[1], e.cover[2], f.se[0], f.se[1],
                                  f.se[2]);
                    out << line;
                }
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
