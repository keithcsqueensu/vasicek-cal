// SPDX-License-Identifier: Apache-2.0
//
// Recovery for the Vasicek-rate MLE (M3): the rate model on its own correctly specified data.
//
// The panels are the n -> infinity limit of the recovery panels: for each (PD, rho, T) the factor
// draws z_t of the n = 10^4 recovery scenario (same seed, scenario id and replicate), mapped to rates
// r_t = Phi((Phi^-1(PD) - sqrt(rho) z_t) / sqrt(1 - rho)) with the DGP's own deterministic functions
// (dgp::det_probit, det_ncdf), so the panels are identical on every platform (D-053). Every rate lies
// in (0, 1), so the parity treatment (refuse) never refuses. Each panel is fitted on the recovery box
// and grid, with the profile intervals for PD, rho and the 99.9% conditional PD q (S-23's interval,
// generic in the objective). Estimates are the polished maxima; the closed-form MLE is checked
// against them.
//
// One summary row per (PD, rho, T): bias and RMSE, profile coverage over all replicates (not
// computed = not covering, D-131) with its raw band class, flagged fits, and the worst
// polished-vs-closed-form distance. Verdicts are reviewed in unit_vasicek_rate against the committed
// tests/golden/vasicek_rate/recovery_summary.csv.
//
// Without --write it compares with that summary: counts and classes exactly, doubles to
// TOL_PROFILE_CROSS_PLATFORM_REL (the engine's libm calls may differ in the last bits across
// platforms, as for recovery_replay).
//
// --replay-out FILE writes the panels and the engine's results that validation/scipy/vasicek_rate_mle.py
// replicates: replicates 0 and 1 of every cell (refuse; rates from the DGP), and replicate 0 of the
// nine binomial recovery scenarios with PD 1% and n = 100 under the censored likelihood (count data
// with zero-default periods, detection limit 1/(2n)).
//
//   rate_recovery_harness [--replicates R] [--write FILE] [--replay-out FILE]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "core/objectives/vasicek_rate.hpp"
#include "engine/conditional_pd.hpp"
#include "tests/harness/golden.hpp"
#include "tests/recovery/recovery.hpp"
#include "tests/tolerances.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
namespace o = vcal::objectives;
namespace tol = vcal::tol;
using vcal::test::to_hex;
using Objective = o::VasicekRate<PrecisionF64, o::ZeroRates::Refuse>;
using Censored = o::VasicekRate<PrecisionF64, o::ZeroRates::Censor>;

struct RateFit {
    double value[3];  // polished PD, rho, and q there
    std::int8_t cover[3];
    std::uint32_t flags;
    double closed_form_gap;  // max over PD, rho of |logit polished - logit closed form|; NaN if outside the box
    double residual;         // the profile code's endpoint residual (PD, rho and q)
};

double logit(double v) { return std::log(v) - std::log1p(-v); }

// The rate panel: the recovery scenario's factor draws mapped through the model (see above).
std::vector<o::RateObs> rate_panel(const rc::Scenario& s, std::uint32_t replicate) {
    std::vector<std::int64_t> n(static_cast<std::size_t>(s.periods), s.obligors), d(n.size());
    std::vector<double> z(n.size());
    const dgp::PanelSpec spec{rc::kSeed, s.id, replicate, s.pd, s.rho, n.data(), s.periods};
    if (dgp::simulate_panel(spec, d.data(), z.data()) != dgp::Status::Ok) std::abort();
    const double c = dgp::det_probit(s.pd), sr = std::sqrt(s.rho), s1 = std::sqrt(1.0 - s.rho);
    std::vector<o::RateObs> obs(n.size());
    for (std::size_t t = 0; t < obs.size(); ++t) obs[t] = {dgp::det_ncdf((c - sr * z[t]) / s1), 0.25};
    return obs;
}

// One fit for the replay file: the polished maximum and the PD, rho and q profile intervals.
struct Replay {
    double value[2];
    double lo[3], hi[3];
    std::uint32_t flags, iflags[3];
    double loglik_max;
};

template <class O>
Replay replay_fit(const std::vector<o::RateObs>& obs) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    const auto T = static_cast<std::int64_t>(obs.size());
    std::vector<double> L;
    engine::Estimate2 est{};
    if (engine::calibrate(backends::CpuBackend{1}, O{}, primary, check, obs.data(), T, g, L, est) != engine::Status::Ok) {
        std::abort();
    }
    const auto prof = engine::profile_intervals(O{}, primary, obs.data(), T, g, L, est);
    const auto qi = engine::conditional_pd_interval(O{}, primary, obs.data(), T, g, L, est, prof);
    return {{prof.max_at[0], prof.max_at[1]},
            {prof.lo[0], prof.lo[1], qi.lo},
            {prof.hi[0], prof.hi[1], qi.hi},
            est.flags,
            {prof.flags[0], prof.flags[1], qi.flags},
            prof.loglik_max};
}

RateFit rate_fit(const rc::Scenario& s, std::uint32_t replicate) {
    const auto obs = rate_panel(s, replicate);
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    std::vector<double> L;
    engine::Estimate2 est{};
    if (engine::calibrate(backends::CpuBackend{1}, Objective{}, primary, check, obs.data(), s.periods, g, L, est) !=
        engine::Status::Ok) {
        std::abort();
    }
    const auto prof = engine::profile_intervals(Objective{}, primary, obs.data(), s.periods, g, L, est);
    const auto qi = engine::conditional_pd_interval(Objective{}, primary, obs.data(), s.periods, g, L, est, prof);
    RateFit out{};
    out.value[0] = prof.max_at[0];
    out.value[1] = prof.max_at[1];
    out.value[2] = engine::conditional_pd(prof.max_at[0], prof.max_at[1]);
    out.flags = est.flags;
    out.residual = std::fmax(prof.residual_max, qi.residual_max);
    const double truth[3] = {s.pd, s.rho, engine::conditional_pd(s.pd, s.rho)};
    for (int a = 0; a < 3; ++a) {
        const double lo = a < 2 ? prof.lo[a] : qi.lo, hi = a < 2 ? prof.hi[a] : qi.hi;
        const std::uint32_t f = a < 2 ? prof.flags[a] : qi.flags;
        out.cover[a] = static_cast<std::int8_t>((f & engine::kIntervalNotComputed) ? -1 : (lo <= truth[a] && truth[a] <= hi));
    }
    o::RateClosedForm cf{};
    out.closed_form_gap = std::nan("");
    if (o::vasicek_rate_closed_form(obs.data(), s.periods, cf) && cf.pd > g.axis[0].lo && cf.pd < g.axis[0].hi &&
        cf.rho > g.axis[1].lo && cf.rho < g.axis[1].hi) {
        out.closed_form_gap =
            std::fmax(std::fabs(logit(out.value[0]) - logit(cf.pd)), std::fabs(logit(out.value[1]) - logit(cf.rho)));
    }
    return out;
}

const char* band_class(double coverage, double lo, double hi) {
    return coverage < lo ? "BELOW" : coverage > hi ? "ABOVE" : "PASS";
}

// The data lines of a CSV text (comments dropped), each split at commas.
std::vector<std::vector<std::string>> split_csv(const std::string& text) {
    std::vector<std::vector<std::string>> out;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> cells;
        std::size_t start = 0;
        for (std::size_t i = 0; i <= line.size(); ++i) {
            if (i == line.size() || line[i] == ',') {
                cells.push_back(line.substr(start, i - start));
                start = i + 1;
            }
        }
        out.push_back(cells);
    }
    return out;
}

// Compares the summary with the committed one: every non-hex field exactly, every *_hex double to
// TOL_PROFILE_CROSS_PLATFORM_REL. Returns the number of fields that differ.
int compare(const std::string& got_text, const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::stringstream buf;
    buf << in.rdbuf();
    const auto want = split_csv(buf.str());
    const auto got = split_csv(got_text);
    if (want.empty() || want.size() != got.size() || want[0] != got[0]) {
        std::fprintf(stderr, "the summary's shape differs from %s\n", path.c_str());
        return 1;
    }
    const auto& header = want[0];
    int bad = 0;
    double worst = 0.0;
    for (std::size_t i = 1; i < want.size(); ++i) {
        for (std::size_t j = 0; j < header.size(); ++j) {
            const std::string& a = want[i][j];
            const std::string& b = got[i][j];
            const bool is_hex = header[j].size() > 4 && header[j].compare(header[j].size() - 4, 4, "_hex") == 0;
            bool ok = a == b;
            if (is_hex && !ok) {
                const double x = vcal::test::parse_double(a), y = vcal::test::parse_double(b);
                const double rel = std::fabs(y - x) / std::fmax(std::fabs(x), 1e-300);
                worst = std::fmax(worst, rel);
                ok = rel <= tol::TOL_PROFILE_CROSS_PLATFORM_REL;
            }
            if (!ok) {
                std::fprintf(stderr, "row %zu %s: %s, committed %s\n", i, header[j].c_str(), b.c_str(), a.c_str());
                ++bad;
            }
        }
    }
    std::printf("worst relative difference in the doubles: %.3g\n", worst);
    return bad;
}

}  // namespace

int main(int argc, char** argv) {
    std::uint32_t R = rc::kReplicates;
    std::string write, replay_out;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--replicates") == 0 && i + 1 < argc) {
            R = static_cast<std::uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "--write") == 0 && i + 1 < argc) {
            write = argv[++i];
        } else if (std::strcmp(argv[i], "--replay-out") == 0 && i + 1 < argc) {
            replay_out = argv[++i];
        } else {
            std::fprintf(stderr, "usage: rate_recovery_harness [--replicates R] [--write FILE] [--replay-out FILE]\n");
            return 2;
        }
    }
    // One scenario per (PD, rho, T): the recovery scenario with n = 10^4 supplies the factor draws.
    std::vector<rc::Scenario> scen;
    for (std::uint32_t id = 0; id < rc::kScenarios; ++id) {
        if (id % 3 == 2) scen.push_back(rc::scenario(id));
    }
    const auto per = static_cast<std::int64_t>(scen.size());
    std::vector<RateFit> fits(static_cast<std::size_t>(R) * scen.size());
    const auto start = std::chrono::steady_clock::now();
    backends::CpuBackend{}.parallel_for(static_cast<std::int64_t>(R) * per, [&](std::int64_t j) {
        const auto r = static_cast<std::uint32_t>(j / per);
        const auto k = static_cast<std::size_t>(j % per);
        fits[static_cast<std::size_t>(r) * scen.size() + k] = rate_fit(scen[k], r);
    });
    double lo = 0.0, hi = 0.0;
    rc::coverage_band(static_cast<std::int64_t>(R), lo, hi);
    std::ostringstream csv;
    csv << "# Recovery of the Vasicek-rate MLE on its own model (the n -> infinity limit of the recovery\n"
           "# panels; factor draws of the n = 10^4 scenario). Written by rate_recovery_harness --write.\n"
           "# *_class is the raw position of the profile coverage against the band; reviewed in unit_vasicek_rate.\n"
           "z_source_scenario,pd,rho,periods,replicates,flagged,flagged_edge,flagged_flat,flagged_near_bound,pd_bias_hex,pd_rmse_hex,rho_bias_hex,rho_rmse_hex,"
           "q_median_rel_err_hex,pd_covered,rho_covered,q_covered,not_computed,pd_class,rho_class,q_class,"
           "closed_form_gap_max_hex,residual_max_hex\n";
    for (std::size_t k = 0; k < scen.size(); ++k) {
        const rc::Scenario& s = scen[k];
        rc::CompensatedSum b[2], e[2];
        std::int64_t cov[3] = {}, flagged = 0, nc = 0, edge = 0, flat = 0, near = 0;
        double gap = 0.0, res = 0.0;
        std::vector<double> qrel;
        const double qt = engine::conditional_pd(s.pd, s.rho);
        for (std::uint32_t r = 0; r < R; ++r) {
            const RateFit& f = fits[static_cast<std::size_t>(r) * scen.size() + k];
            const double truth[2] = {s.pd, s.rho};
            for (int a = 0; a < 2; ++a) {
                b[a].add(f.value[a] - truth[a]);
                e[a].add((f.value[a] - truth[a]) * (f.value[a] - truth[a]));
            }
            for (int a = 0; a < 3; ++a) {
                cov[a] += f.cover[a] == 1;
                nc += f.cover[a] < 0;
            }
            flagged += (f.flags & rc::kNoReliableInterval) ? 1 : 0;
            edge += (f.flags & engine::kFlagGridEdge) ? 1 : 0;
            flat += (f.flags & engine::kFlagFlatSurface) ? 1 : 0;
            near += (f.flags & engine::kFlagNearBound) ? 1 : 0;
            if (std::isfinite(f.closed_form_gap)) gap = std::fmax(gap, f.closed_form_gap);
            res = std::fmax(res, f.residual);
            qrel.push_back(f.value[2] / qt - 1.0);
        }
        std::sort(qrel.begin(), qrel.end());
        const double qmed = R % 2 ? qrel[R / 2] : 0.5 * (qrel[R / 2 - 1] + qrel[R / 2]);
        const double Rd = static_cast<double>(R);
        csv << s.id << ',' << s.pd << ',' << s.rho << ',' << s.periods << ',' << R << ',' << flagged << ',' << edge << ',' << flat << ',' << near << ','
            << to_hex(b[0].value() / Rd) << ',' << to_hex(std::sqrt(e[0].value() / Rd)) << ','
            << to_hex(b[1].value() / Rd) << ',' << to_hex(std::sqrt(e[1].value() / Rd)) << ',' << to_hex(qmed) << ','
            << cov[0] << ',' << cov[1] << ',' << cov[2] << ',' << nc << ','
            << band_class(static_cast<double>(cov[0]) / Rd, lo, hi) << ','
            << band_class(static_cast<double>(cov[1]) / Rd, lo, hi) << ','
            << band_class(static_cast<double>(cov[2]) / Rd, lo, hi) << ',' << to_hex(gap) << ',' << to_hex(res) << '\n';
        std::printf("PD %-6g rho %-5g T %-3lld | bias PD %+.2e rho %+.4f | cover PD %.3f rho %.3f q %.3f | flagged %lld | "
                    "closed-form gap %.2g\n",
                    s.pd, s.rho, static_cast<long long>(s.periods), b[0].value() / Rd, b[1].value() / Rd,
                    static_cast<double>(cov[0]) / Rd, static_cast<double>(cov[1]) / Rd, static_cast<double>(cov[2]) / Rd,
                    static_cast<long long>(flagged), gap);
    }
    std::printf("%u replicates x %zu scenarios in %.0f s\n", R, scen.size(),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
    if (!replay_out.empty()) {
        std::ostringstream rp;
        rp << "# Panels and engine results for validation/scipy/vasicek_rate_mle.py. Written by\n"
              "# rate_recovery_harness --replay-out. rates_hex: the period rates, space-separated (a count panel's\n"
              "# are d/n); detect_hex: the detection limit (count data 1/(2n)); ends and estimates on the natural scale.\n"
              "kind,source_scenario,replicate,treatment,detect_hex,rates_hex,pd_hex,rho_hex,loglik_max_hex,pd_lo_hex,"
              "pd_hi_hex,rho_lo_hex,rho_hi_hex,q_lo_hex,q_hi_hex,flags,interval_flags\n";
        const auto emit = [&](const char* kind, std::uint32_t id, std::uint32_t r, const char* treat,
                              const std::vector<o::RateObs>& obs, const Replay& f) {
            rp << kind << ',' << id << ',' << r << ',' << treat << ',' << to_hex(obs[0].detect) << ',';
            for (std::size_t t = 0; t < obs.size(); ++t) rp << (t ? " " : "") << to_hex(obs[t].rate);
            rp << ',' << to_hex(f.value[0]) << ',' << to_hex(f.value[1]) << ',' << to_hex(f.loglik_max);
            for (int a = 0; a < 3; ++a) rp << ',' << to_hex(f.lo[a]) << ',' << to_hex(f.hi[a]);
            rp << ',' << f.flags << ',' << (f.iflags[0] | (f.iflags[1] << 8) | (f.iflags[2] << 16)) << '\n';
        };
        for (const auto& s : scen) {
            for (std::uint32_t r = 0; r < 2; ++r) {
                const auto obs = rate_panel(s, r);
                emit("rates", s.id, r, "refuse", obs, replay_fit<Objective>(obs));
            }
        }
        for (std::uint32_t id = 27; id < 54; id += 3) {  // PD 1%, n = 100: zero-default periods are common
            const rc::Scenario s = rc::scenario(id);
            const auto d = rc::panel(s, 0);
            std::vector<o::RateObs> obs(d.size());
            for (std::size_t t = 0; t < d.size(); ++t) obs[t] = o::rate_obs(s.obligors, d[t]);
            emit("counts", s.id, 0, "censor", obs, replay_fit<Censored>(obs));
        }
        std::ofstream rf(replay_out, std::ios::binary);
        rf << rp.str();
        if (!rf) return 1;
        std::printf("wrote %s\n", replay_out.c_str());
    }
    if (!write.empty()) {
        std::ofstream f(write, std::ios::binary);
        f << csv.str();
        if (!f) return 1;
        std::printf("wrote %s\n", write.c_str());
        return 0;
    }
    const std::string path = vcal::test::golden_path("vasicek_rate/recovery_summary.csv");
    const int bad = compare(csv.str(), path);
    if (bad) {
        std::fprintf(stderr, "%d fields differ from %s\n", bad, path.c_str());
        return 1;
    }
    std::printf("summary matches %s\n", path.c_str());
    return 0;
}
