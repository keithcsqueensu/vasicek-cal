// SPDX-License-Identifier: Apache-2.0
//
// Study S-1, Z-sign invariance, calibration half (studies/z-sign-invariance/PREDICTION.md, D-300).
//
// The engine writes p(z) = Phi((c - sqrt(rho) z) / sqrt(1 - rho)). A test-only objective writes
// the mirrored convention p-(z) = Phi((c + sqrt(rho) z) / sqrt(1 - rho)) = p(-z), with its own
// quadrature hint derived for that convention, and is integrated by the parity rule unchanged.
// Z is integrated out against a symmetric density, so every surface cell, estimate, SE, interval
// and flag must agree with parity within rounding. Checks C1-C12 of the prediction, on its 34
// fixed panels; a wrong-hint control (the mirrored integrand with the parity hint) must be seen.
//
//   study_z_sign_invariance   prints the census and the check table (Markdown); exit status 0
//                             only if the control fires and every check holds, apart from the
//                             reviewed finding of D-301 (panel 26, where rho is not identified).
//
// Never parity: nothing here is reachable from the library.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "tests/recovery/recovery.hpp"

namespace {

using namespace vcal;
namespace rc = vcal::recovery;
using Parity = objectives::BinomialMixture<PrecisionF64>;
using Obs = Parity::Obs;

// --- the mirrored convention (test-only) ------------------------------------------------------

struct MirroredVasicek1F {
    double c, sqrt_rho, sqrt_one_minus_rho;
    double threshold(double z) const { return (c + sqrt_rho * z) / sqrt_one_minus_rho; }  // dx/dz = +beta
    double beta() const { return sqrt_rho / sqrt_one_minus_rho; }
};

struct MirroredLogIntegrand {
    MirroredVasicek1F model;
    double defaults, survivors;
    double operator()(double z) const {
        const double x = model.threshold(z);
        double v = 0.0;
        if (defaults > 0.0) v += defaults * special::log_phi(x);
        if (survivors > 0.0) v += survivors * special::log_phi(-x);
        return v;
    }
};

MirroredLogIntegrand mirrored_integrand(double pd, double rho, std::int64_t n, std::int64_t d) {
    const model::Vasicek1F m = model::make_vasicek1f(pd, rho);  // c, sqrt(rho), sqrt(1 - rho)
    return {{m.c, m.sqrt_rho, m.sqrt_one_minus_rho}, static_cast<double>(d), static_cast<double>(n - d)};
}

struct Slopes {
    double d1, d2;
};

// h-(z) = g-(z) - z^2/2 with x = x-(z): h-'(z) = +beta (d lambda(x) - s lambda(-x)) - z,
// h-''(z) = beta^2 (d lambda'(x) + s lambda'(-x)) - 1.
Slopes mirrored_slopes(const MirroredLogIntegrand& f, double z) {
    const double x = f.model.threshold(z);
    const double b = f.model.beta();
    double g1 = 0.0, g2 = 0.0;
    if (f.defaults > 0.0) {
        g1 += f.defaults * special::inverse_mills(x);
        g2 += f.defaults * model::detail::inverse_mills_slope(x);
    }
    if (f.survivors > 0.0) {
        g1 -= f.survivors * special::inverse_mills(-x);
        g2 += f.survivors * model::detail::inverse_mills_slope(-x);
    }
    return {b * g1 - z, b * b * g2 - 1.0};
}

// The hint derived for the mirrored convention: the same algorithm as parity (large-n start
// shrunk to the prior, bracket, Newton with bisection), on h-; for d in {0, n} the survival
// factor's half-point x_c mapped back through z = (sqrt(1 - rho) x - c) / sqrt(rho).
quadrature::IntegrandHint mirrored_hint(const MirroredLogIntegrand& f) {
    const double n = f.defaults + f.survivors;
    if (!(n > 0.0)) return {0.0, 1.0, 0.0, 1.0, false};
    const double x_l = special::probit((f.defaults + 0.5) / (n + 1.0));
    const double z_l = (f.model.sqrt_one_minus_rho * x_l - f.model.c) / f.model.sqrt_rho;
    const double tau_l = -(mirrored_slopes(f, z_l).d2 + 1.0);
    double z = tau_l * z_l / (tau_l + 1.0);
    const double step0 = 1.0 / std::sqrt(tau_l + 1.0);
    Slopes s = mirrored_slopes(f, z);
    double lo = z, hi = z, step = step0;
    for (int i = 0; i < model::kHintMaxIterations && s.d1 != 0.0; ++i) {
        if (s.d1 > 0.0) {
            lo = hi;
            hi = z + step;
            if (mirrored_slopes(f, hi).d1 < 0.0) break;
        } else {
            hi = lo;
            lo = z - step;
            if (mirrored_slopes(f, lo).d1 > 0.0) break;
        }
        step *= 2.0;
    }
    for (int i = 0; i < model::kHintMaxIterations && s.d1 != 0.0; ++i) {
        if (s.d1 > 0.0) lo = z; else hi = z;
        double next = z - s.d1 / s.d2;
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        const bool done = std::fabs(next - z) <= model::kHintRelativeTolerance * (1.0 + std::fabs(z));
        z = next;
        s = mirrored_slopes(f, z);
        if (done) break;
    }
    const double scale = 1.0 / std::sqrt(-s.d2);
    if (f.defaults > 0.0 && f.survivors > 0.0) return {z, scale, z, scale, false};
    const double q = -std::expm1(-special::constants::kLn2 / n);
    const double b = f.model.beta();
    double x_c, width;
    if (f.defaults == 0.0) {
        x_c = special::probit(q);
        width = 1.0 / (n * b * special::inverse_mills(-x_c));
    } else {
        x_c = special::probit_upper(q);
        width = 1.0 / (n * b * special::inverse_mills(x_c));
    }
    const double z_c = (f.model.sqrt_one_minus_rho * x_c - f.model.c) / f.model.sqrt_rho;
    return {z, scale, z_c, width, true};
}

// The mirrored objective; `wrong_hint` is the control: the mirrored integrand with the parity
// hint, centred on the unmirrored mode and half-point.
template <bool wrong_hint>
struct MirroredObjective : Parity {
    template <class I>
    double log_contrib(const Obs& y, const Theta& th, const I& integrator) const {
        const auto f = mirrored_integrand(th.pd, th.rho, y.n, y.d);
        quadrature::IntegrandHint h{};
        if constexpr (wrong_hint) {
            h = model::binomial_mixture_hint(
                model::make_binomial_log_integrand(model::make_vasicek1f(th.pd, th.rho), y.n, y.d));
        } else {
            h = mirrored_hint(f);
        }
        return special::lbinom(y.n, y.d) + integrator.log_integrate(f, h);
    }
};
using Mirrored = MirroredObjective<false>;
using Control = MirroredObjective<true>;
static_assert(objectives::check_objective<Mirrored>());

// --- the panels (PREDICTION.md) ---------------------------------------------------------------

constexpr std::uint64_t kStudySeed = 0x53315A5349474E53ull;  // "S1ZSIGNS"

struct Panel {
    int id;
    std::string what;
    std::vector<Obs> obs;
};

std::vector<Obs> dgp_panel(std::uint64_t seed, std::uint32_t scenario, double pd, double rho,
                           const std::vector<std::int64_t>& n) {
    std::vector<std::int64_t> d(n.size());
    const dgp::PanelSpec spec{seed, scenario, 0, pd, rho, n.data(), static_cast<std::int64_t>(n.size())};
    if (dgp::simulate_panel(spec, d.data(), nullptr) != dgp::Status::Ok) std::abort();
    std::vector<Obs> o(n.size());
    for (std::size_t t = 0; t < n.size(); ++t) o[t] = {n[t], d[t]};
    return o;
}

std::vector<Obs> cycle(std::int64_t T, std::int64_t n, const std::vector<std::int64_t>& pattern) {
    std::vector<Obs> o(static_cast<std::size_t>(T));
    for (std::int64_t t = 0; t < T; ++t) o[static_cast<std::size_t>(t)] = {n, pattern[static_cast<std::size_t>(t) % pattern.size()]};
    return o;
}

std::vector<Panel> panels() {
    std::vector<Panel> p;
    const std::uint32_t subset[9] = {29, 37, 72, 4, 68, 49, 7, 43, 51};
    for (int i = 0; i < 9; ++i) {
        const rc::Scenario s = rc::scenario(subset[i]);
        p.push_back({i + 1, "recovery " + std::to_string(subset[i]) + " r0",
                     dgp_panel(rc::kSeed, s.id, s.pd, s.rho, std::vector<std::int64_t>(static_cast<std::size_t>(s.periods), s.obligors))});
    }
    struct D {
        int id;
        double pd, rho;
        std::int64_t T, n;
    };
    const D ds[] = {{10, 0.001, 0.02, 40, 100},    {11, 0.0005, 0.30, 20, 50},    {12, 0.15, 0.40, 20, 5},
                    {13, 0.18, 0.45, 20, 2},       {14, 0.01, 0.45, 20, 1000},    {15, 0.05, 0.49, 40, 10000},
                    {16, 0.01, 0.0015, 40, 10000}, {17, 0.0002, 0.10, 40, 10000}, {18, 0.01, 0.12, 20, 1000000},
                    {19, 0.001, 0.05, 20, 1000000}, {20, 0.05, 0.24, 20, 100000}, {22, 0.19, 0.20, 20, 1000}};
    for (const D& d : ds) {
        char buf[96];
        std::snprintf(buf, sizeof buf, "DGP PD %g rho %g T %lld n %lld", d.pd, d.rho, static_cast<long long>(d.T),
                      static_cast<long long>(d.n));
        p.push_back({d.id, buf,
                     dgp_panel(kStudySeed, static_cast<std::uint32_t>(d.id), d.pd, d.rho,
                               std::vector<std::int64_t>(static_cast<std::size_t>(d.T), d.n))});
    }
    std::vector<std::int64_t> n21(20);
    for (int t = 0; t < 20; ++t) n21[static_cast<std::size_t>(t)] = std::llround(std::pow(10.0, 1.0 + 4.0 * t / 19.0));
    p.push_back({21, "DGP PD 0.02 rho 0.15 T 20 n 10..1e5", dgp_panel(kStudySeed, 21, 0.02, 0.15, n21)});
    p.push_back({23, "every d = 0 (n 1000)", cycle(20, 1000, {0})});
    p.push_back({24, "every d = n (n 3)", cycle(10, 3, {3})});
    auto one = cycle(20, 1000, {0});
    one[0].d = 1;
    p.push_back({25, "one default, period 0", one});
    p.push_back({26, "n = 1, d 0,1", cycle(20, 1, {0, 1})});
    p.push_back({27, "d = 10 every period", cycle(20, 1000, {10})});
    p.push_back({28, "d 0,500", cycle(20, 1000, {0, 500})});
    p.push_back({29, "n = 4, d 0..4", cycle(20, 4, {0, 1, 2, 3, 4})});
    p.push_back({30, "n = 1e6, d 0,0,1,3,0", cycle(20, 1000000, {0, 0, 1, 3, 0})});
    p.push_back({31, "d = 2000 of 10000", cycle(20, 10000, {2000})});
    p.push_back({32, "n = 100, d 0,100", cycle(20, 100, {0, 100})});
    p.push_back({33, "n = 1e4, d 1,1,1,1,30", cycle(20, 10000, {1, 1, 1, 1, 30})});
    auto late = cycle(40, 1000, {0});
    late[39].d = 3;
    p.push_back({34, "one late period, d = 3", late});
    std::sort(p.begin(), p.end(), [](const Panel& a, const Panel& b) { return a.id < b.id; });
    return p;
}

// --- one fit, as a user's (and the recovery harness's) --------------------------------------------

struct FitOut {
    engine::Status status;
    std::vector<double> L;
    engine::Estimate2 est;
    engine::ProfileIntervals2 prof;
    std::vector<resample::Replicate2> reps;
    resample::PercentileInterval boot[2];
    std::uint32_t boot_edge;
};

template <class Objective>
FitOut fit(const Panel& p) {
    static const auto primary = quadrature::parity_rule();
    static const auto check = quadrature::parity_rule(true);
    const Grid<2> g = rc::grid();
    const auto T = static_cast<std::int64_t>(p.obs.size());
    FitOut o{};
    o.status = engine::calibrate(backends::CpuBackend{1}, Objective{}, primary, check, p.obs.data(), T, g, o.L, o.est);
    if (o.status != engine::Status::Ok) return o;
    o.prof = engine::profile_intervals(Objective{}, primary, p.obs.data(), T, g, o.L, o.est);
    const std::int64_t B = rc::kBootstrapReplicates;
    const auto idx = resample::bootstrap_indices(kStudySeed ^ (static_cast<std::uint64_t>(p.id) << 32),
                                                 resample::Scheme::IidBootstrap, B, T);
    const auto W = resample::weights_from_indices(idx.data(), B, T, T);
    o.reps.resize(static_cast<std::size_t>(B));
    resample::replicate_estimates(backends::CpuBackend{1}, g, o.L.data(), T, W.data(), B, o.reps.data());
    for (const auto& r : o.reps) o.boot_edge += (r.flags & engine::kFlagGridEdge) ? 1u : 0u;
    o.boot[0] = resample::percentile_interval(o.reps.data(), B, 0);
    o.boot[1] = resample::percentile_interval(o.reps.data(), B, 1);
    return o;
}

// --- comparison -----------------------------------------------------------------------------------

constexpr double kEps = DBL_EPSILON;

double rel(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return 0.0;
    if (a == b) return 0.0;  // includes equal infinities
    if (!std::isfinite(a) || !std::isfinite(b)) return HUGE_VAL;
    return std::fabs(a - b) / std::fmax(std::fabs(a), std::fabs(b));
}
double absdiff(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return 0.0;
    if (a == b) return 0.0;
    if (!std::isfinite(a) || !std::isfinite(b)) return HUGE_VAL;
    return std::fabs(a - b);
}
double logit_diff(double a, double b) {
    return absdiff(grid::to_scaled(AxisScale::Logit, a), grid::to_scaled(AxisScale::Logit, b));
}

// One check: the worst ratio of observed to bound (<= 1 holds), the worst raw value, and where.
struct Check {
    Check(const char* id_, const char* what_, const char* bound_) : id(id_), what(what_), bound(bound_) {}
    const char* id;
    const char* what;
    const char* bound;
    double worst_ratio = 0.0;
    double worst_value = 0.0;
    std::string where = "-";
    std::int64_t compared = 0;
    std::int64_t mismatches = 0;  // comparisons over the bound
    std::string failing = "";     // panels with a comparison over the bound
    double worst_elsewhere = 0.0; // the largest value outside panel 26 (the reviewed finding)
    void see(double value, double bound_value, const std::string& at) {
        ++compared;
        if (at.rfind("panel 26 ", 0) != 0 && !(value <= worst_elsewhere)) worst_elsewhere = value;
        const double r = bound_value > 0.0 ? value / bound_value : (value == 0.0 ? 0.0 : HUGE_VAL);
        if (r > worst_ratio || (std::isnan(r) && !std::isnan(worst_ratio))) {
            worst_ratio = std::isnan(r) ? HUGE_VAL : r;
            worst_value = value;
            where = at;
        }
        if (!(value <= bound_value)) {
            ++mismatches;
            const std::string panel = at.substr(0, at.find(' ', 6));
            if (failing.find(panel + ";") == std::string::npos) failing += panel + ";";
        }
    }
    bool holds() const { return mismatches == 0; }
};

struct CellStats {
    std::int64_t cells = 0, bitwise = 0;
    double worst_interior = 0.0, worst_one_sided = 0.0;  // |diff| in units of eps * scale
};

}  // namespace

int main() {
    const auto all = panels();
    const std::size_t P = all.size();
    std::vector<FitOut> par(P), mir(P), ctl(P);
    backends::CpuBackend{0}.parallel_for(static_cast<std::int64_t>(3 * P), [&](std::int64_t job) {
        const auto i = static_cast<std::size_t>(job % static_cast<std::int64_t>(P));
        switch (job / static_cast<std::int64_t>(P)) {
            case 0: par[i] = fit<Parity>(all[i]); break;
            case 1: mir[i] = fit<Mirrored>(all[i]); break;
            default: ctl[i] = fit<Control>(all[i]); break;
        }
    });

    const Grid<2> g = rc::grid();
    const std::int64_t K = g.size();
    Check c1{"C1", "hint: mode and centre negated, scale and width equal", "bitwise"};
    Check c2{"C2", "surface cells, 0 < d < n: |diff|", "256 eps + 4 eps scale"};
    Check c3{"C3", "surface cells, d in {0, n}: |diff|", "1e-12 + 4 eps scale"};
    Check c4{"C4", "panel surface sum_t l_t: |diff|", "sum of per-period bounds"};
    Check c5{"C5", "grid argmax k*, and every bootstrap replicate's", "identical"};
    Check c6{"C6", "flags: fit, quad-check count, profile, bootstrap edge/excluded/per-replicate", "identical"};
    Check c7{"C7", "PD, rho, loglik at the estimate, relative", "1e-12"};
    Check c8{"C8", "SEs and correlation, relative", "5e-10"};
    Check c9a{"C9", "profile end points, logit units", "2e-9"};
    Check c9b{"C9", "polished maximum l_max, relative", "1e-12"};
    Check c10{"C10", "bootstrap percentile end points, relative", "2e-11"};
    Check c11{"C11", "quadrature check max |l(rule) - l(doubled)|, absolute diff", "1e-12"};
    std::vector<CellStats> cell(P);
    double control_worst_large_n = 0.0;  // largest |diff| - bound over panels 18, 19, 20, 30
    int control_flagged = 0;
    std::string control_flagged_at;

    for (std::size_t i = 0; i < P; ++i) {
        const Panel& p = all[i];
        const auto T = static_cast<std::int64_t>(p.obs.size());
        const std::string pid = "panel " + std::to_string(p.id);
        const FitOut& a = par[i];
        const FitOut& b = mir[i];
        if (a.status != engine::Status::Ok || b.status != engine::Status::Ok) {
            c6.see(1.0, 0.0, pid + ": calibrate status");
            continue;
        }
        // C1: the hints, at every distinct observation and grid point.
        for (std::int64_t t = 0; t < T; ++t) {
            bool seen = false;
            for (std::int64_t u = 0; u < t; ++u) seen = seen || p.obs[static_cast<std::size_t>(u)] == p.obs[static_cast<std::size_t>(t)];
            if (seen) continue;
            const Obs& y = p.obs[static_cast<std::size_t>(t)];
            for (std::int64_t k = 0; k < K; ++k) {
                double v[2];
                g.values(k, v);
                const auto hp = model::binomial_mixture_hint(
                    model::make_binomial_log_integrand(model::make_vasicek1f(v[0], v[1]), y.n, y.d));
                const auto hm = mirrored_hint(mirrored_integrand(v[0], v[1], y.n, y.d));
                const bool same = hm.mode == -hp.mode && hm.scale == hp.scale && hm.center == -hp.center &&
                                  hm.width == hp.width && hm.one_sided == hp.one_sided;
                c1.see(same ? 0.0 : 1.0, 0.0, pid + " t " + std::to_string(t) + " k " + std::to_string(k));
            }
        }
        // C2-C4: surface cells and the panel surface; the control's cells for C12.
        const bool large_n = p.id == 18 || p.id == 19 || p.id == 20 || p.id == 30;
        const std::vector<double> ones(static_cast<std::size_t>(T), 1.0);
        for (std::int64_t k = 0; k < K; ++k) {
            double sum_bound = 0.0;
            for (std::int64_t t = 0; t < T; ++t) {
                const Obs& y = p.obs[static_cast<std::size_t>(t)];
                const double lp = a.L[static_cast<std::size_t>(t * K + k)];
                const double lm = b.L[static_cast<std::size_t>(t * K + k)];
                const double scale = Parity::rounding_scale(y, lp);
                const bool one_sided = y.d == 0 || y.d == y.n;
                const double bound = one_sided ? 1e-12 + 4.0 * kEps * scale : 256.0 * kEps + 4.0 * kEps * scale;
                sum_bound += bound;
                const double diff = absdiff(lp, lm);
                const std::string at = pid + " t " + std::to_string(t) + " k " + std::to_string(k);
                (one_sided ? c3 : c2).see(diff, bound, at);
                ++cell[i].cells;
                if (std::memcmp(&lp, &lm, sizeof lp) == 0) ++cell[i].bitwise;
                double& w = one_sided ? cell[i].worst_one_sided : cell[i].worst_interior;
                w = std::fmax(w, diff / (kEps * std::fmax(scale, 1.0)));
                if (large_n && ctl[i].L.size() == a.L.size()) {
                    const double dc = absdiff(lp, ctl[i].L[static_cast<std::size_t>(t * K + k)]);
                    control_worst_large_n = std::fmax(control_worst_large_n, std::isfinite(dc) ? dc - bound : HUGE_VAL);
                }
            }
            c4.see(absdiff(engine::weighted_sum(a.L.data(), T, K, ones.data(), k),
                           engine::weighted_sum(b.L.data(), T, K, ones.data(), k)),
                   sum_bound, pid + " k " + std::to_string(k));
        }
        // C5, C6: identities.
        c5.see(a.est.grid_index == b.est.grid_index ? 0.0 : 1.0, 0.0, pid + " k*");
        for (std::size_t r = 0; r < a.reps.size(); ++r) {
            c5.see(a.reps[r].grid_index == b.reps[r].grid_index ? 0.0 : 1.0, 0.0, pid + " boot " + std::to_string(r));
            c6.see(a.reps[r].flags == b.reps[r].flags ? 0.0 : 1.0, 0.0, pid + " boot flags " + std::to_string(r));
        }
        c6.see(a.est.flags == b.est.flags ? 0.0 : 1.0, 0.0, pid + " fit flags");
        c6.see(a.est.quad_check_flagged == b.est.quad_check_flagged ? 0.0 : 1.0, 0.0, pid + " quad-check count");
        c6.see(a.prof.flags[0] == b.prof.flags[0] && a.prof.flags[1] == b.prof.flags[1] ? 0.0 : 1.0, 0.0, pid + " profile flags");
        c6.see(a.boot_edge == b.boot_edge && a.boot[0].excluded == b.boot[0].excluded &&
                       a.boot[1].excluded == b.boot[1].excluded ? 0.0 : 1.0,
               0.0, pid + " bootstrap counts");
        // C7-C11.
        const char* nm[2] = {"PD", "rho"};
        for (int q = 0; q < 2; ++q) {
            c7.see(rel(a.est.value[q], b.est.value[q]), 1e-12, pid + " " + nm[q]);
            c8.see(rel(a.est.se[q], b.est.se[q]), 5e-10, pid + " se " + nm[q]);
            c9a.see(logit_diff(a.prof.lo[q], b.prof.lo[q]), 2e-9, pid + " lo " + nm[q]);
            c9a.see(logit_diff(a.prof.hi[q], b.prof.hi[q]), 2e-9, pid + " hi " + nm[q]);
        }
        c7.see(rel(a.est.loglik, b.est.loglik), 1e-12, pid + " loglik");
        c8.see(rel(a.est.corr, b.est.corr), 5e-10, pid + " corr");
        c9b.see(rel(a.prof.loglik_max, b.prof.loglik_max), 1e-12, pid + " l_max");
        c10.see(rel(a.boot[0].lo, b.boot[0].lo), 2e-11, pid + " boot lo PD");
        c10.see(rel(a.boot[0].hi, b.boot[0].hi), 2e-11, pid + " boot hi PD");
        c10.see(rel(a.boot[1].lo, b.boot[1].lo), 2e-11, pid + " boot lo rho");
        c10.see(rel(a.boot[1].hi, b.boot[1].hi), 2e-11, pid + " boot hi rho");
        c11.see(absdiff(a.est.quad_check_max, b.est.quad_check_max), 1e-12, pid + " quad-check max");
        // C12: the control's own quadrature check.
        if (ctl[i].status == engine::Status::Ok && (ctl[i].est.flags & engine::kFlagQuadratureUnconverged)) {
            ++control_flagged;
            control_flagged_at += (control_flagged_at.empty() ? "" : ", ") + std::to_string(p.id);
        }
    }

    // --- report -----------------------------------------------------------------------------------
    std::printf("## Census\n\n");
    std::printf("| # | Panel | T | n | d = 0 | d = n | PD hat | rho hat | fit flags | profile truncation (PD, rho) | rho spread of S at the PD column of k* | rho at k* (parity, mirrored) | bitwise cells | worst C2 (eps x scale) | worst C3 (eps x scale) |\n");
    std::printf("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
    for (std::size_t i = 0; i < P; ++i) {
        const Panel& p = all[i];
        std::int64_t nmin = p.obs[0].n, nmax = p.obs[0].n, zero = 0, full = 0;
        for (const Obs& y : p.obs) {
            nmin = std::min(nmin, y.n);
            nmax = std::max(nmax, y.n);
            zero += y.d == 0;
            full += y.d == y.n;
        }
        const auto& e = par[i].est;
        std::string fl;
        if (e.flags & engine::kFlagGridEdge) fl += "edge ";
        if (e.flags & engine::kFlagFlatSurface) fl += "flat ";
        if (e.flags & engine::kFlagNearBound) fl += "near-bound ";
        if (e.flags & engine::kFlagRefinementRejected) fl += "rejected ";
        if (e.flags & engine::kFlagQuadratureUnconverged) fl += "quad ";
        if (e.flags & engine::kFlagNumeric) fl += "numeric ";
        if (fl.empty()) fl = "-";
        const auto tr = [](std::uint32_t f) {
            std::string s;
            if (f & engine::kIntervalNotComputed) return std::string("n/c");
            if (f & engine::kIntervalLowerTruncated) s += "lo";
            if (f & engine::kIntervalUpperTruncated) s += s.empty() ? "hi" : "+hi";
            return s.empty() ? std::string("-") : s;
        };
        // How far the panel surface moves along rho at the estimate's PD column: ~0 means rho is not
        // identified, and the argmax along rho is decided by rounding.
        double spread = 0.0, rho_p = HUGE_VAL, rho_m = HUGE_VAL;
        if (par[i].status == engine::Status::Ok && mir[i].status == engine::Status::Ok) {
            const auto T = static_cast<std::int64_t>(p.obs.size());
            const std::vector<double> ones(static_cast<std::size_t>(T), 1.0);
            std::int32_t ki[2], km[2];
            g.unflatten(par[i].est.grid_index, ki);
            g.unflatten(mir[i].est.grid_index, km);
            double lo = HUGE_VAL, hi = -HUGE_VAL;
            for (std::int32_t j = 0; j < g.axis[1].n; ++j) {
                const std::int32_t idx[2] = {ki[0], j};
                const double v = engine::weighted_sum(par[i].L.data(), T, g.size(), ones.data(), g.flatten(idx));
                lo = std::fmin(lo, v);
                hi = std::fmax(hi, v);
            }
            spread = hi - lo;
            rho_p = g.axis[1].value_at(ki[1]);
            rho_m = g.axis[1].value_at(km[1]);
        }
        std::printf("| %d | %s | %zu | %lld%s%lld | %lld | %lld | %.6g | %.6g | %s| %s, %s | %.3g | %.4g, %.4g | %.3f | %.2f | %.2f |\n", p.id,
                    p.what.c_str(), p.obs.size(), static_cast<long long>(nmin), nmin == nmax ? "" : "..",
                    nmin == nmax ? 0LL : static_cast<long long>(nmax), static_cast<long long>(zero),
                    static_cast<long long>(full), e.value[0], e.value[1], fl.c_str(), tr(par[i].prof.flags[0]).c_str(),
                    tr(par[i].prof.flags[1]).c_str(), spread, rho_p, rho_m,
                    static_cast<double>(cell[i].bitwise) / static_cast<double>(std::max<std::int64_t>(cell[i].cells, 1)),
                    cell[i].worst_interior, cell[i].worst_one_sided);
    }

    const bool c12 = control_worst_large_n >= 1e-6 && control_flagged >= 1;
    std::printf("\n## Checks\n\n| # | Quantity | Bound | Compared | Over the bound | Largest | Largest / bound | Where | Largest outside panel 26 | Held |\n|---|---|---|---|---|---|---|---|---|---|\n");
    bool ok = true;
    for (const Check* c : {&c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9a, &c9b, &c10, &c11}) {
        std::printf("| %s | %s | %s | %lld | %lld%s%s | %.3g | %.3g | %s | %.3g | %s |\n", c->id, c->what, c->bound,
                    static_cast<long long>(c->compared), static_cast<long long>(c->mismatches),
                    c->failing.empty() ? "" : " in ", c->failing.c_str(), c->worst_value, c->worst_ratio, c->where.c_str(),
                    c->worst_elsewhere, c->holds() ? "yes" : "NO");
        ok = ok && c->holds();
    }
    std::printf("| C12 | control (wrong hint): excess over the C2/C3 bound in panels 18, 19, 20, 30; fits flagged by the quadrature check | >= 1e-6; >= 1 | 34 | | %.3g; %d (panels %s) | | | | %s |\n",
                control_worst_large_n, control_flagged, control_flagged_at.empty() ? "none" : control_flagged_at.c_str(),
                c12 ? "yes" : "NO");
    // The reviewed finding (D-301): panel 26 (n = 1 in every period) does not identify rho, so the
    // panel surface is flat along rho to rounding and the argmax along rho, and the bootstrap's, is
    // decided by rounding. It is pinned: C5, C7 and C10 may fail there and nowhere else, and the
    // exit status stays 0 only while every other comparison holds.
    bool unreviewed = false;
    for (const Check* c : {&c1, &c2, &c3, &c4, &c6, &c8, &c9a, &c9b, &c11}) unreviewed = unreviewed || !c->holds();
    for (const Check* c : {&c5, &c7, &c10}) unreviewed = unreviewed || !(c->failing.empty() || c->failing == "panel 26;");
    std::printf("\nVerdict: %s\n", !c12 ? "NOT TRUSTED (control did not fire)"
                                   : ok ? "PASS"
                                   : unreviewed ? "FINDING, UNREVIEWED"
                                                : "FINDING, reviewed (D-301): panel 26 only, rho not identified at n = 1");
    return c12 && !unreviewed ? 0 : 1;
}
