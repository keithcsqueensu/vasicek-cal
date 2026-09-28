// SPDX-License-Identifier: Apache-2.0
// engine/moments.hpp (M3): the method of moments' joint default probability against mpmath goldens,
// its inversion, its edge cases and weights, large-T consistency on a panel of its own seed, and the
// committed reference results (a cross-platform regression; the scipy replication reads the same file).
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "engine/moments.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/recovery/mom_panels.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace e = vcal::engine;
namespace tol = vcal::tol;
using Obs = vcal::objectives::BinomialMixture<vcal::PrecisionF64>::Obs;
constexpr double kPdLo = 1e-4, kPdHi = 0.2, kRhoLo = 1e-3, kRhoHi = 0.5;

double hex(const std::string& s) { return vcal::test::parse_double(s); }
double logit(double v) { return std::log(v) - std::log1p(-v); }

std::string sci(double v) {
    char b[32];
    std::snprintf(b, sizeof b, "%.3g", v);
    return b;
}

}  // namespace

VCAL_TEST(mom_joint_default_matches_mpmath) {
    const auto g = vcal::test::read_golden_csv("moments/joint_default.csv");
    VCAL_REQUIRE(g.rows.size() == 25u);
    const auto integrator = vcal::quadrature::parity_rule();
    double worst = 0.0;
    for (const auto& r : g.rows) {
        const double got = e::log_joint_default(integrator, hex(r[g.column("pd_hex")]), hex(r[g.column("rho_hex")]));
        const double err = std::fabs(got - hex(r[g.column("log_pd2_hex")]));
        worst = std::fmax(worst, err);
        VCAL_CHECK(err <= tol::TOL_MOM_JOINT_DEFAULT_LOG_ABS);
    }
    vcal::test::note("worst |log PD_2 - mpmath|: " + sci(worst));
}

// Inverting the golden joint default probability gives back rho.
VCAL_TEST(mom_inversion_recovers_rho) {
    const auto g = vcal::test::read_golden_csv("moments/joint_default.csv");
    const auto integrator = vcal::quadrature::parity_rule();
    double worst = 0.0;
    for (const auto& r : g.rows) {
        const double pd = hex(r[g.column("pd_hex")]), rho = hex(r[g.column("rho_hex")]);
        if (!(rho > kRhoLo && rho < kRhoHi)) continue;  // the ends are the box's bounds: flagged, below
        const auto m = e::mom_from_moments(integrator, pd, std::exp(hex(r[g.column("log_pd2_hex")])), kPdLo, kPdHi, kRhoLo, kRhoHi);
        VCAL_CHECK_EQ(m.flags, 0u);
        const double du = std::fabs(logit(m.rho) - logit(rho));
        worst = std::fmax(worst, du);
        VCAL_CHECK(du <= tol::TOL_MOM_INVERSION_U);
    }
    vcal::test::note("worst |logit rho-hat - logit rho|: " + sci(worst));
}

VCAL_TEST(mom_edge_cases_are_flagged_not_extrapolated) {
    const auto integrator = vcal::quadrature::parity_rule();
    // No defaults: refused.
    const std::vector<Obs> none = {{100, 0}, {100, 0}, {100, 0}};
    const auto m0 = e::mom_from_counts(integrator, none.data(), 3, kPdLo, kPdHi, kRhoLo, kRhoHi);
    VCAL_CHECK(m0.flags & e::kMomRefused);
    VCAL_CHECK(std::isnan(m0.rho));
    // Defaults but never two in a period: PD_2-hat = 0, rho at the floor.
    const std::vector<Obs> single = {{100, 1}, {100, 0}, {100, 1}};
    const auto m1 = e::mom_from_counts(integrator, single.data(), 3, kPdLo, kPdHi, kRhoLo, kRhoHi);
    VCAL_CHECK(m1.flags & e::kMomRhoAtFloor);
    VCAL_CHECK(m1.rho == kRhoLo);
    // Clustering beyond the cap: every default in one period.
    const std::vector<Obs> clustered = {{100, 60}, {100, 0}, {100, 0}, {100, 0}};
    const auto m2 = e::mom_from_counts(integrator, clustered.data(), 4, kPdLo, kPdHi, kRhoLo, kRhoHi);
    VCAL_CHECK(m2.flags & e::kMomRhoAtCap);
    VCAL_CHECK(m2.rho == kRhoHi);
    // PD-hat outside the PD box is flagged, and rho is still solved at PD-hat.
    const std::vector<Obs> high = {{100, 30}, {100, 20}, {100, 25}};
    const auto m3 = e::mom_from_counts(integrator, high.data(), 3, kPdLo, kPdHi, kRhoLo, kRhoHi);
    VCAL_CHECK(m3.flags & e::kMomPdOutsideBox);
    VCAL_CHECK(m3.pd == 75.0 / 300.0);
}

// A weight row is a resampling of periods (D-042): weight 2 on a period equals the period twice.
VCAL_TEST(mom_weights_reweight_the_sums) {
    const auto integrator = vcal::quadrature::parity_rule();
    const std::vector<Obs> panel = {{500, 7}, {500, 2}, {500, 12}, {500, 4}};
    const std::vector<Obs> doubled = {{500, 7}, {500, 2}, {500, 2}, {500, 12}, {500, 4}};
    const double w[4] = {1.0, 2.0, 1.0, 1.0};
    const auto a = e::mom_from_counts(integrator, panel.data(), 4, kPdLo, kPdHi, kRhoLo, kRhoHi, w);
    const auto b = e::mom_from_counts(integrator, doubled.data(), 5, kPdLo, kPdHi, kRhoLo, kRhoHi);
    VCAL_CHECK(a.pd == b.pd && a.pd2 == b.pd2 && a.rho == b.rho && a.flags == b.flags);
}

// Consistency: on one long panel of its own seed (T = 20,000, n = 1,000, PD 5%, rho 0.12), rho-hat
// lies within TOL_MOM_CONSISTENCY_Z of its large-T standard deviation from the truth, and PD-hat
// within the same number of its binomial-mixture standard deviations.
VCAL_TEST(mom_is_consistent_on_a_long_panel) {
    const std::int64_t T = 20000, n = 1000;
    const double pd = 0.05, rho = 0.12;
    std::vector<std::int64_t> ns(static_cast<std::size_t>(T), n), d(ns.size());
    const vcal::dgp::PanelSpec spec{vcal::momref::kMomSeed, 9999, 0, pd, rho, ns.data(), T};
    VCAL_REQUIRE(vcal::dgp::simulate_panel(spec, d.data(), nullptr) == vcal::dgp::Status::Ok);
    std::vector<Obs> obs(ns.size());
    for (std::size_t t = 0; t < obs.size(); ++t) obs[t] = {n, d[t]};
    const auto integrator = vcal::quadrature::parity_rule();
    const auto m = e::mom_from_counts(integrator, obs.data(), T, kPdLo, kPdHi, kRhoLo, kRhoHi);
    VCAL_REQUIRE(m.flags == 0u);
    // Standard deviations from the model, not fitted: PD-hat's is sqrt(var(p(Z)) / T) (the binomial
    // part is negligible at n = 1,000); rho-hat's, by the delta method through PD_2, is bounded by
    // 2 x the MLE's large-T sd rho (1 - rho) sqrt(2 / T) (the MoM efficiency of S-8's reference).
    const double pd2 = std::exp(e::log_joint_default(integrator, pd, rho));
    const double sd_pd = std::sqrt((pd2 - pd * pd) / static_cast<double>(T));
    const double sd_rho = 2.0 * rho * (1.0 - rho) * std::sqrt(2.0 / static_cast<double>(T));
    const double zp = (m.pd - pd) / sd_pd, zr = (m.rho - rho) / sd_rho;
    vcal::test::note("PD-hat " + sci(m.pd) + " (z " + sci(zp) + "), rho-hat " + sci(m.rho) + " (z " + sci(zr) + ")");
    VCAL_CHECK(std::fabs(zp) <= tol::TOL_MOM_CONSISTENCY_Z);
    VCAL_CHECK(std::fabs(zr) <= tol::TOL_MOM_CONSISTENCY_Z);
}

// The committed reference results (mom_reference --write): recomputed, flags exactly, doubles to
// TOL_MOM_REFERENCE_REL (a root found to 1e-12 in logit, and libm's last digits across platforms).
VCAL_TEST(mom_reference_results_reproduce) {
    const auto want = vcal::test::read_golden_csv("moments/reference.csv");
    const auto got_text = vcal::momref::reference_csv();
    VCAL_REQUIRE(want.rows.size() == 108u);
    std::vector<std::vector<std::string>> got;
    {
        std::size_t start = 0;
        bool header = true;
        while (start < got_text.size()) {
            const std::size_t end = got_text.find('\n', start);
            const std::string line = got_text.substr(start, end - start);
            start = end + 1;
            if (line.empty() || line[0] == '#') continue;
            if (header) {
                header = false;
                continue;
            }
            std::vector<std::string> cells;
            std::size_t a = 0;
            for (std::size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == ',') {
                    cells.push_back(line.substr(a, i - a));
                    a = i + 1;
                }
            }
            got.push_back(cells);
        }
    }
    VCAL_REQUIRE(got.size() == want.rows.size());
    double worst = 0.0;
    int bad = 0;
    for (std::size_t i = 0; i < got.size(); ++i) {
        const auto& w = want.rows[i];
        const auto& g = got[i];
        for (const char* col : {"kind", "cell_scenario", "periods", "replicate", "flags"}) {
            if (g[want.column(col)] != w[want.column(col)]) ++bad;
        }
        // The panel itself, exactly: counts as text; rates by value, since printf's %a text for the
        // same double differs between C runtimes (MSVC keeps trailing zeros, glibc drops them).
        if (w[want.column("kind")] == "counts") {
            if (g[want.column("data")] != w[want.column("data")]) ++bad;
        } else {
            const auto split = [](const std::string& s) {
                std::vector<double> v;
                std::size_t a = 0;
                for (std::size_t k = 0; k <= s.size(); ++k) {
                    if (k == s.size() || s[k] == ' ') {
                        v.push_back(vcal::test::parse_double(s.substr(a, k - a)));
                        a = k + 1;
                    }
                }
                return v;
            };
            if (split(g[want.column("data")]) != split(w[want.column("data")])) ++bad;
        }
        for (const char* col : {"pd_hat_hex", "pd2_hat_hex", "rho_hat_hex"}) {
            const double x = hex(w[want.column(col)]), y = hex(g[want.column(col)]);
            if (std::isnan(x) && std::isnan(y)) continue;
            const double rel = x == y ? 0.0 : std::fabs(y / x - 1.0);
            worst = std::fmax(worst, rel);
            if (!(rel <= tol::TOL_MOM_REFERENCE_REL)) ++bad;
        }
    }
    vcal::test::note("worst relative difference: " + sci(worst));
    VCAL_CHECK_EQ(bad, 0);
}
