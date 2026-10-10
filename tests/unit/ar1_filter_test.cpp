// SPDX-License-Identifier: Apache-2.0
// The AR(1)-factor likelihood (M6; D-195; benchmark_models/ar1/filter.hpp): phi = 0 against the per-period
// likelihood bit for bit; the dense reference (ref/, independent) on panels chosen to be hard for the node
// placement; the placement check; smoothness in theta; and the score against Richardson-extrapolated differences.
#include "benchmark_models/ar1/filter.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "benchmark_models/ar1/transition.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/surface.hpp"
#include "ref/ref.hpp"
#include "tests/harness/golden.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace ar1 = vcal::benchmark_models::ar1;
namespace tol = vcal::tol;
using vcal::test::describe;

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

struct Panel {
    std::vector<ar1::Obs> obs;
    std::vector<vcalref::Period> ref;
    std::string name;
};

Panel make_panel(std::int64_t n, const std::vector<std::int64_t>& d, const std::string& name) {
    Panel p;
    for (const auto x : d) {
        p.obs.push_back({n, x});
        p.ref.push_back({n, x});
    }
    p.name = name + " (n " + std::to_string(n) + ")";
    return p;
}

// Adjacent periods with very different modes (owner): a zero-default period next to one at about the 99.9%
// conditional quantile of PD 1%, rho 0.12 (p ~ 0.09), in either order and sandwiched.
std::vector<Panel> adversarial_panels(const std::vector<std::int64_t>& ns) {
    std::vector<Panel> out;
    for (const auto n : ns) {
        const std::int64_t hi = std::max<std::int64_t>(1, std::llround(0.09 * static_cast<double>(n)));
        out.push_back(make_panel(n, {0, hi}, "0 then high"));
        out.push_back(make_panel(n, {hi, 0}, "high then 0"));
        out.push_back(make_panel(n, {0, hi, 0}, "0, high, 0"));
        out.push_back(make_panel(n, {hi, 0, hi}, "high, 0, high"));
    }
    return out;
}

double ll(const Panel& p, const ar1::Theta& th, ar1::FilterWork& w) {
    return ar1::log_likelihood(p.obs.data(), static_cast<std::int64_t>(p.obs.size()), th, vcal::quadrature::parity_rule(), w)
        .value;
}

}  // namespace

// phi = 0 is the per-period likelihood bit for bit, under every placement, with no special case in the code:
// each period's term is BinomialMixture::log_contrib's, and the total is engine::weighted_sum's, as calibrate
// sums it. Panels with zero- and all-default periods (the composite rule) and repeated observations.
VCAL_TEST(ar1_phi_zero_is_the_per_period_likelihood_bit_for_bit) {
    const vcal::objectives::BinomialMixture<vcal::PrecisionF64> objective;
    const auto rule = vcal::quadrature::parity_rule();
    const std::vector<std::vector<ar1::Obs>> panels = {
        {{100, 3}, {1000, 0}, {10000, 250}, {20, 20}, {500, 7}, {100, 3}},
        {{1, 0}, {1, 1}, {2, 1}},
        {{10000, 0}, {10000, 0}, {10000, 1}, {10000, 400}}};
    ar1::FilterWork w;
    int compared = 0;
    for (const auto& obs : panels) {
        const auto T = static_cast<std::int64_t>(obs.size());
        for (const double pd : {0.001, 0.01, 0.05, 0.3}) {
            for (const double rho : {0.02, 0.12, 0.5, 0.9}) {
                std::vector<double> want;
                for (const auto& o : obs) want.push_back(objective.log_contrib({o.n, o.d}, {pd, rho}, rule));
                const std::vector<double> ones(obs.size(), 1.0);
                const double total = vcal::engine::weighted_sum(want.data(), T, 1, ones.data(), 0);
                for (const auto placement : {ar1::Placement::Prior, ar1::Placement::Predictive}) {
                    const double got = ar1::detail::forward(obs.data(), T, {pd, rho, 0.0}, rule, placement, w);
                    VCAL_CHECK(bits_equal(got, total));
                }
                const auto est = ar1::log_likelihood(obs.data(), T, {pd, rho, 0.0}, rule, w);
                VCAL_CHECK(bits_equal(est.value, total));
                VCAL_CHECK_EQ(est.placement, 0.0);
                for (std::size_t t = 0; t < obs.size(); ++t) VCAL_CHECK(bits_equal(w.contrib[t], want[t]));
                ++compared;
            }
        }
    }
    vcal::test::note(std::to_string(compared) + " (panel, PD, rho) points bit for bit");
}

// The estimate (Neighbours placement) against ref's dense filter, which shares no node placement or quadrature
// with it, on the adversarial panels at high, low and negative phi.
VCAL_TEST(ar1_matches_the_dense_reference_on_adversarial_panels) {
    ar1::FilterWork w;
    double worst = 0.0;
    for (const auto& p : adversarial_panels({20, 100, 1000, 10000})) {
        for (const double phi : {-0.6, 0.6, 0.95}) {
            const ar1::Theta th{0.01, 0.12, phi};
            const double got = ll(p, th, w);
            const double want = vcalref::ar1_loglik(th.pd, th.rho, th.phi, p.ref);
            const double rel = std::fabs(got - want) / std::fmax(1.0, std::fabs(want));
            worst = std::fmax(worst, rel);
            if (!(rel <= tol::TOL_AR1_REFERENCE_REL)) {
                vcal::test::report_failure(__FILE__, __LINE__,
                                           p.name + ", phi " + describe(phi) + ": " + describe(got) + " vs ref " +
                                               describe(want));
            }
        }
    }
    vcal::test::note("worst relative difference from ref: " + describe(worst));
}

// D-195's finding, kept as a test: with nodes placed for the prior, the doubled-node check agrees to 1e-11 while
// both rules miss the mass by most of a nat, so a self-consistency check cannot detect placement error. The
// placement check (Neighbours against Predictive) is what sees it; the estimate is right.
VCAL_TEST(ar1_doubled_nodes_cannot_see_placement_error_and_the_placement_check_does) {
    const auto p = make_panel(10000, {900, 0, 900}, "high, 0, high");
    const ar1::Theta th{0.01, 0.12, 0.95};
    const auto T = static_cast<std::int64_t>(p.obs.size());
    ar1::FilterWork w;
    const double want = vcalref::ar1_loglik(th.pd, th.rho, th.phi, p.ref);
    const double prior_n = ar1::detail::forward(p.obs.data(), T, th, vcal::quadrature::parity_rule(false),
                                                ar1::Placement::Prior, w);
    const double prior_2n = ar1::detail::forward(p.obs.data(), T, th, vcal::quadrature::parity_rule(true),
                                                 ar1::Placement::Prior, w);
    vcal::test::note("prior placement: N vs 2N " + describe(prior_n - prior_2n) + ", vs ref " + describe(prior_n - want));
    VCAL_CHECK(std::fabs(prior_n - prior_2n) < 1e-9);  // the self-consistency check passes ...
    VCAL_CHECK(std::fabs(prior_n - want) > 0.5);       // ... while the value is wrong

    // A zero-default period before a crisis: the predictive (one forward pass) misses the earlier period's tail.
    const auto q = make_panel(10000, {0, 900}, "0 then high");
    const auto est = ar1::log_likelihood(q.obs.data(), 2, th, vcal::quadrature::parity_rule(), w);
    const double ref_q = vcalref::ar1_loglik(th.pd, th.rho, th.phi, q.ref);
    vcal::test::note("0 then high: placement check " + describe(est.placement) + ", estimate vs ref " +
                     describe(est.value - ref_q));
    VCAL_CHECK(est.placement > 1e-6);
    VCAL_CHECK(std::fabs(est.value - ref_q) <= tol::TOL_AR1_REFERENCE_REL * std::fabs(ref_q));
}

// The score against Richardson-extrapolated central differences of the log-likelihood, on ordinary and
// adversarial panels, at negative, zero and high phi.
VCAL_TEST(ar1_score_matches_richardson_differences) {
    struct Case {
        std::int64_t n;
        std::vector<std::int64_t> d;
        ar1::Theta th;
    };
    const std::vector<Case> cases = {{100, {0, 3, 1, 7, 0, 2}, {0.01, 0.12, 0.5}},
                                     {1000, {90, 0, 90}, {0.01, 0.12, 0.9}},
                                     {10000, {0, 900, 0}, {0.01, 0.12, 0.95}},
                                     {10000, {50, 120, 80, 0, 300}, {0.02, 0.3, -0.4}},
                                     {20, {0, 0, 1, 0}, {0.05, 0.05, 0.0}},
                                     {1000, {5, 12, 30, 8}, {0.01, 0.12, 0.0}}};
    const auto rule = vcal::quadrature::parity_rule();
    ar1::FilterWork w;
    double worst = 0.0;
    for (const auto& c : cases) {
        std::vector<ar1::Obs> obs;
        for (const auto x : c.d) obs.push_back({c.n, x});
        const auto T = static_cast<std::int64_t>(obs.size());
        const auto L = [&](const double (&v)[3]) {
            return ar1::log_likelihood(obs.data(), T, {v[0], v[1], v[2]}, rule, w).value;
        };
        (void)ar1::log_likelihood(obs.data(), T, c.th, rule, w);
        double g[3];
        ar1::score(obs.data(), T, c.th, w, g);
        const double base[3] = {c.th.pd, c.th.rho, c.th.phi};
        for (int i = 0; i < 3; ++i) {
            const double h0 = i == 2 ? 1e-3 : 1e-3 * base[i];
            const auto central = [&](double h) {
                double p[3] = {base[0], base[1], base[2]}, m[3] = {base[0], base[1], base[2]};
                p[i] += h;
                m[i] -= h;
                return (L(p) - L(m)) / (2.0 * h);
            };
            const double d1 = central(h0), d2 = central(h0 / 2), d3 = central(h0 / 4);
            const double r1 = (4.0 * d2 - d1) / 3.0, r2 = (4.0 * d3 - d2) / 3.0;
            const double fd = (16.0 * r2 - r1) / 15.0;
            const double rel = std::fabs(g[i] - fd) / std::fmax(1.0, std::fabs(fd));
            worst = std::fmax(worst, rel);
            if (!(rel <= tol::TOL_AR1_SCORE_REL)) {
                vcal::test::report_failure(__FILE__, __LINE__,
                                           "n " + std::to_string(c.n) + ", parameter " + std::to_string(i) + ": " +
                                               describe(g[i]) + " vs " + describe(fd));
            }
        }
    }
    vcal::test::note("worst relative score error: " + describe(worst));
}

// The per-kernel CTest entries' guard: a kernel VCAL_AR1_KERNEL names that this CPU or build lacks skips the entry
// visibly (exit 77, SKIP_RETURN_CODE), never passing without having run that path.
VCAL_TEST(ar1_forced_kernel_is_available) {
    namespace k = vcal::benchmark_models::ar1::kernel;
    const std::string env = vcal::test::env_or_empty("VCAL_AR1_KERNEL");
    if (!env.empty() && k::selected() == k::Kernel::Unavailable) {
        std::printf("SKIPPED: VCAL_AR1_KERNEL=%s is not a kernel this CPU or build has\n", env.c_str());
        std::fflush(stdout);
        std::exit(77);
    }
    vcal::test::note(std::string("transition kernel: ") + k::name(k::selected()));
}

// P-18 (D-199): the filter's transition sums run on whichever kernel is selected (scalar, SSE2, AVX2), and the
// log-likelihood, the placement check, the per-period terms and the score are the same bit for bit under each, on
// the adversarial panels at negative, zero and high phi. CTest also runs the whole filter suite with each kernel
// forced (unit_ar1_filter_kernel_<name>).
VCAL_TEST(ar1_filter_is_the_same_under_every_kernel) {
    namespace k = vcal::benchmark_models::ar1::kernel;
    struct Restore {
        k::Kernel saved = k::selected();
        ~Restore() { k::select(saved); }
    } restore;
    std::vector<double> ref;
    bool have_ref = false;
    int kernels = 0;
    for (const auto kern : {k::Kernel::Scalar, k::Kernel::Sse2, k::Kernel::Avx2}) {
        if (!k::select(kern)) {
            vcal::test::note(std::string("kernel ") + k::name(kern) + " not available here");
            continue;
        }
        ++kernels;
        std::vector<double> got;
        ar1::FilterWork w;
        for (const auto& p : adversarial_panels({20, 1000, 10000})) {
            for (const double phi : {-0.6, 0.0, 0.6, 0.95}) {
                const ar1::Theta th{0.01, 0.12, phi};
                const auto T = static_cast<std::int64_t>(p.obs.size());
                const auto l = ar1::log_likelihood(p.obs.data(), T, th, vcal::quadrature::parity_rule(), w);
                double g[3];
                ar1::score(p.obs.data(), T, th, w, g);
                got.push_back(l.value);
                got.push_back(l.placement);
                got.insert(got.end(), w.contrib.begin(), w.contrib.end());
                got.insert(got.end(), g, g + 3);
            }
        }
        if (!have_ref) {
            ref = got;
            have_ref = true;
            continue;
        }
        VCAL_REQUIRE(got.size() == ref.size());
        for (std::size_t i = 0; i < got.size(); ++i) VCAL_CHECK(bits_equal(got[i], ref[i]));
    }
    vcal::test::note(std::to_string(kernels) + " kernels, " + std::to_string(ref.size()) + " numbers bit for bit");
}

// The committed goldens (tests/golden/ar1/likelihood.csv, tests/studies/ar1_goldens.cpp): the engine reproduces its
// own values to TOL_AR1_REFERENCE_REL on every platform, and they agree with ref's to the same tolerance. The scipy
// replication (validation/scipy/ar1_binomial_mle.py) checks the same rows three ways, in two formulations of the
// factor's law (transition densities, and the explicit joint normal density).
VCAL_TEST(ar1_likelihood_goldens_reproduce) {
    const auto t = vcal::test::read_golden_csv("ar1/likelihood.csv");
    const auto rule = vcal::quadrature::parity_rule();
    ar1::FilterWork w;
    double worst_engine = 0.0, worst_ref = 0.0;
    for (const auto& r : t.rows) {
        std::vector<ar1::Obs> obs;
        std::istringstream ns(r[t.column("n")]), ds(r[t.column("d")]);
        std::int64_t n = 0, d = 0;
        while (ns >> n && ds >> d) obs.push_back({n, d});
        const ar1::Theta th{vcal::test::parse_double(r[t.column("pd_hex")]), vcal::test::parse_double(r[t.column("rho_hex")]),
                            vcal::test::parse_double(r[t.column("phi_hex")])};
        const double want = vcal::test::parse_double(r[t.column("loglik_hex")]);
        const double ref = vcal::test::parse_double(r[t.column("ref_loglik_hex")]);
        const double got = ar1::log_likelihood(obs.data(), static_cast<std::int64_t>(obs.size()), th, rule, w).value;
        const double scale = std::fmax(1.0, std::fabs(want));
        worst_engine = std::fmax(worst_engine, std::fabs(got - want) / scale);
        worst_ref = std::fmax(worst_ref, std::fabs(want - ref) / scale);
        VCAL_CHECK(std::fabs(got - want) <= tol::TOL_AR1_REFERENCE_REL * scale);
        VCAL_CHECK(std::fabs(want - ref) <= tol::TOL_AR1_REFERENCE_REL * scale);
    }
    vcal::test::note(std::to_string(t.rows.size()) + " rows; engine vs golden " + describe(worst_engine) +
                     ", golden vs ref " + describe(worst_ref));
}

// D-195's smoothness condition: the node placement moves with theta, so a clamp, branch or floor in it would put
// a kink into the likelihood. Second differences of core minus ref (smooth by construction: one fixed grid)
// along lines in phi and rho, on the adversarial panels, against core's own second differences. Slow: ref is
// O(T G^2) per point.
VCAL_TEST(ar1_likelihood_is_smooth_in_phi_and_rho) {
    std::vector<Panel> panels = {make_panel(1000, {90, 0, 90}, "high, 0, high"),
                                 make_panel(10000, {0, 900, 0}, "0, high, 0"),
                                 make_panel(100, {0, 9, 0, 3}, "mixed"),
                                 make_panel(10000, {50, 120, 80, 0, 300}, "five periods")};
    ar1::FilterWork w;
    double worst = 0.0;
    for (const auto& p : panels) {
        for (int line = 0; line < 2; ++line) {
            constexpr std::size_t K = 41;
            std::vector<double> core(K), diff(K);
            for (std::size_t k = 0; k < K; ++k) {
                const double f = static_cast<double>(k) / static_cast<double>(K - 1);
                ar1::Theta th{0.01, 0.12, 0.6};
                if (line == 0) {
                    th.phi = -ar1::kPhiBound + 2.0 * ar1::kPhiBound * f;
                } else {
                    th.rho = 0.02 + 0.48 * f;
                }
                core[k] = ll(p, th, w);
                diff[k] = core[k] - vcalref::ar1_loglik(th.pd, th.rho, th.phi, p.ref);
            }
            double max_core = 0.0, max_diff = 0.0;
            for (std::size_t k = 1; k + 1 < K; ++k) {
                max_core = std::fmax(max_core, std::fabs(core[k + 1] - 2.0 * core[k] + core[k - 1]));
                max_diff = std::fmax(max_diff, std::fabs(diff[k + 1] - 2.0 * diff[k] + diff[k - 1]));
            }
            const double ratio = max_diff / max_core;
            worst = std::fmax(worst, ratio);
            if (!(ratio <= tol::TOL_AR1_SMOOTH_D2_REL)) {
                vcal::test::report_failure(__FILE__, __LINE__,
                                           p.name + (line == 0 ? " along phi" : " along rho") +
                                               ": second differences of core - ref at " + describe(ratio) +
                                               " of core's");
            }
        }
    }
    vcal::test::note("worst second-difference ratio: " + describe(worst));
}
