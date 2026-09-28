// SPDX-License-Identifier: Apache-2.0
// engine/surface_cache.hpp (D-167): rows cached across panels are identical to evaluate_surface's; a
// change in any one setting that determines a row is a cache miss, never a stale hit; rows are filled
// on demand; and a recovery fit through the cache equals the fit without it, bit for bit.
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/grid.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/objectives/vasicek_rate.hpp"
#include "core/quadrature/composite_legendre.hpp"
#include "core/quadrature/gauss_hermite.hpp"
#include "core/quadrature/parity.hpp"
#include "engine/surface.hpp"
#include "engine/surface_cache.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/recovery/recovery.hpp"

namespace {

namespace e = vcal::engine;
namespace q = vcal::quadrature;
using Backend = vcal::backends::CpuBackend;
using Binomial = vcal::objectives::BinomialMixture<vcal::PrecisionF64>;
using BinomialMixed = vcal::objectives::BinomialMixture<vcal::PrecisionMixed>;

vcal::Grid<2> small_grid() { return {{{1e-3, 0.2, 9, vcal::AxisScale::Logit}, {1e-2, 0.5, 7, vcal::AxisScale::Logit}}}; }

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

// The row evaluate_surface computes for one observation under a configuration.
template <class Objective, class Integrator>
std::vector<double> direct_row(const Objective& o, const Integrator& in, const vcal::Grid<2>& g,
                               const typename Objective::Obs& y) {
    std::vector<double> L(static_cast<std::size_t>(g.size()));
    e::evaluate_surface(Backend{1}, o, in, &y, 1, g, L.data());
    return L;
}

}  // namespace

VCAL_TEST(cached_surfaces_equal_evaluate_surface_bit_for_bit) {
    const auto in = q::parity_rule();
    const auto g = small_grid();
    const std::vector<std::vector<Binomial::Obs>> panels = {
        {{1000, 3}, {1000, 12}, {1000, 3}, {1000, 0}}, {{1000, 12}, {1000, 7}, {1000, 0}}, {{500, 12}, {1000, 12}}};
    e::SurfaceRowCache cache;
    for (const auto& p : panels) {
        const auto T = static_cast<std::int64_t>(p.size());
        std::vector<double> a(static_cast<std::size_t>(T * g.size())), b(a.size());
        e::evaluate_surface(Backend{1}, Binomial{}, in, p.data(), T, g, a.data());
        e::evaluate_surface_cached(Backend{1}, cache, Binomial{}, in, p.data(), T, g, b.data());
        VCAL_CHECK(same_bits(a, b));
    }
    // Rows filled on demand, once each: (1000, 3), (1000, 12), (1000, 0), (1000, 7), (500, 12).
    VCAL_CHECK_EQ(cache.size(), 5u);
    VCAL_CHECK_EQ(cache.misses(), 5);
    VCAL_CHECK_EQ(cache.hits(), 3);  // (1000, 12) and (1000, 0) in panel 2, (1000, 12) in panel 3
}

// Every setting that determines a row is part of its key: changing any one of them is a miss, and the
// row returned is the one that configuration computes.
VCAL_TEST(changing_any_row_setting_is_a_cache_miss) {
    const auto base_in = q::parity_rule();
    const auto base_g = small_grid();
    const Binomial::Obs y{1000, 12};
    e::SurfaceRowCache cache;
    const Backend b{1};
    const auto request = [&](const auto& objective, const auto& integrator, const vcal::Grid<2>& g, const auto& obs) {
        const auto config = e::surface_configuration(objective, integrator, g);
        const auto before = cache.misses();
        const auto& row = cache.row(b, config, objective, integrator, g, obs);
        return std::make_pair(cache.misses() - before, same_bits(row, direct_row(objective, integrator, g, obs)));
    };
    VCAL_CHECK(request(Binomial{}, base_in, base_g, y) == std::make_pair(std::int64_t{1}, true));
    VCAL_CHECK(request(Binomial{}, base_in, base_g, y) == std::make_pair(std::int64_t{0}, true));  // the same: a hit

    std::vector<std::pair<std::string, vcal::Grid<2>>> grids;
    for (int a = 0; a < 2; ++a) {
        auto g = base_g;
        g.axis[a].lo *= 1.5;
        grids.push_back({"lower bound " + std::to_string(a), g});
        g = base_g;
        g.axis[a].hi *= 0.9;
        grids.push_back({"upper bound " + std::to_string(a), g});
        g = base_g;
        g.axis[a].n += 2;
        grids.push_back({"spacing " + std::to_string(a), g});
        g = base_g;
        g.axis[a].scale = vcal::AxisScale::Probit;
        grids.push_back({"scale " + std::to_string(a), g});
    }
    for (const auto& [what, g] : grids) {
        const auto r = request(Binomial{}, base_in, g, y);
        if (r.first != 1 || !r.second) vcal::test::note("grid change not a clean miss: " + what);
        VCAL_CHECK(r.first == 1 && r.second);
    }
    // The integrator: Gauss-Hermite N, the composite rule's panels and its Legendre points, and the
    // doubled check rule.
    auto gh = base_in;
    gh.interior.rule = q::gauss_hermite_rule(64);
    auto panels = base_in;
    panels.one_sided.panels = 32;
    auto points = base_in;
    points.one_sided.rule = q::gauss_legendre_rule(12);
    for (const auto& in : {gh, panels, points, q::parity_rule(true)}) {
        const auto r = request(Binomial{}, in, base_g, y);
        VCAL_CHECK(r.first == 1 && r.second);
    }
    // The precision policy (a different objective type).
    const auto rm = request(BinomialMixed{}, q::SplitRule<vcal::PrecisionMixed>{
                                                 q::GaussHermiteAdaptive<vcal::PrecisionMixed>{base_in.interior.rule},
                                                 q::CompositeLegendre<vcal::PrecisionMixed>{base_in.one_sided.rule,
                                                                                            base_in.one_sided.panels}},
                            base_g, BinomialMixed::Obs{1000, 12});
    VCAL_CHECK(rm.first == 1);
    // The zero/all-default rule: the rate objective's treatments, on the same observation.
    using namespace vcal::objectives;
    const RateObs z{0.0, 5e-3};
    VCAL_CHECK(request(VasicekRate<vcal::PrecisionF64, ZeroRates::Censor>{}, base_in, base_g, z).first == 1);
    VCAL_CHECK(request(VasicekRate<vcal::PrecisionF64, ZeroRates::Substitute>{}, base_in, base_g, z).first == 1);
    VCAL_CHECK(request(VasicekRate<vcal::PrecisionF64, ZeroRates::Censor>{}, base_in, base_g, z).first == 0);
    // ... and a different detection limit is a different observation.
    VCAL_CHECK(request(VasicekRate<vcal::PrecisionF64, ZeroRates::Censor>{}, base_in, base_g, RateObs{0.0, 5e-4}).first == 1);
}

// A recovery fit through the cache equals the fit without it, bit for bit, field by field.
VCAL_TEST(recovery_fits_through_the_cache_are_identical) {
    namespace rc = vcal::recovery;
    e::SurfaceRowCache cache;
    int differing = 0;
    for (const std::uint32_t id : {29u, 55u}) {  // T = 20, n = 10^4 and 10^3
        for (std::uint32_t r = 0; r < 2; ++r) {
            const auto a = rc::fit(rc::scenario(id), r);
            const auto c = rc::fit(rc::scenario(id), r, &cache);
            const auto eq = [](const auto& x, const auto& y) { return std::memcmp(&x, &y, sizeof x) == 0; };
            const bool same = eq(a.value, c.value) && eq(a.se, c.se) && eq(a.loglik, c.loglik) && a.flags == c.flags &&
                              eq(a.prof_lo, c.prof_lo) && eq(a.prof_hi, c.prof_hi) && eq(a.prof_flags, c.prof_flags) &&
                              eq(a.boot_lo, c.boot_lo) && eq(a.boot_hi, c.boot_hi) && a.boot_edge == c.boot_edge &&
                              eq(a.q_hat, c.q_hat) && eq(a.q_se_s, c.q_se_s) && eq(a.q_prof_lo, c.q_prof_lo) &&
                              eq(a.q_prof_hi, c.q_prof_hi) && a.q_prof_flags == c.q_prof_flags &&
                              eq(a.q_boot_lo, c.q_boot_lo) && eq(a.q_boot_hi, c.q_boot_hi);
            differing += same ? 0 : 1;
        }
    }
    VCAL_CHECK_EQ(differing, 0);
    VCAL_CHECK(cache.hits() > 0);  // replicates of one scenario share rows
}
