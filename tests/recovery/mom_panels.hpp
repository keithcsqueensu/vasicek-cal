// SPDX-License-Identifier: Apache-2.0
//
// The method of moments' reference panels and results (see mom_reference.cpp): built by one function,
// so the writer and the regression test (unit_moments) produce them identically.
#pragma once

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "engine/moments.hpp"
#include "tests/harness/golden.hpp"
#include "tests/recovery/recovery.hpp"

namespace vcal::momref {

inline constexpr std::uint64_t kMomSeed = 0x4D334D4F4D56414Cull;  // "M3MOMVAL"; not the recovery seed
inline constexpr std::int64_t kN = 1000;

struct Panel {
    recovery::Scenario cell;  // PD, rho, T of the recovery cell (n is kN here)
    std::uint32_t replicate;
    std::vector<std::int64_t> d;
    std::vector<double> rates;  // the same factor draws mapped through the model (dgp's det_ncdf)
};

inline Panel make_panel(const recovery::Scenario& s, std::uint32_t replicate) {
    Panel p{s, replicate, std::vector<std::int64_t>(static_cast<std::size_t>(s.periods)), {}};
    std::vector<std::int64_t> n(p.d.size(), kN);
    std::vector<double> z(p.d.size());
    const dgp::PanelSpec spec{kMomSeed, s.id, replicate, s.pd, s.rho, n.data(), s.periods};
    if (dgp::simulate_panel(spec, p.d.data(), z.data()) != dgp::Status::Ok) std::abort();
    const double c = dgp::det_probit(s.pd), sr = std::sqrt(s.rho), s1 = std::sqrt(1.0 - s.rho);
    for (const double zt : z) p.rates.push_back(dgp::det_ncdf((c - sr * zt) / s1));
    return p;
}

// Every reference panel: the 27 (PD, rho, T) cells, replicates 0 and 1.
inline std::vector<Panel> reference_panels() {
    std::vector<Panel> out;
    for (std::uint32_t id = 0; id < recovery::kScenarios; id += 3) {
        for (std::uint32_t r = 0; r < 2; ++r) out.push_back(make_panel(recovery::scenario(id), r));
    }
    return out;
}

inline std::string reference_csv() {
    using Obs = objectives::BinomialMixture<PrecisionF64>::Obs;
    const auto integrator = quadrature::parity_rule();
    const Grid<2> g = recovery::grid();
    std::ostringstream csv;
    csv << "# Method-of-moments reference results (engine/moments.hpp). Written by mom_reference --write.\n"
           "# Panels from the recovery DGP under the validation seed M3MOMVAL (not the recovery panels), n = 1000;\n"
           "# kind counts: the default counts; kind rates: the same factor draws as rates. Box: the recovery grid's.\n"
           "kind,cell_scenario,pd,rho,periods,replicate,n,data,pd_hat_hex,pd2_hat_hex,rho_hat_hex,flags\n";
    for (const auto& p : reference_panels()) {
        std::vector<Obs> obs(p.d.size());
        for (std::size_t t = 0; t < obs.size(); ++t) obs[t] = {kN, p.d[t]};
        const auto mc = engine::mom_from_counts(integrator, obs.data(), static_cast<std::int64_t>(obs.size()),
                                                g.axis[0].lo, g.axis[0].hi, g.axis[1].lo, g.axis[1].hi);
        const auto mr = engine::mom_from_rates(integrator, p.rates.data(), static_cast<std::int64_t>(p.rates.size()),
                                               g.axis[0].lo, g.axis[0].hi, g.axis[1].lo, g.axis[1].hi);
        for (int k = 0; k < 2; ++k) {
            const auto& m = k == 0 ? mc : mr;
            csv << (k == 0 ? "counts" : "rates") << ',' << p.cell.id << ',' << p.cell.pd << ',' << p.cell.rho << ','
                << p.cell.periods << ',' << p.replicate << ',' << kN << ',';
            for (std::size_t t = 0; t < p.d.size(); ++t) {
                csv << (t ? " " : "");
                if (k == 0) {
                    csv << p.d[t];
                } else {
                    csv << test::to_hex(p.rates[t]);
                }
            }
            csv << ',' << test::to_hex(m.pd) << ',' << test::to_hex(m.pd2) << ',' << test::to_hex(m.rho) << ','
                << m.flags << '\n';
        }
    }
    return csv.str();
}

}  // namespace vcal::momref
