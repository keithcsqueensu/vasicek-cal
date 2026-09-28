// SPDX-License-Identifier: Apache-2.0
//
// The grid-Bayesian estimator's reference panels and results (see posterior_reference.cpp), built by
// one function each so the writer and the regression test produce them identically.
#pragma once

#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "backends/cpu/cpu_backend.hpp"
#include "core/objectives/binomial_mixture.hpp"
#include "core/quadrature/parity.hpp"
#include "dgp/dgp.hpp"
#include "engine/calibrate.hpp"
#include "engine/posterior.hpp"
#include "tests/harness/golden.hpp"
#include "tests/recovery/recovery.hpp"

namespace vcal::postref {

inline constexpr std::uint64_t kBayesSeed = 0x4D3342415956414Cull;  // "M3BAYVAL"; not the recovery seed
inline constexpr std::int64_t kN = 100;
inline constexpr std::int64_t kT = 20;
using Objective = objectives::BinomialMixture<PrecisionF64>;

inline Grid<2> jeffreys_check_grid() {
    return {{{1e-4, 0.2, 13, AxisScale::Logit}, {1e-3, 0.5, 9, AxisScale::Logit}}};
}

inline std::string jeffreys_csv() {
    const auto integrator = quadrature::parity_rule();
    const Grid<2> g = jeffreys_check_grid();
    const auto jt = engine::jeffreys_table(backends::CpuBackend{}, Objective{}, integrator, kN, g);
    std::ostringstream csv;
    csv << "# Jeffreys table (engine::jeffreys_table) for n = 100, 13 x 9 logit grid over the box. Written by\n"
           "# posterior_reference --write. log sqrt det I_u, I_u the Fisher information of one period in logit u.\n"
           "n,pd_hex,rho_hex,log_prior_hex\n";
    for (std::int64_t k = 0; k < g.size(); ++k) {
        double v[2];
        g.values(k, v);
        csv << kN << ',' << test::to_hex(v[0]) << ',' << test::to_hex(v[1]) << ','
            << test::to_hex(jt.log_prior[static_cast<std::size_t>(k)]) << '\n';
    }
    return csv.str();
}

inline std::string posterior_csv() {
    const auto primary = quadrature::parity_rule();
    const auto check = quadrature::parity_rule(true);
    const Grid<2> g = recovery::grid();
    const backends::CpuBackend backend{};
    const auto jt = engine::jeffreys_table(backend, Objective{}, primary, kN, g);
    std::ostringstream csv;
    csv << "# Grid-Bayesian reference results (engine::grid_posterior). Written by posterior_reference --write.\n"
           "# Panels: recovery DGP under the validation seed M3BAYVAL (not the recovery panels), T = 20, n = 100.\n"
           "# The Jeffreys prior uses the parity grid's table for n = 100. Ends on the natural scale.\n"
           "pd,rho,replicate,prior,counts,flags,refinements,pd_lo_bound_hex,pd_hi_bound_hex,pd_points,rho_lo_bound_hex,"
           "rho_hi_bound_hex,rho_points,pd_et_lo_hex,pd_et_hi_hex,pd_hpd_lo_hex,pd_hpd_hi_hex,rho_et_lo_hex,"
           "rho_et_hi_hex,rho_hpd_lo_hex,rho_hpd_hi_hex\n";
    std::uint32_t scenario = 0;
    for (const double pd : {0.01, 0.05}) {
        for (const double rho : {0.12, 0.24}) {
            ++scenario;
            for (std::uint32_t r = 0; r < 2; ++r) {
                std::vector<std::int64_t> n(static_cast<std::size_t>(kT), kN), d(n.size());
                const dgp::PanelSpec spec{kBayesSeed, scenario, r, pd, rho, n.data(), kT};
                if (dgp::simulate_panel(spec, d.data(), nullptr) != dgp::Status::Ok) std::abort();
                std::vector<Objective::Obs> obs(n.size());
                for (std::size_t t = 0; t < obs.size(); ++t) obs[t] = {kN, d[t]};
                std::vector<double> L;
                engine::Estimate2 est{};
                if (engine::calibrate(backend, Objective{}, primary, check, obs.data(), kT, g, L, est) != engine::Status::Ok) {
                    std::abort();
                }
                double se_u[2];
                for (int a = 0; a < 2; ++a) {
                    se_u[a] = est.se[a] / grid::dvalue_dscaled(g.axis[a].scale, grid::to_scaled(g.axis[a].scale, est.value[a]));
                }
                for (const auto prior : {engine::Prior::Flat, engine::Prior::Jeffreys}) {
                    const auto p = engine::grid_posterior(backend, Objective{}, primary, obs.data(), kT, g, L, prior, &jt, se_u);
                    csv << pd << ',' << rho << ',' << r << ',' << (prior == engine::Prior::Flat ? "flat" : "jeffreys") << ',';
                    for (std::size_t t = 0; t < d.size(); ++t) csv << (t ? " " : "") << d[t];
                    csv << ',' << p.flags << ',' << p.refinements;
                    for (int a = 0; a < 2; ++a) {
                        csv << ',' << test::to_hex(p.grid.axis[a].lo) << ',' << test::to_hex(p.grid.axis[a].hi) << ','
                            << p.grid.axis[a].n;
                    }
                    for (int a = 0; a < 2; ++a) {
                        csv << ',' << test::to_hex(p.et_lo[a]) << ',' << test::to_hex(p.et_hi[a]) << ','
                            << test::to_hex(p.hpd_lo[a]) << ',' << test::to_hex(p.hpd_hi[a]);
                    }
                    csv << '\n';
                }
            }
        }
    }
    return csv.str();
}

}  // namespace vcal::postref
