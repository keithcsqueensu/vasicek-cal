// SPDX-License-Identifier: Apache-2.0
//
// Reference results of the grid-Bayesian estimator (engine/posterior.hpp; M3) for its scipy
// replication (validation/scipy/grid_posterior.py) and its regression test (unit_posterior).
// Panels come from the recovery DGP under a validation seed of their own (kBayesSeed), not the
// recovery seed, so nothing here computes study S-9's or S-15's result before their registered runs
// (D-161).
//
// Two files:
//   jeffreys.csv   the Jeffreys table (log sqrt det I_u) for n = 100 on a 13 x 9 logit grid over the box;
//   posterior.csv  for (PD, rho) in {1%, 5%} x {0.12, 0.24}, T = 20, n = 100, replicates 0 and 1, both
//                  priors: the panel, and the estimator's flags, refinements, final grid, and
//                  equal-tailed and HPD ends (natural scale) for PD and rho. The Jeffreys prior of
//                  these fits uses the parity grid's table for n = 100.
//
//   posterior_reference --write DIR
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>

#include "tests/recovery/posterior_panels.hpp"

int main(int argc, char** argv) {
    if (argc != 3 || std::strcmp(argv[1], "--write") != 0) {
        std::fprintf(stderr, "usage: posterior_reference --write DIR\n");
        return 2;
    }
    const std::string dir = argv[2];
    for (const auto& [name, text] : {std::pair<std::string, std::string>{"jeffreys.csv", vcal::postref::jeffreys_csv()},
                                     {"posterior.csv", vcal::postref::posterior_csv()}}) {
        std::ofstream f(dir + "/" + name, std::ios::binary);
        f << text;
        if (!f) return 1;
        std::printf("wrote %s/%s\n", dir.c_str(), name.c_str());
    }
    return 0;
}
