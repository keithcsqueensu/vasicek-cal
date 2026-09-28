// SPDX-License-Identifier: Apache-2.0
//
// Reference results of the method of moments (engine/moments.hpp; M3) for its scipy replication and
// its regression test. Panels come from the recovery DGP under a validation seed of their own
// (kMomSeed, not the recovery seed), so that nothing here computes study S-8's result on the recovery
// panels before S-8's own run (D-161): for each (PD, rho, T) of the recovery matrix, replicates 0 and
// 1 with n = 1,000 (count panels), and the same factor draws mapped to rates (rate series).
//
// One row per panel: the panel, and the engine's PD-hat, PD_2-hat, rho-hat and flags. The box is the
// recovery grid's. unit_moments recomputes every row (a cross-platform regression);
// validation/scipy/method_of_moments.py replicates it with scipy.
//
//   mom_reference [--write FILE]      without --write, print the rows
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "tests/recovery/mom_panels.hpp"

int main(int argc, char** argv) {
    std::string write;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--write") == 0 && i + 1 < argc) {
            write = argv[++i];
        } else {
            std::fprintf(stderr, "usage: mom_reference [--write FILE]\n");
            return 2;
        }
    }
    const std::string csv = vcal::momref::reference_csv();
    if (write.empty()) {
        std::fputs(csv.c_str(), stdout);
        return 0;
    }
    std::ofstream f(write, std::ios::binary);
    f << csv;
    if (!f) return 1;
    std::printf("wrote %s\n", write.c_str());
    return 0;
}
