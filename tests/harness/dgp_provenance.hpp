// SPDX-License-Identifier: Apache-2.0
//
// Which Bernoulli kernel the DGP uses in this run (P-17, D-194). Every kernel gives the same panels bit for
// bit; the tools record which one ran, so a run's log and timing sample say how it was generated. The
// kernel is the fastest available unless VCAL_DGP_KERNEL (scalar, sse2 or avx2) names one.
#pragma once

#include <cstdio>
#include <string>

#include "dgp/dgp.hpp"

namespace vcal::test {

// "scalar", "sse2", "avx2", or "unavailable".
inline std::string dgp_kernel() {
    return vcal::dgp::bernoulli_kernel_name(vcal::dgp::bernoulli_kernel_selected());
}

// Prints "<tool>: dgp_kernel=<name>" to stderr. False, after saying why, if VCAL_DGP_KERNEL names a kernel
// this CPU or build lacks: the tool should then exit rather than let every simulation fail.
inline bool announce_dgp_kernel(const char* tool) {
    const std::string k = dgp_kernel();
    std::fprintf(stderr, "%s: dgp_kernel=%s\n", tool, k.c_str());
    if (vcal::dgp::bernoulli_kernel_available(vcal::dgp::bernoulli_kernel_selected())) return true;
    std::fprintf(stderr, "%s: VCAL_DGP_KERNEL names a kernel this CPU or build does not have (scalar, sse2 or avx2)\n",
                 tool);
    return false;
}

}  // namespace vcal::test
