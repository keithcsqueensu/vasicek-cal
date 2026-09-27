// SPDX-License-Identifier: Apache-2.0
#include "abi/build_info.hpp"

#include <cstdio>
#include <string>

#include "vcal/vcal.h"
#include "vcal_git_info.h"  // generated at build time: VCAL_GIT_COMMIT, VCAL_GIT_DIRTY

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifndef VCAL_VERSION
#error "VCAL_VERSION must be defined by the build"
#endif
#ifndef VCAL_BUILD_TYPE
#error "VCAL_BUILD_TYPE must be defined by the build"
#endif

namespace vcal::abi {

namespace {

std::string compiler() {
    char b[128];
#if defined(__clang__)
    std::snprintf(b, sizeof b, "Clang %s", __clang_version__);
#elif defined(__GNUC__)
    std::snprintf(b, sizeof b, "GCC %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif defined(_MSC_FULL_VER)
    std::snprintf(b, sizeof b, "MSVC %d.%d.%d", _MSC_FULL_VER / 10000000, (_MSC_FULL_VER / 100000) % 100,
                  _MSC_FULL_VER % 100000);
#else
    std::snprintf(b, sizeof b, "unknown");
#endif
    return b;
}

std::string make() {
    std::string s;
    s += "vcal_version=" VCAL_VERSION "\n";
    s += "abi_version=" + std::to_string(VCAL_ABI_VERSION_MAJOR) + "." + std::to_string(VCAL_ABI_VERSION_MINOR) + "\n";
    s += "git_commit=" VCAL_GIT_COMMIT "\n";
    s += "git_dirty=" VCAL_GIT_DIRTY "\n";
    s += "compiler=" + compiler() + "\n";
    s += "cxx_standard=" + std::to_string(__cplusplus) + "\n";
    s += "build_type=" VCAL_BUILD_TYPE "\n";
#if defined(_OPENMP)
    s += "openmp=" + std::to_string(_OPENMP) + "\n";
#else
    s += "openmp=off\n";
#endif
    s += "cuda=off\n";
    s += "backends=cpu\n";
    s += "profiles=parity\n";
    s += "fp_contract=off\n";
    return s;
}

}  // namespace

const char* build_info() {
    static const std::string s = make();
    return s.c_str();
}

}  // namespace vcal::abi
