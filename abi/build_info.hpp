// SPDX-License-Identifier: Apache-2.0
//
// Build provenance for vcal_build_info (D-142): "key=value" lines, fixed for the life of the
// library. The git commit is read at build time (cmake/GitInfo.cmake), the compiler and language
// settings from the compiler's own macros.
#pragma once

namespace vcal::abi {

const char* build_info();

}  // namespace vcal::abi
