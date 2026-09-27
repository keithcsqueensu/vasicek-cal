// SPDX-License-Identifier: Apache-2.0
//
// SHA-256 (FIPS 180-4) for test-side content hashes, e.g. the DGP's cross-platform
// reference-panel check (D-057). Verified against Python hashlib vectors
// (tests/golden/dgp/sha256_vectors.csv).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace vcal::test {

std::string sha256_hex(const std::vector<std::uint8_t>& message);

}  // namespace vcal::test
