// SPDX-License-Identifier: Apache-2.0
//
// Reader for the golden CSV tables under tests/golden/ (D-033). Lines starting with '#'
// are comments; the first other line is the header. Doubles are read from the *_hex
// columns, which are exact; decimal columns are for human readers only.
// Every function throws std::runtime_error on malformed input, which the runner reports
// as a test failure.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace vcal::test {

struct CsvTable {
    std::string path;
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;

    std::size_t column(const std::string& name) const;
};

// `relative_path` is relative to tests/golden/.
CsvTable read_golden_csv(const std::string& relative_path);
std::string golden_path(const std::string& relative_path);  // for non-CSV golden files

double parse_double(const std::string& text);  // hex float, "inf", "-inf", "nan"
std::string to_hex(double v);                   // exact hex float (printf %a), parseable by parse_double
std::int64_t parse_int(const std::string& text);

}  // namespace vcal::test
