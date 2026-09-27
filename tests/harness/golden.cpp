// SPDX-License-Identifier: Apache-2.0
#include "tests/harness/golden.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifndef VCAL_GOLDEN_DIR
#error "VCAL_GOLDEN_DIR must be defined by the build"
#endif

namespace vcal::test {
namespace {

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> cells;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) cells.push_back(cell);
    return cells;
}

}  // namespace

std::size_t CsvTable::column(const std::string& name) const {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) return i;
    }
    throw std::runtime_error(path + ": no column '" + name + "'");
}

std::string golden_path(const std::string& relative_path) { return std::string(VCAL_GOLDEN_DIR) + "/" + relative_path; }

CsvTable read_golden_csv(const std::string& relative_path) {
    CsvTable table;
    table.path = std::string(VCAL_GOLDEN_DIR) + "/" + relative_path;
    std::ifstream in(table.path);
    if (!in) throw std::runtime_error("cannot open " + table.path);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        if (table.header.empty()) {
            table.header = split(line);
            continue;
        }
        auto cells = split(line);
        if (cells.size() != table.header.size()) {
            throw std::runtime_error(table.path + ": row has " + std::to_string(cells.size()) +
                                     " cells, header has " + std::to_string(table.header.size()));
        }
        table.rows.push_back(std::move(cells));
    }
    if (table.rows.empty()) throw std::runtime_error(table.path + ": no data rows");
    return table;
}

double parse_double(const std::string& text) {
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(text.c_str(), &end);
    // ERANGE is expected for exact subnormal hex values; only a partial parse is an error.
    if (end == text.c_str() || *end != '\0') throw std::runtime_error("bad double: '" + text + "'");
    return v;
}

std::string to_hex(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%a", v);
    return buf;
}

std::int64_t parse_int(const std::string& text) {
    char* end = nullptr;
    errno = 0;
    const long long v = std::strtoll(text.c_str(), &end, 10);
    if (end == text.c_str() || *end != '\0' || errno == ERANGE) {
        throw std::runtime_error("bad integer: '" + text + "'");
    }
    return static_cast<std::int64_t>(v);
}

}  // namespace vcal::test
