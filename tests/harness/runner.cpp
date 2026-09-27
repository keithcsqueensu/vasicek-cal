// SPDX-License-Identifier: Apache-2.0
#include "tests/harness/vcal_test.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <exception>
#include <string_view>
#include <vector>

namespace vcal::test {
namespace {

struct Case {
    const char* name;
    TestFn fn;
    const char* file;
    int line;
};

std::vector<Case>& registry() {
    static std::vector<Case> cases;
    return cases;
}

int g_failed_checks = 0;

bool run_case(const Case& c) {
    const int before = g_failed_checks;
    try {
        c.fn();
    } catch (const RequireFailed&) {
        // already reported
    } catch (const std::exception& e) {
        report_failure(c.file, c.line, std::string("unhandled exception: ") + e.what());
    } catch (...) {
        report_failure(c.file, c.line, "unhandled non-standard exception");
    }
    const bool passed = g_failed_checks == before;
    std::printf("[ %s ] %s\n", passed ? "PASS" : "FAIL", c.name);
    std::fflush(stdout);
    return passed;
}

}  // namespace

Registrar::Registrar(const char* name, TestFn fn, const char* file, int line) {
    registry().push_back({name, fn, file, line});
}

void report_failure(const char* file, int line, const std::string& message) {
    ++g_failed_checks;
    std::fflush(stdout);
    std::fprintf(stderr, "%s:%d: FAILED: %s\n", file, line, message.c_str());
    std::fflush(stderr);
}

bool near_abs(double actual, double expected, double abs_tol) {
    if (std::isnan(actual) || std::isnan(expected)) return false;
    if (actual == expected) return true;  // also covers equal infinities
    return std::fabs(actual - expected) <= abs_tol;
}

bool near_rel(double actual, double expected, double rel_tol) {
    if (std::isnan(actual) || std::isnan(expected)) return false;
    if (actual == expected) return true;
    return std::fabs(actual - expected) <= rel_tol * std::fabs(expected);
}

namespace {
// Maps the bit pattern to an integer that is monotone in the double's value,
// with -0 and +0 both at 0.
std::int64_t ordered_bits(double x) {
    std::int64_t i = 0;
    std::memcpy(&i, &x, sizeof i);
    return i < 0 ? std::numeric_limits<std::int64_t>::min() - i : i;
}
}  // namespace

std::uint64_t ulp_distance(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return std::numeric_limits<std::uint64_t>::max();
    const std::int64_t ia = ordered_bits(a);
    const std::int64_t ib = ordered_bits(b);
    return ia >= ib ? static_cast<std::uint64_t>(ia) - static_cast<std::uint64_t>(ib)
                    : static_cast<std::uint64_t>(ib) - static_cast<std::uint64_t>(ia);
}

void note(const std::string& message) { std::printf("  note: %s\n", message.c_str()); }

std::string env_or_empty(const char* name) {
#if defined(_MSC_VER)
    char* value = nullptr;
    std::size_t length = 0;
    if (_dupenv_s(&value, &length, name) != 0 || value == nullptr) return {};
    std::string out(value);
    std::free(value);
    return out;
#else
    const char* value = std::getenv(name);
    return value ? std::string(value) : std::string();
#endif
}

void UlpStats::add(double input, double actual, double expected) {
    ++count;
    const std::uint64_t d = ulp_distance(actual, expected);
    if (d > max_ulps || count == 1) {
        max_ulps = d;
        worst_input = input;
        worst_actual = actual;
        worst_expected = expected;
    }
}

std::string UlpStats::summary(const std::string& label) const {
    return label + ": max " + std::to_string(max_ulps) + " ulp over " + std::to_string(count) +
           " points (input " + describe(worst_input) + ", got " + describe(worst_actual) +
           ", expected " + describe(worst_expected) + ")";
}

}  // namespace vcal::test

int main(int argc, char** argv) {
    using vcal::test::Case;
    auto& cases = vcal::test::registry();

    std::sort(cases.begin(), cases.end(), [](const Case& a, const Case& b) {
        return std::strcmp(a.name, b.name) < 0;
    });
    for (size_t i = 1; i < cases.size(); ++i) {
        if (std::strcmp(cases[i - 1].name, cases[i].name) == 0) {
            std::fprintf(stderr, "duplicate test name: %s\n", cases[i].name);
            return 2;
        }
    }
    if (cases.empty()) {
        std::fprintf(stderr, "no tests registered\n");
        return 2;
    }

    if (argc == 2 && std::string_view(argv[1]) == "--list") {
        for (const Case& c : cases) std::printf("%s\n", c.name);
        return 0;
    }

    std::vector<const Case*> selected;
    if (argc == 1) {
        for (const Case& c : cases) selected.push_back(&c);
    } else {
        for (int i = 1; i < argc; ++i) {
            auto it = std::find_if(cases.begin(), cases.end(), [&](const Case& c) {
                return std::string_view(c.name) == argv[i];
            });
            if (it == cases.end()) {
                std::fprintf(stderr, "unknown test: %s\n", argv[i]);
                return 2;
            }
            selected.push_back(&*it);
        }
    }

    size_t failed = 0;
    for (const Case* c : selected) {
        if (!vcal::test::run_case(*c)) ++failed;
    }
    std::printf("%zu of %zu tests failed (%d failed checks)\n", failed, selected.size(),
                vcal::test::g_failed_checks);
    return failed == 0 ? 0 : 1;
}
