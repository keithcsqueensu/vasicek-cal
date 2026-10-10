// SPDX-License-Identifier: Apache-2.0
// dgp/ (M1.6; D-052..D-057, D-061, D-109, D-110): known answers, bitwise identity with the
// Python mirror on every platform, accuracy against mpmath, stream addressing, and moments.
#include "dgp/dgp.hpp"

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "tests/harness/golden.hpp"
#include "tests/harness/sha256.hpp"
#include "tests/harness/vcal_test.hpp"
#include "tests/tolerances.hpp"

namespace {

namespace dgp = vcal::dgp;
namespace tol = vcal::tol;
using vcal::test::describe;
using vcal::test::parse_double;
using vcal::test::parse_int;
using vcal::test::UlpStats;

bool bits_equal(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

std::uint32_t hex32(const std::string& s) { return static_cast<std::uint32_t>(std::stoul(s, nullptr, 16)); }

std::vector<std::uint8_t> hex_bytes(const std::string& s) {
    std::vector<std::uint8_t> out;
    for (std::size_t i = 0; i + 1 < s.size(); i += 2) out.push_back(static_cast<std::uint8_t>(std::stoul(s.substr(i, 2), nullptr, 16)));
    return out;
}

void append_le(std::vector<std::uint8_t>& out, std::uint64_t v) {
    for (int i = 0; i < 8; ++i) out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}

constexpr dgp::BernoulliKernel kAllKernels[] = {dgp::BernoulliKernel::Scalar, dgp::BernoulliKernel::Sse2,
                                                dgp::BernoulliKernel::Avx2};

// The selection as the process started, read during static initialisation, before any test changes it.
const dgp::BernoulliKernel kInitialSelection = dgp::bernoulli_kernel_selected();

// Runs fn(kernel) with each Bernoulli kernel this CPU and build have selected in turn (P-17, D-194), so every
// code path runs on one machine, not only the one CPUID picks; the selection is restored afterwards. On an
// x86-64 CPU with AVX2 (every CI runner) that is all three. A kernel this machine lacks is noted, not failed.
template <class F>
void for_each_kernel(F fn) {
    struct Restore {
        dgp::BernoulliKernel saved = dgp::bernoulli_kernel_selected();
        ~Restore() { dgp::select_bernoulli_kernel(saved); }
    } restore;
    for (const auto kernel : kAllKernels) {
        if (!dgp::select_bernoulli_kernel(kernel)) {
            vcal::test::note(std::string("Bernoulli kernel ") + dgp::bernoulli_kernel_name(kernel) + " not available here");
            continue;
        }
        VCAL_REQUIRE(dgp::bernoulli_kernel_selected() == kernel);
        fn(kernel);
    }
}

}  // namespace

// --- Philox and the building blocks ------------------------------------------------------------

VCAL_TEST(philox4x32_matches_random123_known_answers) {
    std::ifstream in(vcal::test::golden_path("dgp/random123_philox4x32_kat.txt"));
    VCAL_REQUIRE(static_cast<bool>(in));
    std::string line;
    int checked = 0;
    while (std::getline(in, line)) {
        if (line.rfind("philox4x32 ", 0) != 0) continue;
        std::istringstream ss(line);
        std::string name, rounds, f[10];
        ss >> name >> rounds;
        for (auto& x : f) ss >> x;
        const auto out = dgp::philox4x32({hex32(f[0]), hex32(f[1]), hex32(f[2]), hex32(f[3])}, {hex32(f[4]), hex32(f[5])},
                                         std::stoi(rounds));
        VCAL_CHECK(out[0] == hex32(f[6]) && out[1] == hex32(f[7]) && out[2] == hex32(f[8]) && out[3] == hex32(f[9]));
        ++checked;
    }
    VCAL_CHECK_EQ(checked, 6);
}

// C++ and the Python mirror must agree bit for bit, function by function (D-057).
VCAL_TEST(det_math_matches_the_python_mirror_bitwise) {
    const auto t = vcal::test::read_golden_csv("dgp/mirror_values.csv");
    int checked = 0;
    for (const auto& r : t.rows) {
        const std::string& fn = r[t.column("function")];
        const double want = parse_double(r[t.column("value_hex")]);
        double got = 0.0;
        if (fn == "uniform") {
            const std::string w = r[t.column("x_hex")];
            got = dgp::uniform_from_words(hex32(w.substr(0, 8)), hex32(w.substr(8, 8)));
        } else {
            const double x = parse_double(r[t.column("x_hex")]);
            got = fn == "log" ? dgp::det_log(x) : fn == "exp" ? dgp::det_exp(x) : fn == "ncdf" ? dgp::det_ncdf(x) : dgp::det_probit(x);
        }
        if (!bits_equal(got, want)) {
            vcal::test::report_failure(__FILE__, __LINE__, fn + "(" + r[t.column("x")] + "): C++ " + describe(got) +
                                                               ", mirror " + describe(want));
        }
        ++checked;
    }
    VCAL_CHECK(checked > 40);
}

VCAL_TEST(uniform_is_strictly_inside_the_unit_interval) {
    VCAL_CHECK_EQ(dgp::uniform_from_words(0u, 0u), 0x1p-53);
    VCAL_CHECK_EQ(dgp::uniform_from_words(0xFFFFFFFFu, 0xFFFFFFFFu), 1.0 - 0x1p-53);
}

VCAL_TEST(fp_environment_keeps_subnormals) { VCAL_CHECK(dgp::fp_environment_ok()); }

// --- the Bernoulli kernels (P-17, D-194) -------------------------------------------------------


// The threshold K(p) at its edges (docs/methodology/dgp.md, "An equivalent integer form of step 3").
// U(k) = (k + 1/2) 2^-52 is key k's uniform and K(p) counts the keys with U(k) < p, so K steps from k
// to k + 1 just above p = U(k): K(U(k)) = k and K(the next double above U(k)) = k + 1. Confusing "<"
// with "<=" anywhere would move every one of these by one.
VCAL_TEST(bernoulli_threshold_at_its_edges) {
    const std::uint64_t kAll = std::uint64_t{1} << 52;
    const auto K = [](double p) { return dgp::bernoulli_threshold(p); };
    // Nothing defaults at or below the smallest uniform, U(0) = 2^-53.
    for (const double p : {0.0, -0.0, std::numeric_limits<double>::quiet_NaN(), DBL_TRUE_MIN, DBL_MIN, 1e-300, 0x1p-54,
                           0x1p-53}) {
        VCAL_CHECK_EQ(K(p), std::uint64_t{0});
    }
    // The smallest p that changes K: the next double above U(0).
    VCAL_CHECK_EQ(K(std::nextafter(0x1p-53, 1.0)), std::uint64_t{1});
    // The largest p that changes K is the next double above U(2^52 - 1) = 1 - 2^-53, which is 1.
    VCAL_CHECK_EQ(K(1.0 - 0x1p-53), kAll - 1);
    VCAL_CHECK_EQ(std::nextafter(1.0 - 0x1p-53, 2.0), 1.0);
    for (const double p : {1.0, std::nextafter(1.0, 2.0), 2.0, std::numeric_limits<double>::infinity()}) {
        VCAL_CHECK_EQ(K(p), kAll);
    }
    // Steps: chosen keys (the ends, and either side of the 2^26 boundary between the two words) and
    // pseudo-random ones. For each, p exactly on U(k) and one ulp either side, and the comparison itself
    // on the words that give k (with junk in the 6 discarded low bits): U < p iff k < K(p).
    std::vector<std::uint64_t> keys = {0, 1, 2, 3, (1u << 26) - 1, 1u << 26, (1u << 26) + 1, (kAll >> 1) - 1, kAll >> 1,
                                       kAll - 2, kAll - 1};
    std::uint64_t x = 0x9E3779B97F4A7C15ull;
    for (int i = 0; i < 20000; ++i) {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        keys.push_back(x >> 12);  // uniform over [0, 2^52)
    }
    for (const std::uint64_t k : keys) {
        const double u = (static_cast<double>(k) + 0.5) * 0x1p-52;  // exact: k + 1/2 has at most 53 bits
        const std::uint32_t junk = static_cast<std::uint32_t>(k * 0x2545F491u) & 63u;
        const auto a = static_cast<std::uint32_t>((k >> 26) << 6) | junk;
        const auto b = static_cast<std::uint32_t>((k & ((1u << 26) - 1)) << 6) | (63u - junk);
        VCAL_CHECK(bits_equal(dgp::uniform_from_words(a, b), u));
        const double below = std::nextafter(u, 0.0);
        const double above = std::nextafter(u, 1.0);
        VCAL_CHECK_EQ(K(u), k);
        VCAL_CHECK_EQ(K(below), k);
        VCAL_CHECK_EQ(K(above), k + 1);
        for (const double p : {below, u, above}) VCAL_CHECK_EQ(u < p, k < K(p));
    }
    vcal::test::note(std::to_string(keys.size()) + " threshold steps checked");
}

// Every vector kernel this CPU and build have returns the scalar kernel's count, which is the specification's
// loop. The cases are chosen to be hard for the vector kernels: every n from 0 to 80 (each tail length
// after the 16-block vector steps, and odd n), p at the threshold's edges (0, 1, the smallest and largest
// p that change it, beyond either end), and p exactly on obligors' uniforms, which are threshold steps,
// and one ulp either side.
VCAL_TEST(bernoulli_kernels_equal_the_specification_bit_for_bit) {
    vcal::test::note(std::string("selected Bernoulli kernel ") + dgp::bernoulli_kernel_name(dgp::bernoulli_kernel_selected()));
    std::vector<std::int64_t> ns;
    for (std::int64_t n = 0; n <= 80; ++n) ns.push_back(n);
    ns.insert(ns.end(), {255, 256, 257, 1000, 10001, (std::int64_t{1} << 20) + 3});
    const std::uint64_t seed = 0xC0FFEE1234567890ull;
    long long compared = 0;
    for (const auto kernel : kAllKernels) {
        if (kernel == dgp::BernoulliKernel::Scalar) continue;  // the reference
        if (!dgp::bernoulli_kernel_available(kernel)) {
            vcal::test::note(std::string("Bernoulli kernel ") + dgp::bernoulli_kernel_name(kernel) + " not available here");
            continue;
        }
        for (std::size_t c = 0; c < ns.size(); ++c) {
            const std::int64_t n = ns[c];
            const auto scenario = static_cast<std::uint32_t>(c);
            const std::uint32_t replicate = 0xFFFFFFF0u + static_cast<std::uint32_t>(c % 16);
            const std::uint32_t period = 7u * static_cast<std::uint32_t>(c);
            std::vector<double> ps = {0.0,
                                      -0.0,
                                      1.0,
                                      2.0,
                                      std::numeric_limits<double>::infinity(),
                                      std::numeric_limits<double>::quiet_NaN(),
                                      DBL_TRUE_MIN,
                                      1e-300,
                                      0x1p-53,
                                      std::nextafter(0x1p-53, 1.0),
                                      0x1p-52,
                                      0x3p-53,
                                      1e-4,
                                      0.02,
                                      0.5,
                                      0.97,
                                      std::nextafter(1.0 - 0x1p-53, 0.0),
                                      1.0 - 0x1p-53};
            // p on the uniforms of the period's first, middle and last obligors, and one ulp either side.
            const dgp::PhiloxKey key{static_cast<std::uint32_t>(seed), static_cast<std::uint32_t>(seed >> 32)};
            for (std::int64_t i : {std::int64_t{0}, n / 2, n - 1}) {
                if (i < 0 || i >= n) continue;
                const auto block = static_cast<std::uint32_t>(dgp::kFactorBlocks + static_cast<std::uint64_t>(i / 2));
                const auto w = dgp::philox4x32({scenario, replicate, period, block}, key);
                const double u = i % 2 == 0 ? dgp::uniform_from_words(w[0], w[1]) : dgp::uniform_from_words(w[2], w[3]);
                ps.insert(ps.end(), {u, std::nextafter(u, 0.0), std::nextafter(u, 1.0)});
            }
            if (n > (std::int64_t{1} << 16)) ps = {0.02, 0.5};  // the large n: two values keep it quick
            for (const double p : ps) {
                const std::int64_t want =
                    dgp::count_defaults(seed, scenario, replicate, period, n, p, dgp::BernoulliKernel::Scalar);
                const std::int64_t got = dgp::count_defaults(seed, scenario, replicate, period, n, p, kernel);
                ++compared;
                if (got != want) {
                    vcal::test::report_failure(__FILE__, __LINE__,
                                               std::string("kernel ") + dgp::bernoulli_kernel_name(kernel) + ", n " +
                                                   std::to_string(n) + ", p " + describe(p) + ": " +
                                                   std::to_string(got) + " against " + std::to_string(want));
                }
            }
        }
    }
    vcal::test::note(std::to_string(compared) + " counts compared");
    // The specification's own boundary behaviour, which the comparisons above inherit.
    VCAL_CHECK_EQ(dgp::count_defaults(seed, 1, 2, 3, 1001, 1.0, dgp::BernoulliKernel::Scalar), 1001);
    VCAL_CHECK_EQ(dgp::count_defaults(seed, 1, 2, 3, 1001, 0.0, dgp::BernoulliKernel::Scalar), 0);
}

// The selection at start-up follows VCAL_DGP_KERNEL (CTest runs this executable under each value, and
// under an invalid one); unset, it is the fastest available. An unavailable selection makes the
// simulations refuse rather than fall back to another kernel.
VCAL_TEST(bernoulli_kernel_selection_follows_vcal_dgp_kernel) {
    const std::string env = vcal::test::env_or_empty("VCAL_DGP_KERNEL");
    dgp::BernoulliKernel want = dgp::BernoulliKernel::Scalar;
    if (env.empty()) {
        if (dgp::bernoulli_kernel_available(dgp::BernoulliKernel::Sse2)) want = dgp::BernoulliKernel::Sse2;
        if (dgp::bernoulli_kernel_available(dgp::BernoulliKernel::Avx2)) want = dgp::BernoulliKernel::Avx2;
    } else {
        want = dgp::BernoulliKernel::Unavailable;
        for (const auto k : kAllKernels) {
            if (env == dgp::bernoulli_kernel_name(k) && dgp::bernoulli_kernel_available(k)) want = k;
        }
    }
    vcal::test::note("VCAL_DGP_KERNEL '" + env + "': selected at start-up " +
                     dgp::bernoulli_kernel_name(kInitialSelection));
    VCAL_CHECK(kInitialSelection == want);
#if defined(__x86_64__) || defined(_M_X64)
    VCAL_CHECK(dgp::bernoulli_kernel_available(dgp::BernoulliKernel::Sse2));  // SSE2 is part of x86-64
#endif

    struct Restore {
        dgp::BernoulliKernel saved = dgp::bernoulli_kernel_selected();
        ~Restore() { dgp::select_bernoulli_kernel(saved); }
    } restore;
    VCAL_REQUIRE(dgp::select_bernoulli_kernel(dgp::BernoulliKernel::Unavailable));
    const std::int64_t n[] = {100, 100};
    std::int64_t d[2];
    VCAL_CHECK(dgp::simulate_panel({1, 0, 0, 0.02, 0.1, n, 2}, d, nullptr) == dgp::Status::KernelUnavailable);
    VCAL_CHECK(dgp::simulate_panel_ar1({1, 0, 0, 0.02, 0.1, n, 2}, 0.5, d, nullptr) == dgp::Status::KernelUnavailable);
}

VCAL_TEST(bernoulli_kernels_refuse_bad_input) {
    VCAL_CHECK(dgp::bernoulli_kernel_available(dgp::BernoulliKernel::Scalar));
    VCAL_CHECK(!dgp::bernoulli_kernel_available(dgp::BernoulliKernel::Unavailable));
    VCAL_CHECK(!dgp::select_bernoulli_kernel(static_cast<dgp::BernoulliKernel>(7)));
    VCAL_CHECK_EQ(dgp::count_defaults(1, 0, 0, 0, -1, 0.5, dgp::BernoulliKernel::Scalar), std::int64_t{-1});
    VCAL_CHECK_EQ(dgp::count_defaults(1, 0, 0, 0, dgp::kMaxObligors + 1, 0.5, dgp::BernoulliKernel::Scalar),
                  std::int64_t{-1});
    VCAL_CHECK_EQ(dgp::count_defaults(1, 0, 0, 0, 10, 0.5, dgp::BernoulliKernel::Unavailable), std::int64_t{-1});
    VCAL_CHECK_EQ(dgp::count_defaults(1, 0, 0, 0, 10, 0.5, static_cast<dgp::BernoulliKernel>(7)), std::int64_t{-1});
}

// --- accuracy against mpmath (target 1e-12 relative, D-052) ------------------------------------

VCAL_TEST(det_log_exp_ncdf_probit_accuracy) {
    struct Case {
        const char* file;
        const char* xcol;
        double (*fn)(double);
        std::uint64_t limit;
        const char* label;
    };
    const Case cases[] = {
        {"dgp/log.csv", "x_hex", dgp::det_log, tol::TOL_DGP_LOG_ULP, "det_log"},
        {"dgp/exp.csv", "x_hex", dgp::det_exp, tol::TOL_DGP_EXP_ULP, "det_exp (normal results)"},
        {"dgp/ncdf.csv", "x_hex", dgp::det_ncdf, tol::TOL_DGP_NCDF_ULP, "det_ncdf (normal results)"},
        {"special/probit.csv", "p_hex", dgp::det_probit, tol::TOL_DGP_PROBIT_ULP, "det_probit"},
    };
    for (const auto& c : cases) {
        const auto t = vcal::test::read_golden_csv(c.file);
        UlpStats s;
        double subnormal_abs = 0.0;
        for (const auto& r : t.rows) {
            const double x = parse_double(r[t.column(c.xcol)]);
            const double e = parse_double(r[t.column("expected_hex")]);
            const double a = c.fn(x);
            if (e != 0.0 && std::fabs(e) < DBL_MIN) {
                subnormal_abs = std::fmax(subnormal_abs, std::fabs(a - e));
            } else {
                s.add(x, a, e);
            }
        }
        VCAL_CHECK_ULP_STATS(s, c.limit, c.label);
        VCAL_CHECK(subnormal_abs <= DBL_MIN);
    }
}

// --- panels ------------------------------------------------------------------------------------

// Q14's goal, tested directly: every platform reproduces the mirror's reference panels bit for
// bit, and the SHA-256 of their canonical bytes equals the single checked-in value.
VCAL_TEST(reference_panels_are_bitwise_identical_to_the_mirror) {
    const auto t = vcal::test::read_golden_csv("dgp/reference_panels.csv");
    struct Rows {
        std::uint64_t seed;
        std::uint32_t scenario, replicate;
        double pd, rho;
        std::vector<std::int64_t> n, d;
        std::vector<double> z;
    };
    std::vector<Rows> panels;
    for (const auto& r : t.rows) {
        const auto id = static_cast<std::size_t>(parse_int(r[t.column("panel")]));
        if (panels.size() <= id) panels.resize(id + 1);
        auto& p = panels[id];
        p.seed = std::stoull(r[t.column("seed_hex")], nullptr, 16);
        p.scenario = static_cast<std::uint32_t>(parse_int(r[t.column("scenario")]));
        p.replicate = static_cast<std::uint32_t>(parse_int(r[t.column("replicate")]));
        p.pd = parse_double(r[t.column("pd_hex")]);
        p.rho = parse_double(r[t.column("rho_hex")]);
        p.n.push_back(parse_int(r[t.column("n")]));
        p.d.push_back(parse_int(r[t.column("d")]));
        p.z.push_back(parse_double(r[t.column("z_hex")]));
    }
    std::ifstream in(vcal::test::golden_path("dgp/reference_panels.sha256"));
    std::string want;
    in >> want;
    for_each_kernel([&](dgp::BernoulliKernel kernel) {
        std::vector<std::uint8_t> blob;
        for (const auto& p : panels) {
            const auto T = static_cast<std::int64_t>(p.n.size());
            std::vector<std::int64_t> d(p.n.size());
            std::vector<double> z(p.n.size());
            const dgp::PanelSpec spec{p.seed, p.scenario, p.replicate, p.pd, p.rho, p.n.data(), T};
            VCAL_REQUIRE(dgp::simulate_panel(spec, d.data(), z.data()) == dgp::Status::Ok);
            for (std::size_t i = 0; i < p.n.size(); ++i) {
                VCAL_CHECK_EQ(d[i], p.d[i]);
                VCAL_CHECK(bits_equal(z[i], p.z[i]));
                std::uint64_t zb = 0;
                std::memcpy(&zb, &z[i], sizeof zb);
                append_le(blob, static_cast<std::uint64_t>(d[i]));
                append_le(blob, zb);
            }
        }
        const std::string got = vcal::test::sha256_hex(blob);
        vcal::test::note(std::string(dgp::bernoulli_kernel_name(kernel)) + ": reference-panel SHA-256 " + got);
        VCAL_CHECK_EQ(got, want);
    });
}

// The AR(1)-factor variant (D-190): every platform reproduces the mirror's AR(1) reference panels bit for
// bit, with their own checked-in SHA-256 (same canonical bytes as the base panels).
VCAL_TEST(ar1_reference_panels_are_bitwise_identical_to_the_mirror) {
    const auto t = vcal::test::read_golden_csv("dgp/reference_panels_ar1.csv");
    struct Rows {
        std::uint64_t seed;
        std::uint32_t scenario, replicate;
        double pd, rho, phi;
        std::vector<std::int64_t> n, d;
        std::vector<double> z;
    };
    std::vector<Rows> panels;
    for (const auto& r : t.rows) {
        const auto id = static_cast<std::size_t>(parse_int(r[t.column("panel")]));
        if (panels.size() <= id) panels.resize(id + 1);
        auto& p = panels[id];
        p.seed = std::stoull(r[t.column("seed_hex")], nullptr, 16);
        p.scenario = static_cast<std::uint32_t>(parse_int(r[t.column("scenario")]));
        p.replicate = static_cast<std::uint32_t>(parse_int(r[t.column("replicate")]));
        p.pd = parse_double(r[t.column("pd_hex")]);
        p.rho = parse_double(r[t.column("rho_hex")]);
        p.phi = parse_double(r[t.column("phi_hex")]);
        p.n.push_back(parse_int(r[t.column("n")]));
        p.d.push_back(parse_int(r[t.column("d")]));
        p.z.push_back(parse_double(r[t.column("z_hex")]));
    }
    VCAL_REQUIRE(panels.size() == 6);
    std::ifstream in(vcal::test::golden_path("dgp/reference_panels_ar1.sha256"));
    std::string want;
    in >> want;
    for_each_kernel([&](dgp::BernoulliKernel kernel) {
        std::vector<std::uint8_t> blob;
        for (const auto& p : panels) {
            const auto T = static_cast<std::int64_t>(p.n.size());
            std::vector<std::int64_t> d(p.n.size());
            std::vector<double> z(p.n.size());
            const dgp::PanelSpec spec{p.seed, p.scenario, p.replicate, p.pd, p.rho, p.n.data(), T};
            VCAL_REQUIRE(dgp::simulate_panel_ar1(spec, p.phi, d.data(), z.data()) == dgp::Status::Ok);
            for (std::size_t i = 0; i < p.n.size(); ++i) {
                VCAL_CHECK_EQ(d[i], p.d[i]);
                VCAL_CHECK(bits_equal(z[i], p.z[i]));
                std::uint64_t zb = 0;
                std::memcpy(&zb, &z[i], sizeof zb);
                append_le(blob, static_cast<std::uint64_t>(d[i]));
                append_le(blob, zb);
            }
        }
        const std::string got = vcal::test::sha256_hex(blob);
        vcal::test::note(std::string(dgp::bernoulli_kernel_name(kernel)) + ": AR(1) reference-panel SHA-256 " + got);
        VCAL_CHECK_EQ(got, want);
    });
}

// The variant's definition, checked against the base DGP directly rather than through the mirror: phi = 0
// is the base panel bit for bit; for any phi, Z_0 is the base factor, Z_t = phi Z_{t-1} + sqrt(1 - phi^2)
// e_t with e_t the base factor (computed in that order), and d_t is the base Bernoulli count at Z_t.
VCAL_TEST(ar1_variant_is_the_base_draws_through_the_recursion) {
    const std::int64_t n[] = {100, 1000, 0, 5000, 300, 1, 2000, 100};
    constexpr std::int64_t T = 8;
    for (const std::uint32_t rep : {0u, 1u, 77u}) {
        const dgp::PanelSpec spec{0x1234ABCDull, 9, rep, 0.03, 0.15, n, T};
        std::int64_t d0[T], d[T];
        double e[T], z[T];
        VCAL_REQUIRE(dgp::simulate_panel(spec, d0, e) == dgp::Status::Ok);
        VCAL_REQUIRE(dgp::simulate_panel_ar1(spec, 0.0, d, z) == dgp::Status::Ok);
        for (std::int64_t t = 0; t < T; ++t) {
            VCAL_CHECK_EQ(d[t], d0[t]);
            VCAL_CHECK(bits_equal(z[t], e[t]));
        }
        const dgp::PhiloxKey key{static_cast<std::uint32_t>(spec.seed), static_cast<std::uint32_t>(spec.seed >> 32)};
        for (const double phi : {0.3, 0.6, -0.4}) {
            VCAL_REQUIRE(dgp::simulate_panel_ar1(spec, phi, d, z) == dgp::Status::Ok);
            const double se = std::sqrt(1.0 - phi * phi);
            for (std::int64_t t = 0; t < T; ++t) {
                const double want = t == 0 ? e[0] : phi * z[t - 1] + se * e[t];
                VCAL_CHECK(bits_equal(z[t], want));
                // d_t: the base DGP's Bernoulli uniforms for this period, compared with p(Z_t).
                const double p = dgp::det_ncdf((dgp::det_probit(spec.pd) - std::sqrt(spec.rho) * z[t]) /
                                               std::sqrt(1.0 - spec.rho));
                std::int64_t count = 0;
                for (std::int64_t i = 0; i < n[t]; i += 2) {
                    const auto w = dgp::philox4x32({spec.scenario, rep, static_cast<std::uint32_t>(t),
                                                    dgp::kFactorBlocks + static_cast<std::uint32_t>(i / 2)},
                                                   key);
                    count += dgp::uniform_from_words(w[0], w[1]) < p ? 1 : 0;
                    if (i + 1 < n[t]) count += dgp::uniform_from_words(w[2], w[3]) < p ? 1 : 0;
                }
                VCAL_CHECK_EQ(d[t], count);
            }
        }
    }
}

VCAL_TEST(ar1_rejects_phi_outside_the_open_unit_interval) {
    const std::int64_t n[] = {10, 10};
    const dgp::PanelSpec spec{1, 0, 0, 0.01, 0.1, n, 2};
    std::int64_t d[2];
    for (const double phi : {1.0, -1.0, 1.5, std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity()}) {
        VCAL_CHECK(dgp::simulate_panel_ar1(spec, phi, d, nullptr) == dgp::Status::InvalidSpec);
    }
    VCAL_CHECK(dgp::simulate_panel_ar1(spec, 0.999, d, nullptr) == dgp::Status::Ok);
}

VCAL_TEST(sha256_matches_hashlib) {
    const auto t = vcal::test::read_golden_csv("dgp/sha256_vectors.csv");
    for (const auto& r : t.rows) {
        VCAL_CHECK_EQ(vcal::test::sha256_hex(hex_bytes(r[t.column("message_hex")])), r[t.column("sha256_hex")]);
    }
}

// Stream addressing (D-109): a panel's first periods do not depend on how many periods follow,
// and changing n in one period changes no other period.
VCAL_TEST(stream_addressing_isolates_periods) {
    const std::int64_t n_short[] = {100, 200, 300};
    const std::int64_t n_long[] = {100, 200, 300, 400, 500};
    const std::int64_t n_changed[] = {100, 5000, 300};
    std::int64_t d1[3], d2[5], d3[3];
    double z1[3], z2[5], z3[3];
    VCAL_REQUIRE(dgp::simulate_panel({42, 3, 9, 0.02, 0.2, n_short, 3}, d1, z1) == dgp::Status::Ok);
    VCAL_REQUIRE(dgp::simulate_panel({42, 3, 9, 0.02, 0.2, n_long, 5}, d2, z2) == dgp::Status::Ok);
    VCAL_REQUIRE(dgp::simulate_panel({42, 3, 9, 0.02, 0.2, n_changed, 3}, d3, z3) == dgp::Status::Ok);
    for (int t = 0; t < 3; ++t) {
        VCAL_CHECK(d1[t] == d2[t] && bits_equal(z1[t], z2[t]));
        VCAL_CHECK(bits_equal(z1[t], z3[t]));  // the factor never depends on n
    }
    VCAL_CHECK(d1[0] == d3[0] && d1[2] == d3[2]);
}

VCAL_TEST(factor_stream_exhaustion_and_invalid_specs_are_refused) {
    double z = 0.0;
    VCAL_CHECK(dgp::factor_draw(1, 0, 0, 0, 0, &z) == dgp::Status::FactorStreamExhausted);
    const std::int64_t n[] = {10};
    std::int64_t d[1];
    VCAL_CHECK(dgp::simulate_panel({1, 0, 0, 0.0, 0.1, n, 1}, d, nullptr) == dgp::Status::InvalidSpec);
    VCAL_CHECK(dgp::simulate_panel({1, 0, 0, 0.1, 1.0, n, 1}, d, nullptr) == dgp::Status::InvalidSpec);
    VCAL_CHECK(dgp::simulate_panel({1, 0, 0, 0.1, 0.1, n, 0}, d, nullptr) == dgp::Status::InvalidSpec);
    const std::int64_t negative[] = {-1};
    VCAL_CHECK(dgp::simulate_panel({1, 0, 0, 0.1, 0.1, negative, 1}, d, nullptr) == dgp::Status::InvalidSpec);
}

// Moments over many replicates: E[d/n] = PD and Var(d/n) = (E[p^2] - PD^2) + (PD - E[p^2]) / n,
// with E[p^2] = Phi2(c, c; rho) from the mpmath golden table. Standard errors come from batch
// means (40 batches of 1000 replicates), which stay honest for the right-skewed, heavy-tailed
// default rate where a fourth-moment SE is biased low. The seed is fixed, so the outcome is
// deterministic; the bound is in standard errors.
VCAL_TEST(default_rate_moments_match_the_model) {
    const double pd = 0.05;
    const double rho = 0.12;
    const std::int64_t n = 200;
    const int batches = 40;
    const int per_batch = 1000;
    const auto t = vcal::test::read_golden_csv("quadrature/mixture.csv");
    double ep2 = std::nan("");
    for (const auto& r : t.rows) {
        if (r[t.column("kind")] == "closed" && parse_double(r[t.column("pd_hex")]) == pd &&
            parse_double(r[t.column("rho_hex")]) == rho && parse_int(r[t.column("n")]) == 2 && parse_int(r[t.column("d")]) == 2) {
            ep2 = std::exp(parse_double(r[t.column("log_integral_hex")]));
        }
    }
    VCAL_REQUIRE(!std::isnan(ep2));
    std::vector<double> batch_mean(batches);
    std::vector<double> batch_var(batches);
    for (int b = 0; b < batches; ++b) {
        std::vector<double> x;
        for (int i = 0; i < per_batch; ++i) {
            std::int64_t d = 0;
            const auto rep = static_cast<std::uint32_t>(b * per_batch + i);
            VCAL_REQUIRE(dgp::simulate_panel({7, 0, rep, pd, rho, &n, 1}, &d, nullptr) == dgp::Status::Ok);
            x.push_back(static_cast<double>(d) / static_cast<double>(n));
        }
        double m = 0.0;
        for (const double v : x) m += v;
        m /= per_batch;
        double ss = 0.0;
        for (const double v : x) ss += (v - m) * (v - m);
        batch_mean[static_cast<std::size_t>(b)] = m;
        batch_var[static_cast<std::size_t>(b)] = ss / (per_batch - 1);
    }
    const auto mean_and_se = [&](const std::vector<double>& v) {
        double m = 0.0;
        for (const double x : v) m += x;
        m /= batches;
        double ss = 0.0;
        for (const double x : v) ss += (x - m) * (x - m);
        return std::pair<double, double>{m, std::sqrt(ss / (batches - 1) / batches)};
    };
    const auto [mean, se_mean] = mean_and_se(batch_mean);
    const auto [var, se_var] = mean_and_se(batch_var);
    const double var_model = (ep2 - pd * pd) + (pd - ep2) / static_cast<double>(n);
    const double z_mean = (mean - pd) / se_mean;
    const double z_var = (var - var_model) / se_var;
    vcal::test::note("mean " + describe(mean) + " vs " + describe(pd) + " (z " + describe(z_mean) + "); var " +
                     describe(var) + " vs " + describe(var_model) + " (z " + describe(z_var) + ")");
    VCAL_CHECK(std::fabs(z_mean) <= tol::TOL_DGP_MOMENT_Z);
    VCAL_CHECK(std::fabs(z_var) <= tol::TOL_DGP_MOMENT_Z);
}
