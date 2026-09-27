// SPDX-License-Identifier: Apache-2.0
// dgp/ (M1.6; D-052..D-057, D-061, D-109, D-110): known answers, bitwise identity with the
// Python mirror on every platform, accuracy against mpmath, stream addressing, and moments.
#include "dgp/dgp.hpp"

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
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
    std::ifstream in(vcal::test::golden_path("dgp/reference_panels.sha256"));
    std::string want;
    in >> want;
    const std::string got = vcal::test::sha256_hex(blob);
    vcal::test::note("reference-panel SHA-256 " + got);
    VCAL_CHECK_EQ(got, want);
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
