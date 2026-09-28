/* SPDX-License-Identifier: Apache-2.0
 *
 * The C ABI from plain C (M2c, D-144): compiled as C11, never as C++, and linked against the
 * shared library like any outside consumer. It calibrates, resamples and simulates through
 * vcal.h only, and checks the conventions: versioning, struct_size handling, size queries,
 * error codes and messages, and layouts without padding. Bitwise agreement with the engine is
 * abi_engine_equivalence's job; this test checks what a C caller can see.
 *
 * Usage: abi_c_test <tests/golden/dgp/reference_panels.csv>
 */
#define _CRT_SECURE_NO_WARNINGS /* fopen and strtok on MSVC */

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vcal/vcal.h"

/* ---- layouts: sizeof is the sum of the fields, and nothing trails the last one ---- */

#define VCAL_LAYOUT(type, last, size)                                                        \
    _Static_assert(sizeof(type) == (size), #type " is not " #size " bytes");                 \
    _Static_assert(offsetof(type, last) + sizeof(((type*)0)->last) == sizeof(type), #type " has trailing padding")

VCAL_LAYOUT(vcal_context_options, reserved, 16);
VCAL_LAYOUT(vcal_panel, n_defaults, 32);
VCAL_LAYOUT(vcal_grid, rho_scale, 56);
VCAL_LAYOUT(vcal_estimate, nan_count, 96);
VCAL_LAYOUT(vcal_profile_intervals, evaluations, 88);
VCAL_LAYOUT(vcal_resample_spec, weights, 88);
VCAL_LAYOUT(vcal_replicates, flags, 56);
VCAL_LAYOUT(vcal_percentile_intervals, grid_edge, 80);
VCAL_LAYOUT(vcal_dgp_spec, n_obligors, 56);
VCAL_LAYOUT(vcal_rate_series, detection_limits, 32);
VCAL_LAYOUT(vcal_moments_estimate, pd2, 32);
VCAL_LAYOUT(vcal_posterior, rho_sd_logit, 112);

/* ---- a minimal runner ---- */

static int g_checks = 0;
static int g_failures = 0;

static const char* last_error(void) {
    static char buffer[1024];
    int64_t required = 0;
    if (vcal_last_error(buffer, (int64_t)sizeof buffer, &required) != VCAL_OK) return "(message longer than 1024 bytes)";
    return buffer;
}

#define CHECK(cond)                                                                   \
    do {                                                                              \
        ++g_checks;                                                                   \
        if (!(cond)) {                                                                \
            ++g_failures;                                                             \
            fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
        }                                                                             \
    } while (0)

#define CHECK_STATUS(call, expected)                                                                   \
    do {                                                                                               \
        const vcal_status st_ = (call);                                                                \
        ++g_checks;                                                                                    \
        if (st_ != (expected)) {                                                                       \
            ++g_failures;                                                                              \
            fprintf(stderr, "%s:%d: %s returned %s, expected %s; last error: \"%s\"\n", __FILE__,    \
                    __LINE__, #call, vcal_status_string(st_), vcal_status_string(expected), last_error()); \
        }                                                                                              \
    } while (0)

#define CHECK_ERROR_MENTIONS(text) CHECK(strstr(last_error(), (text)) != NULL)

static int bits_equal(double a, double b) { return memcmp(&a, &b, sizeof a) == 0; }

/* ---- fixtures ---- */

enum { T = 40, B = 199 };
static int64_t g_n[T];
static int64_t g_d[T];

static vcal_panel make_panel(void) {
    vcal_panel p;
    memset(&p, 0, sizeof p);
    p.struct_size = (uint32_t)sizeof p;
    p.n_periods = T;
    p.n_obligors = g_n;
    p.n_defaults = g_d;
    return p;
}

static vcal_grid default_grid(void) {
    vcal_grid g;
    memset(&g, 0, sizeof g);
    g.struct_size = (uint32_t)sizeof g;
    CHECK_STATUS(vcal_grid_default(&g), VCAL_OK);
    return g;
}

static vcal_resample_spec make_spec(int32_t scheme) {
    vcal_resample_spec s;
    memset(&s, 0, sizeof s);
    s.struct_size = (uint32_t)sizeof s;
    s.scheme = scheme;
    return s;
}

/* ---- tests ---- */

static void test_version_and_build_info(void) {
    char small[4];
    char info[2048];
    int64_t required = 0;
    int64_t again = 0;
    CHECK(vcal_abi_version() == VCAL_ABI_VERSION);
    CHECK((vcal_abi_version() >> 16) == VCAL_ABI_VERSION_MAJOR);

    CHECK_STATUS(vcal_build_info(NULL, 0, &required), VCAL_OK); /* size query */
    CHECK(required > 1 && required <= (int64_t)sizeof info);
    CHECK_STATUS(vcal_build_info(small, (int64_t)sizeof small, &again), VCAL_E_BUFFER_TOO_SMALL);
    CHECK(again == required);
    CHECK_ERROR_MENTIONS("capacity 4");
    CHECK_STATUS(vcal_build_info(info, (int64_t)sizeof info, &again), VCAL_OK);
    CHECK((int64_t)strlen(info) + 1 == required);
    CHECK(strstr(info, "abi_version=0.3\n") != NULL);
    CHECK(strstr(info, "git_commit=") != NULL);
    CHECK(strstr(info, "profiles=parity,native\n") != NULL);
    printf("%s", info);
}

static void test_errors_and_context(void) {
    vcal_context* ctx = NULL;
    vcal_context_options o;
    struct {
        vcal_context_options options;
        uint64_t newer_field; /* what a caller built against a later minor version would pass */
    } newer;
    int64_t required = 0;
    char tiny[2];

    CHECK(strcmp(vcal_status_string(VCAL_OK), "VCAL_OK") == 0);
    CHECK(strcmp(vcal_status_string(VCAL_E_INVALID_ARGUMENT), "VCAL_E_INVALID_ARGUMENT") == 0);
    CHECK(strcmp(vcal_status_string(VCAL_E_BUFFER_TOO_SMALL), "VCAL_E_BUFFER_TOO_SMALL") == 0);
    CHECK(strcmp(vcal_status_string(VCAL_E_UNSUPPORTED), "VCAL_E_UNSUPPORTED") == 0);
    CHECK(strcmp(vcal_status_string(VCAL_E_NUMERIC), "VCAL_E_NUMERIC") == 0);
    CHECK(strcmp(vcal_status_string(VCAL_E_OUT_OF_MEMORY), "VCAL_E_OUT_OF_MEMORY") == 0);
    CHECK(strcmp(vcal_status_string(VCAL_E_FP_ENVIRONMENT), "VCAL_E_FP_ENVIRONMENT") == 0);
    CHECK(strcmp(vcal_status_string(VCAL_E_INTERNAL), "VCAL_E_INTERNAL") == 0);
    CHECK(vcal_status_string(12345) != NULL);

    CHECK_STATUS(vcal_context_create(NULL, NULL), VCAL_E_INVALID_ARGUMENT);
    CHECK_ERROR_MENTIONS("vcal_context_create: context is NULL");
    /* A successful call clears the message; vcal_status_string and vcal_last_error do not. */
    CHECK_STATUS(vcal_context_create(NULL, &ctx), VCAL_OK);
    CHECK(ctx != NULL);
    CHECK(strcmp(last_error(), "") == 0);
    CHECK_STATUS(vcal_context_destroy(ctx), VCAL_OK);
    CHECK_STATUS(vcal_context_destroy(NULL), VCAL_OK);

    memset(&o, 0, sizeof o);
    o.struct_size = (uint32_t)sizeof o;
    o.n_threads = -1;
    CHECK_STATUS(vcal_context_create(&o, &ctx), VCAL_E_INVALID_ARGUMENT);
    CHECK(ctx == NULL);
    (void)vcal_status_string(VCAL_OK);
    CHECK_ERROR_MENTIONS("options.n_threads = -1");
    CHECK_STATUS(vcal_last_error(tiny, (int64_t)sizeof tiny, &required), VCAL_E_BUFFER_TOO_SMALL);
    CHECK(required == (int64_t)strlen(last_error()) + 1);
    CHECK_STATUS(vcal_last_error(NULL, 0, NULL), VCAL_E_INVALID_ARGUMENT);

    o.n_threads = 0;
    o.profile = 7;
    CHECK_STATUS(vcal_context_create(&o, &ctx), VCAL_E_INVALID_ARGUMENT);
    CHECK_ERROR_MENTIONS("options.profile = 7");
    o.profile = VCAL_PROFILE_PARITY;
    o.struct_size = 8;
    CHECK_STATUS(vcal_context_create(&o, &ctx), VCAL_E_INVALID_ARGUMENT);
    CHECK_ERROR_MENTIONS("options.struct_size = 8");
    o.struct_size = (uint32_t)sizeof o;
    o.reserved = 1;
    CHECK_STATUS(vcal_context_create(&o, &ctx), VCAL_E_UNSUPPORTED);
    o.reserved = 0;

    /* A larger struct from a newer caller: fine while the fields this library lacks are zero. */
    memset(&newer, 0, sizeof newer);
    newer.options.struct_size = (uint32_t)sizeof newer;
    newer.options.n_threads = 2;
    CHECK_STATUS(vcal_context_create(&newer.options, &ctx), VCAL_OK);
    CHECK_STATUS(vcal_context_destroy(ctx), VCAL_OK);
    newer.newer_field = 1;
    CHECK_STATUS(vcal_context_create(&newer.options, &ctx), VCAL_E_UNSUPPORTED);
    CHECK_ERROR_MENTIONS("nonzero byte at offset 16");
}

/* Every reference panel of the Python mirror (tools/gen_dgp_tables.py), bit for bit. */
static void test_dgp_reference_panels(vcal_context* ctx, const char* path) {
    enum { kMaxRows = 256 };
    char line[512];
    long panel_of[kMaxRows], period[kMaxRows];
    vcal_dgp_spec spec[kMaxRows];
    int64_t n[kMaxRows], d[kMaxRows];
    double z[kMaxRows];
    int rows = 0, panels = 0, r0 = 0;
    FILE* f = fopen(path, "r");
    CHECK(f != NULL);
    if (f == NULL) return;
    while (fgets(line, (int)sizeof line, f) != NULL && rows < kMaxRows) {
        char* field[10];
        int k = 0;
        if (line[0] == '#' || strncmp(line, "panel,", 6) == 0) continue;
        for (char* tok = strtok(line, ",\r\n"); tok != NULL && k < 10; tok = strtok(NULL, ",\r\n")) field[k++] = tok;
        CHECK(k == 10);
        if (k != 10) break;
        /* panel,seed_hex,scenario,replicate,pd_hex,rho_hex,period,n,d,z_hex */
        memset(&spec[rows], 0, sizeof spec[rows]);
        spec[rows].struct_size = (uint32_t)sizeof spec[rows];
        panel_of[rows] = strtol(field[0], NULL, 10);
        spec[rows].seed = (uint64_t)strtoull(field[1], NULL, 16);
        spec[rows].scenario = (uint32_t)strtoul(field[2], NULL, 10);
        spec[rows].replicate = (uint32_t)strtoul(field[3], NULL, 10);
        spec[rows].pd = strtod(field[4], NULL); /* C99 hex floats: exact */
        spec[rows].rho = strtod(field[5], NULL);
        period[rows] = strtol(field[6], NULL, 10);
        n[rows] = strtoll(field[7], NULL, 10);
        d[rows] = strtoll(field[8], NULL, 10);
        z[rows] = strtod(field[9], NULL);
        ++rows;
    }
    fclose(f);
    CHECK(rows > 0);
    for (int r = 1; r <= rows; ++r) {
        if (r < rows && panel_of[r] == panel_of[r0]) continue;
        {
            const int periods = r - r0;
            int64_t got_d[kMaxRows];
            double got_z[kMaxRows];
            int64_t required = 0;
            vcal_dgp_spec s = spec[r0];
            s.n_periods = periods;
            s.n_obligors = &n[r0];
            CHECK_STATUS(vcal_dgp_simulate(ctx, &s, NULL, NULL, 0, &required), VCAL_OK);
            CHECK(required == periods);
            CHECK_STATUS(vcal_dgp_simulate(ctx, &s, got_d, got_z, periods, NULL), VCAL_OK);
            for (int t = 0; t < periods; ++t) {
                CHECK(period[r0 + t] == t);
                CHECK(got_d[t] == d[r0 + t]);
                CHECK(bits_equal(got_z[t], z[r0 + t]));
            }
            ++panels;
            r0 = r;
        }
    }
    printf("DGP: %d reference panels (%d periods) reproduced bit for bit\n", panels, rows);
}

static void test_calibrate_resample(vcal_context* ctx) {
    const vcal_panel panel = make_panel();
    const vcal_grid grid = default_grid();
    const int64_t K = (int64_t)grid.pd_points * grid.rho_points;
    vcal_estimate est;
    vcal_profile_intervals prof;
    vcal_dgp_spec dgp;
    int64_t required = 0;

    /* The panel: T = 40 periods of 1,000 obligors, PD 1%, rho 0.12. */
    memset(&dgp, 0, sizeof dgp);
    dgp.struct_size = (uint32_t)sizeof dgp;
    dgp.seed = 0x4D32434142494341ull;
    dgp.pd = 0.01;
    dgp.rho = 0.12;
    dgp.n_periods = T;
    dgp.n_obligors = g_n;
    for (int t = 0; t < T; ++t) g_n[t] = 1000;
    CHECK_STATUS(vcal_dgp_simulate(ctx, &dgp, g_d, NULL, T, NULL), VCAL_OK);

    /* Calibrate, with profile intervals. */
    memset(&est, 0, sizeof est);
    memset(&prof, 0, sizeof prof);
    est.struct_size = (uint32_t)sizeof est;
    prof.struct_size = (uint32_t)sizeof prof;
    CHECK_STATUS(vcal_calibrate(ctx, &panel, &grid, &est, &prof), VCAL_OK);
    printf("estimate: pd %.17g (se %.3g), rho %.17g (se %.3g), loglik %.17g, flags %u\n", est.pd, est.se_pd,
           est.rho, est.se_rho, est.loglik, est.flags);
    printf("profile:  pd [%.17g, %.17g], rho [%.17g, %.17g], flags %u/%u, residual %.3g\n", prof.pd_lo, prof.pd_hi,
           prof.rho_lo, prof.rho_hi, prof.pd_flags, prof.rho_flags, prof.residual_max);
    CHECK(est.struct_size == (uint32_t)sizeof est);
    CHECK((est.flags & (VCAL_FLAG_GRID_EDGE | VCAL_FLAG_FLAT_SURFACE | VCAL_FLAG_NUMERIC)) == 0);
    CHECK(est.se_pd > 0.0 && est.se_rho > 0.0 && fabs(est.corr_pd_rho) < 1.0);
    CHECK(isfinite(est.loglik) && est.loglik < 0.0);
    CHECK(est.grid_index >= 0 && est.grid_index < K);
    CHECK(prof.pd_lo < est.pd && est.pd < prof.pd_hi);
    CHECK(prof.rho_lo < est.rho && est.rho < prof.rho_hi);
    CHECK(prof.pd_flags == 0 && prof.rho_flags == 0 && prof.reserved == 0);
    CHECK(prof.evaluations > 0);
    CHECK((est.flags & VCAL_FLAG_RHO_NOT_IDENTIFIED) == 0);

    /* D-302: single-obligor periods carry no information about rho. The numbers are still
     * reported, the flag says rho is meaningless, and rho's profile interval is the whole box. */
    {
        int64_t n1[20], d1[20];
        vcal_panel single;
        vcal_estimate e1;
        vcal_profile_intervals p1;
        for (int t = 0; t < 20; ++t) {
            n1[t] = 1;
            d1[t] = t % 4 == 0 ? 1 : 0;
        }
        memset(&single, 0, sizeof single);
        single.struct_size = (uint32_t)sizeof single;
        single.n_periods = 20;
        single.n_obligors = n1;
        single.n_defaults = d1;
        memset(&e1, 0, sizeof e1);
        memset(&p1, 0, sizeof p1);
        e1.struct_size = (uint32_t)sizeof e1;
        p1.struct_size = (uint32_t)sizeof p1;
        CHECK_STATUS(vcal_calibrate(ctx, &single, &grid, &e1, &p1), VCAL_OK);
        CHECK((e1.flags & VCAL_FLAG_RHO_NOT_IDENTIFIED) != 0);
        CHECK(isfinite(e1.pd) && isfinite(e1.rho) && isfinite(e1.loglik));
        CHECK(p1.rho_lo == grid.rho_lo && p1.rho_hi == grid.rho_hi);
        CHECK(p1.rho_flags == (VCAL_INTERVAL_LOWER_TRUNCATED | VCAL_INTERVAL_UPPER_TRUNCATED));
        n1[7] = 2; /* one period with two obligors identifies rho */
        CHECK_STATUS(vcal_calibrate(ctx, &single, &grid, &e1, NULL), VCAL_OK);
        CHECK((e1.flags & VCAL_FLAG_RHO_NOT_IDENTIFIED) == 0);
    }

    /* Profile is optional, and results do not depend on the thread count. */
    {
        vcal_context* one = NULL;
        vcal_context_options o;
        vcal_estimate again = est;
        memset(&o, 0, sizeof o);
        o.struct_size = (uint32_t)sizeof o;
        o.n_threads = 1;
        CHECK_STATUS(vcal_context_create(&o, &one), VCAL_OK);
        CHECK_STATUS(vcal_calibrate(one, &panel, &grid, &again, NULL), VCAL_OK);
        CHECK(memcmp(&again, &est, sizeof est) == 0);
        CHECK_STATUS(vcal_context_destroy(one), VCAL_OK);
    }

    /* An output struct from a newer caller: known fields written, the rest untouched. */
    {
        struct {
            vcal_estimate estimate;
            double newer_field;
        } newer;
        memset(&newer, 0, sizeof newer);
        newer.estimate.struct_size = (uint32_t)sizeof newer;
        newer.newer_field = 12345.0;
        CHECK_STATUS(vcal_calibrate(ctx, &panel, &grid, &newer.estimate, NULL), VCAL_OK);
        CHECK(newer.estimate.struct_size == (uint32_t)sizeof newer);
        CHECK(bits_equal(newer.estimate.pd, est.pd) && newer.newer_field == 12345.0);
    }

    /* Grid values: the endpoints are exactly the bounds. */
    {
        double pd_values[61];
        CHECK_STATUS(vcal_grid_values(&grid, VCAL_AXIS_PD, NULL, 0, &required), VCAL_OK);
        CHECK(required == 61);
        CHECK_STATUS(vcal_grid_values(&grid, VCAL_AXIS_PD, pd_values, 61, NULL), VCAL_OK);
        CHECK(pd_values[0] == grid.pd_lo && pd_values[60] == grid.pd_hi);
        CHECK_STATUS(vcal_grid_values(&grid, 2, pd_values, 61, NULL), VCAL_E_INVALID_ARGUMENT);
    }

    /* The surface: its argmax column is the estimate's grid index. */
    {
        double* L;
        int64_t best = -1;
        double best_sum = -INFINITY;
        CHECK_STATUS(vcal_surface(ctx, &panel, &grid, NULL, 0, &required), VCAL_OK);
        CHECK(required == T * K);
        L = (double*)malloc((size_t)required * sizeof(double));
        CHECK(L != NULL);
        if (L == NULL) return;
        CHECK_STATUS(vcal_surface(ctx, &panel, &grid, L, required - 1, NULL), VCAL_E_BUFFER_TOO_SMALL);
        CHECK_STATUS(vcal_surface(ctx, &panel, &grid, L, required, NULL), VCAL_OK);
        for (int64_t k = 0; k < K; ++k) {
            double s = 0.0;
            for (int t = 0; t < T; ++t) s += L[t * K + k];
            if (s > best_sum) {
                best_sum = s;
                best = k;
            }
        }
        CHECK(best == est.grid_index);
        free(L);
    }

    /* iid bootstrap: size query, replicates, intervals, determinism. */
    {
        vcal_resample_spec s = make_spec(VCAL_RESAMPLE_IID_BOOTSTRAP);
        vcal_replicates reps;
        vcal_percentile_intervals iv;
        double pd[B], rho[B], ll[B], pd2[B];
        int64_t gi[B];
        uint32_t flags[B];
        s.seed = 20260927u;
        s.replicates = B;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, NULL, NULL, &required), VCAL_OK);
        CHECK(required == B);
        memset(&reps, 0, sizeof reps);
        memset(&iv, 0, sizeof iv);
        reps.struct_size = (uint32_t)sizeof reps;
        iv.struct_size = (uint32_t)sizeof iv;
        reps.capacity = B - 1;
        reps.pd = pd;
        reps.rho = rho;
        reps.grid_loglik = ll;
        reps.grid_index = gi;
        reps.flags = flags;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, &reps, &iv, &required), VCAL_E_BUFFER_TOO_SMALL);
        reps.capacity = B;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, &reps, &iv, NULL), VCAL_OK);
        CHECK(iv.replicates == B && iv.level == 0.95 && iv.pd_excluded == 0 && iv.rho_excluded == 0);
        CHECK(iv.pd_lo < est.pd && est.pd < iv.pd_hi);
        CHECK(iv.rho_lo <= iv.rho_hi);
        for (int b = 0; b < B; ++b) {
            CHECK(pd[b] >= grid.pd_lo && pd[b] <= grid.pd_hi && rho[b] >= grid.rho_lo && rho[b] <= grid.rho_hi);
            CHECK(isfinite(ll[b]) && gi[b] >= 0 && gi[b] < K);
        }
        printf("iid bootstrap (B = %d): pd [%.6g, %.6g], rho [%.6g, %.6g], %lld on the grid edge\n", B, iv.pd_lo,
               iv.pd_hi, iv.rho_lo, iv.rho_hi, (long long)iv.grid_edge);
        reps.pd = pd2;
        reps.rho = NULL;
        reps.grid_loglik = NULL;
        reps.grid_index = NULL;
        reps.flags = NULL;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, &reps, NULL, NULL), VCAL_OK);
        CHECK(memcmp(pd, pd2, sizeof pd) == 0);

        /* Its weights are counts: each row sums to T. */
        {
            double* W;
            CHECK_STATUS(vcal_resample_weights(ctx, &s, T, NULL, 0, &required), VCAL_OK);
            CHECK(required == (int64_t)B * T);
            W = (double*)malloc((size_t)required * sizeof(double));
            CHECK(W != NULL);
            if (W != NULL) {
                CHECK_STATUS(vcal_resample_weights(ctx, &s, T, W, required, NULL), VCAL_OK);
                for (int b = 0; b < B; ++b) {
                    double sum = 0.0;
                    for (int t = 0; t < T; ++t) sum += W[b * T + t];
                    CHECK(sum == (double)T);
                }
                free(W);
            }
        }
    }

    /* Supplied indices 0..T-1 are the original panel: the replicate is the estimate. */
    {
        vcal_resample_spec s = make_spec(VCAL_RESAMPLE_INDICES);
        vcal_replicates reps;
        int64_t idx[T];
        double pd, rho;
        for (int t = 0; t < T; ++t) idx[t] = t;
        s.replicates = 1;
        s.draws = T;
        s.indices = idx;
        memset(&reps, 0, sizeof reps);
        reps.struct_size = (uint32_t)sizeof reps;
        reps.capacity = 1;
        reps.pd = &pd;
        reps.rho = &rho;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, &reps, NULL, NULL), VCAL_OK);
        CHECK(bits_equal(pd, est.pd) && bits_equal(rho, est.rho));
        idx[5] = T;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, &reps, NULL, NULL), VCAL_E_INVALID_ARGUMENT);
    }

    /* Replicate counts of the deterministic schemes, and fields a scheme does not use. */
    {
        vcal_resample_spec s = make_spec(VCAL_RESAMPLE_JACKKNIFE);
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, NULL, NULL, &required), VCAL_OK);
        CHECK(required == T);
        s.seed = 1;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, NULL, NULL, &required), VCAL_E_INVALID_ARGUMENT);
        CHECK_ERROR_MENTIONS("spec.seed is not used by VCAL_RESAMPLE_JACKKNIFE");
        s = make_spec(VCAL_RESAMPLE_WALK_FORWARD);
        s.first_end = 20;
        s.step = 5;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, NULL, NULL, &required), VCAL_OK);
        CHECK(required == 5);
        s = make_spec(VCAL_RESAMPLE_BLOCK_BOOTSTRAP);
        s.replicates = 3;
        s.block_length = T + 1;
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, NULL, NULL, &required), VCAL_E_INVALID_ARGUMENT);
        s = make_spec(99);
        CHECK_STATUS(vcal_resample(ctx, &panel, &grid, &s, NULL, NULL, &required), VCAL_E_INVALID_ARGUMENT);
    }

    /* Invalid inputs name the field. */
    {
        vcal_panel bad = panel;
        vcal_grid g = grid;
        const int64_t saved = g_d[3];
        g_d[3] = 1001;
        CHECK_STATUS(vcal_calibrate(ctx, &bad, &grid, &est, NULL), VCAL_E_INVALID_ARGUMENT);
        CHECK_ERROR_MENTIONS("vcal_calibrate: panel.n_defaults[3] = 1001");
        g_d[3] = saved;
        bad.n_defaults = NULL;
        CHECK_STATUS(vcal_calibrate(ctx, &bad, &grid, &est, NULL), VCAL_E_INVALID_ARGUMENT);
        g.pd_points = 2;
        CHECK_STATUS(vcal_calibrate(ctx, &panel, &g, &est, NULL), VCAL_E_INVALID_ARGUMENT);
        CHECK_ERROR_MENTIONS("at least 3 points");
        g = grid;
        g.rho_hi = 1.0;
        CHECK_STATUS(vcal_calibrate(ctx, &panel, &g, &est, NULL), VCAL_E_INVALID_ARGUMENT);
        CHECK_STATUS(vcal_calibrate(NULL, &panel, &grid, &est, NULL), VCAL_E_INVALID_ARGUMENT);
        CHECK_STATUS(vcal_calibrate(ctx, &panel, &grid, NULL, NULL), VCAL_E_INVALID_ARGUMENT);
        dgp.rho = 1.0;
        CHECK_STATUS(vcal_dgp_simulate(ctx, &dgp, g_d, NULL, T, NULL), VCAL_E_INVALID_ARGUMENT);
        CHECK_ERROR_MENTIONS("spec.rho = 1");
    }
}

/* ABI 0.3 (D-169): the M3 estimators from C, on the test panel; the numbers are checked against
 * the engine in abi_engine_equivalence, so here: they run, fill their outputs, and refuse as
 * documented. */
static void test_m3_estimators(vcal_context* parity) {
    vcal_panel p = make_panel();
    vcal_grid g = default_grid();
    vcal_estimate est;
    vcal_moments_estimate mom;
    vcal_posterior post;
    vcal_context_options o;
    vcal_context* native = NULL;
    double rates[3] = {0.01, 0.0, 0.03};
    double limits[3] = {0.005, 0.005, 0.005};
    vcal_rate_series rs;
    memset(&est, 0, sizeof est);
    est.struct_size = sizeof est;
    memset(&mom, 0, sizeof mom);
    mom.struct_size = sizeof mom;
    memset(&post, 0, sizeof post);
    post.struct_size = sizeof post;
    memset(&rs, 0, sizeof rs);
    rs.struct_size = sizeof rs;
    rs.n_periods = 3;
    rs.rates = rates;

    CHECK_STATUS(vcal_calibrate_moments(parity, &p, NULL, &g, &mom), VCAL_OK);
    CHECK(mom.pd > 0.0 && mom.pd2 > 0.0);
    CHECK_STATUS(vcal_calibrate_moments(parity, &p, &rs, &g, &mom), VCAL_E_INVALID_ARGUMENT);
    CHECK_ERROR_MENTIONS("exactly one of counts and rates");
    CHECK_STATUS(vcal_calibrate_posterior(parity, &p, &g, VCAL_PRIOR_FLAT, &post), VCAL_OK);
    CHECK((post.flags & VCAL_POSTERIOR_REFUSED) == 0 && post.pd_et_lo < post.pd_et_hi && post.rho_et_lo < post.rho_et_hi);
    /* Parity refuses a zero rate, and names the period; the native treatments need a native context. */
    CHECK_STATUS(vcal_calibrate_rate(parity, NULL, &rs, &g, VCAL_ZERO_RATES_REFUSE, &est, NULL), VCAL_E_INVALID_ARGUMENT);
    CHECK_ERROR_MENTIONS("periods 1 have a default rate of 0 or 1");
    CHECK_STATUS(vcal_calibrate_rate(parity, NULL, &rs, &g, VCAL_ZERO_RATES_CENSOR, &est, NULL), VCAL_E_INVALID_ARGUMENT);
    CHECK_ERROR_MENTIONS("VCAL_PROFILE_NATIVE");
    memset(&o, 0, sizeof o);
    o.struct_size = sizeof o;
    o.profile = VCAL_PROFILE_NATIVE;
    CHECK_STATUS(vcal_context_create(&o, &native), VCAL_OK);
    if (native != NULL) {
        CHECK_STATUS(vcal_calibrate_rate(native, NULL, &rs, &g, VCAL_ZERO_RATES_CENSOR, &est, NULL), VCAL_E_INVALID_ARGUMENT);
        CHECK_ERROR_MENTIONS("detection_limits is needed");
        rs.detection_limits = limits;
        CHECK_STATUS(vcal_calibrate_rate(native, NULL, &rs, &g, VCAL_ZERO_RATES_CENSOR, &est, NULL), VCAL_OK);
        CHECK(est.pd > 0.0 && est.rho > 0.0);
        CHECK_STATUS(vcal_calibrate_rate(native, NULL, &rs, &g, VCAL_ZERO_RATES_DROP, &est, NULL), VCAL_E_INVALID_ARGUMENT);
        CHECK_ERROR_MENTIONS("need at least 3");
        CHECK_STATUS(vcal_context_destroy(native), VCAL_OK);
    }
}

int main(int argc, char** argv) {
    vcal_context* ctx = NULL;
    if (argc != 2) {
        fprintf(stderr, "usage: abi_c_test <reference_panels.csv>\n");
        return 2;
    }
    test_version_and_build_info();
    test_errors_and_context();
    CHECK_STATUS(vcal_context_create(NULL, &ctx), VCAL_OK);
    if (ctx != NULL) {
        test_dgp_reference_panels(ctx, argv[1]);
        test_calibrate_resample(ctx);
        test_m3_estimators(ctx);
        CHECK_STATUS(vcal_context_destroy(ctx), VCAL_OK);
    }
    printf("abi_c_test: %d checks, %d failed\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
