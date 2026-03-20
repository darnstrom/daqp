/*
 * bench_qp.c — standalone C benchmark for DAQP.
 *
 * Measures solve time for dense QPs at several problem sizes using a simple
 * LCG random number generator (no external dependencies beyond libdaqpstat
 * and libm, which are part of the standard C environment).
 *
 * Build:
 *   gcc -O3 -I<include_dir> bench_qp.c <path/to/libdaqpstat.a> -lm -o bench_qp
 *
 * Example (from repo root after a CMake build in ./build):
 *   gcc -O3 -I include -I codegen bench_qp.c build/libdaqpstat.a -lm -o bench_qp
 *   ./bench_qp
 *
 * To compare two builds (e.g. v0.8.0 vs current):
 *   ./bench_qp_v080   | tee v080.txt
 *   ./bench_qp_head   | tee head.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "api.h"

/* -------------------------------------------------------------------------
 * Minimal LCG — no stdlib rand() so results are identical across platforms.
 * ------------------------------------------------------------------------- */
static unsigned long long lcg_state;

static void lcg_seed(unsigned long long s) { lcg_state = s; }

static double lcg_next(void) {
    lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(lcg_state >> 32) * (1.0 / 4294967296.0);
}

/* -------------------------------------------------------------------------
 * Wall-clock time in seconds (POSIX monotonic clock).
 * ------------------------------------------------------------------------- */
static double wall_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

/* -------------------------------------------------------------------------
 * Generate a symmetric positive-definite Hessian H (n×n, upper triangular,
 * column-major) that is diagonally dominant so the solver stays well-conditioned.
 * The same H is reused across repetitions to isolate solve-time cost.
 * ------------------------------------------------------------------------- */
static void make_hessian(double *H, int n) {
    int i, j;
    lcg_seed(0xDEADBEEFCAFEBABEULL);
    /* Fill upper triangle */
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            double v = lcg_next() - 0.5;
            H[i * n + j] = v;
            H[j * n + i] = v;
        }
        /* Make diagonally dominant */
        H[i * n + i] = (double)(2 * n) + 1.0;
    }
}

/* -------------------------------------------------------------------------
 * Single benchmark pass: solve N_REP QPs and return elapsed wall-clock seconds.
 * Every repetition generates a fresh f, A, b to avoid dead-code elimination.
 * ------------------------------------------------------------------------- */
#define N_PASSES 7        /* best-of N_PASSES to reduce OS jitter          */
#define N_REP    3000     /* solves per pass                                 */

static double bench_once(double *H, double *f, double *A,
                         double *bupper, double *blower, int *sense,
                         double *x, double *lam,
                         int n, int m, int ms,
                         DAQPSettings *settings) {
    int rep, i, k;
    double t0, t1;
    DAQPResult res;
    res.x   = x;
    res.lam = lam;

    lcg_seed(0x0123456789ABCDEFULL);
    t0 = wall_time();
    for (rep = 0; rep < N_REP; rep++) {
        /* Fresh cost vector */
        for (i = 0; i < n; i++) f[i] = lcg_next() - 0.5;
        /* Fresh constraint matrix */
        for (k = 0; k < (m - ms) * n; k++) A[k] = lcg_next() - 0.5;
        /* Fresh bounds */
        for (i = 0; i < ms; i++) {
            bupper[i] =  3.0 + lcg_next();
            blower[i] = -(3.0 + lcg_next());
        }
        for (i = ms; i < m; i++) {
            bupper[i] =  2.0 + lcg_next();
            blower[i] = -(2.0 + lcg_next());
        }

        DAQPProblem qp;
        qp.n       = n;
        qp.m       = m;
        qp.ms      = ms;
        qp.H       = H;
        qp.f       = f;
        qp.A       = A;
        qp.bupper  = bupper;
        qp.blower  = blower;
        qp.sense   = sense;
        qp.nh      = 0;
        qp.problem_type = 0;

        daqp_quadprog(&res, &qp, settings);
    }
    t1 = wall_time();
    return t1 - t0;
}

int main(void) {
    /* Problem sizes: n = decision variables, m = total constraints,
     * ms = simple (box) constraints.                                         */
    static const int ns[]  = { 10,  20,  40,  60,  80 };
    static const int ms_[] = {  5,  10,  20,  30,  40 };
    static const int mg[]  = { 10,  20,  40,  60,  80 }; /* general constraints */
    const int N_SIZES = (int)(sizeof(ns) / sizeof(ns[0]));

    DAQPSettings settings;
    daqp_default_settings(&settings);

    printf("DAQP benchmark  (%d problems × %d passes, best-of-%d)\n\n",
           N_SIZES, N_REP, N_PASSES);
    printf("%-6s  %-8s  %-10s  %-12s  %s\n",
           "n", "ms", "m", "avg_us", "total_ms");
    printf("%-6s  %-8s  %-10s  %-12s  %s\n",
           "------", "--------", "----------", "------------", "----------");

    double grand_total = 0.0;
    int sz;
    for (sz = 0; sz < N_SIZES; sz++) {
        const int n  = ns[sz];
        const int ms = ms_[sz];
        const int m  = ms + mg[sz];

        double *H      = (double*)malloc((size_t)(n * n)   * sizeof(double));
        double *f      = (double*)malloc((size_t)n          * sizeof(double));
        double *A      = (double*)malloc((size_t)(mg[sz]*n) * sizeof(double));
        double *bupper = (double*)malloc((size_t)m          * sizeof(double));
        double *blower = (double*)malloc((size_t)m          * sizeof(double));
        int    *sense  = (int*)  calloc( (size_t)m,           sizeof(int));
        double *x      = (double*)malloc((size_t)n          * sizeof(double));
        double *lam    = (double*)malloc((size_t)m          * sizeof(double));

        if (!H || !f || !A || !bupper || !blower || !sense || !x || !lam) {
            fprintf(stderr, "malloc failed for n=%d\n", n);
            return 1;
        }

        make_hessian(H, n);

        double best = 1e30;
        int pass;
        for (pass = 0; pass < N_PASSES; pass++) {
            double elapsed = bench_once(H, f, A, bupper, blower, sense,
                                        x, lam, n, m, ms, &settings);
            if (elapsed < best) best = elapsed;
        }

        grand_total += best;
        printf("%-6d  %-8d  %-10d  %-12.2f  %.2f\n",
               n, ms, m,
               best / N_REP * 1.0e6,   /* µs per solve */
               best * 1.0e3);           /* ms total for N_REP solves */

        free(H); free(f); free(A); free(bupper); free(blower);
        free(sense); free(x); free(lam);
    }

    printf("%-6s  %-8s  %-10s  %-12s  %s\n",
           "------", "--------", "----------", "------------", "----------");
    printf("Grand total (best-of-%d, %d solves each): %.3f s\n",
           N_PASSES, N_REP, grand_total);
    return 0;
}
