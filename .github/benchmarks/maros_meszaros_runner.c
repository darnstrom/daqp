/* Benchmark DAQP on the exported Maros-Meszaros dense subset. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "api.h"
#include "constants.h"

typedef struct {
    int n, m, ms;
    double *H, *f, *A, *bupper, *blower;
    int *sense;
} Problem;

typedef struct {
    int exitflag, iterations, solved;
    double solve_time, primal, dual, gap;
} Result;

typedef struct {
    const char *name;
    double solver_tol;
    double check_tol;
    int use_defaults;
} Tolerance;

static int read_problem(const char *path, Problem *p) {
    FILE *file = fopen(path, "rb");
    if (file == NULL) return 0;
    if (fread(&p->n, sizeof(int), 1, file) != 1 ||
        fread(&p->m, sizeof(int), 1, file) != 1 ||
        fread(&p->ms, sizeof(int), 1, file) != 1) {
        fclose(file);
        return 0;
    }

    const int n = p->n, m = p->m, general = m - p->ms;
    p->H = malloc(sizeof(double) * (size_t)n * n);
    p->f = malloc(sizeof(double) * (size_t)n);
    p->A = general > 0 ? malloc(sizeof(double) * (size_t)general * n) : NULL;
    p->bupper = malloc(sizeof(double) * (size_t)m);
    p->blower = malloc(sizeof(double) * (size_t)m);
    p->sense = malloc(sizeof(int) * (size_t)m);

    int ok = p->H != NULL && p->f != NULL && p->bupper != NULL &&
             p->blower != NULL && p->sense != NULL &&
             (general == 0 || p->A != NULL);
    ok = ok && fread(p->H, sizeof(double), (size_t)n * n, file) == (size_t)n * n;
    ok = ok && fread(p->f, sizeof(double), (size_t)n, file) == (size_t)n;
    if (general > 0)
        ok = ok && fread(p->A, sizeof(double), (size_t)general * n, file) ==
                   (size_t)general * n;
    ok = ok && fread(p->bupper, sizeof(double), (size_t)m, file) == (size_t)m;
    ok = ok && fread(p->blower, sizeof(double), (size_t)m, file) == (size_t)m;
    ok = ok && fread(p->sense, sizeof(int), (size_t)m, file) == (size_t)m;
    fclose(file);
    return ok;
}

static void free_problem(Problem *p) {
    free(p->H);
    free(p->f);
    free(p->A);
    free(p->bupper);
    free(p->blower);
    free(p->sense);
}

static double row_dot(const Problem *p, int row, const double *x) {
    if (row < p->ms) return x[row];
    const double *a = p->A + (size_t)(row - p->ms) * p->n;
    double value = 0.0;
    for (int j = 0; j < p->n; ++j) value += a[j] * x[j];
    return value;
}

/* Scaled KKT residuals following the supplied Maros-Meszaros native runner. */
static void kkt_residuals(const Problem *p, const double *x, const double *lam,
                          double *primal, double *dual, double *gap) {
    double residual = 0.0, scale = 0.0;
    for (int i = 0; i < p->m; ++i) {
        const double value = row_dot(p, i, x);
        double violation = 0.0;
        if (p->bupper[i] < 1e29 && value > p->bupper[i])
            violation = value - p->bupper[i];
        if (p->blower[i] > -1e29 && p->blower[i] - value > violation)
            violation = p->blower[i] - value;
        if (violation > residual) residual = violation;
        if (fabs(value) > scale) scale = fabs(value);
    }
    *primal = residual / (1.0 + scale);

    double *gradient = calloc((size_t)p->n, sizeof(double));
    double dual_scale = 0.0;
    for (int i = 0; i < p->n; ++i) {
        double hx = 0.0;
        for (int j = 0; j < p->n; ++j)
            hx += p->H[(size_t)i * p->n + j] * x[j];
        gradient[i] = hx + p->f[i];
        if (fabs(hx) > dual_scale) dual_scale = fabs(hx);
        if (fabs(p->f[i]) > dual_scale) dual_scale = fabs(p->f[i]);
    }
    for (int i = 0; i < p->m; ++i) {
        if (lam[i] == 0.0) continue;
        if (fabs(lam[i]) > dual_scale) dual_scale = fabs(lam[i]);
        if (i < p->ms) {
            gradient[i] += lam[i];
        } else {
            const double *a = p->A + (size_t)(i - p->ms) * p->n;
            for (int j = 0; j < p->n; ++j) gradient[j] += a[j] * lam[i];
        }
    }
    residual = 0.0;
    for (int i = 0; i < p->n; ++i)
        if (fabs(gradient[i]) > residual) residual = fabs(gradient[i]);
    *dual = residual / (1.0 + dual_scale);
    free(gradient);

    double gap_value = 0.0, gap_scale = 0.0;
    for (int i = 0; i < p->m; ++i) {
        if (lam[i] == 0.0) continue;
        const double bound = lam[i] > 0.0 ? p->bupper[i] : p->blower[i];
        if (fabs(bound) > 1e29) continue;
        const double term = lam[i] * (row_dot(p, i, x) - bound);
        gap_value += term;
        gap_scale += fabs(term);
    }
    *gap = fabs(gap_value) / (1.0 + gap_scale);
}

static Result solve_problem(Problem *p, const Tolerance *tolerance, int repeats) {
    Result output;
    memset(&output, 0, sizeof(output));
    output.solve_time = HUGE_VAL;
    double *x = malloc(sizeof(double) * (size_t)p->n);
    double *lam = malloc(sizeof(double) * (size_t)p->m);
    int *sense = malloc(sizeof(int) * (size_t)p->m);
    memcpy(sense, p->sense, sizeof(int) * (size_t)p->m);

    for (int repeat = 0; repeat < repeats; ++repeat) {
        DAQPProblem qp = {
            p->n, p->m, p->ms, p->H, p->f, p->A, p->bupper, p->blower,
            p->sense, NULL, 0, 0
        };
        DAQPSettings settings;
        daqp_default_settings(&settings);
        if (!tolerance->use_defaults) {
            settings.primal_tol = tolerance->solver_tol;
            settings.dual_tol = tolerance->solver_tol;
        }
        settings.time_limit = 1000.0;

        DAQPResult result;
        memset(&result, 0, sizeof(result));
        result.x = x;
        result.lam = lam;
        daqp_quadprog(&result, &qp, &settings);
        memcpy(p->sense, sense, sizeof(int) * (size_t)p->m);

        output.exitflag = result.exitflag;
        output.iterations = result.iter;
        if (result.solve_time < output.solve_time) output.solve_time = result.solve_time;
    }

    if (output.exitflag > 0) {
        kkt_residuals(p, x, lam, &output.primal, &output.dual, &output.gap);
        output.solved = output.primal <= tolerance->check_tol &&
                        output.dual <= tolerance->check_tol &&
                        output.gap <= tolerance->check_tol;
    } else {
        output.primal = output.dual = output.gap = HUGE_VAL;
    }

    free(x);
    free(lam);
    free(sense);
    return output;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s <export-dir> <output.csv> <repeats>\n", argv[0]);
        return 2;
    }
    const char *directory = argv[1];
    const int repeats = atoi(argv[3]);
    if (repeats < 1) return 2;

    char index_path[4096];
    snprintf(index_path, sizeof(index_path), "%s/index.txt", directory);
    FILE *index = fopen(index_path, "r");
    FILE *csv = fopen(argv[2], "w");
    if (index == NULL || csv == NULL) {
        fprintf(stderr, "could not open benchmark input or output\n");
        return 1;
    }

    const Tolerance tolerances[] = {
        {"default", 0.0, 1.0, 1},
        {"low", 1e-3, 1e-3, 0},
        {"med", 1e-6, 1e-6, 0},
        {"high", 1e-9, 1e-9, 0},
    };
    fprintf(csv, "problem,n,m,posdef,tolerance,solved,solve_time_s,exitflag,"
                 "iterations,primal_residual,dual_residual,duality_gap\n");

    char name[256], path[4096];
    int n, m, posdef, problem_count = 0;
    while (fscanf(index, "%255s %d %d %d", name, &n, &m, &posdef) == 4) {
        snprintf(path, sizeof(path), "%s/%s.bin", directory, name);
        Problem problem;
        memset(&problem, 0, sizeof(problem));
        if (!read_problem(path, &problem)) {
            fprintf(stderr, "could not read %s\n", path);
            return 1;
        }
        ++problem_count;
        fprintf(stderr, "[%d] %s (n=%d, m=%d)\n", problem_count, name, n, m);
        for (size_t i = 0; i < sizeof(tolerances) / sizeof(tolerances[0]); ++i) {
            const Result result = solve_problem(&problem, &tolerances[i], repeats);
            fprintf(csv, "%s,%d,%d,%d,%s,%d,%.17g,%d,%d,%.17g,%.17g,%.17g\n",
                    name, n, m, posdef, tolerances[i].name, result.solved,
                    result.solve_time, result.exitflag, result.iterations,
                    result.primal, result.dual, result.gap);
            fflush(csv);
        }
        free_problem(&problem);
    }

    fclose(index);
    fclose(csv);
    return problem_count == 62 ? 0 : 1;
}
