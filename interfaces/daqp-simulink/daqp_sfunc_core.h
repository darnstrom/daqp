#ifndef DAQP_SFUNC_CORE_H
#define DAQP_SFUNC_CORE_H
/*
 * The solver logic of the DAQP S-function, free of Simulink dependencies
 * (so that it can be tested without MATLAB).
 *
 * Usage: create a solver with daqp_sfunc_new, write the problem data into
 * its qp (H, f, A, bupper, blower, sense, all row-major) before each call to
 * daqp_sfunc_solve, and read the solution from res.
 */
#include "api.h"
#include "utils.h"

// Warm start policies
#define DAQP_SFUNC_COLD 0       // Every solve starts from the sense input
#define DAQP_SFUNC_WARM 1       // Start from the working set of the previous solve
#define DAQP_SFUNC_WARM_FIXED 2 // As WARM, but H and A are only read at the first solve

typedef struct{
    DAQPProblem qp;     // The problem data (owned, written by the caller)
    DAQPResult res;     // The result of the latest solve
    DAQPWorkspace work;
    DAQPSettings settings;
    int warm_start;
    int is_setup;       // work is setup for qp
    int has_solved;     // The latest solve succeeded (=> its working set can be reused)
    int *sense_prev;    // The sense of the previous solve (to detect changes)
}DAQPSfunc;

// break_points (nh of them) are copied; nh <= 1 means a normal (non-hierarchical) QP
DAQPSfunc* daqp_sfunc_new(int n, int m, int ms, const int* break_points, int nh,
        int warm_start, int iter_limit);
// Solve the problem in s->qp; returns the exit flag (also in s->res.exitflag)
int daqp_sfunc_solve(DAQPSfunc* s);
void daqp_sfunc_free(DAQPSfunc* s);

#endif // DAQP_SFUNC_CORE_H
