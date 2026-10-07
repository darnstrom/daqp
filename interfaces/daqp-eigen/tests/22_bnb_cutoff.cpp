#include "api.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>

// If fval_bound prunes nodes of the branch and bound and no integer-feasible
// solution with a lower objective is found, the exit flag is DAQP_EXIT_CUTOFF.
// A problem without an integer-feasible solution remains DAQP_EXIT_INFEASIBLE,
// also with a finite fval_bound.
//
// The problems have f = 0, so that fval_bound refers to the objective
// 0.5*x'*H*x itself.

namespace {

// Four binary variables with H = I and the general constraint
// lower <= x1 + x2 + x3 + x4 <= upper
DAQPResult solve(c_float lower, c_float upper, c_float fval_bound, c_float* x) {
    c_float H[16] = {1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1};
    c_float f[4] = {0, 0, 0, 0};
    c_float A[4] = {1, 1, 1, 1};
    c_float bu[5] = {1, 1, 1, 1, upper};
    c_float bl[5] = {0, 0, 0, 0, lower};
    int sense[5] = {DAQP_BINARY, DAQP_BINARY, DAQP_BINARY, DAQP_BINARY, 0};
    DAQPProblem qp = {4, 5, 4, H, f, A, bu, bl, sense, nullptr, 0, 0};
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.fval_bound = fval_bound;
    DAQPResult res{};
    res.x = x;
    res.lam = nullptr;
    daqp_quadprog(&res, &qp, &settings);
    return res;
}

} // namespace

int main() {
    c_float x[4];

    // At least two of the binary variables are one, so the optimal value is 1
    DAQPResult res = solve(1.5, DAQP_INF, DAQP_INF, x);
    assert(res.exitflag == DAQP_EXIT_OPTIMAL);
    assert(std::abs(res.fval - 1) < 1e-9);

    res = solve(1.5, DAQP_INF, 1.1, x);
    assert(res.exitflag == DAQP_EXIT_OPTIMAL);
    assert(std::abs(res.fval - 1) < 1e-9);

    // The relaxation at the root (objective 0.28125) is below the bound, but
    // every integer-feasible solution is above it
    res = solve(1.5, DAQP_INF, 0.9, x);
    assert(res.exitflag == DAQP_EXIT_CUTOFF);
    assert(res.nodes > 1);

    // The sum of the binary variables cannot lie in [0.2, 0.8]
    res = solve(0.2, 0.8, DAQP_INF, x);
    assert(res.exitflag == DAQP_EXIT_INFEASIBLE);
    res = solve(0.2, 0.8, 10, x);
    assert(res.exitflag == DAQP_EXIT_INFEASIBLE);

    return 0;
}
