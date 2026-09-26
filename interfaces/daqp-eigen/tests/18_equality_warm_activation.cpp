#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>

// An equality that is linearly dependent on inequalities that are active from
// a warm start must not be dropped when the working set is formed: it is only
// redundant with respect to the other equalities. Here x1 <= 0.5, x2 <= 0.5
// are marked active, and x1 + x2 = 1 is dependent on them (and consistent at
// these bounds). Once the bounds are relaxed, the equality has to hold still.
//
// With soft constraints, more than n+1 constraints can be active. The second
// case marks four violated soft constraints active in a problem with n = 2 (so
// that they are parked beyond the first n+1 entries of the working set while
// the equality is activated), and compares the solution with the one from a
// cold start. The third case marks more constraints without a soft slack
// active than the working set can hold (four inequalities and an equality
// with n = 2), so that some of them are left out.

namespace {

void solve_soft(bool warm, c_float* x) {
    c_float H[4] = {1, 0, 0, 1};
    c_float f[2] = {0, 0};
    c_float A[10] = { 1,  0,   // x1 <= -1  (soft)
                      0,  1,   // x2 <= -1  (soft)
                     -1,  0,   // x1 >= 2   (soft)
                      0, -1,   // x2 >= 2   (soft)
                      1,  2};  // x1 + 2*x2 = 1
    c_float bu[5] = {-1, -1, -2, -2, 1};
    c_float bl[5] = {-DAQP_INF, -DAQP_INF, -DAQP_INF, -DAQP_INF, 1};
    const int active = warm ? DAQP_ACTIVE : 0;
    int sense[5] = {DAQP_SOFT | active, DAQP_SOFT | active,
                    DAQP_SOFT | active, DAQP_SOFT | active,
                    DAQP_ACTIVE | DAQP_IMMUTABLE};
    DAQPProblem qp = {2, 5, 0, H, f, A, bu, bl, sense, nullptr, 0, 0};
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.rho_soft = 1e-2; // Keep the violations moderate
    c_float lam[5];
    DAQPResult res{};
    res.x = x;
    res.lam = lam;
    daqp_quadprog(&res, &qp, &settings);
    assert(res.exitflag > 0);
    assert(std::abs(x[0] + 2 * x[1] - 1) < 1e-9); // The equality holds
}

void solve_hard(bool warm, c_float* x) {
    c_float H[4] = {1, 0, 0, 1};
    c_float f[2] = {-3, -3};
    c_float A[10] = {1, 0,   // x1 <= 1
                     0, 1,   // x2 <= 1
                     1, 1,   // x1 + x2 <= 1.5
                     1, 0,   // x1 <= 5
                     1, -1}; // x1 - x2 = 0
    c_float bu[5] = {1, 1, 1.5, 5, 0};
    c_float bl[5] = {-DAQP_INF, -DAQP_INF, -DAQP_INF, -DAQP_INF, 0};
    const int active = warm ? DAQP_ACTIVE : 0;
    int sense[5] = {active, active, active, active, DAQP_ACTIVE | DAQP_IMMUTABLE};
    DAQPProblem qp = {2, 5, 0, H, f, A, bu, bl, sense, nullptr, 0, 0};
    c_float lam[5];
    DAQPResult res{};
    res.x = x;
    res.lam = lam;
    daqp_quadprog(&res, &qp, nullptr);
    assert(res.exitflag == DAQP_EXIT_OPTIMAL);
    assert(std::abs(x[0] - x[1]) < 1e-9); // The equality holds
}

} // namespace

int main() {
    c_float H[4] = {1, 0, 0, 1};
    c_float f[2] = {-1, -1};
    c_float A[6] = {1, 0,
                    0, 1,
                    1, 1};
    c_float bu[3] = {0.5, 0.5, 1};
    c_float bl[3] = {-DAQP_INF, -DAQP_INF, 1};
    int sense[3] = {DAQP_ACTIVE, DAQP_ACTIVE, DAQP_ACTIVE | DAQP_IMMUTABLE};
    DAQPProblem qp = {2, 3, 0, H, f, A, bu, bl, sense, nullptr, 0, 0};

    DAQPWorkspace work{};
    assert(setup_daqp(&qp, &work, nullptr) >= 0);
    c_float x[2], lam[3];
    DAQPResult res{};
    res.x = x;
    res.lam = lam;
    daqp_solve(&res, &work);
    assert(res.exitflag == DAQP_EXIT_OPTIMAL);
    assert(std::abs(x[0] - 0.5) < 1e-9 && std::abs(x[1] - 0.5) < 1e-9);

    // Relax the inequalities: the solution is still (0.5, 0.5), on x1 + x2 = 1
    bu[0] = 0.7;
    bu[1] = 0.7;
    assert(daqp_update_ldp(DAQP_UPDATE_d, &work, &qp) >= 0);
    daqp_solve(&res, &work);
    assert(res.exitflag == DAQP_EXIT_OPTIMAL);
    assert(std::abs(x[0] + x[1] - 1) < 1e-9);
    assert(std::abs(x[0] - 0.5) < 1e-9 && std::abs(x[1] - 0.5) < 1e-9);

    free_daqp_workspace(&work);
    free_daqp_ldp(&work);

    c_float x_cold[2], x_warm[2];
    solve_soft(false, x_cold);
    solve_soft(true, x_warm);
    assert(std::abs(x_cold[0] - x_warm[0]) < 1e-8);
    assert(std::abs(x_cold[1] - x_warm[1]) < 1e-8);

    solve_hard(false, x_cold);
    solve_hard(true, x_warm);
    assert(std::abs(x_cold[0] - 0.75) < 1e-9 && std::abs(x_cold[1] - 0.75) < 1e-9);
    assert(std::abs(x_warm[0] - 0.75) < 1e-9 && std::abs(x_warm[1] - 0.75) < 1e-9);
}
