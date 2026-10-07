#include "api.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

// If branch and bound reaches the time limit after an integer-feasible
// solution has been found, the best such solution is returned, converted to
// the variables of the QP, with the exit flag DAQP_EXIT_TIMELIMIT_FEASIBLE.
// Without such a solution, the exit flag remains DAQP_EXIT_TIMELIMIT.
//
// A time limit of 1 ns is reached at the first check in the tree, after 32
// nodes, since the relaxations below take fewer than 32 iterations each.

namespace {

struct Problem {
    int n, m, nb;
    std::vector<c_float> H, f, A, bu, bl;
    std::vector<int> sense;
};

// n variables, of which the first nb are binary and the others are bounded by
// [-1, 1]. The Hessian is H = I + c*ones and the unconstrained minimizer is
// x = 0.5, so that every binary variable is fractional in the relaxation at the
// root. The general constraint sum(x) <= 0.4*n is active at the root.
Problem make_problem(int n, int nb, c_float c) {
    Problem p;
    p.n = n;
    p.m = n + 1;
    p.nb = nb;
    p.H.assign(n * n, c);
    for (int i = 0; i < n; i++) p.H[i * n + i] += 1;
    p.f.assign(n, 0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) p.f[i] -= 0.5 * p.H[i * n + j];
    p.A.assign(n, 1);
    p.bu.assign(p.m, 1);
    p.bl.assign(p.m, -1);
    for (int i = 0; i < nb; i++) p.bl[i] = 0;
    p.bu[n] = 0.4 * n;
    p.bl[n] = -DAQP_INF;
    p.sense.assign(p.m, 0);
    for (int i = 0; i < nb; i++) p.sense[i] = DAQP_BINARY;
    return p;
}

DAQPResult solve(const Problem& p, c_float time_limit, std::vector<c_float>& x) {
    std::vector<c_float> H = p.H, f = p.f, A = p.A, bu = p.bu, bl = p.bl;
    std::vector<int> sense = p.sense;
    DAQPProblem qp = {p.n, p.m, p.n, H.data(), f.data(), A.data(),
                      bu.data(), bl.data(), sense.data(), nullptr, 0, 0};
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.time_limit = time_limit;
    x.assign(p.n, 0);
    DAQPResult res{};
    res.x = x.data();
    res.lam = nullptr;
    daqp_quadprog(&res, &qp, &settings);
    return res;
}

c_float objective(const Problem& p, const std::vector<c_float>& x) {
    c_float J = 0;
    for (int i = 0; i < p.n; i++) {
        J += p.f[i] * x[i];
        for (int j = 0; j < p.n; j++) J += 0.5 * x[i] * p.H[i * p.n + j] * x[j];
    }
    return J;
}

// Binary variables at 0 or 1, and all constraints satisfied
bool is_integer_feasible(const Problem& p, const std::vector<c_float>& x, c_float tol) {
    c_float sum = 0;
    for (int i = 0; i < p.n; i++) {
        if (x[i] > p.bu[i] + tol || x[i] < p.bl[i] - tol) return false;
        if (i < p.nb && std::min(std::abs(x[i]), std::abs(x[i] - 1)) > tol) return false;
        sum += x[i];
    }
    return sum <= p.bu[p.n] + tol;
}

} // namespace

int main() {
#ifndef PROFILING
    return 0; // The time limit is only enforced with PROFILING
#else
    const c_float tol = 1e-6;
    std::vector<c_float> x;

    // An integer-feasible solution is found before the time limit
    Problem p = make_problem(12, 8, 0.2);
    DAQPResult ref = solve(p, 0, x);
    assert(ref.exitflag == DAQP_EXIT_OPTIMAL);
    assert(ref.nodes > 32); // The time limit below ends the exploration early

    DAQPResult res = solve(p, 1e-9, x);
    assert(res.exitflag == DAQP_EXIT_TIMELIMIT_FEASIBLE);
    assert(res.nodes <= 32);
    assert(is_integer_feasible(p, x, tol));
    assert(std::abs(res.fval - objective(p, x)) < tol * (1 + std::abs(res.fval)));
    assert(res.fval >= ref.fval - tol);

    // No integer-feasible solution is found before the time limit (the first
    // leaf of the tree is at depth 40)
    Problem q = make_problem(40, 40, 0);
    res = solve(q, 1e-9, x);
    assert(res.exitflag == DAQP_EXIT_TIMELIMIT);
    assert(res.nodes <= 32);

    return 0;
#endif
}
