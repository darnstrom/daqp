#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

// Equality reduction of a problem whose Hessian is given by its Cholesky
// factor (problem_type 2) gives the same reduction and solution as for H = R'R

namespace {

struct Problem {
    int n, m, neq;
    std::vector<c_float> R, H, f, A, bu, bl; // R: upper factor of H, packed by rows
    std::vector<int> sense;

    // n bounded variables, neq equalities and nineq inequalities, feasible at a
    // random point
    Problem(int n_, int neq_, int nineq, unsigned seed, bool diagonal = false)
        : n(n_), m(n_ + neq_ + nineq), neq(neq_), R(n_ * (n_ + 1) / 2, 0.0),
          H(n_ * n_, 0.0), f(n_), A((neq_ + nineq) * n_), bu(m), bl(m), sense(m, 0) {
        std::mt19937 rng(seed);
        std::normal_distribution<double> N(0, 1);
        auto r = [&](int i, int j) -> c_float& { return R[DAQP_R_OFFSET(i, n) + j]; };
        for (int i = 0; i < n; ++i) {
            r(i, i) = 1.0 + std::abs(N(rng));
            for (int j = i + 1; j < n; ++j) r(i, j) = diagonal ? 0.0 : 0.3 * N(rng);
        }
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k <= std::min(i, j); ++k) H[i * n + j] += r(k, i) * r(k, j);
        for (c_float& a : A) a = N(rng);
        for (c_float& fi : f) fi = N(rng);
        std::vector<c_float> x0(n);
        for (int i = 0; i < n; ++i) {
            x0[i] = N(rng);
            bu[i] = 2.0 + std::abs(x0[i]);
            bl[i] = -bu[i];
        }
        for (int k = 0; k < neq + nineq; ++k) {
            const c_float ax = std::inner_product(&A[k * n], &A[k * n] + n, x0.begin(), 0.0);
            bu[n + k] = (k < neq) ? ax : ax + 0.5;
            bl[n + k] = (k < neq) ? ax : -DAQP_INF;
            if (k < neq) sense[n + k] = DAQP_ACTIVE | DAQP_IMMUTABLE;
        }
    }

    DAQPProblem qp(bool factored) {
        return {n, m, n, factored ? R.data() : H.data(), f.data(), A.data(),
                bu.data(), bl.data(), sense.data(), nullptr, 0, factored ? 2 : 0};
    }
};

struct Solution {
    std::vector<c_float> x;
    c_float fval = 0;
    int exitflag = 0, path = -1, nz = -1, metric = -1;
};

Solution solve(DAQPWorkspace& work, int n) {
    Solution s;
    s.x.resize(n);
    DAQPResult res{};
    res.x = s.x.data();
    daqp_solve(&res, &work);
    s.exitflag = res.exitflag;
    s.fval = res.fval;
    if (DAQP_IS_REDUCED(&work)) {
        s.path = work.eq->path;
        s.nz = work.eq->nz;
        s.metric = work.eq->metric;
    }
    return s;
}

void setup(DAQPWorkspace& work, DAQPProblem& qp, int policy) {
    allocate_daqp_settings(&work);
    work.settings->eq_reduction = policy;
    assert(setup_daqp_main(&qp, &work, nullptr, 0) > 0);
}

void cleanup(DAQPWorkspace& work) {
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
}

Solution setup_solve(Problem& p, bool factored, int policy) {
    DAQPWorkspace work{};
    DAQPProblem qp = p.qp(factored);
    setup(work, qp, policy);
    Solution s = solve(work, p.n);
    cleanup(work);
    return s;
}

c_float max_diff(const std::vector<c_float>& a, const std::vector<c_float>& b) {
    c_float d = 0;
    for (size_t i = 0; i < a.size(); ++i) d = std::max(d, std::abs(a[i] - b[i]));
    return d;
}

// The factored reduction agrees with the reduction of H and with no reduction
void compare(Problem& p) {
    const Solution full = setup_solve(p, false, DAQP_EQ_REDUCTION_OFF);
    const Solution unf = setup_solve(p, false, DAQP_EQ_REDUCTION_ON);
    const Solution fac = setup_solve(p, true, DAQP_EQ_REDUCTION_ON);
    assert(full.exitflag > 0 && unf.exitflag > 0 && fac.exitflag > 0);
    assert(fac.path == unf.path && fac.path >= 0);
    assert(fac.nz == unf.nz && fac.metric == unf.metric);
    assert(max_diff(fac.x, unf.x) < 1e-8);
    assert(max_diff(fac.x, full.x) < 1e-8);
    assert(std::abs(fac.fval - full.fval) < 1e-8 * std::max<c_float>(1, std::abs(full.fval)));
}

} // namespace

int main() {
    Problem dense(40, 16, 20, 1);
    compare(dense);

    Problem diagonal(40, 16, 20, 2, true); // The factor is used as a metric
    compare(diagonal);
    assert(setup_solve(diagonal, true, DAQP_EQ_REDUCTION_ON).metric == 1);

    Problem determined(20, 20, 10, 3); // Nothing is left of the reduced problem
    compare(determined);
    assert(setup_solve(determined, true, DAQP_EQ_REDUCTION_ON).nz == 0);

    // Update the equalities and f of a reduced workspace (H x as R'(R x))
    Problem p(40, 16, 20, 4);
    DAQPWorkspace work{};
    DAQPProblem qp = p.qp(true);
    setup(work, qp, DAQP_EQ_REDUCTION_ON);
    std::mt19937 rng(5);
    std::normal_distribution<double> N(0, 1);
    for (int t = 0; t < 3; ++t) {
        std::vector<c_float> x1(p.n);
        for (c_float& x : x1) x = 0.5 * N(rng);
        for (int k = 0; k < p.m - p.n; ++k) { // x1 is feasible
            const c_float ax = std::inner_product(&p.A[k * p.n], &p.A[k * p.n] + p.n,
                                                  x1.begin(), 0.0);
            if (k < p.neq) p.bu[p.n + k] = p.bl[p.n + k] = ax;
            else p.bu[p.n + k] = std::max(p.bu[p.n + k], ax + 0.1);
        }
        for (c_float& fi : p.f) fi = N(rng);
        assert(daqp_update_ldp(DAQP_UPDATE_d | DAQP_UPDATE_v, &work, &qp) >= 0);
        const Solution s = solve(work, p.n);
        const Solution full = setup_solve(p, false, DAQP_EQ_REDUCTION_OFF);
        assert(s.exitflag > 0 && full.exitflag > 0 && s.path >= 0);
        assert(max_diff(s.x, full.x) < 1e-8);
    }
    cleanup(work);
}
