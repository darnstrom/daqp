#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

// Equality reduction of a problem whose Hessian is given by its Cholesky
// factor (problem_type 2): the reduction is formed from the factor and gives
// the solution of the reduction of the Hessian itself, with the same path

namespace {

struct Problem {
    int n, m, ms, neq;
    std::vector<c_float> R, H, f, A, bu, bl; // R: upper factor of H = R'R, packed by rows
    std::vector<int> sense;

    c_float& r(int i, int j) { return R[DAQP_R_OFFSET(i, n) + j]; }

    /*
     * n variables with bounds, neq equalities (of which ndep are dependent on
     * the ones before) and nineq inequalities, feasible at a random point. With
     * diagonal, R is diagonal; with tiny_pivot, R_{n-1,n-1} = tiny_pivot and
     * the near-null direction R^{-1} e_{n-1} of R is in the null space of the
     * equalities.
     */
    Problem(int n_, int neq_, int nineq, unsigned seed, bool diagonal = false,
            int ndep = 0, c_float tiny_pivot = 0, int nbin = 0)
        : n(n_), m(n_ + neq_ + nineq), ms(n_), neq(neq_),
          R(n_ * (n_ + 1) / 2, 0.0), H(n_ * n_, 0.0), f(n_), A((neq_ + nineq) * n_),
          bu(m), bl(m), sense(m, 0) {
        std::mt19937 rng(seed);
        std::normal_distribution<double> N(0, 1);
        for (int i = 0; i < n; ++i) {
            r(i, i) = 1.0 + std::abs(N(rng));
            for (int j = i + 1; j < n; ++j) r(i, j) = diagonal ? 0.0 : 0.3 * N(rng);
        }
        if (tiny_pivot > 0) r(n - 1, n - 1) = tiny_pivot;
        for (int i = 0; i < n; ++i) // H = R'R
            for (int j = 0; j < n; ++j) {
                c_float s = 0;
                for (int k = 0; k <= std::min(i, j); ++k) s += r(k, i) * r(k, j);
                H[i * n + j] = s;
            }
        // Near-null direction v = R^{-1} e_{n-1} (back substitution)
        std::vector<c_float> v(n, 0.0);
        if (tiny_pivot > 0)
            for (int i = n - 1; i >= 0; --i) {
                c_float s = (i == n - 1) ? 1.0 : 0.0;
                for (int j = i + 1; j < n; ++j) s -= r(i, j) * v[j];
                v[i] = s / r(i, i);
            }
        const c_float vv = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
        for (c_float& a : A) a = N(rng);
        for (int k = 0; k < neq; ++k) {
            c_float* a = &A[k * n];
            if (k >= neq - ndep) { // A duplicate or a combination of two rows before
                const c_float* a1 = &A[(k - neq + ndep) * n];
                const c_float* a2 = a1 + n;
                for (int j = 0; j < n; ++j) a[j] = (k % 2) ? a1[j] : a1[j] + 2 * a2[j];
            }
            else if (vv > 0) { // v in the null space
                const c_float av = std::inner_product(a, a + n, v.begin(), 0.0);
                for (int j = 0; j < n; ++j) a[j] -= av / vv * v[j];
            }
        }
        std::vector<c_float> x0(n);
        for (c_float& x : x0) x = N(rng);
        for (int i = 0; i < nbin; ++i) x0[i] = (x0[i] > 0) ? 1.0 : 0.0;
        for (int i = 0; i < n; ++i) {
            bu[i] = 2.0 + std::abs(x0[i]);
            bl[i] = -bu[i];
        }
        for (int k = 0; k < neq + nineq; ++k) {
            const c_float ax = std::inner_product(&A[k * n], &A[k * n] + n, x0.begin(), 0.0);
            bu[n + k] = (k < neq) ? ax : ax + 0.5 * std::abs(N(rng));
            bl[n + k] = (k < neq) ? ax : -DAQP_INF;
            if (k < neq) sense[n + k] = DAQP_ACTIVE | DAQP_IMMUTABLE;
        }
        for (int i = 0; i < nbin; ++i) {
            sense[i] = DAQP_BINARY;
            bu[i] = 1.0;
            bl[i] = 0.0;
        }
        for (int i = 0; i < n; ++i) f[i] = N(rng) * r(i, i) * r(i, i);
    }

    DAQPProblem qp(bool factored) {
        return {n, m, ms, factored ? R.data() : H.data(), f.data(), A.data(),
                bu.data(), bl.data(), sense.data(), nullptr, 0, factored ? 2 : 0};
    }
};

struct Solution {
    std::vector<c_float> x, lam;
    c_float fval = 0;
    int exitflag = 0;
    bool reduced = false;
    int path = -1, neq = 0, nz = 0, metric = 0, n_prox = -1;
};

void cleanup(DAQPWorkspace& work) {
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
}

Solution read_reduction(DAQPWorkspace& work, Solution s) {
    s.reduced = DAQP_IS_REDUCED(&work);
    if (s.reduced) {
        s.path = work.eq->path;
        s.neq = work.eq->neq;
        s.nz = work.eq->nz;
        s.metric = work.eq->metric;
        s.n_prox = work.eq->other.n_prox; // Of the reduced problem
    }
    return s;
}

Solution solve(DAQPWorkspace& work, Problem& p) {
    Solution s;
    s.x.resize(p.n);
    s.lam.resize(p.m);
    DAQPResult res{};
    res.x = s.x.data();
    res.lam = s.lam.data();
    daqp_solve(&res, &work);
    s.exitflag = res.exitflag;
    s.fval = res.fval;
    return read_reduction(work, s);
}

Solution setup_solve(Problem& p, bool factored, int policy, int mask = 0,
                     c_float eps_prox = DAQP_DEFAULT_EPS_PROX) {
    DAQPWorkspace work{};
    allocate_daqp_settings(&work);
    work.settings->eq_reduction = policy;
    work.settings->eps_prox = eps_prox;
    DAQPProblem qp = p.qp(factored);
    assert(setup_daqp_main(&qp, &work, nullptr, mask) > 0);
    Solution s = solve(work, p);
    cleanup(work);
    return s;
}

c_float max_diff(const std::vector<c_float>& a, const std::vector<c_float>& b) {
    c_float d = 0;
    for (size_t i = 0; i < a.size(); ++i) d = std::max(d, std::abs(a[i] - b[i]));
    return d;
}

// The factored reduction against the reduction of H, and for a well-posed
// problem also against the solve without the reduction
void compare(Problem& p, c_float tol, int path, bool well_posed = true) {
    const Solution unf = setup_solve(p, false, DAQP_EQ_REDUCTION_ON);
    const Solution fac = setup_solve(p, true, DAQP_EQ_REDUCTION_ON);
    assert(unf.exitflag > 0 && fac.exitflag > 0);
    assert(fac.reduced && unf.reduced);
    assert(fac.path == path && unf.path == path);
    assert(fac.neq == unf.neq && fac.nz == unf.nz && fac.metric == unf.metric);
    const c_float scale = std::max<c_float>(1, std::abs(unf.fval));
    assert(std::abs(fac.fval - unf.fval) < tol * scale);
    if (well_posed) {
        const Solution full = setup_solve(p, false, DAQP_EQ_REDUCTION_OFF);
        assert(full.exitflag > 0);
        assert(std::abs(fac.fval - full.fval) < tol * scale);
        assert(max_diff(fac.x, unf.x) < tol);
        assert(max_diff(fac.x, full.x) < tol);
        assert(max_diff(fac.lam, unf.lam) < 1e3 * tol);
    }
    // The equality constraints hold
    for (int k = 0; k < p.neq; ++k) {
        const c_float ax = std::inner_product(&p.A[k * p.n], &p.A[k * p.n] + p.n,
                                              fac.x.begin(), 0.0);
        assert(std::abs(ax - p.bu[p.n + k]) < 1e-9);
    }
}

} // namespace

int main() {
    // A dense factor: PATH_LDP, with the same solution and multipliers
    Problem dense(40, 16, 20, 1);
    compare(dense, 1e-8, DAQP_EQ_PATH_LDP);

    // AUTO reduces a factored one-shot problem where it reduces the Hessian
    // (and the heuristics for a diagonal Hessian apply to a diagonal factor)
    Problem diag_few(40, 6, 20, 2, true);
    for (Problem* p : {&dense, &diag_few}) {
        const bool unf = setup_solve(*p, false, DAQP_EQ_REDUCTION_AUTO, DAQP_UPDATE_eliminate).reduced;
        const bool fac = setup_solve(*p, true, DAQP_EQ_REDUCTION_AUTO, DAQP_UPDATE_eliminate).reduced;
        assert(fac == unf);
        assert(fac == (p == &dense));
        // A workspace that is to be updated is only reduced with ON
        assert(!setup_solve(*p, true, DAQP_EQ_REDUCTION_AUTO).reduced);
        assert(!setup_solve(*p, true, DAQP_EQ_REDUCTION_OFF, DAQP_UPDATE_eliminate).reduced);
    }

    // Dependent equalities are kept as constraints
    Problem dependent(40, 16, 20, 3, false, 4);
    compare(dependent, 1e-8, DAQP_EQ_PATH_LDP);
    assert(setup_solve(dependent, true, DAQP_EQ_REDUCTION_ON).neq == 12);

    // Equalities that determine x (nothing is left of the reduced problem)
    Problem determined(20, 20, 10, 4);
    compare(determined, 1e-8, DAQP_EQ_PATH_LDP);
    assert(setup_solve(determined, true, DAQP_EQ_REDUCTION_ON).nz == 0);

    // A diagonal factor is used as a metric
    Problem diagonal(40, 16, 20, 5, true);
    compare(diagonal, 1e-8, DAQP_EQ_PATH_LDP);
    assert(setup_solve(diagonal, true, DAQP_EQ_REDUCTION_ON).metric == 1);

    // A nearly singular factor whose near-null direction is in the null space
    // of the equalities: Z'HZ is singular in the sense of eq_chol, and the
    // reduced problem (PATH_QP) is solved by the proximal method. The minimizer
    // is not well determined, so only the objective function values of the two
    // reductions are compared
    Problem singular(30, 10, 20, 6, false, 0, 1e-9);
    compare(singular, 1e-6, DAQP_EQ_PATH_QP, false);

    // Binary variables (branch and bound in the reduced problem)
    Problem binaries(30, 10, 20, 7, false, 0, 0, 6);
    compare(binaries, 1e-8, DAQP_EQ_PATH_LDP);

    // A supplied factor is not regularized: also the reduced problem is posed
    // with a factor, so eps_prox > 0 does not select the proximal method
    {
        const Solution full = setup_solve(dense, true, DAQP_EQ_REDUCTION_OFF);
        const Solution fac = setup_solve(dense, true, DAQP_EQ_REDUCTION_ON, 0, 1e-3);
        assert(fac.reduced && fac.n_prox == 0 && fac.exitflag > 0);
        assert(max_diff(fac.x, full.x) < 1e-8);
    }

    // Updates of the bounds and of f of a reduced workspace: the right-hand
    // side is reduced with the factor (H x is formed as R'(R x))
    for (int mask : {int(DAQP_UPDATE_d), int(DAQP_UPDATE_d | DAQP_UPDATE_v)}) {
        Problem p(40, 16, 20, 8);
        DAQPWorkspace work[2] = {};
        DAQPProblem qp[2] = {p.qp(false), p.qp(true)};
        for (int k = 0; k < 2; ++k) {
            allocate_daqp_settings(&work[k]);
            work[k].settings->eq_reduction = DAQP_EQ_REDUCTION_ON;
            assert(setup_daqp_main(&qp[k], &work[k], nullptr, 0) > 0);
            assert(solve(work[k], p).exitflag > 0);
        }
        std::mt19937 rng(9);
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
            if (mask & DAQP_UPDATE_v)
                for (c_float& fi : p.f) fi = N(rng);
            Solution s[2];
            for (int k = 0; k < 2; ++k) {
                assert(daqp_update_ldp(mask, &work[k], &qp[k]) >= 0);
                s[k] = solve(work[k], p);
                assert(s[k].exitflag > 0 && s[k].reduced);
            }
            const Solution full = setup_solve(p, false, DAQP_EQ_REDUCTION_OFF);
            assert(max_diff(s[1].x, s[0].x) < 1e-8);
            assert(max_diff(s[1].x, full.x) < 1e-8);
            assert(std::abs(s[1].fval - full.fval) < 1e-8 * std::max<c_float>(1, std::abs(full.fval)));
        }
        for (int k = 0; k < 2; ++k) cleanup(work[k]);
    }
}
