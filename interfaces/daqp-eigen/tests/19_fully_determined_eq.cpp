#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <vector>

// Equality constraints that determine all variables (neq == n): with the
// reduction, x is found from the QR of the equalities, and the other
// constraints only have to be consistent with it

namespace {

// QP dumped by PlaCo (dynamics/differential.py): 6 variables determined by 6
// ill-conditioned equalities (the last rows), and 8 inequalities. The LDL of
// A_E H^{-1} A_E' found it overdetermined (DAQP_EXIT_OVERDETERMINED_INITIAL)
const c_float trace_H[] = {1.0001229408273133e-08, 1.322387163385806e-12, 0.0, 0.0, 2.9410036845470272e-08, 1.9090783270497926e-08, 1.322387163385806e-12, 1.0001953143576155e-08, 0.0, 0.0, 1.9090783270497926e-08, 3.98583187085612e-08, 0.0, 0.0, 1.0000000100002615, 0.0, -8.086680895858504e-09, 1.617336179171701e-08, 0.0, 0.0, 0.0, 1.0000000100002615, -8.086680895858504e-09, -1.617336179171701e-08, 2.9410036845470272e-08, 1.9090783270497926e-08, -8.086680895858504e-09, -8.086680895858504e-09, 0.00150001, 0.0, 1.9090783270497926e-08, 3.98583187085612e-08, 1.617336179171701e-08, -1.617336179171701e-08, 0.0, 0.00300001};
const c_float trace_f[] = {4.016078732501674e-11, 8.406977455174777e-11, 1.5434634016451068e-10, -5.000000000230088, -5.219741072709436e-09, 2.111715351385761e-06};
const c_float trace_A[] = {-2.9410036827391805e-05, -1.9090783258762745e-05, 0.0, 0.0, -0.9999999993852959, 0.0, -1.9090783251854405e-05, -3.9858318669636696e-05, 0.0, 0.0, 0.0, -0.9999999990234283, 0.0, 0.0, -1.4465894554877226e-05, 0.0, 0.44721359545316547, -0.8944271909063309, 0.0, 0.0, 0.0, -1.4465894554877226e-05, 0.44721359545316547, 0.8944271909063309, 2.9410036827391805e-05, 1.9090783258762745e-05, 0.0, 0.0, 0.9999999993852959, 0.0, 1.9090783251854405e-05, 3.9858318669636696e-05, 0.0, 0.0, 0.0, 0.9999999990234283, 0.0, 0.0, 1.4465894554877226e-05, 0.0, -0.44721359545316547, 0.8944271909063309, 0.0, 0.0, 0.0, 1.4465894554877226e-05, -0.44721359545316547, -0.8944271909063309, -0.8164965809277261, 0.0, 0.4082482904638631, 0.4082482904638631, 0.0, 0.0, 0.0, -0.5773502691896258, -0.5773502691896258, 0.5773502691896258, 0.0, 0.0, 2.9410036827391805e-05, 1.9090783258762745e-05, 0.0, 0.0, 0.9999999993852959, 0.0, 1.9090783251854405e-05, 3.9858318669636696e-05, 0.0, 0.0, 0.0, 0.9999999990234283, 0.0, 0.0, 1.4465894554877226e-05, 0.0, -0.44721359545316547, 0.8944271909063309, 0.0, 0.0, 0.0, 1.4465894554877226e-05, -0.44721359545316547, -0.8944271909063309};
const c_float trace_bupper[] = {1e+30, 1e+30, 1e+30, 1e+30, 1e+30, 1e+30, 1e+30, 1e+30, 6.445223966399031e-13, 6.366882492203749e-13, 5.521974103878115e-05, -0.0020117153494210727, 0.0, -7.944854452274644e-35};
const c_float trace_blower[] = {-1.0000052191263655, -0.9978882836741048, -0.8944271909063309, -0.8944271909063309, -0.9999947796442265, -1.0021117143727516, -0.8944271909063309, -0.8944271909063309, 6.445223966399031e-13, 6.366882492203749e-13, 5.521974103878115e-05, -0.0020117153494210727, 0.0, -7.944854452274644e-35};
const c_float trace_qpmad_solution[] = {15.810979766546584, -48.25444485345284, 39.93820219327324, -8.316242660178489, 0.0005114333921021341, -0.0003902182973334334};

struct Problem {
    int n, m, ms;
    std::vector<c_float> H, f, A, bu, bl;
    std::vector<int> sense;
    DAQPProblem qp{};

    Problem(int n_, int m_, int ms_)
        : n(n_), m(m_), ms(ms_), H(n * n, 0.0), f(n, 0.0),
          A((m - ms) * n, 0.0), bu(m, DAQP_INF), bl(m, -DAQP_INF), sense(m, 0) {}

    void finalize() {
        qp = {n, m, ms, H.data(), f.data(), A.data(), bu.data(), bl.data(),
              sense.data(), nullptr, 0, 0};
    }
    // Row of constraint i (a unit vector for a simple bound)
    c_float row_times(int i, const c_float* x) const {
        if (i < ms) return x[i];
        c_float s = 0;
        for (int j = 0; j < n; ++j) s += A[(i - ms) * n + j] * x[j];
        return s;
    }
    // ||H x + f + [I 0; A]' lam||_inf
    c_float stationarity(const c_float* x, const c_float* lam) const {
        c_float res = 0;
        for (int j = 0; j < n; ++j) {
            c_float g = f[j] + (j < ms ? lam[j] : 0);
            for (int k = 0; k < n; ++k) g += H[j * n + k] * x[k];
            for (int i = ms; i < m; ++i) g += A[(i - ms) * n + j] * lam[i];
            res = std::fmax(res, std::fabs(g));
        }
        return res;
    }
};

// n variables, determined by the equalities E x = b (E dense and well
// conditioned), and a band of inequalities -10 <= x_i + x_{i+1} <= 10.
// With diagonal, H is diagonal (the QR is then formed in the metric of H).
Problem determined(int n, bool diagonal, int ms = 0) {
    const int neq = n, nineq = n - 1;
    Problem p(n, ms + neq + nineq, ms);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            p.H[i * n + j] = i == j ? 2.0 + i : (diagonal ? 0.0 : 0.1);
        p.f[i] = 1.0 - 0.3 * i;
    }
    for (int i = 0; i < neq; ++i) {
        for (int j = 0; j < n; ++j)
            p.A[i * n + j] = i == j ? 2.0 : 0.1 * ((i + 2 * j) % 5 - 2);
        p.bu[ms + i] = p.bl[ms + i] = 0.5 * (i % 3) - 0.4;
    }
    for (int i = 0; i < nineq; ++i) {
        p.A[(neq + i) * n + i] = p.A[(neq + i) * n + i + 1] = 1.0;
        p.bu[ms + neq + i] = 10.0;
        p.bl[ms + neq + i] = -10.0;
    }
    for (int i = 0; i < ms; ++i) {
        p.bu[i] = 5.0;
        p.bl[i] = -5.0;
    }
    p.finalize();
    return p;
}

struct Solution {
    std::vector<c_float> x, lam;
    DAQPResult res{};
    Solution(const Problem& p) : x(p.n), lam(p.m) {
        res.x = x.data();
        res.lam = lam.data();
    }
};

Solution quadprog(Problem& p, int policy) {
    Solution s(p);
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.eq_reduction = policy;
    daqp_quadprog(&s.res, &p.qp, &settings);
    return s;
}

// The solution of a well-conditioned problem is the same with and without the
// reduction, and satisfies the KKT conditions
void check_against_full(Problem& p) {
    Solution red = quadprog(p, DAQP_EQ_REDUCTION_ON);
    Solution full = quadprog(p, DAQP_EQ_REDUCTION_OFF);
    assert(red.res.exitflag == DAQP_EXIT_OPTIMAL);
    assert(full.res.exitflag == DAQP_EXIT_OPTIMAL);
    assert(std::fabs(red.res.fval - full.res.fval) < 1e-9 * (1 + std::fabs(full.res.fval)));
    for (int j = 0; j < p.n; ++j) assert(std::fabs(red.x[j] - full.x[j]) < 1e-9);
    for (int i = 0; i < p.m; ++i) assert(std::fabs(red.lam[i] - full.lam[i]) < 1e-8);
    assert(p.stationarity(red.x.data(), red.lam.data()) < 1e-9);
}

} // namespace

int main() {
    // The PlaCo trace: solved, with the solution that qpmad found
    {
        Problem p(6, 14, 0);
        p.H.assign(trace_H, trace_H + 36);
        p.f.assign(trace_f, trace_f + 6);
        p.A.assign(trace_A, trace_A + 14 * 6);
        p.bu.assign(trace_bupper, trace_bupper + 14);
        p.bl.assign(trace_blower, trace_blower + 14);
        p.finalize();
        Solution s = quadprog(p, DAQP_EQ_REDUCTION_ON);
        assert(s.res.exitflag == DAQP_EXIT_OPTIMAL);
        for (int j = 0; j < 6; ++j)
            assert(std::fabs(s.x[j] - trace_qpmad_solution[j]) < 1e-8 * (1 + std::fabs(trace_qpmad_solution[j])));
        assert(std::fabs(s.res.fval - 873.6911777271131) < 1e-8 * 873.7);
        for (int i = 0; i < 8; ++i) assert(s.lam[i] == 0); // Inactive inequalities
        assert(p.stationarity(s.x.data(), s.lam.data()) < 1e-8);
    }

    // Dense and diagonal Hessians, with and without simple bounds
    for (bool diagonal : {false, true})
        for (int ms : {0, 8}) {
            Problem p = determined(8, diagonal, ms);
            check_against_full(p);
        }

    // An inequality that the determined x violates: infeasible
    for (bool diagonal : {false, true}) {
        Problem p = determined(8, diagonal);
        Solution s = quadprog(p, DAQP_EQ_REDUCTION_ON);
        const c_float v = p.row_times(p.m - 1, s.x.data());
        p.bu[p.m - 1] = v - 1e-3;
        p.bl[p.m - 1] = v - 2e-3;
        assert(quadprog(p, DAQP_EQ_REDUCTION_ON).res.exitflag == DAQP_EXIT_INFEASIBLE);
        assert(quadprog(p, DAQP_EQ_REDUCTION_OFF).res.exitflag == DAQP_EXIT_INFEASIBLE);
    }

    // A simple bound that the determined x violates: infeasible
    {
        Problem p = determined(8, false, 8);
        Solution s = quadprog(p, DAQP_EQ_REDUCTION_ON);
        p.bu[3] = s.x[3] - 1e-3;
        assert(quadprog(p, DAQP_EQ_REDUCTION_ON).res.exitflag == DAQP_EXIT_INFEASIBLE);
    }

    // More equalities than variables: a consistent extra equality is kept as
    // a constraint (and holds), an inconsistent one is overdetermined
    for (c_float offset : {0.0, 1e-2}) {
        Problem base = determined(8, false);
        Solution s = quadprog(base, DAQP_EQ_REDUCTION_ON);
        Problem p(8, base.m + 1, 0);
        p.H = base.H;
        p.f = base.f;
        p.A = base.A;
        p.bu = base.bu;
        p.bl = base.bl;
        // x_0 + x_1 = its value at the solution (+ offset)
        p.A.insert(p.A.begin() + 8 * 8, 8, 0.0);
        p.A[8 * 8] = p.A[8 * 8 + 1] = 1.0;
        p.bu.insert(p.bu.begin() + 8, s.x[0] + s.x[1] + offset);
        p.bl.insert(p.bl.begin() + 8, s.x[0] + s.x[1] + offset);
        p.finalize();
        Solution r = quadprog(p, DAQP_EQ_REDUCTION_ON);
        if (offset == 0) {
            assert(r.res.exitflag == DAQP_EXIT_OPTIMAL);
            for (int j = 0; j < 8; ++j) assert(std::fabs(r.x[j] - s.x[j]) < 1e-9);
            assert(p.stationarity(r.x.data(), r.lam.data()) < 1e-9);
        }
        else assert(r.res.exitflag == DAQP_EXIT_OVERDETERMINED_INITIAL);
    }

    // A soft constraint that the equalities determine (all of them for neq ==
    // n, one in the span of the equalities for neq < n) is not reduced away,
    // which would make it hard: a violation gives a soft optimum, as without
    // the reduction
    for (int neq : {3, 2}) {
        Problem p(3, neq + 1, 0);
        for (int i = 0; i < 3; ++i) p.H[i * 3 + i] = 1;
        for (int i = 0; i < neq; ++i) {
            p.A[i * 3 + i] = 1;
            p.bu[i] = p.bl[i] = 1;
        }
        p.A[neq * 3] = p.A[neq * 3 + 1] = 1; // x_0 + x_1 <= 1, violated by 1
        p.bu[neq] = 1;
        p.sense[neq] = DAQP_SOFT;
        p.finalize();
        Solution red = quadprog(p, DAQP_EQ_REDUCTION_ON);
        Solution full = quadprog(p, DAQP_EQ_REDUCTION_OFF);
        assert(full.res.exitflag == DAQP_EXIT_SOFT_OPTIMAL);
        assert(red.res.exitflag == DAQP_EXIT_SOFT_OPTIMAL);
        assert(std::fabs(red.res.soft_slack - full.res.soft_slack) < 1e-9);
        for (int j = 0; j < 3; ++j) assert(std::fabs(red.x[j] - full.x[j]) < 1e-9);
    }

    // A workspace that is updated and solved repeatedly: the right-hand side
    // of the equalities changes (only one of them nonzero, which forms the
    // responses to b_E), and the solution follows
    for (int mask : {int(DAQP_UPDATE_d),
                     int(DAQP_UPDATE_Rinv | DAQP_UPDATE_M | DAQP_UPDATE_v |
                         DAQP_UPDATE_d | DAQP_UPDATE_sense)}) {
        Problem p = determined(8, false);
        for (int i = 0; i < 8; ++i) p.bu[i] = p.bl[i] = 0;
        DAQPWorkspace work{};
        allocate_daqp_settings(&work);
        work.settings->eq_reduction = DAQP_EQ_REDUCTION_ON;
        assert(setup_daqp_main(&p.qp, &work, nullptr, 0) > 0);
        assert(DAQP_IS_REDUCED(&work));
        Solution s(p);
        for (int k = 0; k < 4; ++k) {
            p.bu[2] = p.bl[2] = 0.5 * k;
            assert(daqp_update_ldp(mask, &work, &p.qp) >= 0);
            assert(DAQP_IS_REDUCED(&work));
            daqp_solve(&s.res, &work);
            Solution full = quadprog(p, DAQP_EQ_REDUCTION_OFF);
            assert(s.res.exitflag == DAQP_EXIT_OPTIMAL);
            assert(std::fabs(s.res.fval - full.res.fval) < 1e-9 * (1 + std::fabs(full.res.fval)));
            for (int j = 0; j < 8; ++j) assert(std::fabs(s.x[j] - full.x[j]) < 1e-9);
            assert(p.stationarity(s.x.data(), s.lam.data()) < 1e-9);
        }
        free_daqp_workspace(&work);
        free_daqp_ldp(&work);
    }
}
