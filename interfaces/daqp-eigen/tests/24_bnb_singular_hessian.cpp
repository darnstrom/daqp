#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <vector>

// Branch and bound with a singular Hessian (regularized, see daqp_bnb_prox)

namespace {

const int policies[3] = {DAQP_EQ_REDUCTION_OFF, DAQP_EQ_REDUCTION_AUTO, DAQP_EQ_REDUCTION_ON};

struct Problem {
    int n, m, ms;
    std::vector<c_float> H, f, A, bu, bl;
    std::vector<int> sense;
    bool quadratic;
    Problem(int n_, int mg, bool quad)
        : n(n_), m(n_ + mg), ms(n_), H(quad ? n_ * n_ : 0, 0.0), f(n_, 0.0),
          A(mg * n_, 0.0), bu(m, DAQP_INF), bl(m, -DAQP_INF), sense(m, 0), quadratic(quad) {}
    DAQPProblem qp() {
        return {n, m, ms, quadratic ? H.data() : nullptr, f.data(), A.data(),
                bu.data(), bl.data(), sense.data(), nullptr, 0, 0};
    }
    c_float a(int r, int j) const { return A[r * n + j]; }
    c_float& a(int r, int j) { return A[r * n + j]; }
    c_float objective(const std::vector<c_float>& x) const {
        c_float J = 0;
        for (int i = 0; i < n; ++i) {
            J += f[i] * x[i];
            if (quadratic) for (int j = 0; j < n; ++j) J += 0.5 * x[i] * H[i * n + j] * x[j];
        }
        return J;
    }
    c_float violation(const std::vector<c_float>& x) const {
        c_float v = 0;
        for (int i = 0; i < m; ++i) {
            c_float ax = 0;
            if (i < ms) ax = x[i];
            else for (int j = 0; j < n; ++j) ax += a(i - ms, j) * x[j];
            v = std::fmax(v, std::fmax(ax - bu[i], bl[i] - ax));
            if (sense[i] & DAQP_BINARY)
                v = std::fmax(v, std::fmin(std::fabs(ax - bl[i]), std::fabs(ax - bu[i])));
        }
        return v;
    }
};

struct Solution {
    int exitflag;
    c_float fval;
    std::vector<c_float> x;
};

Solution solve(Problem p, int policy, c_float rel_subopt = 0) {
    Solution s{0, 0, std::vector<c_float>(p.n, 0)};
    std::vector<c_float> lam(p.m);
    DAQPResult res{};
    res.x = s.x.data();
    res.lam = lam.data();
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.eq_reduction = policy;
    settings.rel_subopt = rel_subopt;
    DAQPProblem qp = p.qp();
    daqp_quadprog(&res, &qp, &settings);
    s.exitflag = res.exitflag;
    s.fval = res.fval;
    return s;
}

// Fixed charge MILP: max profit of y_i <= M_i b_i (y_i in [0,U]) with cost F_i b_i and
// sum y >= D. Since the profit is positive, y_i = min(M_i,U) for open i
Problem fixed_charge(int k, std::vector<c_float>& cap, std::vector<c_float>& fixed, c_float& D) {
    Problem p(2 * k, k + 1, false);
    const c_float U = 100;
    cap.assign(k, 0);
    fixed.assign(k, 0);
    for (int i = 0; i < k; ++i) {
        const c_float M = U * (0.6 + 0.2 * i), profit = 0.4 + 0.15 * ((3 * i) % k);
        fixed[i] = U * (0.5 + 0.1 * ((5 * i) % k));
        cap[i] = std::fmin(M, U);
        p.f[i] = fixed[i];
        p.f[k + i] = -profit;
        p.sense[i] = DAQP_BINARY;
        p.bl[i] = 0; p.bu[i] = 1;
        p.bl[k + i] = 0; p.bu[k + i] = U;
        p.a(i, k + i) = 1; p.a(i, i) = -M; p.bu[p.ms + i] = 0; // y_i - M_i b_i <= 0
        p.a(k, k + i) = -1;                                    // -sum y <= -D
    }
    D = 1.5 * U;
    p.bu[p.ms + k] = -D;
    return p;
}

// Hybrid MPC of a scalar piecewise-affine system in MLD form (mode d_t = 1 iff x_t >= 0),
// with loose big-M constants: x_{t+1} = a1 x_t + b1 u_t + z_t, z_t = d_t*(da x_t + db u_t).
// Variables per step: d_t (binary), u_t, z_t, x_{t+1}. The binaries (and z unless costed) have
// no quadratic cost
Problem hybrid_mpc(int N, c_float x0, c_float cost_d) {
    Problem p(4 * N, 7 * N, true);
    const c_float a1 = 0.75, b1 = 1, a2 = 1.25, b2 = 0.5, da = a2 - a1, db = b2 - b1;
    const c_float X = 10, Umax = 5, Mz = 10 * (std::fabs(da) * X + std::fabs(db) * Umax), Mx = 10 * X;
    for (int t = 0; t < N; ++t) {
        const int d = t, u = N + t, z = 2 * N + t, x = 3 * N + t, xp = t > 0 ? 3 * N + t - 1 : -1;
        p.sense[d] = DAQP_BINARY;
        p.bl[d] = 0; p.bu[d] = 1;
        p.bl[u] = -Umax; p.bu[u] = Umax;
        p.bl[z] = -Mz; p.bu[z] = Mz;
        p.bl[x] = -X; p.bu[x] = X;
        p.H[x * p.n + x] = 1; p.H[u * p.n + u] = 0.1; p.H[d * p.n + d] = cost_d;
        const int r = 7 * t, R = p.ms + r;
        const c_float c = xp >= 0 ? 0 : x0; // Constant part of x_t
        p.a(r, x) = 1; p.a(r, u) = -b1; p.a(r, z) = -1; // Dynamics (equality)
        if (xp >= 0) p.a(r, xp) = -a1;
        p.bu[R] = p.bl[R] = xp >= 0 ? 0 : a1 * x0;
        p.a(r + 1, d) = -Mx; if (xp >= 0) p.a(r + 1, xp) = 1; p.bu[R + 1] = -c;     // x_t <= Mx d
        p.a(r + 2, d) = Mx;  if (xp >= 0) p.a(r + 2, xp) = -1; p.bu[R + 2] = Mx + c; // -x_t <= Mx (1-d)
        p.a(r + 3, z) = 1;  p.a(r + 3, d) = -Mz; p.bu[R + 3] = 0;                   // |z| <= Mz d
        p.a(r + 4, z) = -1; p.a(r + 4, d) = -Mz; p.bu[R + 4] = 0;
        p.a(r + 5, z) = 1;  p.a(r + 5, u) = -db; p.a(r + 5, d) = Mz;                // |z - (da x + db u)| <= Mz (1-d)
        if (xp >= 0) p.a(r + 5, xp) = -da;
        p.bu[R + 5] = Mz + (xp >= 0 ? 0 : da * x0);
        p.a(r + 6, z) = -1; p.a(r + 6, u) = db; p.a(r + 6, d) = Mz;
        if (xp >= 0) p.a(r + 6, xp) = da;
        p.bu[R + 6] = Mz - (xp >= 0 ? 0 : da * x0);
    }
    return p;
}

bool close(c_float a, c_float b, c_float tol) { return std::fabs(a - b) <= tol * (1 + std::fabs(b)); }

} // namespace

int main() {
    // MILP: correct assignment, exitflag 4 unless the bias is within rel_subopt
    {
        const int k = 5;
        std::vector<c_float> cap, fixed;
        c_float D;
        Problem p = fixed_charge(k, cap, fixed, D);
        c_float J_best = DAQP_INF;
        for (int S = 0; S < (1 << k); ++S) { // Enumerate the open facilities
            c_float J = 0, supply = 0;
            for (int i = 0; i < k; ++i)
                if (S >> i & 1) { J += p.f[i] + p.f[k + i] * cap[i]; supply += cap[i]; }
            if (supply >= D && J < J_best) J_best = J;
        }
        for (int policy : policies) {
            Solution s = solve(p, policy);
            assert(s.exitflag == DAQP_EXIT_OPTIMAL_INEXACT);
            assert(close(s.fval, J_best, 1e-6) && close(p.objective(s.x), J_best, 1e-6));
            assert(p.violation(s.x) < 1e-6);
            assert(solve(p, policy, 1e-3).exitflag == DAQP_EXIT_OPTIMAL);
        }
    }

    // Hybrid MPC with zero-cost binaries and equalities (reducing them gave false infeasibility)
    for (c_float x0 : {-6.0, -2.5, 0.5, 3.0, 7.5}) {
        Problem p = hybrid_mpc(6, x0, 0);
        Solution s_off = solve(p, DAQP_EQ_REDUCTION_OFF);
        assert(s_off.exitflag == DAQP_EXIT_OPTIMAL || s_off.exitflag == DAQP_EXIT_OPTIMAL_INEXACT);
        for (int policy : policies) {
            Solution s = solve(p, policy);
            assert(s.exitflag == s_off.exitflag);
            assert(p.violation(s.x) < 1e-5);
            assert(close(s.fval, s_off.fval, 1e-6) && close(p.objective(s.x), s.fval, 1e-6));
        }
    }

    // Positive definite after eliminating the equalities (only z is singular): solved exactly
    for (c_float x0 : {-6.0, 0.5, 7.5}) {
        Problem p = hybrid_mpc(6, x0, 1e-2);
        Solution s_off = solve(p, DAQP_EQ_REDUCTION_OFF);
        for (int policy : {DAQP_EQ_REDUCTION_AUTO, DAQP_EQ_REDUCTION_ON}) {
            Solution s = solve(p, policy);
            assert(s.exitflag == DAQP_EXIT_OPTIMAL);
            assert(close(s.fval, s_off.fval, 1e-5) && p.violation(s.x) < 1e-5);
        }
    }

    // Repeated updates of a workspace agree with cold solves (the reduction that is tried
    // with DAQP_EQ_REDUCTION_ON must not change the full problem)
    for (c_float cost_d : {0.0, 1e-2}) {
        for (int policy : policies) {
            Problem p = hybrid_mpc(6, 0, cost_d);
            DAQPProblem qp = p.qp();
            DAQPWorkspace work{};
            DAQPSettings settings;
            daqp_default_settings(&settings);
            settings.eq_reduction = policy;
            work.settings = &settings;
            assert(setup_daqp(&qp, &work, nullptr) >= 0);
            std::vector<c_float> x(p.n), lam(p.m);
            DAQPResult res{};
            res.x = x.data();
            res.lam = lam.data();
            for (int step = 0; step < 30; ++step) {
                Problem q = hybrid_mpc(6, 7.5 * std::sin(0.4 * step), cost_d);
                p.bu = q.bu;
                p.bl = q.bl;
                qp.bupper = p.bu.data();
                qp.blower = p.bl.data();
                daqp_update_ldp(DAQP_UPDATE_d, &work, &qp);
                daqp_solve(&res, &work);
                Solution s = solve(q, policy);
                assert(res.exitflag == s.exitflag);
                assert(close(res.fval, s.fval, 1e-6));
            }
            work.settings = nullptr;
            free_daqp_workspace(&work);
            free_daqp_ldp(&work);
        }
    }
}
