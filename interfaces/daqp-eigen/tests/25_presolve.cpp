#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>

// settings.presolve: tighter mixed-integer problems with the same solution

namespace {

struct Problem {
    int n, m, ms;
    std::vector<c_float> H, f, A, bu, bl;
    std::vector<int> sense;
    Problem(int n_, int mg)
        : n(n_), m(n_ + mg), ms(n_), H(n_ * n_, 0.0), f(n_, 0.0), A(mg * n_, 0.0),
          bu(n_ + mg, DAQP_INF), bl(n_ + mg, -DAQP_INF), sense(n_ + mg, 0) {}
    DAQPProblem qp() {
        return {n, m, ms, H.data(), f.data(), A.data(), bu.data(), bl.data(), sense.data(), nullptr, 0, 0};
    }
    c_float& a(int r, int j) { return A[r * n + j]; }
};

struct Solution {
    int exitflag, nodes;
    c_float fval;
};

// By value, since daqp_quadprog may change sense. presolve sets settings.presolve
Solution solve(Problem p, int presolve = 0) {
    std::vector<c_float> x(p.n), lam(p.m);
    DAQPResult res{};
    res.x = x.data();
    res.lam = lam.data();
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.presolve = presolve;
    DAQPProblem qp = p.qp();
    daqp_quadprog(&res, &qp, &settings);
    return {res.exitflag, res.nodes, res.fval};
}

bool same(const Solution& s, const Solution& t) {
    return s.exitflag == t.exitflag && (s.exitflag <= 0 || std::fabs(s.fval - t.fval) <= 1e-6 * (1 + std::fabs(s.fval)));
}

// The same exit flag and objective with and without presolve
Solution check_same(const Problem& p) {
    Solution s = solve(p), t = solve(p, 1);
    // Cycling is numerical and can come and go with the tightened data
    if (s.exitflag != DAQP_EXIT_CYCLE && t.exitflag != DAQP_EXIT_CYCLE) assert(same(s, t));
    return t;
}

uint64_t state = 7;
c_float urand() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (state >> 11) * (1.0 / 9007199254740992.0);
}

} // namespace

int main() {
    // Fixed charge with loose big-M constants: y_i <= M_i b_i, y_i in [0,U], sum y >= D
    {
        const int k = 12;
        const c_float U = 100;
        Problem p(2 * k, k + 1);
        for (int i = 0; i < k; ++i) {
            p.sense[i] = DAQP_BINARY;
            p.bl[i] = 0; p.bu[i] = 1;
            p.bl[k + i] = 0; p.bu[k + i] = U;
            p.f[i] = U * (0.5 + 0.08 * ((7 * i) % k));
            p.f[k + i] = -(0.3 + 0.07 * ((5 * i) % k));
            p.H[i * p.n + i] = p.H[(k + i) * p.n + k + i] = 1e-2;
            p.a(i, k + i) = 1; p.a(i, i) = -U * (1 + 0.5 * (i % 3)); p.bu[p.ms + i] = 0;
            p.a(k, k + i) = -1;
        }
        p.bu[p.ms + k] = -4.5 * U;
        Solution s = solve(p), t = check_same(p);
        assert(s.exitflag == DAQP_EXIT_OPTIMAL && t.nodes < s.nodes);
    }

    // Infeasible: binaries with 0.2 <= b_0 + b_1 + b_2 <= 0.4
    {
        Problem p(3, 1);
        for (int i = 0; i < 3; ++i) {
            p.sense[i] = DAQP_BINARY; p.bl[i] = 0; p.bu[i] = 1;
            p.H[i * 4] = 1; p.a(0, i) = 1;
        }
        p.bl[3] = 0.2; p.bu[3] = 0.4;
        assert(solve(p, 1).exitflag == DAQP_EXIT_INFEASIBLE);
    }

    // Soft constraints are not relied on (the soft x_1 <= -10 would make it infeasible),
    // and immutable bounds/constraints are kept
    {
        Problem p(3, 3);
        for (int i = 0; i < 3; ++i) { p.H[i * 4] = 1; p.f[i] = -1; p.bl[i] = -5; p.bu[i] = 5; }
        p.sense[0] = DAQP_BINARY; p.bl[0] = 0; p.bu[0] = 1;
        p.sense[2] = DAQP_IMMUTABLE; // x_2 <= 5 must not change
        p.a(0, 1) = 1; p.a(0, 2) = 1; p.bu[3] = 1;                       // x_1 + x_2 <= 1
        p.a(1, 1) = 1; p.a(1, 0) = -50; p.bu[4] = 0;                     // x_1 <= 50 b (big-M)
        p.a(2, 1) = 1; p.bu[5] = -10; p.sense[5] = DAQP_SOFT;             // soft x_1 <= -10
        assert(check_same(p).exitflag == DAQP_EXIT_OPTIMAL);
    }

    // The problem data is not changed, an infeasibility is reported by daqp_quadprog,
    // and an equality constraint does not keep the presolve from applying
    {
        const int k = 8;
        const c_float U = 100;
        Problem p(2 * k, k + 1);
        for (int i = 0; i < k; ++i) {
            p.sense[i] = DAQP_BINARY;
            p.bl[i] = 0; p.bu[i] = 1;
            p.bl[k + i] = 0; p.bu[k + i] = U;
            p.f[i] = U * (0.6 + 0.05 * i); p.f[k + i] = -(0.4 + 0.05 * ((3 * i) % k));
            p.H[i * p.n + i] = p.H[(k + i) * p.n + k + i] = 1e-2;
            p.a(i, k + i) = 1; p.a(i, i) = -2 * U; p.bu[p.ms + i] = 0;
            p.a(k, k + i) = -1;
        }
        p.bu[p.ms + k] = -3.5 * U;
        Problem q = p;
        std::vector<c_float> x(q.n), lam(q.m);
        DAQPResult res{};
        res.x = x.data();
        res.lam = lam.data();
        DAQPSettings settings;
        daqp_default_settings(&settings);
        settings.presolve = 1;
        DAQPProblem qp = q.qp();
        daqp_quadprog(&res, &qp, &settings);
        assert(res.exitflag == DAQP_EXIT_OPTIMAL);
        assert(q.A == p.A && q.bu == p.bu && q.bl == p.bl);
        Solution s = solve(p), t = solve(p, 1);
        assert(same(s, t) && t.nodes < s.nodes);

        Problem e = p; // sum y = D as an equality (eliminated without presolve)
        e.bl[e.ms + k] = e.bu[e.ms + k];
        assert(same(solve(e), solve(e, 1)));

        Problem g(4, 2); // Infinite bounds that an equality fixes (w = 3), and fval_bound
        const c_float hd[4] = {1, 1e-4, 1, 1}, fd[4] = {-0.6, 1, 1, 0};
        for (int i = 0; i < 4; ++i) { g.H[i * 5] = hd[i]; g.f[i] = fd[i]; }
        g.sense[0] = DAQP_BINARY; g.bl[0] = 0; g.bu[0] = 1; g.bl[1] = 0;
        g.a(0, 0) = 1; g.a(0, 1) = -1; g.bu[4] = 0.6;
        g.a(1, 3) = 1; g.bu[5] = g.bl[5] = 3; g.sense[5] = DAQP_ACTIVE + DAQP_IMMUTABLE;
        for (c_float fb : {4.01, 3.99}) {
            std::vector<c_float> xg(4), lg(6);
            DAQPResult rg{};
            rg.x = xg.data(); rg.lam = lg.data();
            DAQPSettings sg;
            daqp_default_settings(&sg);
            sg.presolve = 1; sg.fval_bound = fb;
            Problem h = g;
            DAQPProblem hq = h.qp();
            daqp_quadprog(&rg, &hq, &sg);
            assert(rg.exitflag == (fb > 4 ? DAQP_EXIT_OPTIMAL : DAQP_EXIT_INFEASIBLE));
            if (fb > 4) assert(std::fabs(rg.fval - 4) < 1e-6);
        }

        Problem c(3, 2); // Strictly convex only after the equality x_1 = x_0 + b is eliminated
        c.H[0] = 1; c.H[8] = 1e-2;
        c.f[0] = -1; c.f[1] = -1; c.f[2] = 2 - 0.5e-2;
        c.sense[2] = DAQP_BINARY;
        for (int i = 0; i < 3; ++i) { c.bl[i] = 0; c.bu[i] = i < 2 ? 20 : 1; }
        c.a(0, 0) = -1; c.a(0, 1) = 1; c.a(0, 2) = -1; c.bu[3] = c.bl[3] = 0;
        c.a(1, 0) = 1; c.a(1, 2) = -10; c.bu[4] = 0; // x_0 <= 10 b (tightened by the presolve)
        {
            std::vector<c_float> xs(3), ls(5);
            DAQPResult rs{};
            rs.x = xs.data(); rs.lam = ls.data();
            DAQPSettings ss;
            daqp_default_settings(&ss);
            ss.presolve = 1; ss.eq_reduction = DAQP_EQ_REDUCTION_ON;
            DAQPProblem sq = c.qp();
            daqp_quadprog(&rs, &sq, &ss);
            assert(rs.exitflag == DAQP_EXIT_OPTIMAL);
            assert(std::fabs(xs[0] - 2) < 1e-9 && std::fabs(xs[1] - 3) < 1e-9 && xs[2] == 1);
        }

        Problem f = p; // Infeasible: the demand exceeds the capacity
        f.bu[f.ms + k] = -(k + 1) * U;
        assert(solve(f).exitflag == DAQP_EXIT_INFEASIBLE && solve(f, 1).exitflag == DAQP_EXIT_INFEASIBLE);

        // Updates of a workspace: the presolve is redone for the new bounds (also for
        // loosened bounds), and switching it off restores the problem
        DAQPWorkspace work{};
        work.settings = &settings;
        Problem w = p;
        DAQPProblem wq = w.qp();
        assert(setup_daqp(&wq, &work, nullptr) >= 0);
        for (int step = 0; step < 12; ++step) {
            if (step == 8) settings.presolve = 0;
            w.bu[w.ms + k] = -(1.5 + 0.5 * (step % 6)) * U;  // The demand
            for (int i = 0; i < k; ++i) w.bu[k + i] = U * (0.5 + 0.1 * ((step + i) % 6)); // Capacities
            daqp_update_ldp(DAQP_UPDATE_d, &work, &wq);
            daqp_solve(&res, &work);
            Solution c = solve(w);
            assert(res.exitflag == c.exitflag);
            if (c.exitflag > 0) assert(std::fabs(res.fval - c.fval) <= 1e-6 * (1 + std::fabs(c.fval)));
        }
        work.settings = nullptr;
        free_daqp_workspace(&work);
        free_daqp_ldp(&work);
    }

    // Random MIQPs: binaries in {l,u}, binary constraints, big-M, one- and two-sided rows
    for (int t = 0; t < 300; ++t) {
        const int n = 6 + t % 4, nb = 2 + t % 3, mg = 5 + t % 5;
        Problem p(n, mg);
        std::vector<c_float> x0(n), G(n * n);
        for (int i = 0; i < n; ++i) {
            if (i < nb) {
                p.sense[i] = DAQP_BINARY;
                p.bl[i] = urand() < 0.5 ? 0 : -1;
                p.bu[i] = p.bl[i] + (urand() < 0.5 ? 1 : 2);
                x0[i] = urand() < 0.5 ? p.bl[i] : p.bu[i];
            }
            else { p.bl[i] = -5; p.bu[i] = 5; x0[i] = 4 * urand() - 2; }
            p.f[i] = 4 * urand() - 2;
        }
        for (auto& g : G) g = 2 * urand() - 1;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                c_float s = 0.1 * (i == j);
                for (int k = 0; k < n; ++k) s += 0.3 * G[k * n + i] * G[k * n + j];
                p.H[i * n + j] = s;
            }
        for (int r = 0; r < mg; ++r) {
            c_float ax = 0;
            for (int j = 0; j < n; ++j) { p.a(r, j) = urand() < 0.5 ? 4 * urand() - 2 : 0; ax += p.a(r, j) * x0[j]; }
            if (r == 0) { // Big-M
                for (int j = 0; j < n; ++j) p.a(r, j) = 0;
                p.a(r, nb) = 1; p.a(r, 0) = -20;
                ax = x0[nb] - 20 * x0[0];
            }
            const c_float kind = urand();
            const int id = p.ms + r;
            if (kind < 0.15) { p.sense[id] = DAQP_BINARY; p.bl[id] = ax - 1; p.bu[id] = ax; }
            else if (kind < 0.4) { p.bl[id] = ax - urand(); p.bu[id] = ax + urand(); }
            else if (kind < 0.55) p.bl[id] = ax - urand();
            else p.bu[id] = ax + urand();
        }
        check_same(p);
    }
}
