#include "api.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

// Binary constraints can be relaxed (DAQP_BINARY cleared) and restored by an
// update of the senses (DAQP_UPDATE_sense): branch and bound then branches over
// the binary constraints of the updated senses. An updated workspace must give
// the same exit flag, number of nodes and solution as a workspace that is set
// up for the updated senses, with and without equality reduction. This includes
// making constraints binary in a workspace that has room for fewer binary
// constraints, or none, and an update of the Hessian (which forms the reduction
// anew) while binary constraints are relaxed.

namespace {

constexpr int N = 12; // Variables
constexpr int NB = 6; // The first NB variables can be binary
constexpr int M = N + 3; // Simple bounds and three general constraints

// H = I + c*ones (c = 0.3 by default), and the unconstrained minimizer is 0.5 for the first NB
// variables and 0.2 for the others. The general constraints are
// sum(x[0:NB]) <= 2.6 (active at the root, so that branching is needed) and the
// equality constraints x[6] - x[7] = 0 and x[8] + x[9] = 0.1.
struct Problem {
    c_float H[N*N], f[N], A[3*N], bu[M], bl[M];
    int sense[M];
    DAQPProblem qp{};
    Problem(c_float c = 0.3) {
        for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++) H[i*N+j] = (i == j) + c;
        for(int i = 0; i < N; i++){
            f[i] = 0;
            for(int j = 0; j < N; j++) f[i] -= H[i*N+j]*(j < NB ? 0.5 : 0.2);
        }
        for(int k = 0; k < 3*N; k++) A[k] = 0;
        for(int j = 0; j < NB; j++) A[j] = 1;
        A[N+6] = 1; A[N+7] = -1;
        A[2*N+8] = 1; A[2*N+9] = 1;
        for(int i = 0; i < N; i++){
            bu[i] = i < NB ? 1 : 2;
            bl[i] = i < NB ? 0 : -2;
            sense[i] = 0;
        }
        bu[N] = 2.6; bl[N] = -DAQP_INF; sense[N] = 0;
        bu[N+1] = bl[N+1] = 0; sense[N+1] = DAQP_ACTIVE | DAQP_IMMUTABLE;
        bu[N+2] = bl[N+2] = 0.1; sense[N+2] = DAQP_ACTIVE | DAQP_IMMUTABLE;
        qp = {N, M, N, H, f, A, bu, bl, sense, nullptr, 1, 0};
    }
    // Binary constraints on the variables in [first, last)
    void set_binary(int first, int last) {
        for(int i = 0; i < NB; i++)
            sense[i] = (i >= first && i < last) ? DAQP_BINARY : 0;
    }
};

struct Solution {
    int exitflag, nodes;
    c_float x[N];
};

Solution solve(DAQPWorkspace* work){
    Solution s;
    DAQPResult res{};
    res.x = s.x;
    res.lam = nullptr;
    daqp_solve(&res, work);
    s.exitflag = res.exitflag;
    s.nodes = res.nodes;
    return s;
}

void free_work(DAQPWorkspace* work){
    work->settings = nullptr;
    free_daqp_workspace(work);
    free_daqp_ldp(work);
}

// The solution of a workspace that is set up for p
Solution reference(Problem& p, DAQPSettings* settings){
    DAQPWorkspace ref{};
    ref.settings = settings;
    setup_daqp(&p.qp, &ref, nullptr);
    Solution s = solve(&ref);
    free_work(&ref);
    return s;
}

c_float dist_to_integer(c_float x){ return std::min(std::fabs(x), std::fabs(x - 1)); }

// Compare the solution of the updated workspace with that of a workspace set up
// for the same problem, and check which of the first NB variables are integer
bool check(const char* name, DAQPWorkspace* work, Problem& p, DAQPSettings* settings,
        int first, int last, bool reduced){
    Solution s = solve(work), r = reference(p, settings);
    bool pass = s.exitflag == DAQP_EXIT_OPTIMAL && s.exitflag == r.exitflag && s.nodes == r.nodes;
    pass = pass && (DAQP_IS_REDUCED(work) != 0) == reduced;
    for(int i = 0; i < N; i++) pass = pass && std::fabs(s.x[i] - r.x[i]) < 1e-8;
    bool fractional = false;
    for(int i = 0; i < NB; i++){
        if(i >= first && i < last) pass = pass && dist_to_integer(s.x[i]) < 1e-8;
        else fractional = fractional || dist_to_integer(s.x[i]) > 1e-3;
    }
    // A relaxed binary constraint is not integer in the solutions below
    if(last-first < NB) pass = pass && fractional;
    std::printf("%-62s nodes %3d (set up: %3d) %s\n", name, s.nodes, r.nodes, pass ? "PASS" : "FAIL");
    if(!pass){
        std::printf("  exitflag %d (set up: %d), x =", s.exitflag, r.exitflag);
        for(int i = 0; i < NB; i++) std::printf(" %.4f", s.x[i]);
        std::printf("\n");
    }
    return pass;
}

bool run(int eq_reduction){
    const bool reduced = eq_reduction == DAQP_EQ_REDUCTION_ON;
    bool all_pass = true;
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.eq_reduction = eq_reduction;
    std::printf("eq_reduction = %d\n", eq_reduction);

    // Relax and restore binary constraints
    {
        Problem p;
        p.set_binary(0, NB);
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&p.qp, &work, nullptr);
        all_pass &= check("All binary", &work, p, &settings, 0, NB, reduced);
        const int nodes_all = solve(&work).nodes;
        p.set_binary(2, NB);
        daqp_update_ldp(DAQP_UPDATE_sense, &work, &p.qp);
        all_pass &= check("Two binary constraints relaxed", &work, p, &settings, 2, NB, reduced);
        const int nodes_relaxed = solve(&work).nodes;
        p.set_binary(0, 0);
        daqp_update_ldp(DAQP_UPDATE_sense, &work, &p.qp);
        Solution s = solve(&work);
        bool pass = s.exitflag == DAQP_EXIT_OPTIMAL && s.nodes == 1 && nodes_relaxed < nodes_all;
        std::printf("%-62s nodes %3d %s\n", "All binary constraints relaxed (root only)", s.nodes,
                pass ? "PASS" : "FAIL");
        all_pass &= pass;
        p.set_binary(0, NB);
        daqp_update_ldp(DAQP_UPDATE_sense, &work, &p.qp);
        all_pass &= check("Binary constraints restored", &work, p, &settings, 0, NB, reduced);
        free_work(&work);
    }

    // More binary constraints than the tree has room for
    {
        Problem p;
        p.set_binary(0, 1);
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&p.qp, &work, nullptr);
        p.set_binary(0, NB);
        daqp_update_ldp(DAQP_UPDATE_sense, &work, &p.qp);
        all_pass &= check("Set up with one binary constraint, updated to all", &work, p, &settings, 0, NB, reduced);
        free_work(&work);
    }

    // No binary constraint at setup (no branch-and-bound workspace)
    {
        Problem p;
        p.set_binary(0, 0);
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&p.qp, &work, nullptr);
        bool pass = work.bnb == nullptr;
        p.set_binary(2, NB);
        daqp_update_ldp(DAQP_UPDATE_sense, &work, &p.qp);
        pass = pass && work.bnb != nullptr;
        std::printf("%-62s %s\n", "Branch-and-bound workspace formed by the update", pass ? "PASS" : "FAIL");
        all_pass &= pass;
        all_pass &= check("Set up without binary constraints, updated", &work, p, &settings, 2, NB, reduced);
        free_work(&work);
    }

    // A new Hessian while binary constraints are relaxed (which forms the
    // reduction anew), then restored
    {
        Problem p;
        p.set_binary(0, NB);
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&p.qp, &work, nullptr);
        p.set_binary(4, NB);
        daqp_update_ldp(DAQP_UPDATE_sense, &work, &p.qp);
        Problem q(0.5);
        for(int i = 0; i < M; i++) q.sense[i] = p.sense[i];
        daqp_update_ldp(DAQP_UPDATE_Rinv | DAQP_UPDATE_v, &work, &q.qp);
        all_pass &= check("New Hessian while relaxed", &work, q, &settings, 4, NB, reduced);
        q.set_binary(0, NB);
        daqp_update_ldp(DAQP_UPDATE_sense, &work, &q.qp);
        all_pass &= check("Binary constraints restored after new Hessian", &work, q, &settings, 0, NB, reduced);
        free_work(&work);
    }
    return all_pass;
}

} // namespace

int main(){
    bool all_pass = run(DAQP_EQ_REDUCTION_OFF);
    all_pass &= run(DAQP_EQ_REDUCTION_ON);
    return all_pass ? 0 : 1;
}
