#include "api.h"
#include "utils.h"

#include <cmath>
#include <cstdio>

// Updating a workspace must give the same solution as setting up the updated
// problem from scratch, also when a previous update exited early because the
// unconstrained optimum was feasible (which skips forming parts of the LDP),
// or when a previous update failed.

namespace {

// Two variables with a non-diagonal Hessian, two simple bounds, and one
// general constraint x0 + x1 >= -0.8
struct Problem {
    c_float H[4] = {4, 1, 1, 2};
    c_float f[2] = {0.5, 0.5}; // x_unc = [-1/14, -3/14] is feasible
    c_float A[2] = {1, 1};
    c_float bu[3] = {5, 5, 2};
    c_float bl[3] = {-5, -5, -0.8};
    int sense[3] = {0, 0, 0};
    DAQPProblem qp{};
    Problem(int problem_type = 0) {
        qp = {2, 3, 2, H, f, A, bu, bl, sense, nullptr, 1, problem_type};
    }
    // x_unc = [-2/7, -6/7] violates the general constraint
    void make_constrained() { f[0] = 2; f[1] = 2; }
    void make_unconstrained() { f[0] = 0.5; f[1] = 0.5; }
};

constexpr int vec_mask = DAQP_UPDATE_v | DAQP_UPDATE_d;
constexpr int unc_mask = vec_mask | DAQP_UPDATE_unconstrained;

bool check(const char* name, DAQPWorkspace* work, DAQPProblem* qp){
    c_float x[2], lam[3], x_ref[2], lam_ref[3];
    DAQPResult res{}, res_ref{};
    res.x = x; res.lam = lam;
    res_ref.x = x_ref; res_ref.lam = lam_ref;
    daqp_solve(&res, work);

    DAQPSettings settings;
    daqp_default_settings(&settings);
    DAQPWorkspace ref{};
    ref.settings = &settings;
    setup_daqp(qp, &ref, nullptr);
    daqp_solve(&res_ref, &ref);
    ref.settings = nullptr;
    free_daqp_workspace(&ref);
    free_daqp_ldp(&ref);

    bool pass = res.exitflag == res_ref.exitflag;
    for(int i = 0; i < 2; i++) pass = pass && std::fabs(x[i] - x_ref[i]) < 1e-8;
    for(int i = 0; i < qp->m; i++) pass = pass && std::fabs(lam[i] - lam_ref[i]) < 1e-8;
    std::printf("%-58s %s\n", name, pass ? "PASS" : "FAIL");
    if(!pass)
        std::printf("  got x = [%g, %g] (exitflag %d), expected x = [%g, %g] (exitflag %d)\n",
                    x[0], x[1], res.exitflag, x_ref[0], x_ref[1], res_ref.exitflag);
    return pass;
}

void free_work(DAQPWorkspace* work){
    work->settings = nullptr;
    free_daqp_workspace(work);
    free_daqp_ldp(work);
}

} // namespace

int main(){
    bool all_pass = true;
    DAQPSettings settings;
    daqp_default_settings(&settings);

    // The unconstrained check reuses a normalized Rinv
    {
        Problem p;
        p.make_constrained();
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&p.qp, &work, nullptr);
        daqp_update_ldp(unc_mask, &work, &p.qp);
        all_pass &= check("Unconstrained check after normalization", &work, &p.qp);
        free_work(&work);
    }

    // An unconstrained exit at setup skips forming M
    {
        Problem p;
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp_main(&p.qp, &work, nullptr, DAQP_UPDATE_unconstrained);
        all_pass &= check("Unconstrained exit at setup", &work, &p.qp);
        p.make_constrained();
        daqp_update_ldp(vec_mask, &work, &p.qp);
        all_pass &= check("Vector update after unconstrained setup", &work, &p.qp);
        p.make_unconstrained();
        daqp_update_ldp(unc_mask, &work, &p.qp);
        all_pass &= check("Unconstrained again", &work, &p.qp);
        p.make_constrained();
        daqp_update_ldp(unc_mask, &work, &p.qp);
        all_pass &= check("Constrained again", &work, &p.qp);
        free_work(&work);
    }

    // Consecutive unconstrained exits keep M pending
    {
        Problem p;
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp_main(&p.qp, &work, nullptr, DAQP_UPDATE_unconstrained);
        daqp_update_ldp(unc_mask, &work, &p.qp);
        all_pass &= check("Two unconstrained exits", &work, &p.qp);
        p.make_constrained();
        daqp_update_ldp(unc_mask, &work, &p.qp);
        all_pass &= check("Constrained after two unconstrained exits", &work, &p.qp);
        free_work(&work);
    }

    // An unconstrained exit in a vector update does not make M pending
    {
        Problem p;
        p.make_constrained();
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&p.qp, &work, nullptr);
        p.make_unconstrained();
        daqp_update_ldp(unc_mask, &work, &p.qp);
        bool pass = (work.state & DAQP_STATE_UNCONSTRAINED) && !(work.state & DAQP_UPDATE_M);
        std::printf("%-58s %s\n", "M is kept after unconstrained vector update", pass ? "PASS" : "FAIL");
        all_pass &= pass;
        all_pass &= check("Unconstrained vector update", &work, &p.qp);
        p.make_constrained();
        daqp_update_ldp(unc_mask, &work, &p.qp);
        all_pass &= check("Constrained vector update", &work, &p.qp);
        free_work(&work);
    }

    // A failed update leaves the skipped parts pending
    {
        Problem p;
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp_main(&p.qp, &work, nullptr, DAQP_UPDATE_unconstrained);
        p.bl[0] = 6; // Infeasible bounds
        int flag = daqp_update_ldp(vec_mask, &work, &p.qp);
        bool pass = flag < 0;
        std::printf("%-58s %s\n", "Infeasible bounds are rejected", pass ? "PASS" : "FAIL");
        all_pass &= pass;
        p.bl[0] = -5;
        p.make_constrained();
        daqp_update_ldp(vec_mask, &work, &p.qp);
        all_pass &= check("Vector update after failed update", &work, &p.qp);
        free_work(&work);
    }

    // A new Hessian in a workspace that has already been normalized
    {
        Problem p;
        p.make_constrained();
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&p.qp, &work, nullptr);
        p.make_unconstrained();
        p.H[0] = 3; p.H[1] = 0.5; p.H[2] = 0.5; p.H[3] = 3;
        daqp_update_ldp(DAQP_UPDATE_Rinv | unc_mask, &work, &p.qp);
        all_pass &= check("Unconstrained exit after new Hessian", &work, &p.qp);
        p.make_constrained();
        daqp_update_ldp(vec_mask, &work, &p.qp);
        all_pass &= check("Vector update after new Hessian", &work, &p.qp);
        free_work(&work);
    }

    // Non-symmetric AVI: the unconstrained exit skips the factorization
    {
        Problem p(1);
        p.H[1] = 1.5; p.H[2] = 0.5;
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp_main(&p.qp, &work, nullptr, DAQP_UPDATE_unconstrained);
        all_pass &= check("AVI unconstrained exit at setup", &work, &p.qp);
        p.make_constrained();
        daqp_update_ldp(vec_mask, &work, &p.qp);
        all_pass &= check("AVI vector update after unconstrained setup", &work, &p.qp);
        free_work(&work);
    }

    // Non-symmetric AVI: a vector update has to form d
    {
        Problem p(1);
        p.H[1] = 1.5; p.H[2] = 0.5;
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&p.qp, &work, nullptr);
        p.make_constrained();
        daqp_update_ldp(vec_mask, &work, &p.qp);
        all_pass &= check("AVI vector update", &work, &p.qp);
        free_work(&work);
    }

    // LP: updating Rinv without a Hessian
    {
        c_float f[2] = {1, -1}, bu[2] = {1, 1}, bl[2] = {-1, -1};
        DAQPProblem qp = {2, 2, 2, nullptr, f, nullptr, bu, bl, nullptr, nullptr, 1, 0};
        DAQPWorkspace work{};
        work.settings = &settings;
        setup_daqp(&qp, &work, nullptr);
        daqp_update_ldp(DAQP_UPDATE_Rinv | vec_mask, &work, &qp);
        all_pass &= check("LP update with DAQP_UPDATE_Rinv", &work, &qp);
        free_work(&work);
    }

    return all_pass ? 0 : 1;
}
