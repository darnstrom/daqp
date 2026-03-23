// Test: badly scaled QPs must not be falsely declared infeasible.
//
// For a poorly-scaled QP the LDP transformation produces v = R^{-T}*f that
// can be very large.  The internal LDP objective ||u||^2 therefore equals
// roughly ||v||^2 even at the optimum, which used to exceed the internal
// fval_bound (= 2*settings->fval_bound) and trigger a spurious
// DAQP_EXIT_INFEASIBLE.
//
// The fix adds ||v||^2 to the internal bound so that the check is equivalent
// to testing whether the QP objective 0.5*(||u||^2 - ||v||^2) exceeds
// settings->fval_bound.

#include <Eigen/Dense>
#include <daqp.hpp>
#include <iostream>

// Helper: solve a 1-D bounded QP  min 0.5*H*x^2 + f*x  s.t.  bl <= x <= bu
// and return the exit flag.  The expected optimum is x_expected.
static int solve_1d(double H_val, double f_val, double bl, double bu,
                    double x_expected, double tol, const char* label) {
    Eigen::MatrixXd H = H_val * Eigen::MatrixXd::Identity(1, 1);
    Eigen::VectorXd f(1);
    f(0) = f_val;

    // A has 0 rows and 1 column -> ms = 1 (pure simple bound on x)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(0, 1);
    Eigen::VectorXd buv(1), blv(1);
    buv(0) = bu;
    blv(0) = bl;

    Eigen::VectorXi sense(1);
    sense(0) = 0;
    Eigen::VectorXi bp(0);

    EigenDAQPResult res = daqp_solve(H, f, A, buv, blv, sense, bp);

    bool ok = (res.exitflag == DAQP_EXIT_OPTIMAL ||
               res.exitflag == DAQP_EXIT_SOFT_OPTIMAL);
    if (!ok)
        std::cerr << label << ": exitflag = " << res.exitflag
                  << " (expected OPTIMAL=1)\n";

    double err = std::abs(res.get_primal()(0) - x_expected);
    if (ok && err > tol) {
        std::cerr << label << ": x = " << res.get_primal()(0)
                  << " (expected " << x_expected << ", tol " << tol << ")\n";
        ok = false;
    }
    return ok ? 0 : 1;
}

int main() {
    int failures = 0;

    // ------------------------------------------------------------------
    // Case 1 (non-prox path): H is just above zero_tol so the solver
    // uses the Cholesky factor directly without proximal regularisation.
    // H = 1e-10, f = 1e13 -> v = f/sqrt(H) = 1e18, ||v||^2 = 1e36.
    // Without the fix ||u||^2 = 1e36 > 2*DAQP_INF = 2e30 and the solver
    // would return INFEASIBLE.  With the fix it returns OPTIMAL with x*=0.
    // ------------------------------------------------------------------
    failures += solve_1d(/*H=*/1e-10, /*f=*/1e13,
                         /*bl=*/0.0,  /*bu=*/1.0,
                         /*x_expected=*/0.0, /*tol=*/1e-4,
                         "case1_no_prox");

    // ------------------------------------------------------------------
    // Case 2 (prox path): H is below zero_tol so proximal regularisation
    // is applied.  The effective Hessian is eps_prox = 1e-6, giving
    // Rinv ~ 1e3 and v = f * Rinv = 1e16, ||v||^2 = 1e32.
    // Without the fix the inner daqp_ldp returns INFEASIBLE.
    // ------------------------------------------------------------------
    failures += solve_1d(/*H=*/1e-20, /*f=*/1e13,
                         /*bl=*/0.0,  /*bu=*/1.0,
                         /*x_expected=*/0.0, /*tol=*/1e-4,
                         "case2_with_prox");

    // ------------------------------------------------------------------
    // Case 3: a 2-D problem that mixes a large-eigenvalue direction (well
    // scaled) with a small-eigenvalue direction (badly scaled).  Both
    // constraints must be satisfied at the solution.
    // min  0.5*(1e8*x1^2 + 1e-8*x2^2) + 0*x
    // s.t. 0 <= x1 <= 1,  0 <= x2 <= 1
    // x_unc = 0 -> x* = (0,0).
    // ------------------------------------------------------------------
    {
        Eigen::MatrixXd H(2, 2);
        H << 1e8, 0, 0, 1e-8;
        Eigen::VectorXd f = Eigen::VectorXd::Zero(2);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(0, 2);
        Eigen::VectorXd bu(2), bl(2);
        bu << 1.0, 1.0;
        bl << 0.0, 0.0;

        Eigen::VectorXi sense = Eigen::VectorXi::Zero(2);
        Eigen::VectorXi bp(0);

        EigenDAQPResult res = daqp_solve(H, f, A, bu, bl, sense, bp);

        bool ok = (res.exitflag == DAQP_EXIT_OPTIMAL ||
                   res.exitflag == DAQP_EXIT_SOFT_OPTIMAL);
        if (!ok) {
            std::cerr << "case3_2d: exitflag = " << res.exitflag
                      << " (expected OPTIMAL=1)\n";
            failures++;
        } else {
            double err = res.get_primal().norm();
            if (err > 1e-4) {
                std::cerr << "case3_2d: x = " << res.get_primal().transpose()
                          << " (expected [0,0])\n";
                failures++;
            }
        }
    }

    if (failures == 0)
        std::cout << "All badly-scaled tests passed.\n";
    return failures;
}
