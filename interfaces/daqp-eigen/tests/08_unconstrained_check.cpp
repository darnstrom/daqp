#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>

// Test that the unconstrained-feasibility shortcut is correctly bypassed for
// LP objectives and for problems with equality constraints.

int main() {
    double precision = 1e-5;
    bool all_pass = true;

    // -----------------------------------------------------------------------
    // Test 1: LP (no Hessian) where x_unc = -f satisfies the simple bounds
    // but is NOT the LP optimal.
    //
    // min  x1 + x2        (f = [1, 1], H = empty -> LP)
    // s.t. -10 <= x1 <= 10
    //      -10 <= x2 <= 10
    //
    // Unconstrained "optimum": x_unc = -f = [-1, -1] (would be reported
    // by the buggy shortcut as optimal, yielding fval = -2).
    // True LP optimum:         x* = [-10, -10], fval = -20.
    // -----------------------------------------------------------------------
    {
        int n = 2, m = 2;
        Eigen::MatrixXd H(0, 0);  // Empty Hessian signals an LP to DAQP (H_ptr == nullptr)
        Eigen::VectorXd f  = (Eigen::VectorXd(2) << 1, 1).finished();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(0, n);  // no general constraints
        Eigen::VectorXd bu = (Eigen::VectorXd(2) << 10, 10).finished();
        Eigen::VectorXd bl = (Eigen::VectorXd(2) << -10, -10).finished();
        Eigen::VectorXi sense        = Eigen::VectorXi::Zero(m);
        Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

        EigenDAQPResult result = daqp_solve(H, f, A, bu, bl, sense, break_points);

        Eigen::VectorXd expected = (Eigen::VectorXd(2) << -10, -10).finished();
        bool pass = (result.exitflag == DAQP_EXIT_OPTIMAL) &&
                    result.get_primal().isApprox(expected, precision);
        std::cout << "Test 1 (LP optimal, not unconstrained): "
                  << (pass ? "PASS" : "FAIL") << std::endl;
        if (!pass) {
            std::cout << "  Expected: " << expected.transpose() << std::endl;
            std::cout << "  Got:      " << result.get_primal().transpose() << std::endl;
            std::cout << "  exitflag: " << result.exitflag << std::endl;
        }
        all_pass = all_pass && pass;
    }

    // -----------------------------------------------------------------------
    // Test 2: QP with a general equality constraint (sense = ACTIVE+IMMUTABLE)
    // where the unconstrained optimum x_unc happens to satisfy the equality.
    // The shortcut must not be taken; the solver must run and produce the
    // correct primal AND report the equality as active.
    //
    // min  0.5*(x1^2+x2^2) - x1 + x2   (H=I, f=[-1,1])
    // s.t. x1 - x2 = 2                  (general equality, sense = 5)
    //
    // x_unc = -f = [1, -1]; check: 1 - (-1) = 2 == rhs  => equality satisfied.
    // KKT:  x + f + A'*lam = 0  =>  [x1-1; x2+1] + [1;-1]*lam = 0
    //        x1 = 1-lam, x2 = -1+lam, x1-x2 = 2-2*lam = 2  =>  lam = 0.
    // True optimum: x* = [1, -1], lam = 0.
    // -----------------------------------------------------------------------
    {
        int n = 2, m = 3;  // 2 simple bounds + 1 general equality
        Eigen::MatrixXd H = Eigen::MatrixXd::Identity(2, 2);
        Eigen::VectorXd f = (Eigen::VectorXd(2) << -1, 1).finished();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A =
            (Eigen::MatrixXd(1, 2) << 1, -1).finished();
        // [simple bound x1, simple bound x2, equality x1-x2]
        Eigen::VectorXd bu = (Eigen::VectorXd(3) << 5, 5, 2).finished();
        Eigen::VectorXd bl = (Eigen::VectorXd(3) << -5, -5, 2).finished();
        // Equality constraint has sense = DAQP_ACTIVE (1) + DAQP_IMMUTABLE (4) = 5
        Eigen::VectorXi sense = (Eigen::VectorXi(3) << 0, 0, 5).finished();
        Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

        EigenDAQPResult result = daqp_solve(H, f, A, bu, bl, sense, break_points);

        Eigen::VectorXd expected = (Eigen::VectorXd(2) << 1, -1).finished();
        bool pass = (result.exitflag == DAQP_EXIT_OPTIMAL) &&
                    result.get_primal().isApprox(expected, precision);
        std::cout << "Test 2 (QP with equality, x_unc on constraint): "
                  << (pass ? "PASS" : "FAIL") << std::endl;
        if (!pass) {
            std::cout << "  Expected: " << expected.transpose() << std::endl;
            std::cout << "  Got:      " << result.get_primal().transpose() << std::endl;
            std::cout << "  exitflag: " << result.exitflag << std::endl;
        }
        all_pass = all_pass && pass;
    }

    // -----------------------------------------------------------------------
    // Test 3: QP where blower == bupper (unmarked equality, no explicit sense
    // flag) and the unconstrained optimum satisfies it.  The solver must
    // still run normally and return the correct primal.
    //
    // min  0.5*(x1^2+x2^2) - x1 + x2   (H=I, f=[-1,1])
    // s.t. x1 - x2 = 2                  (general equality via equal bounds)
    //      -5 <= x1 <= 5  (simple bounds)
    //      -5 <= x2 <= 5
    //
    // Same as Test 2 but without the explicit IMMUTABLE sense flag.
    // -----------------------------------------------------------------------
    {
        int n = 2, m = 3;
        Eigen::MatrixXd H = Eigen::MatrixXd::Identity(2, 2);
        Eigen::VectorXd f = (Eigen::VectorXd(2) << -1, 1).finished();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A =
            (Eigen::MatrixXd(1, 2) << 1, -1).finished();
        Eigen::VectorXd bu = (Eigen::VectorXd(3) << 5, 5, 2).finished();
        Eigen::VectorXd bl = (Eigen::VectorXd(3) << -5, -5, 2).finished();  // bl==bu for equality
        Eigen::VectorXi sense        = Eigen::VectorXi::Zero(3);             // no explicit flags
        Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

        EigenDAQPResult result = daqp_solve(H, f, A, bu, bl, sense, break_points);

        Eigen::VectorXd expected = (Eigen::VectorXd(2) << 1, -1).finished();
        bool pass = (result.exitflag == DAQP_EXIT_OPTIMAL) &&
                    result.get_primal().isApprox(expected, precision);
        std::cout << "Test 3 (QP with unmarked equality via equal bounds): "
                  << (pass ? "PASS" : "FAIL") << std::endl;
        if (!pass) {
            std::cout << "  Expected: " << expected.transpose() << std::endl;
            std::cout << "  Got:      " << result.get_primal().transpose() << std::endl;
            std::cout << "  exitflag: " << result.exitflag << std::endl;
        }
        all_pass = all_pass && pass;
    }

    // -----------------------------------------------------------------------
    // Test 4: QP with Hessian but NO linear term (f = empty -> v == NULL).
    //
    // min  0.5*(x1^2 + x2^2)   (H = I, f = empty)
    // s.t. -1 <= x1 <= 1
    //      -1 <= x2 <= 1
    //
    // Unconstrained optimum is x = 0, which is feasible.
    // The shortcut should detect this and return x* = [0, 0] with
    // exitflag == DAQP_EXIT_OPTIMAL without running the iterative solver.
    // Before the fix, work->x pointed to uninitialized memory (the malloc'd
    // xold buffer), so the feasibility check used garbage values and could
    // return incorrect results.
    // -----------------------------------------------------------------------
    {
        int n = 2, m = 2;
        Eigen::MatrixXd H = Eigen::MatrixXd::Identity(2, 2);
        Eigen::VectorXd f(0);  // Empty linear term (f == NULL in C layer)
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(0, n);
        Eigen::VectorXd bu = (Eigen::VectorXd(2) << 1, 1).finished();
        Eigen::VectorXd bl = (Eigen::VectorXd(2) << -1, -1).finished();
        Eigen::VectorXi sense        = Eigen::VectorXi::Zero(m);
        Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

        EigenDAQPResult result = daqp_solve(H, f, A, bu, bl, sense, break_points);

        Eigen::VectorXd expected = Eigen::VectorXd::Zero(2);
        bool pass = (result.exitflag == DAQP_EXIT_OPTIMAL) &&
                    result.get_primal().isApprox(expected, precision);
        std::cout << "Test 4 (QP no linear term, x=0 feasible, v==NULL fix): "
                  << (pass ? "PASS" : "FAIL") << std::endl;
        if (!pass) {
            std::cout << "  Expected: " << expected.transpose() << std::endl;
            std::cout << "  Got:      " << result.get_primal().transpose() << std::endl;
            std::cout << "  exitflag: " << result.exitflag << std::endl;
        }
        all_pass = all_pass && pass;
    }

    // -----------------------------------------------------------------------
    // Test 5: QP with Hessian but NO linear term, where x = 0 violates
    // a constraint. The solver must run normally and find the true optimum.
    //
    // min  0.5*(x1^2 + x2^2)   (H = I, f = empty)
    // s.t. 2 <= x1 <= 5
    //      2 <= x2 <= 5
    //
    // x = 0 violates both lower bounds, so the shortcut must NOT be taken.
    // True optimum: x* = [2, 2].
    // -----------------------------------------------------------------------
    {
        int n = 2, m = 2;
        Eigen::MatrixXd H = Eigen::MatrixXd::Identity(2, 2);
        Eigen::VectorXd f(0);  // Empty linear term
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(0, n);
        Eigen::VectorXd bu = (Eigen::VectorXd(2) << 5, 5).finished();
        Eigen::VectorXd bl = (Eigen::VectorXd(2) << 2, 2).finished();
        Eigen::VectorXi sense        = Eigen::VectorXi::Zero(m);
        Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

        EigenDAQPResult result = daqp_solve(H, f, A, bu, bl, sense, break_points);

        Eigen::VectorXd expected = (Eigen::VectorXd(2) << 2, 2).finished();
        bool pass = (result.exitflag == DAQP_EXIT_OPTIMAL) &&
                    result.get_primal().isApprox(expected, precision);
        std::cout << "Test 5 (QP no linear term, x=0 infeasible, v==NULL fix): "
                  << (pass ? "PASS" : "FAIL") << std::endl;
        if (!pass) {
            std::cout << "  Expected: " << expected.transpose() << std::endl;
            std::cout << "  Got:      " << result.get_primal().transpose() << std::endl;
            std::cout << "  exitflag: " << result.exitflag << std::endl;
        }
        all_pass = all_pass && pass;
    }

    return all_pass ? 0 : 1;
}
