#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>

// Test QP-informed iterative refinement (daqp_refine_active with Rinv != I).
// Setting pivot_tol = DAQP_INF forces refinement to run on every solve so we
// can verify that the new QP-informed path (r_d + eps correction) still
// produces a correct solution for a QP with a non-trivial Hessian.

int main() {
    double precision = 1e-5;

    // min_x  0.5*x'*H*x + f'*x
    //  s.t.  bl <= [I; A]*x <= bu
    // with a non-trivial Hessian H so that Rinv != I.
    Eigen::MatrixXd H = (Eigen::MatrixXd(3, 3) <<
        4, 1, 0,
        1, 3, 0,
        0, 0, 2).finished();
    Eigen::VectorXd f = (Eigen::VectorXd(3) << 1, 2, -1).finished();

    // Simple bounds: -2 <= x <= 2
    // General constraints: A*x in [-3, 3]
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A =
        (Eigen::MatrixXd(2, 3) <<
            1,  1,  1,
            1, -1,  0).finished();
    Eigen::VectorXd bu = (Eigen::VectorXd(5) <<  2,  2,  2,  3,  3).finished();
    Eigen::VectorXd bl = (Eigen::VectorXd(5) << -2, -2, -2, -3, -3).finished();

    Eigen::VectorXi sense        = Eigen::VectorXi::Zero(5);
    Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

    // Reference: solve with default pivot_tol (refinement triggered only for
    // near-singular factorisations).
    EigenDAQPResult ref = daqp_solve(H, f, A, bu, bl, sense, break_points);
    std::cout << "Reference solution: " << ref.get_primal().transpose() << "\n";
    std::cout << "Reference exitflag: " << ref.exitflag << "\n";

    // Solve with pivot_tol = DAQP_INF so that daqp_refine_active is always
    // called, exercising the QP-informed refinement path.
    DAQP solver(3, 10, 10);
    // Setting pivot_tol = DAQP_INF forces daqp_refine_active to run on every
    // solve (normally it only runs when the LDL D diagonal drops below pivot_tol).
    solver.set_pivot_tol(DAQP_INF);
    int status = solver.update(H, f, A, bu, bl, sense, break_points);
    if (status < 0) {
        std::cerr << "update() failed with status " << status << "\n";
        return 1;
    }
    solver.solve();
    std::cout << "Refined solution: " << solver.get_primal().transpose() << "\n";
    std::cout << "Refined exitflag: " << solver.get_status() << "\n";

    if (!solver.get_primal().isApprox(ref.get_primal(), precision)) {
        std::cerr << "FAIL: refined solution differs from reference\n";
        return 1;
    }
    std::cout << "PASS\n";
    return 0;
}
