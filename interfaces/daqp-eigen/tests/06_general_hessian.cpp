#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>

int main() {
    double precision = 1e-5;

    // Solve: min_x 0.5*x'*H*x + f'*x
    //         s.t  bl <= A*x <= bu
    // with a non-identity Hessian using the workspace class
    Eigen::MatrixXd H = (Eigen::MatrixXd(2, 2) << 2, 0, 0, 2).finished();
    Eigen::VectorXd f = (Eigen::VectorXd(2) << 1, 1).finished();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A =
        (Eigen::MatrixXd(2, 2) << 1, 2, 1, -1).finished();
    Eigen::VectorXd bu = (Eigen::VectorXd(4) << 1, 2, 3, 4).finished();
    Eigen::VectorXd bl = (Eigen::VectorXd(4) << -1, -2, -3, -4).finished();

    Eigen::VectorXi sense        = Eigen::VectorXi::Zero(4);
    Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

    // Use the workspace class (DAQP) with a general Hessian
    DAQP solver(2, 10, 10);
    solver.set_time_limit(10.0);
    int status = solver.update(H, f, A, bu, bl, sense, break_points);
    if (status < 0) {
        std::cerr << "update() failed with status " << status << std::endl;
        return 1;
    }
    solver.solve();

    std::cout << "General-Hessian QP solution: \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Status: " << solver.get_status() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds" << std::endl;

    // Expected solution: same as the static daqp_solve result with this H and f
    EigenDAQPResult ref = daqp_solve(H, f, A, bu, bl, sense, break_points);
    std::cout << "Reference solution: \n";
    std::cout << ref.get_primal().transpose() << std::endl;

    return solver.get_primal().isApprox(ref.get_primal(), precision) ? 0 : 1;
}
