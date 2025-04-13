#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>

int main() {
    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(2, 2);
    Eigen::VectorXd f = Eigen::VectorXd::Ones(2);

    Eigen::MatrixXd A  = (Eigen::MatrixXd(2, 2) << 1, 2, 1, -1).finished();
    Eigen::VectorXd bu = (Eigen::VectorXd(4) << 1, 2, 3, 4).finished();
    Eigen::VectorXd bl = (Eigen::VectorXd(4) << -1, -2, -3, -4).finished();

    Eigen::VectorXi sense        = Eigen::VectorXi::Zero(4);
    Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

    // Solve
    EigenDAQPResult result1 = daqp_solve(H, f, A, bu, bl, sense, break_points);
    std::cout << "QP Solution 1: \n";
    std::cout << result1.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << result1.solve_time << " seconds" << std::endl;

    EigenDAQPResult result2 = daqp_solve(H, f, A, bu, bl);
    std::cout << "QP Solution 2: \n";
    std::cout << result2.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << result2.solve_time << " seconds" << std::endl;

    return result1.get_primal().isApprox((Eigen::VectorXd(2) << -1, -1).finished()) &&
               result1.get_primal().isApprox(result2.get_primal())
           ? 0
           : 1;
}
