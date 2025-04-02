#include <iostream>
#include <Eigen/Dense>
#include "daqp.hpp"

// Test basic QP
int main()
{
    // Task 1: -1 <= x <= 1
    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(2, 2);
    Eigen::MatrixXd A = (Eigen::MatrixXd(2, 4) << 1, 2, 1, -1).finished();
    Eigen::VectorXd bu = (Eigen::VectorXd(4) << 1,2,3,4).finished();
    Eigen::VectorXd bl = (Eigen::VectorXd(4) << -1,-2,-3,-4).finished();
    Eigen::VectorXd sense = Eigen::VectorXi::Zeros(4);
    Eigen::VectorXd break_points = Eigen::VectorXi::Zeros(0);

    // Solve
    EigenDAQPResult result = daqp_solve(H, f, A, bu, bl, sense, break_points);
    std::cout << "Solution1: \n";
    std::cout << result.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << result.solve_time << " seconds" << std::endl;

    result = daqp_solve(H, f, A, bu, bl);
    std::cout << "Solution2: \n";
    std::cout << result.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << result.solve_time << " seconds" << std::endl;

    return 0;
}
