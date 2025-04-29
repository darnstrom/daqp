#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>

int main() {
    // Task 1: -1 <= x <= 1
    Eigen::MatrixXd A0  = Eigen::MatrixXd::Identity(3, 3);
    Eigen::VectorXd bu0 = Eigen::VectorXd::Ones(3);
    Eigen::VectorXd bl0 = -Eigen::VectorXd::Ones(3);
    
    // Task 2: x1+x2+x3 <= 1
    Eigen::MatrixXd A1  = (Eigen::MatrixXd(1, 3) << 1, 1, 1).finished();
    Eigen::VectorXd bu1 = Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl1 = Eigen::VectorXd::Constant(1, -DAQP_INF);
    
    // Task 3: x1 - x2 == 0.5
    Eigen::MatrixXd A2  = (Eigen::MatrixXd(1, 3) << 1, -1, 0).finished();
    Eigen::VectorXd bu2 = 0.5 * Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl2 = 0.5 * Eigen::VectorXd::Ones(1);
    
    // Task 4: 10 <= 3*x1+x2-x3 <= 20
    Eigen::MatrixXd A3  = (Eigen::MatrixXd(1, 3) << 3, 1, -1).finished();
    Eigen::VectorXd bu3 = 20 * Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl3 = 10 * Eigen::VectorXd::Ones(1);
    
    // Stack the tasks
    Eigen::MatrixXd A  = (Eigen::MatrixXd(6, 3) << A0, A1, A2, A3).finished();
    Eigen::VectorXd bu = (Eigen::VectorXd(6) << bu0, bu1, bu2, bu3).finished();
    Eigen::VectorXd bl = (Eigen::VectorXd(6) << bl0, bl1, bl2, bl3).finished();
    
    Eigen::VectorXi break_points = (Eigen::VectorXi(4) << 3, 4, 5, 6).finished();
    
    // Solve
    DAQP solver(3, 50, 5);
    solver.solve(A, bu, bl, break_points);
    
    std::cout << "Solution: \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds" << std::endl;

    return solver.get_primal().isApprox((Eigen::VectorXd(3) << 1, 0.5, -1).finished()) ? 0 : 1;
}
