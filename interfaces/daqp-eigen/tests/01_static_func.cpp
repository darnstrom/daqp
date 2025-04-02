#include <iostream>
#include <Eigen/Dense>
#include "daqp.hpp"

int main()
{
    //  DATA
    Eigen::MatrixXd H,A;
    Eigen::VectorXd f,bu,bl;
    Eigen::VectorXi sense,break_points;

    // === BASIC QP ===
    H = Eigen::MatrixXd::Identity(2, 2);
    f = Eigen::VectorXd::Ones(2);
    A = (Eigen::MatrixXd(2, 2) << 1, 2, 1, -1).finished();
    bu = (Eigen::VectorXd(4) << 1,2,3,4).finished();
    bl = (Eigen::VectorXd(4) << -1,-2,-3,-4).finished();
    sense = Eigen::VectorXi::Zero(4);
    break_points = Eigen::VectorXi::Zero(0);

    // Solve
    EigenDAQPResult result1 = daqp_solve(H, f, A, bu, bl, sense, break_points);
    std::cout << "QP Solution 1: \n";
    std::cout << result1.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << result1.solve_time << " seconds" << std::endl;

    EigenDAQPResult result2 = daqp_solve(H, f, A, bu, bl);
    std::cout << "QP Solution 2: \n";
    std::cout << result2.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << result2.solve_time << " seconds" << std::endl;


    // === BASIC HQP ===
    // Task 1: -1 <= x <= 1
    Eigen::MatrixXd A0 = Eigen::MatrixXd::Identity(3, 3);
    Eigen::VectorXd bu0 = Eigen::VectorXd::Ones(3);
    Eigen::VectorXd bl0 = -Eigen::VectorXd::Ones(3);
    auto m0 = bu0.size();

    // Task 2: x1+x2+x3 <= 1
    Eigen::MatrixXd A1 = (Eigen::MatrixXd(1, 3) << 1, 1, 1).finished();
    Eigen::VectorXd bu1 = Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl1 = Eigen::VectorXd::Constant(1, -DAQP_INF);
    auto m1 = bu1.size();

    // Task 3: x1 - x2 == 0.5
    Eigen::MatrixXd A2 = (Eigen::MatrixXd(1, 3) << 1, -1, 0).finished();
    Eigen::VectorXd bu2 = 0.5 * Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl2 = 0.5 * Eigen::VectorXd::Ones(1);
    auto m2 = bu2.size();

    // Task 4: 10 <= 3*x1+x2-x3 <= 20
    Eigen::MatrixXd A3 = (Eigen::MatrixXd(1, 3) << 3, 1, -1).finished();
    Eigen::VectorXd bu3 = 20 * Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl3 = 10 * Eigen::VectorXd::Ones(1);
    auto m3 = bu3.size();

    A = (Eigen::MatrixXd(6, 3) << A0, A1, A2, A3).finished();
    bu = (Eigen::VectorXd(6) << bu0, bu1, bu2, bu3).finished();
    bl = (Eigen::VectorXd(6) << bl0, bl1, bl2, bl3).finished();


    break_points = (Eigen::VectorXi(4) << m0, m0 + m1, m0 + m1 + m2, m0 + m1 + m2 + m3).finished();

    // Solve
    EigenDAQPResult result3 = daqp_solve(A, bu, bl, break_points);
    std::cout << "Solution: \n";
    std::cout << result3.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << result3.solve_time << " seconds" << std::endl;

    return 0;
}
