#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>

int main() {
    // 
    double precision = 1e-5;
    Eigen::MatrixXd A;
    Eigen::VectorXd bu, bl;
    Eigen::VectorXi break_points;

    // Setup solver
    DAQP solver(3, 50, 5);

    // Task 1: -1 <= x <= 1
    Eigen::MatrixXd A0 = Eigen::MatrixXd::Identity(3, 3);
    Eigen::VectorXd bu0 = Eigen::VectorXd::Ones(3);
    Eigen::VectorXd bl0 = -Eigen::VectorXd::Ones(3);

    // Task 2: x1+x2+x3 <= 1
    Eigen::MatrixXd A1 = (Eigen::MatrixXd(1, 3) << 1, 1, 1).finished();
    Eigen::VectorXd bu1 = Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl1 = Eigen::VectorXd::Constant(1, -DAQP_INF);

    // Task 3: x1 - x2 == 0.5
    Eigen::MatrixXd A2 = (Eigen::MatrixXd(1, 3) << 1, -1, 0).finished();
    Eigen::VectorXd bu2 = 0.5 * Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl2 = 0.5 * Eigen::VectorXd::Ones(1);

    // Task 4: 10 <= 3*x1+x2-x3 <= 20
    Eigen::MatrixXd A3 = (Eigen::MatrixXd(1, 3) << 3, 1, -1).finished();
    Eigen::VectorXd bu3 = 20 * Eigen::VectorXd::Ones(1);
    Eigen::VectorXd bl3 = 10 * Eigen::VectorXd::Ones(1);

    // ==== First setup (stack all tasks) ====
    A = (Eigen::MatrixXd(6, 3) << A0, A1, A2, A3).finished();
    bu = (Eigen::VectorXd(6) << bu0, bu1, bu2, bu3).finished();
    bl = (Eigen::VectorXd(6) << bl0, bl1, bl2, bl3).finished();
    break_points = (Eigen::VectorXi(4) << 3, 4, 5, 6).finished();
    solver.solve(A, bu, bl, break_points);

    std::cout << "Solution 1: \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds";
    std::cout << " | Iterations: " << solver.get_iterations() << std::endl;
    if (!solver.get_primal().isApprox((Eigen::VectorXd(3) << 1, 0.5, -1).finished(),precision)) {
        return 1;
    }
    solver.set_warm_start();
    solver.solve(A, bu, bl, break_points);
    std::cout << "Solution 1 (warm start): \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds";
    std::cout << " | Iterations: " << solver.get_iterations() << std::endl;
    solver.set_cold_start();

    // ==== Second setup: remove some tasks ====
    A = (Eigen::MatrixXd(4, 3) << A0, A3).finished();
    bu = (Eigen::VectorXd(4) << bu0, bu3).finished();
    bl = (Eigen::VectorXd(4) << bl0, bl3).finished();
    break_points = (Eigen::VectorXi(2) << 3, 4).finished();
    solver.solve(A, bu, bl, break_points);

    std::cout << "Solution 2: \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds" << std::endl;
    if (!solver.get_primal().isApprox((Eigen::VectorXd(3) << 1, 1, -1).finished(),precision)) {
        return 1;
    }

    // ==== Third setup: remove the third dimension ====
    A = (Eigen::MatrixXd(5, 2) << 1,0,0,1, 1,1,1,-1,3,1).finished();
    bu = (Eigen::VectorXd(5) << 1,1, bu1, bu2, bu3).finished();
    bl = (Eigen::VectorXd(5) << -1,-1, bl1, bl2, bl3).finished();
    break_points = (Eigen::VectorXi(4) << 2, 3, 4, 5).finished();
    solver.solve(A, bu, bl, break_points);

    std::cout << "Solution 3: \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds" << std::endl;
    solver.solve(A, bu, bl, break_points);
    if (!solver.get_primal().isApprox((Eigen::VectorXd(2) << 0.75, 0.25).finished(),precision)) {
        return 1;
    }

    // ==== Fourth setup: reorder some tasks ====
    A = (Eigen::MatrixXd(6, 3) << A1, A3, A0, A2).finished();
    bu = (Eigen::VectorXd(6) << bu1, bu3, bu0, bu2).finished();
    bl = (Eigen::VectorXd(6) << bl1, bl3, bl0, bl2).finished();
    break_points = (Eigen::VectorXi(4) << 1, 2, 5, 6).finished();
    solver.solve(A, bu, bl, break_points);

    std::cout << "Solution 4: \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds" << std::endl;
    if (!solver.get_primal().isApprox((Eigen::VectorXd(3) << 2.25, 1, -2.25).finished(), precision)) {
        return 1;
    }

    return 0;
}
