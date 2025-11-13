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
    solver.set_warm_start();

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
    int cold_iter = solver.get_iterations();

    std::cout << "Solution cold start: \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds";
    std::cout << " | Iterations: " << solver.get_iterations() << std::endl;
    if (!solver.get_primal().isApprox((Eigen::VectorXd(3) << 1, 0.5, -1).finished(),precision)) {
        return 1;
    }
    solver.set_warm_start();
    solver.solve(A, bu, bl, break_points);
    std::cout << "Solution warm start: \n";
    std::cout << solver.get_primal().transpose() << std::endl;
    std::cout << "Solve time: " << solver.get_solve_time() << " seconds";
    std::cout << " | Iterations: " << solver.get_iterations() << std::endl;
    if (!solver.get_primal().isApprox((Eigen::VectorXd(3) << 1, 0.5, -1).finished(),precision)) {
        return 1;
    }

    // Ensure cold start reduced number of iterations
    if(cold_iter <= solver.get_iterations()) return 1;


    return 0;
}
