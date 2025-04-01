#include <Eigen/Dense>
#include "api.h"

DAQPResult daqp_solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A, Eigen::VectorXd& bu, Eigen::VectorXd& bl, Eigen::VectorXi& break_points) {
    int cols = A.cols();
    int rows = A.rows();
    int n_tasks = break_points.size();

    assert(bu.size() == rows);
    assert(bl.size() == rows);
    assert(break_points(Eigen::last) == rows);

    DAQPResult result;
    result.x = new double[cols];
    result.lam = NULL;

    DAQPProblem qp = {cols, rows, 0, NULL, NULL, A.data(), bu.data(), bl.data(), NULL, break_points.data(), n_tasks};
    daqp_quadprog(&result, &qp, NULL);

    return result;
}
