#include <Eigen/Dense>
#include "daqp/api.h"

class EigenDAQPResult : public DAQPResult {
  private:
      Eigen::VectorXd x_;
      Eigen::VectorXd lam_;

  public:
    EigenDAQPResult(int m, int n)
      : x_(Eigen::VectorXd(n))
      , lam_(Eigen::VectorXd(m)) {
          x = x_.data();
          lam = lam_.data();
    }

    Eigen::VectorXd get_x() {
        return x_;
    }

    Eigen::VectorXd get_lam() {
        return lam_;
    }
};

EigenDAQPResult daqp_solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd &bu,
                           Eigen::VectorXd &bl,
                           Eigen::VectorXi &break_points) {

    int cols    = A.cols();
    int rows    = A.rows();
    int n_tasks = break_points.size();

    assert(bu.size() == rows);
    assert(bl.size() == rows);
    assert(break_points(Eigen::last) == rows);

    EigenDAQPResult result(rows, cols);

    DAQPProblem qp = {cols, rows, 0, NULL, NULL, A.data(), bu.data(), bl.data(), NULL, break_points.data(), n_tasks};
    daqp_quadprog(&result, &qp, NULL);

    return result;
}
