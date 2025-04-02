#include <Eigen/Dense>
#include "daqp/api.h"

class EigenDAQPResult : public DAQPResult {
  private:
      Eigen::VectorXd x_;
      Eigen::VectorXd lam_;

  public:
    EigenDAQPResult(int n, int m)
      : x_(Eigen::VectorXd(n))
      , lam_(Eigen::VectorXd(m)) {
          x = x_.data();
          lam = lam_.data();
    }

    Eigen::VectorXd get_primal() {
        return x_;
    }

    Eigen::VectorXd get_dual() {
        return lam_;
    }
};


EigenDAQPResult daqp_solve(Eigen::MatrixXd &H,
                           Eigen::VectorXd &f,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd &bu,
                           Eigen::VectorXd &bl,
                           Eigen::VectorXi &sense,
                           Eigen::VectorXi break_points){
    int n  = A.cols();
    int m  = bu.size();
    int ms = m-A.rows();
    int n_tasks = break_points.size();

    double* H_ptr = H.size() == 0 ? NULL : H.data();
    double* f_ptr = f.size() == 0 ? NULL : f.data();
    double* A_ptr = A.size() == 0 ? NULL : A.data();
    int* sense_ptr = A.size() == 0 ? NULL : sense.data();
    int* bp_ptr = break_points.size() == 0 ? NULL : break_points.data();

    assert(bu.size() == bl.size());
    assert(ms <= n);
    assert(n_tasks == 0 || break_points(Eigen::last) == m);

    EigenDAQPResult result(n, m);

    DAQPProblem qp = {n, m, ms, H_ptr, f_ptr, A_ptr, bu.data(), bl.data(), sense_ptr, bp_ptr, n_tasks};
    daqp_quadprog(&result, &qp, NULL);

    return result;
}

// Solve: min_x 0.5 x'*H*x + f'*x
//        s.t  bl<= A*x <= bu
EigenDAQPResult daqp_solve(Eigen::MatrixXd &H,
                           Eigen::VectorXd &f,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd &bu,
                           Eigen::VectorXd &bl) {
    Eigen::VectorXi sense(0);
    Eigen::VectorXi break_points(0);
    return daqp_solve(H,f,A,bu,bl,sense,break_points);
}


// Solve min_x ||x||
//        s.t  bl <= A*x <= bu
EigenDAQPResult daqp_solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd &bu,
                           Eigen::VectorXd &bl,
                           Eigen::VectorXi &break_points) {
    Eigen::MatrixXd H(0,0);
    Eigen::VectorXd f(0);
    Eigen::VectorXi sense(0);
    return daqp_solve(H,f,A,bu,bl,sense,break_points);
}
