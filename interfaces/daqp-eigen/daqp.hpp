#ifndef DAQP_HPP
#define DAQP_HPP

#include <Eigen/Dense>
#include "api.h"


class EigenDAQPResult : public DAQPResult {
  private:
    Eigen::VectorXd x_;
    Eigen::VectorXd lam_;

  public:
    EigenDAQPResult();
    EigenDAQPResult(int n, int m);
    void resize_primal(int n);
    void resize_dual(int m);
    Eigen::VectorXd get_primal();
    Eigen::VectorXd get_dual();
};


EigenDAQPResult daqp_solve(Eigen::MatrixXd& H,
                           Eigen::VectorXd& f,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd& bu,
                           Eigen::VectorXd& bl,
                           Eigen::VectorXi& sense,
                           Eigen::VectorXi& break_points);

// Solve: min_x 0.5 x'*H*x + f'*x
//         s.t  bl<= A*x <= bu
EigenDAQPResult daqp_solve(Eigen::MatrixXd& H,
                           Eigen::VectorXd& f,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd& bu,
                           Eigen::VectorXd& bl);

// Solve min_x ||x||
//        s.t  bl <= A*x <= bu
EigenDAQPResult daqp_solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd& bu,
                           Eigen::VectorXd& bl,
                           Eigen::VectorXi& break_points);


class DAQP {
  private:
    int max_variables_;
    int max_constraints_;
    int max_constraints_in_level_;
    DAQPWorkspace work_;
    DAQPSettings settings_;
    EigenDAQPResult result_;
    int resize_result(int n, int m, Eigen::VectorXi& break_points);

  public:
    DAQP(int max_variables, int max_constraints, int max_constraints_in_level);
    int update(Eigen::MatrixXd& H,
               Eigen::VectorXd& f,
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
               Eigen::VectorXd& bu,
               Eigen::VectorXd& bl,
               Eigen::VectorXi& sense,
               Eigen::VectorXi& break_points,
               int update_mask = -1);
    const EigenDAQPResult& solve();
    const EigenDAQPResult& solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                                 Eigen::VectorXd& bu,
                                 Eigen::VectorXd& bl,
                                 Eigen::VectorXi& break_points);
};

#endif // DAQP_HPP
