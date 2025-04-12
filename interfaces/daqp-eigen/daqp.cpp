#include "daqp.hpp"
#include "utils.h"


EigenDAQPResult::EigenDAQPResult()
  : DAQPResult() {
}

EigenDAQPResult::EigenDAQPResult(int n, int m)
  : x_{Eigen::VectorXd(n)}
  , lam_{Eigen::VectorXd(m)} {
    x   = x_.data();
    lam = lam_.data();
}

void EigenDAQPResult::resize_primal(int n) {
    x_.resize(n);
    x = x_.data();
}

void EigenDAQPResult::resize_dual(int m) {
    lam_.resize(m);
    lam = lam_.data();
}

Eigen::VectorXd EigenDAQPResult::get_primal() {
    return x_;
}

Eigen::VectorXd EigenDAQPResult::get_dual() {
    return lam_;
}


EigenDAQPResult daqp_solve(Eigen::MatrixXd& H,
                           Eigen::VectorXd& f,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd& bu,
                           Eigen::VectorXd& bl,
                           Eigen::VectorXi& sense,
                           Eigen::VectorXi& break_points) {
    int n       = A.cols();
    int m       = bu.size();
    int ms      = m - A.rows();
    int n_tasks = break_points.size();

    double *H_ptr  = H.size() == 0 ? nullptr : H.data();
    double *f_ptr  = f.size() == 0 ? nullptr : f.data();
    double *A_ptr  = A.size() == 0 ? nullptr : A.data();
    int *sense_ptr = A.size() == 0 ? nullptr : sense.data();
    int *bp_ptr    = break_points.size() == 0 ? nullptr : break_points.data();

    assert(bu.size() == bl.size());
    assert(ms <= n);
    assert(n_tasks == 0 || break_points(Eigen::last) == m);

    EigenDAQPResult result(n, m);

    DAQPProblem qp = {n, m, ms, H_ptr, f_ptr, A_ptr, bu.data(), bl.data(), sense_ptr, bp_ptr, n_tasks};
    daqp_quadprog(&result, &qp, nullptr);

    return result;
}

// Solve: min_x 0.5 x'*H*x + f'*x
//         s.t  bl<= A*x <= bu
EigenDAQPResult daqp_solve(Eigen::MatrixXd& H,
                           Eigen::VectorXd& f,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd& bu,
                           Eigen::VectorXd& bl) {
    Eigen::VectorXi sense(0);
    Eigen::VectorXi break_points(0);
    return daqp_solve(H, f, A, bu, bl, sense, break_points);
}

// Solve min_x ||x||
//        s.t  bl <= A*x <= bu
EigenDAQPResult daqp_solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                           Eigen::VectorXd& bu,
                           Eigen::VectorXd& bl,
                           Eigen::VectorXi& break_points) {
    Eigen::MatrixXd H(0, 0);
    Eigen::VectorXd f(0);
    Eigen::VectorXi sense(0);
    return daqp_solve(H, f, A, bu, bl, sense, break_points);
}


DAQP::DAQP(int max_variables, int max_constraints, int max_constraints_in_level)
  : max_variables_{max_variables}
  , max_constraints_{max_constraints}
  , max_constraints_in_level_{max_constraints_in_level} {
    allocate_daqp_workspace(&work_, max_variables, max_constraints_in_level);
    allocate_daqp_ldp(&work_, max_variables, max_constraints, 0, 0, 0);
    daqp_default_settings(&settings_);
    work_.settings = &settings_;
}

const EigenDAQPResult& DAQP::solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                                   Eigen::VectorXd& bu,
                                   Eigen::VectorXd& bl,
                                   Eigen::VectorXi& break_points) {
    int cols    = A.cols();
    int rows    = A.rows();
    int n_tasks = break_points.size();
    assert(bu.size() == rows);
    assert(bl.size() == rows);
    assert(break_points(Eigen::last) == rows);
    if (result_.x == nullptr) {
        result_.resize_primal(cols);
        result_.resize_dual(rows);
    }
    if (cols != work_.n) {
        result_.resize_primal(cols);
    }
    if (rows != work_.m) {
        result_.resize_dual(rows);
    }
    DAQPProblem qp = {
      cols, rows, 0, nullptr, nullptr, A.data(), bu.data(), bl.data(), nullptr, break_points.data(), n_tasks};
    update_ldp(update_mask_, &work_, &qp);
    int exitflag = daqp_hiqp(&work_);
    daqp_extract_result(&result_, &work_);
    return result_;
}
