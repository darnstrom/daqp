#include "daqp.hpp"
#include "utils.h"


EigenDAQPResult::EigenDAQPResult()
  : DAQPResult() {
    active_set = nullptr;
}

EigenDAQPResult::EigenDAQPResult(int n, int m)
  : x_{Eigen::VectorXd(n)}
  , lam_{Eigen::VectorXd(m)}
  , active_set_{Eigen::VectorXi(m)} {
    x = x_.data();
    lam = lam_.data();
    active_set = active_set_.data();
}

void EigenDAQPResult::resize_primal(int n) {
    x_.resize(n);
    x = x_.data();
}

void EigenDAQPResult::resize_dual(int m) {
    lam_.resize(m);
    lam = lam_.data();
}

void EigenDAQPResult::resize_active_set(int m) {
    active_set_.resize(m);
    active_set = active_set_.data();
}

Eigen::VectorXd EigenDAQPResult::get_primal() const {
    return x_;
}

Eigen::VectorXd EigenDAQPResult::get_dual() const {
    return lam_;
}

Eigen::VectorXi EigenDAQPResult::get_active_set() const {
    if (n_active > 0 && active_set != nullptr) {
        active_set_.conservativeResize(n_active);
    } else {
        active_set_.resize(0);
    }
    return active_set_;
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
    int *sense_ptr = sense.size() == 0 ? nullptr : sense.data();
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

DAQP::~DAQP(){
    work_.settings = nullptr;
    free_daqp_workspace(&work_);
    free_daqp_ldp(&work_);
}

int DAQP::resize_result(const int n, const int m, Eigen::VectorXi& break_points){
    // Ensure that the new dimensions is not too large
    if (n > max_variables_) return 1;
    if (m > max_constraints_) return 2;

    // Check if maximum soft constraints exceeds limit
    int ns = 0;
    for(int i = 0, prev =0; i < break_points.size(); i++){
        ns = (ns  > break_points(i)- prev) ? ns : break_points(i)-prev;
        prev = break_points(i);
    }
    if (ns > max_constraints_in_level_) return 3;

    // Resize primal and dual variables
    if (result_.x == nullptr) {
        result_.resize_primal(n);
        result_.resize_dual(m);
        result_.resize_active_set(m);
    }
    if (n != work_.n) {
        result_.resize_primal(n);
    }
    if (m != work_.m) {
        result_.resize_dual(m);
        result_.resize_active_set(m);
    }
    return 0;
}

// Update the workspace with a new problem
int DAQP::update(Eigen::MatrixXd& H,
                 Eigen::VectorXd& f,
                 Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                 Eigen::VectorXd& bu,
                 Eigen::VectorXd& bl,
                 Eigen::VectorXi& sense,
                 Eigen::VectorXi& break_points,
                 int update_mask) {

    int n       = A.cols();
    int m       = bu.size();
    int ms      = m - A.rows();
    int n_tasks = break_points.size();

    assert(bu.size() == bl.size());
    assert(ms <= n);
    assert(n_tasks == 0 || break_points(Eigen::last) == m);

    // Get correct pointers
    double *H_ptr  = H.size() == 0 ? nullptr : H.data();
    double *f_ptr  = f.size() == 0 ? nullptr : f.data();
    double *A_ptr  = A.size() == 0 ? nullptr : A.data();
    int *sense_ptr = sense.size() == 0 ? nullptr : sense.data();
    int *bp_ptr    = break_points.size() == 0 ? nullptr : break_points.data();

    //
    int resize_status = resize_result(n,m,break_points);
    assert(resize_status == 0);

    DAQPProblem qp = {n, m, ms, H_ptr, f_ptr, A_ptr, bu.data(), bl.data(), sense_ptr, bp_ptr, n_tasks};
    qp_ = qp; // To remember the qp

    if(update_mask < 0){
        // Assume that everythig should be updated
        update_mask = UPDATE_Rinv + UPDATE_M + UPDATE_v + UPDATE_d + UPDATE_sense + UPDATE_hierarchy;
    }

    return update_ldp(update_mask, &work_, &qp_);
}

// Solve the LDP that is in the workspace
const EigenDAQPResult& DAQP::solve(){
    daqp_solve(&result_,&work_);
    return result_;
}

// Solve min_x ||x||
//        s.t  bl <= A*x <= bu
const EigenDAQPResult& DAQP::solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A,
                                   Eigen::VectorXd& bu,
                                   Eigen::VectorXd& bl,
                                   Eigen::VectorXi& break_points) {
    Eigen::MatrixXd H(0, 0);
    Eigen::VectorXd f(0);
    Eigen::VectorXi sense(0);
    update(H,f,A,bu,bl,sense,break_points,UPDATE_M + UPDATE_d + UPDATE_sense + UPDATE_hierarchy);
    return solve();
}

void DAQP::set_primal_tol(double val)  { settings_.primal_tol = val;}
void DAQP::set_dual_tol(double val)    { settings_.dual_tol = val;}
void DAQP::set_zero_tol(double val)    { settings_.zero_tol = val;}
void DAQP::set_pivot_tol(double val)   { settings_.pivot_tol = val;}
void DAQP::set_progress_tol(double val){ settings_.progress_tol = val;}
void DAQP::set_cycle_tol(int val)      { settings_.cycle_tol = val;}
void DAQP::set_iter_limit(int val)     { settings_.iter_limit = val;}
void DAQP::set_fval_bound(double val)  { settings_.fval_bound = val;}
void DAQP::set_eps_prox(double val)    { settings_.eps_prox = val;}
void DAQP::set_eta_prox(double val)    { settings_.eta_prox = val;}
void DAQP::set_rho_soft(double val)    { settings_.rho_soft = val;}
void DAQP::set_rel_subopt(double val)  { settings_.rel_subopt = val;}
void DAQP::set_abs_subopt(double val)  { settings_.abs_subopt = val;}
void DAQP::set_sing_tol(double val)  { settings_.sing_tol = val;}
void DAQP::set_refactor_tol(double val)  { settings_.refactor_tol = val;}

Eigen::VectorXd DAQP::get_primal() {return result_.get_primal();}
Eigen::VectorXd DAQP::get_dual()   {return result_.get_dual();}
Eigen::VectorXi DAQP::get_active_set() {return result_.get_active_set();}
int             DAQP::get_status() {return result_.exitflag;}
int             DAQP::get_iterations() {return result_.iter;}
double          DAQP::get_solve_time() {return result_.solve_time;}




