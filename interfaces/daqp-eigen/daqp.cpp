#include "daqp.hpp"
#include "utils.h"


EigenDAQPResult::EigenDAQPResult()
  : DAQPResult() {
}

EigenDAQPResult::EigenDAQPResult(int n, int m)
  : x_{Eigen::VectorXd(n)}
  , lam_{Eigen::VectorXd(m)}
  , slack_{Eigen::VectorXd(m)} {
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
    slack_.resize(m);
}

Eigen::VectorXd EigenDAQPResult::get_primal() const {
    return x_;
}

Eigen::VectorXd EigenDAQPResult::get_dual() const {
    return lam_;
}

Eigen::VectorXd EigenDAQPResult::get_slack() const {
    return slack_;
}

Eigen::VectorXi EigenDAQPResult::get_active_set() const {
    Eigen::Array<bool, Eigen::Dynamic, 1> mask = (lam_.array().abs() > 0);

    int n_active = mask.cast<int>().sum();
    Eigen::VectorXi active_set(n_active);
    for (int idx = 0, k = 0; k < lam_.size(); ++k) {
        if (mask(k)) {
            active_set(idx++) = k;
        }
    }
    return active_set;
}

void EigenDAQPResult::set_slack(Eigen::VectorXd const& slack) {
    slack_ = slack;
}


Eigen::VectorXd compute_slacks(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& matrix,
                               Eigen::VectorXd const& upper,
                               Eigen::VectorXd const& lower,
                               Eigen::VectorXd const& primal) {
    // For DAQP first ms entries are simple bounds (primal values), rest are A*primal
    int n  = primal.size();
    int m  = upper.size();
    int ms = m - matrix.rows();

    Eigen::VectorXd values(m);
    values.head(ms) = primal.head(ms);
    if (matrix.rows() > 0) {
        values.tail(matrix.rows()) = matrix * primal;
    }

    auto slack_lw         = (values - lower).array();
    auto slack_up         = (values - upper).array();
    Eigen::VectorXd slack = (slack_lw < 0).select(slack_lw, (slack_up > 0).select(slack_up, 0));
    return slack;
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

    double* H_ptr  = H.size() == 0 ? nullptr : H.data();
    double* f_ptr  = f.size() == 0 ? nullptr : f.data();
    double* A_ptr  = A.size() == 0 ? nullptr : A.data();
    double* bu_ptr = bu.size() == 0 ? nullptr : bu.data();
    double* bl_ptr = bl.size() == 0 ? nullptr : bl.data();
    int* sense_ptr = sense.size() == 0 ? nullptr : sense.data();
    int* bp_ptr    = break_points.size() == 0 ? nullptr : break_points.data();

    assert(bu.size() == bl.size());
    assert(ms <= n);
    assert(n_tasks == 0 || break_points(Eigen::last) == m);

    EigenDAQPResult result(n, m);

    DAQPProblem qp = {n, m, ms, H_ptr, f_ptr, A_ptr, bu_ptr, bl_ptr, sense_ptr, bp_ptr, n_tasks};
    daqp_quadprog(&result, &qp, nullptr);

    result.set_slack(compute_slacks(A, bu, bl, result.get_primal()));
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
  : is_solved_{false}
  , is_slack_computed_{false}
  , max_variables_{max_variables}
  , max_constraints_{max_constraints}
  , max_constraints_in_level_{max_constraints_in_level} {
    allocate_daqp_workspace(&work_, max_variables, max_constraints_in_level);
    allocate_daqp_ldp(&work_, max_variables, max_constraints, 0, 0, 0);
    daqp_default_settings(&settings_);
    work_.settings = &settings_;
}

DAQP::~DAQP() {
    work_.settings = nullptr;
    free_daqp_workspace(&work_);
    free_daqp_ldp(&work_);
}

int DAQP::resize_result(int const n, int const m, Eigen::VectorXi const& break_points) {
    // Ensure that the new dimensions is not too large
    if (n > max_variables_) {
        return 1;
    }
    if (m > max_constraints_) {
        return 2;
    }

    // Check if maximum soft constraints exceeds limit
    int ns = 0;
    for (int i = 0, prev = 0; i < break_points.size(); i++) {
        ns   = (ns > break_points(i) - prev) ? ns : break_points(i) - prev;
        prev = break_points(i);
    }
    if (ns > max_constraints_in_level_) {
        return 3;
    }

    // Resize primal and dual variables
    if (result_.x == nullptr) {
        result_.resize_primal(n);
        result_.resize_dual(m);
    }
    if (n != work_.n) {
        result_.resize_primal(n);
    }
    if (m != work_.m) {
        result_.resize_dual(m);
    }
    return 0;
}

// Update the workspace with a new problem
int DAQP::update(Eigen::MatrixXd const& H,
                 Eigen::VectorXd const& f,
                 Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
                 Eigen::VectorXd const& bu,
                 Eigen::VectorXd const& bl,
                 Eigen::VectorXi const& sense,
                 Eigen::VectorXi const& break_points,
                 int update_mask) {
    int n       = A.cols();
    int m       = bu.size();
    int ms      = m - A.rows();
    int n_tasks = break_points.size();

    assert(bu.size() == bl.size());
    assert(ms <= n);
    assert(n_tasks == 0 || break_points(Eigen::last) == m);

    // Get correct pointers (using const_cast because DAQP C API expects non-const pointers but doesn't modify the data)
    double* H_ptr  = H.size() == 0 ? nullptr : const_cast<double*>(H.data());
    double* f_ptr  = f.size() == 0 ? nullptr : const_cast<double*>(f.data());
    double* A_ptr  = A.size() == 0 ? nullptr : const_cast<double*>(A.data());
    double* bu_ptr = bu.size() == 0 ? nullptr : const_cast<double*>(bu.data());
    double* bl_ptr = bl.size() == 0 ? nullptr : const_cast<double*>(bl.data());
    int* sense_ptr = sense.size() == 0 ? nullptr : const_cast<int*>(sense.data());
    int* bp_ptr    = break_points.size() == 0 ? nullptr : const_cast<int*>(break_points.data());

    //
    int resize_status = resize_result(n, m, break_points);
    assert(resize_status == 0);

    qp_ = {n, m, ms, H_ptr, f_ptr, A_ptr, bu_ptr, bl_ptr, sense_ptr, bp_ptr, n_tasks};

    if (update_mask < 0) {
        // Assume that everythig should be updated
        update_mask = UPDATE_Rinv + UPDATE_M + UPDATE_v + UPDATE_d + UPDATE_sense + UPDATE_hierarchy;
    }

    int status = update_ldp(update_mask, &work_, &qp_);
    is_solved_ = is_slack_computed_ = false;
    return status;
}

// Solve the LDP that is in the workspace
EigenDAQPResult const& DAQP::solve() {
    if (!is_solved_) {
        daqp_solve(&result_, &work_);
        is_solved_ = true;
    }
    return result_;
}

// Solve min_x ||x||
//        s.t  bl <= A*x <= bu
EigenDAQPResult const& DAQP::solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
                                   Eigen::VectorXd const& bu,
                                   Eigen::VectorXd const& bl,
                                   Eigen::VectorXi const& break_points) {
    Eigen::MatrixXd H(0, 0);
    Eigen::VectorXd f(0);
    Eigen::VectorXi sense(0);
    update(H, f, A, bu, bl, sense, break_points, UPDATE_M + UPDATE_d + UPDATE_sense + UPDATE_hierarchy);
    return solve();
}


void DAQP::set_primal_tol(double val) {
    settings_.primal_tol = val;
}
void DAQP::set_dual_tol(double val) {
    settings_.dual_tol = val;
}
void DAQP::set_zero_tol(double val) {
    settings_.zero_tol = val;
}
void DAQP::set_pivot_tol(double val) {
    settings_.pivot_tol = val;
}
void DAQP::set_progress_tol(double val) {
    settings_.progress_tol = val;
}
void DAQP::set_cycle_tol(int val) {
    settings_.cycle_tol = val;
}
void DAQP::set_iter_limit(int val) {
    settings_.iter_limit = val;
}
void DAQP::set_fval_bound(double val) {
    settings_.fval_bound = val;
}
void DAQP::set_eps_prox(double val) {
    settings_.eps_prox = val;
}
void DAQP::set_eta_prox(double val) {
    settings_.eta_prox = val;
}
void DAQP::set_rho_soft(double val) {
    settings_.rho_soft = val;
}
void DAQP::set_rel_subopt(double val) {
    settings_.rel_subopt = val;
}
void DAQP::set_abs_subopt(double val) {
    settings_.abs_subopt = val;
}
void DAQP::set_sing_tol(double val) {
    settings_.sing_tol = val;
}
void DAQP::set_refactor_tol(double val) {
    settings_.refactor_tol = val;
}


Eigen::VectorXd DAQP::get_primal() {
    return result_.get_primal();
}
Eigen::VectorXd DAQP::get_dual() {
    return result_.get_dual();
}
Eigen::VectorXd DAQP::get_slack() {
    if (!is_solved_) {
        solve();
    }

    if (!is_slack_computed_) {
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix(
          qp_.A, qp_.m, qp_.n);
        Eigen::Map<const Eigen::VectorXd> upper(qp_.bupper, qp_.m);
        Eigen::Map<const Eigen::VectorXd> lower(qp_.blower, qp_.m);

        result_.set_slack(compute_slacks(matrix, upper, lower, result_.get_primal()));
        is_slack_computed_ = true;
    }

    return result_.get_slack();
}
Eigen::VectorXi DAQP::get_active_set() {
    return result_.get_active_set();
}
int DAQP::get_status() {
    return result_.exitflag;
}
int DAQP::get_iterations() {
    return result_.iter;
}
double DAQP::get_solve_time() {
    return result_.solve_time;
}
