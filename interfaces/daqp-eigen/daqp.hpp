#ifndef DAQP_HPP
#define DAQP_HPP

#include <Eigen/Dense>
#include "api.h"

// Errors returned by DAQP::update (the solver's exit flags are > -100)
constexpr int DAQP_EIGEN_TOO_MANY_VARIABLES   = -101; // n > max_variables
constexpr int DAQP_EIGEN_TOO_MANY_CONSTRAINTS = -102; // m > max_constraints
constexpr int DAQP_EIGEN_TOO_MANY_IN_LEVEL    = -103; // Soft constraints/level > max_constraints_in_level


class EigenDAQPResult : public DAQPResult {
  private:
    Eigen::VectorXd x_;
    Eigen::VectorXd lam_;
    Eigen::VectorXd slack_;

  public:
    EigenDAQPResult();
    EigenDAQPResult(int n, int m);
    void resize_primal(int n);
    void resize_dual(int m);
    Eigen::VectorXd const& get_primal() const;
    Eigen::VectorXd const& get_dual() const;
    Eigen::VectorXd const& get_slack() const;
    Eigen::VectorXi get_active_set() const;
    void set_slack(Eigen::VectorXd const& slack);
    void set_slack(Eigen::VectorXd&& slack);
};


Eigen::VectorXd compute_slacks(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& matrix,
                               Eigen::VectorXd const& upper,
                               Eigen::VectorXd const& lower,
                               Eigen::VectorXd const& primal);


EigenDAQPResult daqp_solve(Eigen::MatrixXd const& H,
                           Eigen::VectorXd const& f,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
                           Eigen::VectorXd const& bu,
                           Eigen::VectorXd const& bl,
                           Eigen::VectorXi const& sense,
                           Eigen::VectorXi const& break_points,
                           DAQPSettings* settings = nullptr);

// Solve: min_x 0.5 x'*H*x + f'*x
//         s.t  bl<= A*x <= bu
EigenDAQPResult daqp_solve(Eigen::MatrixXd const& H,
                           Eigen::VectorXd const& f,
                           Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
                           Eigen::VectorXd const& bu,
                           Eigen::VectorXd const& bl,
                           DAQPSettings* settings = nullptr);

// Solve min_x ||x||
//        s.t  bl <= A*x <= bu
EigenDAQPResult daqp_solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
                           Eigen::VectorXd const& bu,
                           Eigen::VectorXd const& bl,
                           Eigen::VectorXi const& break_points,
                           DAQPSettings* settings = nullptr);


class DAQP {
  private:
    bool is_solved_;
    bool is_slack_computed_;
    bool warm_start_;
    int max_variables_;
    int max_constraints_;
    int max_constraints_in_level_;
    DAQPWorkspace work_;
    DAQPSettings settings_;
    EigenDAQPResult result_;
    DAQPProblem qp_;
    int resize_result(int n, int m, Eigen::VectorXi const& sense, Eigen::VectorXi const& break_points);
    int update_bnb(int const* sense, int n, int m);
    c_float* rinv_buf_;
    c_float* v_buf_;
    // Row-major copy of A if it is given in another form (the workspace refers to it)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A_buf_;

  public:
    DAQP(int max_variables, int max_constraints, int max_constraints_in_level);
    DAQP(DAQP const&) = delete;
    DAQP& operator=(DAQP const&) = delete;
    // Returns < 0 if the update failed (see DAQP_EIGEN_* for errors of the interface).
    // The data is referenced, not copied, so it must stay alive and unchanged until
    // the next update (except for an A that is not row-major, which is copied)
    int update(Eigen::MatrixXd const& H,
               Eigen::VectorXd const& f,
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
               Eigen::VectorXd const& bu,
               Eigen::VectorXd const& bl,
               Eigen::VectorXi const& sense,
               Eigen::VectorXi const& break_points,
               int update_mask = -1);
    // A in any other form (e.g., column-major) is copied into the object
    template <class Derived>
    int update(Eigen::MatrixXd const& H,
               Eigen::VectorXd const& f,
               Eigen::MatrixBase<Derived> const& A,
               Eigen::VectorXd const& bu,
               Eigen::VectorXd const& bl,
               Eigen::VectorXi const& sense,
               Eigen::VectorXi const& break_points,
               int update_mask = -1) {
        A_buf_ = A;
        return update(H, f, A_buf_, bu, bl, sense, break_points, update_mask);
    }
    ~DAQP();
    EigenDAQPResult const& solve();
    EigenDAQPResult const& solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
                                 Eigen::VectorXd const& bu,
                                 Eigen::VectorXd const& bl,
                                 Eigen::VectorXi const& break_points);
    template <class Derived>
    EigenDAQPResult const& solve(Eigen::MatrixBase<Derived> const& A,
                                 Eigen::VectorXd const& bu,
                                 Eigen::VectorXd const& bl,
                                 Eigen::VectorXi const& break_points) {
        A_buf_ = A;
        return solve(A_buf_, bu, bl, break_points);
    }

    // Setters for settings
    void set_warm_start();
    void set_cold_start();
    void set_primal_tol(double val);
    void set_dual_tol(double val);
    void set_zero_tol(double val);
    void set_pivot_tol(double val);
    void set_progress_tol(double val);
    void set_cycle_tol(int val);
    void set_iter_limit(int val);
    void set_fval_bound(double val);
    void set_eps_prox(double val);
    void set_eta_prox(double val);
    void set_rho_soft(double val);
    void set_w_soft(double val);
    void set_eq_reduction(int val);
    void set_presolve(int val);
    // Weights of the soft constraints, one element per constraint (an empty
    // vector leaves that weight at settings.rho_soft/settings.w_soft). False
    // if DAQP was built without support for individual weights.
    bool set_soft_weights(Eigen::VectorXd const& rho_lower = Eigen::VectorXd(),
                          Eigen::VectorXd const& rho_upper = Eigen::VectorXd(),
                          Eigen::VectorXd const& w_lower = Eigen::VectorXd(),
                          Eigen::VectorXd const& w_upper = Eigen::VectorXd());
    void set_rel_subopt(double val);
    void set_abs_subopt(double val);
    void set_sing_tol(double val);
    void set_refactor_tol(double val);
    void set_time_limit(double val);

    // Getters for result
    Eigen::VectorXd const& get_primal();
    Eigen::VectorXd const& get_dual();
    Eigen::VectorXd const& get_slack();
    Eigen::VectorXi get_active_set();
    int get_status();
    int get_iterations();
    int get_nodes();
    double get_objective();
    double get_soft_slack();
    double get_solve_time();
};

#endif // DAQP_HPP
