#ifndef DAQP_HPP
#define DAQP_HPP

#include <Eigen/Dense>
#include "api.h"


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
    Eigen::VectorXd get_primal() const;
    Eigen::VectorXd get_dual() const;
    Eigen::VectorXd get_slack() const;
    Eigen::VectorXi get_active_set() const;
    void set_slack(Eigen::VectorXd const& slack);
};


Eigen::VectorXd compute_slacks(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& matrix,
                               Eigen::VectorXd const& upper,
                               Eigen::VectorXd const& lower,
                               Eigen::VectorXd const& primal);


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
    bool is_solved_;
    bool is_slack_computed_;
    int max_variables_;
    int max_constraints_;
    int max_constraints_in_level_;
    DAQPWorkspace work_;
    DAQPSettings settings_;
    EigenDAQPResult result_;
    DAQPProblem qp_;
    int resize_result(int n, int m, Eigen::VectorXi const& break_points);

  public:
    DAQP(int max_variables, int max_constraints, int max_constraints_in_level);
    int update(Eigen::MatrixXd const& H,
               Eigen::VectorXd const& f,
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
               Eigen::VectorXd const& bu,
               Eigen::VectorXd const& bl,
               Eigen::VectorXi const& sense,
               Eigen::VectorXi const& break_points,
               int update_mask = -1);
    ~DAQP();
    EigenDAQPResult const& solve();
    EigenDAQPResult const& solve(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> const& A,
                                 Eigen::VectorXd const& bu,
                                 Eigen::VectorXd const& bl,
                                 Eigen::VectorXi const& break_points);
    // Setters for settings
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
    void set_rel_subopt(double val);
    void set_abs_subopt(double val);
    void set_sing_tol(double val);
    void set_refactor_tol(double val);

    // Getters for result
    Eigen::VectorXd get_primal();
    Eigen::VectorXd get_dual();
    Eigen::VectorXd get_slack();
    Eigen::VectorXi get_active_set();
    int get_status();
    int get_iterations();
    double get_solve_time();
};

#endif // DAQP_HPP
