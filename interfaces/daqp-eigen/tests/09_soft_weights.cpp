#include <cmath>
#include <Eigen/Dense>
#include <daqp.hpp>

int main() {
    DAQP solver(1, 1, 1);
    solver.set_rho_soft(0.5);
    solver.set_w_soft(2.0);
    // Enabling warm start before the first update must not read an
    // uninitialized previous result or erase the soft classification.
    solver.set_warm_start();

    Eigen::MatrixXd H(1, 1);
    H << 1.0;
    Eigen::VectorXd f(1);
    f << -10.0;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(0, 1);
    Eigen::VectorXd bu(1), bl(1);
    bu << 0.0;
    bl << -DAQP_INF;
    Eigen::VectorXi sense(1), break_points(0);
    sense << DAQP_SOFT;

    if (solver.update(H, f, A, bu, bl, sense, break_points) != 0)
        return 1;
    if (std::abs(solver.solve().get_primal()[0] - 8.0/3.0) > 1e-10)
        return 1;

    // Changing uniform weights after a solve must refresh the active soft row.
    solver.set_w_soft(3.0);
    if (std::abs(solver.solve().get_primal()[0] - 7.0/3.0) > 1e-10)
        return 1;
    solver.set_w_soft(2.0);
    solver.set_rho_soft(0.25);
    if (std::abs(solver.solve().get_primal()[0] - 1.6) > 1e-10)
        return 1;
    solver.set_rho_soft(0.5);
    if (std::abs(solver.solve().get_primal()[0] - 8.0/3.0) > 1e-10)
        return 1;

    // Updating a weight after solving must rebuild the active-set
    // factorization. The optimum becomes rho*(10-w)/(1+rho) = 8/11.
    Eigen::VectorXd rho_upper(1);
    rho_upper << 0.1;
    if (!solver.set_soft_weights(Eigen::VectorXd(), rho_upper))
        return 1;
    if (std::abs(solver.solve().get_primal()[0] - 8.0/11.0) > 1e-10)
        return 1;

    // A subsequent warm update must retain DAQP_SOFT. With f=-9 the mixed
    // L1/L2 optimum is rho*(9-w)/(1+rho) = 7/11, rather than the hard x=0.
    f[0] = -9.0;
    const Eigen::VectorXi reuse_sense(0);
    if (solver.update(H, f, A, bu, bl, reuse_sense, break_points) != 0)
        return 1;
    if (std::abs(solver.solve().get_primal()[0] - 7.0/11.0) > 1e-10)
        return 1;

    // A soft constraint can be active while its L1 slack remains fixed at
    // zero (|lambda| < w). Warm-starting must preserve that state rather than
    // restarting the row on the nonzero-slack branch.
    DAQP fixed_slack_solver(1, 1, 1);
    fixed_slack_solver.set_rho_soft(0.5);
    fixed_slack_solver.set_w_soft(2.0);
    fixed_slack_solver.set_warm_start();
    f[0] = -1.0;
    sense[0] = DAQP_SOFT;
    if (fixed_slack_solver.update(H, f, A, bu, bl, sense, break_points) != 0)
        return 1;
    if (std::abs(fixed_slack_solver.solve().get_primal()[0]) > 1e-10)
        return 1;
    f[0] = -1.1;
    if (fixed_slack_solver.update(H, f, A, bu, bl, reuse_sense, break_points) != 0)
        return 1;
    if (std::abs(fixed_slack_solver.solve().get_primal()[0]) > 1e-10)
        return 1;
    if (fixed_slack_solver.get_iterations() != 1)
        return 1;

    // Explicit sense overrides the stored state; make the row hard again.
    Eigen::VectorXi hard_sense = Eigen::VectorXi::Zero(1);
    f[0] = -10.0;
    if (fixed_slack_solver.update(H, f, A, bu, bl, hard_sense,
                                  break_points) != 0)
        return 1;
    if (std::abs(fixed_slack_solver.solve().get_primal()[0]) > 1e-10)
        return 1;

    // C can change the settings directly and refresh the soft factorization.
    double h_c = 1.0, f_c = -10.0, bu_c = 0.0, bl_c = -DAQP_INF;
    double x_c = 0.0, lam_c = 0.0;
    int sense_c = DAQP_SOFT;
    DAQPProblem qp_c{1, 1, 1, &h_c, &f_c, nullptr, &bu_c, &bl_c,
                     &sense_c, nullptr, 0, 0};
    DAQPSettings settings_c;
    daqp_default_settings(&settings_c);
    settings_c.rho_soft = 0.5;
    DAQPWorkspace work_c{};
    work_c.settings = &settings_c;
    if (setup_daqp(&qp_c, &work_c, nullptr) <= 0) return 1;
    DAQPResult result_c{};
    result_c.x = &x_c;
    result_c.lam = &lam_c;
    daqp_solve(&result_c, &work_c);
    if (std::abs(x_c - 10.0/3.0) > 1e-10) return 1;
    settings_c.w_soft = 2.0;
    daqp_refresh_soft_weights(&work_c);
    daqp_solve(&result_c, &work_c);
    if (std::abs(x_c - 8.0/3.0) > 1e-10) return 1;
    settings_c.rho_soft = 0.25;
    daqp_refresh_soft_weights(&work_c);
    daqp_solve(&result_c, &work_c);
    if (std::abs(x_c - 1.6) > 1e-10) return 1;
    work_c.settings = nullptr;
    free_daqp_workspace(&work_c);
    free_daqp_ldp(&work_c);

    return 0;
}
