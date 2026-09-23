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

    // Updating a weight after solving must rebuild the active-set
    // factorization. The optimum becomes rho*(10-w)/(1+rho) = 8/11.
    Eigen::VectorXd rho_upper(1);
    rho_upper << 0.1;
#ifdef DAQP_NO_SOFT_WEIGHTS
    return solver.set_soft_weights(Eigen::VectorXd(), rho_upper) ? 1 : 0;
#else
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
#endif

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

#ifndef DAQP_NO_SOFT_WEIGHTS
    // Until the first solve, an update leaves the equality-reduced problem in
    // the workspace. The weights must still have one entry per constraint of
    // the original problem, and give the same solution as setting them after
    // a solve.
    const int n_eq = 40, neq = 16, max_m = 60;
    const auto eq_problem = [&](int nsoft, Eigen::MatrixXd& He,
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::RowMajor>& Ae,
            Eigen::VectorXd& fe, Eigen::VectorXd& bue, Eigen::VectorXd& ble,
            Eigen::VectorXi& sensee) {
        const int m = neq + nsoft;
        He = Eigen::MatrixXd::Constant(n_eq, n_eq, 0.01);
        He.diagonal().setConstant(2.0);
        fe = Eigen::VectorXd::LinSpaced(n_eq, -5.0, 5.0);
        Ae = Eigen::MatrixXd::Zero(m, n_eq);
        bue = Eigen::VectorXd::Constant(m, 0.5);
        ble = Eigen::VectorXd::Constant(m, -0.5);
        sensee = Eigen::VectorXi::Constant(m, DAQP_SOFT);
        for (int i = 0; i < m; ++i) Ae(i, i % n_eq) = 1.0;
        for (int i = 0; i < neq; ++i) {
            bue[i] = ble[i] = 0.1;
            sensee[i] = DAQP_ACTIVE | DAQP_IMMUTABLE;
        }
    };
    Eigen::MatrixXd He;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Ae;
    Eigen::VectorXd fe, bue, ble;
    Eigen::VectorXi sensee;
    for (int nsoft : {20, max_m - neq}) {
        eq_problem(nsoft, He, Ae, fe, bue, ble, sensee);
        const int m = neq + nsoft;
        const Eigen::VectorXd rho_eq = Eigen::VectorXd::Constant(m, 0.1);

        DAQP after_solve(n_eq, max_m, max_m);
        if (after_solve.update(He, fe, Ae, bue, ble, sensee, break_points) != 0)
            return 1;
        after_solve.solve();
        if (!after_solve.set_soft_weights(Eigen::VectorXd(), rho_eq))
            return 1;
        const Eigen::VectorXd x_ref = after_solve.solve().get_primal();

        DAQP before_solve(n_eq, max_m, max_m);
        // Allocate the weights for a smaller problem, and then grow it
        if (nsoft != 20) {
            Eigen::MatrixXd Hs;
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::RowMajor> As;
            Eigen::VectorXd fs, bus, bls;
            Eigen::VectorXi senses;
            eq_problem(20, Hs, As, fs, bus, bls, senses);
            if (before_solve.update(Hs, fs, As, bus, bls, senses,
                                    break_points) != 0)
                return 1;
            if (!before_solve.set_soft_weights(Eigen::VectorXd(),
                        Eigen::VectorXd::Constant(neq + 20, 0.1)))
                return 1;
        }
        if (before_solve.update(He, fe, Ae, bue, ble, sensee, break_points) != 0)
            return 1;
        // The reduced number of constraints is rejected
        if (before_solve.set_soft_weights(Eigen::VectorXd(),
                    Eigen::VectorXd::Constant(m - neq, 0.1)))
            return 1;
        if (!before_solve.set_soft_weights(Eigen::VectorXd(), rho_eq))
            return 1;
        if ((before_solve.solve().get_primal() - x_ref).norm() > 1e-8)
            return 1;
    }
#endif

    return 0;
}
