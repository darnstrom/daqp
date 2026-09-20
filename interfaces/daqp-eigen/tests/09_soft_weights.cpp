#include <cmath>
#include <Eigen/Dense>
#include <daqp.hpp>

int main() {
    DAQP solver(1, 1, 1);
    solver.set_rho_soft(0.5);
    solver.set_w_soft(2.0);

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
#endif

    return 0;
}
