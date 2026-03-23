#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>
#include <cmath>

// Compute KKT stationarity residual: ||H*x + f + A'*lam||_inf
// For DAQP with no simple bounds (ms=0), lam is a vector of length nc (general constraints).
// Sign convention: lam[i] > 0 if upper bound active, lam[i] < 0 if lower bound active,
// lam[i] = 0 if inactive.  This matches standard KKT for a minimization problem.
static double kkt_stationarity(Eigen::MatrixXd const& H,
                               Eigen::VectorXd const& f,
                               Eigen::MatrixXd const& A,
                               Eigen::VectorXd const& x,
                               Eigen::VectorXd const& lam)
{
    Eigen::VectorXd res = H * x + f + A.transpose() * lam;
    return res.lpNorm<Eigen::Infinity>();
}

// Compute primal feasibility residual: max violation of A*x <= bu
static double primal_feasibility(Eigen::MatrixXd const& A,
                                 Eigen::VectorXd const& bu,
                                 Eigen::VectorXd const& x)
{
    Eigen::VectorXd viol = (A * x - bu).cwiseMax(0.0);
    return viol.lpNorm<Eigen::Infinity>();
}

int main() {
    // -----------------------------------------------------------------------
    // Test: QP with near-linearly-dependent active constraints
    //
    // The active-set Gram matrix M*M' has one large eigenvalue (~n) and
    // n-1 small eigenvalues O(eps^2).  This stresses the LDL factorisation:
    // the D diagonal shrinks to O(eps^2) for constraints 2..n.  The global-
    // minimum pivot search moves these small-D entries to the tail, keeping
    // numerical errors from propagating into the earlier (larger-D) rows.
    //
    // Problem structure
    //   min  0.5*||x||^2  -  1'*x          (f = -ones, H = I)
    //   s.t. a_i'*x <= 0  for i = 1..n     (n upper-bound constraints)
    //
    // where  a_i = e_1 + eps * e_i  (normalised).
    //
    // Feasibility forces x_0 towards 0.  All n constraints are simultaneously
    // active at the optimum x* = [0, 0, ..., 0], giving a rank-1 active-set
    // matrix perturbed by O(eps) — precisely the ill-conditioned scenario.
    // -----------------------------------------------------------------------

    const int    n      = 5;
    const double eps    = 1e-4;   // perturbation; D_2..n ~ eps^2 ~ 1e-8
    const double kkt_tol = 1e-6; // primal_tol default

    // H = I, f = -ones  (unconstrained min: x* = ones, infeasible due to a_i'x<=0)
    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(n, n);
    Eigen::VectorXd f = -Eigen::VectorXd::Ones(n);

    // Build near-parallel constraint normals: a_i = e_1 + eps*e_i, normalised
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(n, n);
    A.setZero();
    for (int i = 0; i < n; i++) {
        A(i, 0) = 1.0;
        A(i, i) += eps;          // e_1 + eps*e_i  (diagonal add: when i=0, norm = 1+eps)
        A.row(i).normalize();
    }

    // Upper-bound constraints only: bl = -inf, bu = 0
    Eigen::VectorXd bu    = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd bl    = Eigen::VectorXd::Constant(n, -DAQP_INF);
    Eigen::VectorXi sense = Eigen::VectorXi::Zero(n);
    Eigen::VectorXi bp    = Eigen::VectorXi::Zero(0);

    EigenDAQPResult result = daqp_solve(H, f, A, bu, bl, sense, bp);

    if (result.exitflag != DAQP_EXIT_OPTIMAL) {
        std::cerr << "Solver did not return optimal: exitflag = "
                  << result.exitflag << std::endl;
        return 1;
    }

    Eigen::VectorXd x   = result.get_primal();
    Eigen::VectorXd lam = result.get_dual();   // length n, DAQP sign convention

    double kkt   = kkt_stationarity(H, f, A, x, lam);
    double pfeas = primal_feasibility(A, bu, x);

    std::cout << "Near-linearly-dependent active-constraint QP (eps=" << eps << ", n=" << n << ")\n";
    std::cout << "  x*:              " << x.transpose()   << "\n";
    std::cout << "  lam*:            " << lam.transpose() << "\n";
    std::cout << "  KKT stationarity (inf-norm): " << kkt   << "\n";
    std::cout << "  Primal feasibility (inf-norm): " << pfeas << "\n";

    if (kkt > kkt_tol) {
        std::cerr << "FAIL: KKT stationarity residual " << kkt
                  << " exceeds tolerance " << kkt_tol << std::endl;
        return 1;
    }
    if (pfeas > kkt_tol) {
        std::cerr << "FAIL: Primal feasibility violation " << pfeas
                  << " exceeds tolerance " << kkt_tol << std::endl;
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
