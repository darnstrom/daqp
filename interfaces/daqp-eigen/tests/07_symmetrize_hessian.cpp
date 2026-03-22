#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>

// Verify that DAQP symmetrizes H before factorization.
//
// Key properties tested:
//   1. A fully symmetric H gives the correct reference solution (regression).
//   2. A lower-triangular H is no longer silently treated as a diagonal matrix
//      (the pre-fix bug: DAQP read only the upper triangle which was all zero
//      except the diagonal, so it ignored the off-diagonal structure).
//   3. A nearly-symmetric H (tiny numerical noise in one triangle) gives the
//      same solution as the perfectly symmetric reference.
int main() {
    double precision = 1e-6;

    // Symmetric positive-definite Hessian (full storage)
    Eigen::MatrixXd H_full(2, 2);
    H_full << 4, 2,
              2, 3;

    Eigen::VectorXd f = (Eigen::VectorXd(2) << 1, -1).finished();

    // Simple bounds only (no general constraints), so ms == m.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(0, 2);
    Eigen::VectorXd bu = (Eigen::VectorXd(2) << 5, 5).finished();
    Eigen::VectorXd bl = (Eigen::VectorXd(2) << -5, -5).finished();
    Eigen::VectorXi sense        = Eigen::VectorXi::Zero(2);
    Eigen::VectorXi break_points = Eigen::VectorXi::Zero(0);

    // 1. Reference: solve with full symmetric H.
    //    Unconstrained optimum: x* = -H^{-1}*f = [-0.625, 0.75]
    EigenDAQPResult ref = daqp_solve(H_full, f, A, bu, bl, sense, break_points);
    if (ref.exitflag <= 0) {
        std::cerr << "Reference solve failed with exitflag " << ref.exitflag << std::endl;
        return 1;
    }

    // 2. Diagonal-only H – used as the "wrong" baseline.
    //    Before the fix, providing a lower-triangular H would produce exactly
    //    the same result as this because the upper triangle was read as zeros.
    Eigen::MatrixXd H_diag = Eigen::MatrixXd::Zero(2, 2);
    H_diag(0, 0) = H_full(0, 0);
    H_diag(1, 1) = H_full(1, 1);
    EigenDAQPResult res_diag = daqp_solve(H_diag, f, A, bu, bl, sense, break_points);
    if (res_diag.exitflag <= 0) {
        std::cerr << "Diagonal solve failed with exitflag " << res_diag.exitflag << std::endl;
        return 1;
    }

    // 3. Lower-triangular H (upper triangle set to zero).
    //    Before the fix: DAQP reads only the upper triangle, sees all-zero
    //    off-diagonals, and treats the matrix as diagonal → wrong solution.
    //    After the fix: DAQP symmetrizes first, so the off-diagonal structure
    //    is preserved (at 0.5× the original value) → solution differs from diagonal.
    Eigen::MatrixXd H_lower = H_full.triangularView<Eigen::Lower>();
    EigenDAQPResult res_lower = daqp_solve(H_lower, f, A, bu, bl, sense, break_points);
    if (res_lower.exitflag <= 0) {
        std::cerr << "Lower-triangle solve failed with exitflag " << res_lower.exitflag
                  << std::endl;
        return 1;
    }

    // 4. Nearly-symmetric H: full H with tiny noise in one off-diagonal entry.
    //    Symmetrization should average the two triangles and give the same
    //    solution as the perfectly symmetric reference.
    Eigen::MatrixXd H_noisy = H_full;
    H_noisy(0, 1) += 1e-8;   // tiny asymmetry
    EigenDAQPResult res_noisy = daqp_solve(H_noisy, f, A, bu, bl, sense, break_points);
    if (res_noisy.exitflag <= 0) {
        std::cerr << "Noisy-H solve failed with exitflag " << res_noisy.exitflag << std::endl;
        return 1;
    }

    // --- Checks ---

    // Full symmetric H gives the reference solution.
    bool ref_ok = ref.get_primal().isApprox(
        (Eigen::VectorXd(2) << -0.625, 0.75).finished(), precision);

    // After symmetrization, lower-triangular H must NOT equal the diagonal result.
    // (Before the fix they would be identical.)
    bool lower_not_diagonal = !res_lower.get_primal().isApprox(
        res_diag.get_primal(), precision);

    // Nearly-symmetric H: symmetrized result matches the reference.
    bool noisy_ok = res_noisy.get_primal().isApprox(ref.get_primal(), 1e-5);

    if (!ref_ok)
        std::cerr << "Reference solution incorrect: "
                  << ref.get_primal().transpose() << std::endl;
    if (!lower_not_diagonal)
        std::cerr << "Lower-triangle H was incorrectly treated as diagonal!\n"
                  << "lower: " << res_lower.get_primal().transpose() << "\n"
                  << "diag:  " << res_diag.get_primal().transpose() << std::endl;
    if (!noisy_ok)
        std::cerr << "Noisy-H result differs from reference:\n"
                  << "noisy: " << res_noisy.get_primal().transpose() << "\n"
                  << "ref:   " << ref.get_primal().transpose() << std::endl;

    return (ref_ok && lower_not_diagonal && noisy_ok) ? 0 : 1;
}
