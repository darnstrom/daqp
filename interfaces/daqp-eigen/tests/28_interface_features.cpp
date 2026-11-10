#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <daqp.hpp>

// The class copies an A that is not row-major (and does not refer to a
// temporary), solves mixed-integer problems, and reports errors; one-shot
// solves take settings

using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

static int failures = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (!(cond)) {                                                           \
            std::cerr << __FILE__ << ":" << __LINE__ << ": " #cond << std::endl; \
            failures++;                                                          \
        }                                                                        \
    } while (0)

static bool near(Eigen::VectorXd const& a, Eigen::VectorXd const& b) {
    return a.size() == b.size() && (a - b).lpNorm<Eigen::Infinity>() <= 1e-8;
}

int main() {
    Eigen::MatrixXd H(3, 3);
    H << 4, 1, 0,
         1, 3, 1,
         0, 1, 2;
    Eigen::VectorXd f(3);
    f << 10, -8, 3;
    Eigen::MatrixXd A_col(2, 3); // Column-major
    A_col << 1, 1, 1,
             1, -1, 2;
    Eigen::VectorXd bu(5), bl(5);
    bu << 1, 1, 1, 2, 0.5;
    bl << -1, -1, -1, -2, -0.5;
    const RowMatrixXd A_row = A_col;
    const Eigen::VectorXi none;

    const EigenDAQPResult ref = daqp_solve(H, f, A_row, bu, bl);
    CHECK(ref.exitflag == DAQP_EXIT_OPTIMAL);

    // Settings for one-shot solves
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.iter_limit = 1;
    CHECK(daqp_solve(H, f, A_row, bu, bl, &settings).exitflag == DAQP_EXIT_ITERLIMIT);

    // The class does not refer to a column-major A (or an expression) after the update
    {
        DAQP solver(3, 5, 0);
        Eigen::MatrixXd A_tmp = A_col;
        CHECK(solver.update(H, f, A_tmp, bu, bl, none, none) >= 0);
        A_tmp.setZero(); // Changed (or freed) by the caller
        CHECK(near(solver.get_primal(), ref.get_primal()));
        CHECK(near(solver.get_slack(), ref.get_slack()));
        CHECK(std::abs(solver.get_objective() - ref.fval) < 1e-8);
        CHECK(solver.update(H, f, 2.0 * A_col - A_col, bu, bl, none, none) >= 0);
        CHECK(near(solver.get_primal(), ref.get_primal()));
    }

    // Errors of the interface
    {
        DAQP solver(2, 5, 0);
        CHECK(solver.update(H, f, A_row, bu, bl, none, none) == DAQP_EIGEN_TOO_MANY_VARIABLES);
        DAQP small(3, 4, 0);
        CHECK(small.update(H, f, A_row, bu, bl, none, none) == DAQP_EIGEN_TOO_MANY_CONSTRAINTS);
        DAQP no_soft(3, 5, 0);
        const Eigen::VectorXi soft = Eigen::VectorXi::Constant(5, DAQP_SOFT);
        CHECK(no_soft.update(H, f, A_row, bu, bl, soft, none) == DAQP_EIGEN_TOO_MANY_IN_LEVEL);
        CHECK(no_soft.update(H, f, A_row, bu, bl, none, none) >= 0);
        CHECK(near(no_soft.get_primal(), ref.get_primal()));
    }

    // Branch and bound in the class: x1, x2 binary
    {
        Eigen::MatrixXd Hb = Eigen::MatrixXd::Identity(3, 3);
        Eigen::VectorXd fb(3);
        fb << -0.4, -0.6, 1;
        RowMatrixXd Ab(1, 3);
        Ab << 1, 1, -1;
        Eigen::VectorXd bub(4), blb(4);
        bub << 1, 1, 10, 1.5;
        blb << 0, 0, -10, -10;
        Eigen::VectorXi sense = Eigen::VectorXi::Zero(4);
        sense(0) = sense(1) = DAQP_BINARY;
        const EigenDAQPResult ref_b = daqp_solve(Hb, fb, Ab, bub, blb, sense, none);
        CHECK(ref_b.exitflag == DAQP_EXIT_OPTIMAL);

        DAQP solver(3, 4, 0);
        CHECK(solver.update(Hb, fb, Ab, bub, blb, sense, none) >= 0);
        CHECK(solver.get_status() == DAQP_EXIT_OPTIMAL);
        CHECK(near(solver.get_primal(), ref_b.get_primal()));
        CHECK(solver.get_nodes() > 0);
        for (int i = 0; i < 2; i++) {
            const double xi = solver.get_primal()(i);
            CHECK(std::abs(xi) < 1e-9 || std::abs(xi - 1) < 1e-9);
        }
        // Without binaries, the relaxation is solved, and with them again the MIQP
        CHECK(solver.update(Hb, fb, Ab, bub, blb, none, none) >= 0);
        CHECK(near(solver.get_primal(), daqp_solve(Hb, fb, Ab, bub, blb).get_primal()));
        CHECK(solver.update(Hb, fb, Ab, bub, blb, sense, none) >= 0);
        CHECK(near(solver.get_primal(), ref_b.get_primal()));
    }

    // Soft constraints: the largest violation is reported
    {
        Eigen::VectorXd bu_s = bu, bl_s = bl;
        bl_s(3) = bu_s(3) = 2.5; // x1+x2+x3 = 2.5
        bl_s(4) = bu_s(4) = -4;  // Infeasible together with the constraints above
        Eigen::VectorXi sense = Eigen::VectorXi::Zero(5);
        sense(4) = DAQP_SOFT;
        DAQP solver(3, 5, 1);
        CHECK(solver.update(H, f, A_row, bu_s, bl_s, sense, none) >= 0);
        CHECK(solver.get_status() > 0);
        CHECK(solver.get_soft_slack() > 0);
    }

    if (failures > 0) {
        std::cerr << failures << " checks failed" << std::endl;
        return 1;
    }
    return 0;
}
