#include "daqp.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

// Automatic equality reduction applies to solves that start from scratch: a
// DAQP object without warm start reduces its problem on every update (as with
// DAQP_EQ_REDUCTION_ON), while a warm-started one solves the full problem (as
// with DAQP_EQ_REDUCTION_OFF). The reduced and the full problem round
// differently, which identifies which of them was solved.

namespace {

constexpr int n = 40, neq = 16, nineq = 20, m = neq + nineq;
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

struct Data {
    Eigen::MatrixXd H;
    Eigen::VectorXd f;
    RowMatrix A;
    Eigen::VectorXd bu, bl;
};

// A sequence of problems where H, A and the bounds all change
std::vector<Data> make_sequence(int length) {
    std::mt19937 rng(17);
    std::normal_distribution<double> N(0, 1);
    std::vector<Data> seq;
    for (int t = 0; t < length; ++t) {
        Data d;
        Eigen::MatrixXd L(n, n);
        for (int i = 0; i < n * n; ++i) L.data()[i] = N(rng);
        d.H = L * L.transpose() / n + Eigen::MatrixXd::Identity(n, n);
        d.f = Eigen::VectorXd::NullaryExpr(n, [&]() { return N(rng); });
        d.A = RowMatrix::NullaryExpr(m, n, [&]() { return N(rng); });
        Eigen::VectorXd x0 = Eigen::VectorXd::NullaryExpr(n, [&]() { return 0.3 * N(rng); });
        Eigen::VectorXd Ax = d.A * x0;
        d.bu = Ax.array() + 0.2;
        d.bl = Eigen::VectorXd::Constant(m, -DAQP_INF);
        d.bu.head(neq) = Ax.head(neq);
        d.bl.head(neq) = Ax.head(neq);
        seq.push_back(d);
    }
    return seq;
}

std::vector<Eigen::VectorXd> run(const std::vector<Data>& seq, bool warm, int policy) {
    DAQP solver(n, m, 0);
    if (warm) solver.set_warm_start();
    solver.set_eq_reduction(policy);
    Eigen::VectorXi sense = Eigen::VectorXi::Zero(m);
    sense.head(neq).setConstant(DAQP_ACTIVE | DAQP_IMMUTABLE);
    const Eigen::VectorXi no_hierarchy(0), reuse_sense(0);
    std::vector<Eigen::VectorXd> xs;
    for (const Data& d : seq) {
        // A warm-started update reuses the previous working set
        const Eigen::VectorXi& s = (warm && !xs.empty()) ? reuse_sense : sense;
        assert(solver.update(d.H, d.f, d.A, d.bu, d.bl, s, no_hierarchy) >= 0);
        const EigenDAQPResult& res = solver.solve();
        assert(res.exitflag == DAQP_EXIT_OPTIMAL);
        const Eigen::VectorXd x = solver.get_primal();
        // The equality constraints hold
        assert(((d.A.topRows(neq) * x) - d.bu.head(neq)).cwiseAbs().maxCoeff() < 1e-8);
        xs.push_back(x);
    }
    return xs;
}

} // namespace

int main() {
    const std::vector<Data> seq = make_sequence(4);
    for (bool warm : {false, true}) {
        const auto x_auto = run(seq, warm, DAQP_EQ_REDUCTION_AUTO);
        const auto x_on = run(seq, warm, DAQP_EQ_REDUCTION_ON);
        const auto x_off = run(seq, warm, DAQP_EQ_REDUCTION_OFF);
        bool on_differs_from_off = false;
        for (size_t t = 0; t < seq.size(); ++t) {
            // Cold: AUTO reduces like ON. Warm: AUTO solves the full problem.
            assert(x_auto[t] == (warm ? x_off[t] : x_on[t]));
            // The reduced and the full problem give the same solution...
            assert((x_on[t] - x_off[t]).cwiseAbs().maxCoeff() < 1e-8);
            // ...but not bit for bit, so the check above tells them apart
            if (x_on[t] != x_off[t]) on_differs_from_off = true;
        }
        assert(on_differs_from_off);
    }
}
