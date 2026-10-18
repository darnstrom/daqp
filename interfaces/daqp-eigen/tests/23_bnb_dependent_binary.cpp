#include "api.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <vector>

// A binary constraint whose value the equalities determine

namespace {

struct Problem {
    int n, m, ms;
    std::vector<c_float> H, f, A, bu, bl;
    std::vector<int> sense;
    DAQPProblem qp{};

    // Binary x_0 in {-1,1}, -5 <= x_i <= 5 otherwise, and E x = b (small H
    // and large f give large LDP offsets)
    Problem(int n_, const std::vector<std::vector<c_float>>& E, const std::vector<c_float>& b)
        : n(n_), m(n_ + (int)E.size()), ms(n_), H(n * n, 0.0), f(n, 0.0),
          A(E.size() * n, 0.0), bu(m), bl(m), sense(m, 0) {
        for (int i = 0; i < n; ++i) {
            H[i * n + i] = 1e-2 * (1 + i);
            f[i] = 10.0 * (i % 2 ? -1 : 1);
            bu[i] = 5;
            bl[i] = -5;
        }
        sense[0] = DAQP_BINARY;
        bu[0] = 1;
        bl[0] = -1;
        for (size_t k = 0; k < E.size(); ++k) {
            for (int j = 0; j < n; ++j) A[k * n + j] = E[k][j];
            bu[ms + k] = bl[ms + k] = b[k];
        }
        qp = {n, m, ms, H.data(), f.data(), A.data(), bu.data(), bl.data(),
              sense.data(), nullptr, 0, 0};
    }
};

DAQPResult quadprog(Problem p, int policy, std::vector<c_float>& x) {
    std::vector<c_float> lam(p.m);
    DAQPResult res{};
    x.assign(p.n, 0);
    res.x = x.data();
    res.lam = lam.data();
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.eq_reduction = policy;
    daqp_quadprog(&res, &p.qp, &settings);
    return res;
}

// Infeasible unless x0 (= x_0 from the equalities) is -1 or 1
void check(int n, const std::vector<std::vector<c_float>>& E, const std::vector<c_float>& b, c_float x0) {
    const bool feasible = std::fabs(std::fabs(x0) - 1) < 1e-12;
    for (int policy : {DAQP_EQ_REDUCTION_OFF, DAQP_EQ_REDUCTION_AUTO, DAQP_EQ_REDUCTION_ON}) {
        std::vector<c_float> x;
        DAQPResult res = quadprog(Problem(n, E, b), policy, x);
        if (feasible) {
            assert(res.exitflag == DAQP_EXIT_OPTIMAL);
            assert(std::fabs(x[0] - x0) < 1e-6);
            for (size_t k = 0; k < E.size(); ++k) {
                c_float r = -b[k];
                for (int j = 0; j < n; ++j) r += E[k][j] * x[j];
                assert(std::fabs(r) < 1e-6);
            }
        }
        else assert(res.exitflag == DAQP_EXIT_INFEASIBLE);
    }
}

} // namespace

int main() {
    // x = (x0, x0) determined (0.999 overflowed the BnB tree)
    for (c_float x0 : {0.65, -0.3, 0.999, 1.0, -1.0})
        check(2, {{1.0, 1.0}, {0.3, -0.7}}, {2 * x0, -0.4 * x0}, x0);

    // Only x_0 determined
    for (c_float x0 : {0.65, 1.0, -1.0})
        check(3, {{2.0, 0.0, 0.0}}, {2 * x0}, x0);

    // x_0 in the span of the equalities
    for (c_float x0 : {0.65, 1.0})
        check(3, {{1.0, 1.0, 0.0}, {1.0, -1.0, 0.0}}, {x0 + 0.5, x0 - 0.5}, x0);
}
