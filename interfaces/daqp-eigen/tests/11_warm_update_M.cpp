#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

// Updating A (or H) without updating the constraint states warm starts from
// the previous working set. The update clears the factorization, so the
// constraints that are still marked active have to be activated again.

namespace {

constexpr int n = 30;
constexpr int m = 45;

struct Problem {
    std::vector<c_float> H = std::vector<c_float>(n * n, 0.0);
    std::vector<c_float> f = std::vector<c_float>(n, 0.0);
    std::vector<c_float> A = std::vector<c_float>(m * n, 0.0);
    std::vector<c_float> bu = std::vector<c_float>(m, 0.0);
    std::vector<c_float> bl = std::vector<c_float>(m, -DAQP_INF);
    std::vector<int> sense = std::vector<int>(m, 0);
    DAQPProblem qp{};

    Problem() {
        qp = {n, m, 0, H.data(), f.data(), A.data(), bu.data(), bl.data(),
              sense.data(), nullptr, 0, 0};
    }

    void randomize(std::mt19937& rng, bool hessian) {
        std::normal_distribution<c_float> N(0, 1);
        std::uniform_real_distribution<c_float> U(0, 1);
        if (hessian) {
            std::vector<c_float> M(n * n);
            for (c_float& v : M) v = N(rng);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) {
                    c_float s = i == j ? 0.1 : 0.0;
                    for (int k = 0; k < n; ++k) s += M[k*n+i] * M[k*n+j] / n;
                    H[i*n+j] = s;
                }
        }
        for (c_float& v : f) v = 5 * N(rng);
        for (c_float& v : A) v = N(rng);
        for (c_float& v : bu) v = U(rng); // x = 0 is feasible
    }

    c_float violation(const std::vector<c_float>& x) const {
        c_float viol = 0;
        for (int i = 0; i < m; ++i) {
            c_float ax = 0;
            for (int j = 0; j < n; ++j) ax += A[i*n+j] * x[j];
            viol = std::max(viol, ax - bu[i]);
        }
        return viol;
    }
};

void check_warm_updates(int mask, bool hessian) {
    std::mt19937 rng(42);
    Problem p;
    p.randomize(rng, true);

    DAQPWorkspace work{};
    allocate_daqp_settings(&work);
    work.settings->eq_reduction = DAQP_EQ_REDUCTION_OFF;
    assert(setup_daqp(&p.qp, &work, nullptr) > 0);

    std::vector<c_float> x(n), lam(m), xref(n), lamref(m);
    DAQPResult res{};
    res.x = x.data();
    res.lam = lam.data();
    DAQPResult ref{};
    ref.x = xref.data();
    ref.lam = lamref.data();

    for (int step = 0; step < 10; ++step) {
        if (step > 0) {
            p.randomize(rng, hessian);
            assert(daqp_update_ldp(mask, &work, &p.qp) >= 0);
        }
        daqp_solve(&res, &work);
        assert(res.exitflag == DAQP_EXIT_OPTIMAL);
        assert(p.violation(x) < 1e-6);

        daqp_quadprog(&ref, &p.qp, nullptr);
        assert(ref.exitflag == DAQP_EXIT_OPTIMAL);
        for (int i = 0; i < n; ++i) assert(std::fabs(x[i] - xref[i]) < 1e-6);
    }
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
}

} // namespace

int main() {
    check_warm_updates(DAQP_UPDATE_M | DAQP_UPDATE_v | DAQP_UPDATE_d, false);
    check_warm_updates(DAQP_UPDATE_Rinv | DAQP_UPDATE_M | DAQP_UPDATE_v |
                       DAQP_UPDATE_d, true);
    return 0;
}
