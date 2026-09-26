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

// Warm-started LPs with consistent, linearly dependent constraints through the
// optimum. The proximal LP solver moves to a vertex by activating the nearest
// constraint, which can be dependent on the active ones without being
// violated; such a constraint must not be taken as a certificate of
// infeasibility, neither in that solve nor in the next warm-started one.

namespace {

constexpr int n = 20;
constexpr int mi = 30; // Independent general constraints
constexpr int md = 20; // Combinations of the independent ones
constexpr int mA = mi + md;
constexpr int m = n + mA;
constexpr int n_problems = 10;
constexpr int n_solves = 20;

c_float row_dot(const std::vector<c_float>& A, int r, const std::vector<c_float>& x) {
    c_float s = 0;
    for (int j = 0; j < n; ++j) s += A[r*n+j] * x[j];
    return s;
}

} // namespace

int main() {
    std::mt19937 rng(3);
    std::normal_distribution<c_float> N(0, 1);
    std::uniform_real_distribution<c_float> U(0, 1);

    for (int p = 0; p < n_problems; ++p) {
        std::vector<c_float> f(n), f0(n), A(mA * n), bu(m), bl(m), bu0(m), x0(n);
        std::vector<int> sense(m, 0), sense_ref(m, 0);
        for (c_float& v : x0) v = 0.5 * (2 * U(rng) - 1);
        for (int i = 0; i < mi * n; ++i) A[i] = N(rng);
        for (int d = 0; d < md; ++d) {
            const int a = static_cast<int>(U(rng) * mi), b = static_cast<int>(U(rng) * mi);
            const c_float ca = N(rng), cb = N(rng);
            for (int j = 0; j < n; ++j) A[(mi+d)*n+j] = ca * A[a*n+j] + cb * A[b*n+j];
        }
        for (int i = 0; i < n; ++i) { bu0[i] = 2; bl[i] = -2; }
        // Many constraints are tight at x0 (degenerate), the rest are slack
        for (int r = 0; r < mA; ++r) {
            c_float scale = 0;
            for (int j = 0; j < n; ++j) scale += std::fabs(A[r*n+j]);
            bu0[n+r] = row_dot(A, r, x0) + (U(rng) < 0.6 ? 0 : 0.1 * scale * U(rng));
            bl[n+r] = -DAQP_INF;
        }
        // x0 is optimal: -f is a nonnegative combination of tight rows
        std::fill(f0.begin(), f0.end(), 0.0);
        for (int r = 0, cnt = 0; r < mA && cnt < n/2; ++r) {
            if (bu0[n+r] != row_dot(A, r, x0) || U(rng) < 0.5) continue;
            const c_float l = U(rng);
            for (int j = 0; j < n; ++j) f0[j] -= l * A[r*n+j];
            ++cnt;
        }
        f = f0;
        bu = bu0;

        DAQPProblem qp = {n, m, n, nullptr, f.data(), A.data(), bu.data(), bl.data(),
                          sense.data(), nullptr, 0, 0};
        DAQPWorkspace work{};
        assert(setup_daqp(&qp, &work, nullptr) > 0);

        std::vector<c_float> x(n), lam(m), xref(n), lamref(m);
        DAQPResult res{};
        res.x = x.data();
        res.lam = lam.data();
        DAQPResult ref{};
        ref.x = xref.data();
        ref.lam = lamref.data();

        for (int k = 0; k < n_solves; ++k) {
            // Perturb the cost and loosen some bounds (the LP stays feasible)
            for (int i = 0; i < n; ++i) f[i] = f0[i] + 0.1 * (1 + std::fabs(f0[i])) * N(rng);
            for (int i = n; i < m; ++i) bu[i] = bu0[i] + (U(rng) < 0.5 ? 0 : 0.01 * U(rng));
            assert(daqp_update_ldp(DAQP_UPDATE_v | DAQP_UPDATE_d, &work, &qp) >= 0);
            daqp_solve(&res, &work);
            assert(res.exitflag == DAQP_EXIT_OPTIMAL);

            DAQPProblem qp_ref = qp;
            std::fill(sense_ref.begin(), sense_ref.end(), 0);
            qp_ref.sense = sense_ref.data();
            daqp_quadprog(&ref, &qp_ref, nullptr);
            assert(ref.exitflag == DAQP_EXIT_OPTIMAL);
            // LP solutions need not be unique, but the optimal value is
            assert(std::fabs(res.fval - ref.fval) < 1e-6 * (1 + std::fabs(ref.fval)));
        }
        free_daqp_workspace(&work);
        free_daqp_ldp(&work);
    }
    return 0;
}
