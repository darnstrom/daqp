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

// Warm starts keep the LDL factorization between solves, so any error in it
// accumulates. A degenerate problem (constraints that are exact combinations
// of others) makes singular working sets frequent; the pivot of a singular
// constraint must then be kept as computed rather than zeroed, since zeroing
// makes L*D*L' deviate from M*M' and the deviation compounds over solves.

namespace {

constexpr int n = 30;
constexpr int mi = 60; // Independent general constraints
constexpr int md = 40; // Combinations of the independent ones
constexpr int mA = mi + md;
constexpr int m = n + mA;
constexpr int n_solves = 300;

} // namespace

int main() {
    std::mt19937 rng(1);
    std::normal_distribution<c_float> N(0, 1);
    std::uniform_real_distribution<c_float> U(0, 1);

    std::vector<c_float> H(n * n, 0.0), f(n), f0(n), A(mA * n), C(md * mi, 0.0);
    std::vector<c_float> bu(m), bl(m), bu0(m);
    std::vector<int> sense(m, 0), sense_ref(m, 0);

    for (int i = 0; i < n; ++i) H[i*n+i] = 2 + U(rng);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j) {
            const c_float v = 0.3 * N(rng) / n;
            H[i*n+j] += v;
            H[j*n+i] += v;
        }
    for (int i = 0; i < mi * n; ++i) A[i] = N(rng);
    for (int d = 0; d < md; ++d) {
        for (int t = 0; t < 3; ++t) C[d*mi + static_cast<int>(U(rng) * mi)] += N(rng);
        for (int j = 0; j < n; ++j) {
            c_float s = 0;
            for (int k = 0; k < mi; ++k) s += C[d*mi+k] * A[k*n+j];
            A[(mi+d)*n+j] = s;
        }
    }
    for (int i = 0; i < n; ++i) { f0[i] = 20 * N(rng); bu0[i] = 1; }
    for (int i = 0; i < mi; ++i) bu0[n+i] = 0.2 * U(rng);
    for (int d = 0; d < md; ++d) { // Tight, so that the combinations become active
        c_float s = 0;
        for (int k = 0; k < mi; ++k) s += std::fabs(C[d*mi+k]) * bu0[n+k];
        bu0[n+mi+d] = 0.3 * s;
    }

    auto perturb = [&]() {
        for (int i = 0; i < n; ++i) f[i] = f0[i] * (1 + 0.3 * N(rng));
        for (int i = 0; i < m; ++i) {
            bu[i] = i < n ? bu0[i] : bu0[i] * (1 + 0.05 * N(rng));
            bl[i] = i < n ? -bu0[i] : -DAQP_INF;
        }
    };
    perturb();

    DAQPProblem qp = {n, m, n, H.data(), f.data(), A.data(), bu.data(), bl.data(),
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
        perturb();
        assert(daqp_update_ldp(DAQP_UPDATE_v | DAQP_UPDATE_d, &work, &qp) >= 0);
        daqp_solve(&res, &work);
        assert(res.exitflag == DAQP_EXIT_OPTIMAL);

        // Reference: cold solve of the same data
        DAQPProblem qp_ref = qp;
        std::fill(sense_ref.begin(), sense_ref.end(), 0);
        qp_ref.sense = sense_ref.data();
        daqp_quadprog(&ref, &qp_ref, nullptr);
        assert(ref.exitflag == DAQP_EXIT_OPTIMAL);
        for (int i = 0; i < n; ++i) assert(std::fabs(x[i] - xref[i]) < 1e-6);
    }
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
    return 0;
}
