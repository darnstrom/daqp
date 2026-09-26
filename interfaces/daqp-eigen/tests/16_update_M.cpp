#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <random>
#include <vector>

// daqp_update_M forms M = A*Rinv four rows and two columns at a time. Compare
// it with a plain row-by-row product for dimensions that leave rows and
// columns over, both for plain A and for A that is unscaled in place in M
// (simple bounds with a normalized Rinv).

namespace {

// M <-- A*Rinv, one row at a time (the formulation before the tiling)
void reference_M(const DAQPWorkspace& work, const c_float* A, c_float* M) {
    const int n = work.n;
    const int mA = work.m - work.ms;
    const int stop_id =
        (work.state & DAQP_STATE_RINV_NORMALIZED) ? n - work.ms : n;
    for (int k = 0, disp2 = n * mA - 1; k < mA; k++, disp2 -= n) {
        int disp = DAQP_ARSUM(n);
        int j = 0;
        for (; j < stop_id; ++j) {
            for (int i = 0; i < j; ++i)
                M[disp2 - i] += work.Rinv[--disp] * A[disp2 - j];
            M[disp2 - j] = work.Rinv[--disp] * A[disp2 - j];
        }
        for (; j < n; ++j) {
            const c_float a = A[disp2 - j] / work.scaling[n - j - 1];
            for (int i = 0; i < j; ++i) M[disp2 - i] += work.Rinv[--disp] * a;
            M[disp2 - j] = work.Rinv[--disp] * a;
        }
    }
}

void check(int n, int mA, bool bounds, std::mt19937& rng) {
    std::normal_distribution<c_float> N(0, 1);
    const int ms = bounds ? n : 0;
    const int m = ms + mA;
    std::vector<c_float> L(n * n), H(n * n, 0.0), f(n), A(mA * n);
    std::vector<c_float> bu(m, 1e3), bl(m, -1e3);
    for (auto& v : L) v = N(rng);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) H[i * n + j] += L[i * n + k] * L[j * n + k];
            if (i == j) H[i * n + j] += n;
        }
    for (auto& v : f) v = N(rng);
    for (auto& v : A) v = N(rng);
    DAQPProblem qp = {n, m, ms, H.data(), f.data(), A.data(), bu.data(),
                      bl.data(), nullptr, nullptr, 0, 0};

    DAQPWorkspace work{};
    assert(setup_daqp(&qp, &work, nullptr) >= 0);
    if (work.Rinv == nullptr) { // A diagonal Hessian (n = 1) forms no product
        free_daqp_workspace(&work);
        free_daqp_ldp(&work);
        return;
    }

    // Reference, normalized as the workspace does
    std::vector<c_float> M_ref(mA * n);
    c_float* M = work.M;
    work.M = M_ref.data();
    reference_M(work, A.data(), M_ref.data());
    assert(daqp_normalize_M(&work) == 0);
    work.M = M;

    // A with the same dimensions (in place in M if Rinv is normalized)
    assert(daqp_update_M(&work, A.data()) == 0);
    for (int i = 0; i < mA * n; i++)
        assert(std::abs(work.M[i] - M_ref[i]) <= 1e-12 * (1 + std::abs(M_ref[i])));

    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
}

} // namespace

int main() {
    std::mt19937 rng(16);
    for (int n : {1, 2, 3, 4, 5, 7, 8, 13})
        for (int mA : {1, 2, 3, 4, 5, 7, 9})
            for (bool bounds : {false, true}) check(n, mA, bounds, rng);
}
