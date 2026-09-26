#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <vector>

// Warm-started degenerate LPs where the nearest blocking constraint of the
// gradient step (which moves the proximal LP iterate toward a vertex) is
// linearly dependent on the active constraints. Such a constraint has to be
// skipped, since it does not bring the iterate closer to a vertex; otherwise
// the solver repeats the same step until the iteration limit.
// Constraints through a common point x0 make the problems degenerate.

namespace {

unsigned long long rng;
double urand() {
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
    return (rng >> 11) * (1.0 / 9007199254740992.0);
}
double nrand() {
    const double u1 = urand() + 1e-300, u2 = urand();
    return std::sqrt(-2 * std::log(u1)) * std::cos(2 * M_PI * u2);
}

struct Problem {
    int n, m, ms;
    std::vector<c_float> f, A, bu, bl;
    std::vector<int> sense;

    Problem(int n_, unsigned long long seed) : n(n_), ms(n_) {
        rng = seed;
        const int mi = n + static_cast<int>(urand() * n);
        const int md = n / 2 + static_cast<int>(urand() * n);
        const int mA = mi + md;
        m = n + mA;
        f.resize(n); A.resize(static_cast<size_t>(mA) * n);
        bu.resize(m); bl.resize(m); sense.assign(m, 0);
        std::vector<double> x0(n);
        for (int i = 0; i < n; ++i) x0[i] = 0.5 * (2 * urand() - 1);
        for (int i = 0; i < n; ++i) { bu[i] = 1; bl[i] = -1; }
        for (int i = 0; i < mi * n; ++i) A[i] = nrand();
        for (int d = 0; d < md; ++d) { // Combinations of two other rows
            const int a = static_cast<int>(urand() * mi), b = static_cast<int>(urand() * mi);
            const double ca = nrand(), cb = nrand();
            for (int j = 0; j < n; ++j) A[(mi+d)*n+j] = ca * A[a*n+j] + cb * A[b*n+j];
        }
        for (int r = 0; r < mA; ++r) { // Many constraints are tight at x0
            double scale = 0;
            for (int j = 0; j < n; ++j) scale += std::fabs(A[r*n+j]);
            const int tight = urand() < 0.6;
            bu[n+r] = row(r, x0) + (tight ? 0 : 0.1 * scale * urand());
            bl[n+r] = -1e30;
        }
        // x0 is optimal: -f is a nonnegative combination of tight rows
        for (int i = 0; i < n; ++i) f[i] = 0;
        for (int r = 0, cnt = 0; r < mA && cnt < n / 2; ++r) {
            if (bu[n+r] != row(r, x0) || urand() < 0.5) continue;
            const double l = urand();
            ++cnt;
            for (int j = 0; j < n; ++j) f[j] -= l * A[r*n+j];
        }
        for (int i = 0; i < n; ++i) { bu[i] = 2; bl[i] = -2; }
    }

    double row(int r, const std::vector<double>& x) const {
        double s = 0;
        for (int j = 0; j < n; ++j) s += A[r*n+j] * x[j];
        return s;
    }

    DAQPProblem qp() {
        return {n, m, ms, nullptr, f.data(), A.data(), bu.data(),
                bl.data(), sense.data(), nullptr, 0, 0};
    }
};

unsigned long long seed_of(int family, int size_ind, int s) {
    return 0x9E3779B97F4A7C15ULL * (1000 * family + 100 * size_ind + s + 1) + 7;
}

} // namespace

int main() {
    struct Case { int n, size_ind, s, steps; };
    for (const Case& c : {Case{30, 1, 29, 15}, Case{80, 2, 28, 4}}) {
        Problem p(c.n, seed_of(6, c.size_ind, c.s));
        DAQPProblem qp = p.qp();
        std::vector<c_float> x(p.n), lam(p.m);
        DAQPResult res{};
        res.x = x.data();
        res.lam = lam.data();
        const std::vector<c_float> f0 = p.f, bu0 = p.bu;
        DAQPWorkspace work{};
        assert(setup_daqp(&qp, &work, nullptr) >= 0);
        for (int k = 0; k < c.steps; ++k) {
            for (int i = 0; i < p.n; ++i) p.f[i] = f0[i] + 0.1 * (1 + std::fabs(f0[i])) * nrand();
            for (int i = p.n; i < p.m; ++i) p.bu[i] = bu0[i] + (urand() < 0.5 ? 0 : 0.01 * urand());
            assert(daqp_update_ldp(DAQP_UPDATE_v | DAQP_UPDATE_d, &work, &qp) >= 0);
            daqp_solve(&res, &work);
            assert(res.exitflag == DAQP_EXIT_OPTIMAL);
        }
        free_daqp_workspace(&work);
        free_daqp_ldp(&work);
    }
    return 0;
}
