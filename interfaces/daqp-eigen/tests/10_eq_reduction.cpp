#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <vector>

namespace {

struct TestProblem {
    static constexpr int n = 40;
    static constexpr int neq = 16;
    static constexpr int nineq = 20;
    static constexpr int m = neq + nineq;

    std::vector<c_float> H = std::vector<c_float>(n * n, 0.0);
    std::vector<c_float> f = std::vector<c_float>(n, 0.0);
    std::vector<c_float> A = std::vector<c_float>(m * n, 0.0);
    std::vector<c_float> bu = std::vector<c_float>(m, 10.0);
    std::vector<c_float> bl = std::vector<c_float>(m, -10.0);
    std::vector<int> sense = std::vector<int>(m, 0);
    DAQPProblem qp{};

    explicit TestProblem(bool diagonal, bool singular = false,
                         bool inequalities = true) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j)
                H[i*n+j] = i == j ? 2.0 : (diagonal ? 0.0 : 0.01);
            A[(i % m)*n+i] = 1.0;
        }
        if (singular) H[(n-1)*n+n-1] = 0.0;
        for (int i = 0; i < neq; ++i) {
            A[i*n+i] = 1.0;
            bu[i] = bl[i] = 1.0;
            sense[i] = DAQP_ACTIVE | DAQP_IMMUTABLE;
        }
        for (int i = 0; i < nineq; ++i)
            A[(neq+i)*n+(i+neq)%n] = 1.0;
        qp = {n, inequalities ? m : neq, 0, H.data(), f.data(), A.data(),
              bu.data(), bl.data(), sense.data(), nullptr, 0, 0};
    }
};

struct SizedProblem {
    int n, neq, nineq, m;
    std::vector<c_float> H, f, A, bu, bl;
    std::vector<int> sense;
    DAQPProblem qp{};

    SizedProblem(int n_, int neq_, int nineq_, bool diagonal)
        : n(n_), neq(neq_), nineq(nineq_), m(neq_ + nineq_),
          H(n * n, 0.0), f(n, 0.0), A(m * n, 0.0),
          bu(m, 10.0), bl(m, -10.0), sense(m, 0) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                H[i*n+j] = i == j ? 2.0 : (diagonal ? 0.0 : 0.01);
        for (int i = 0; i < neq; ++i) {
            A[i*n+i] = 1.0;
            bu[i] = bl[i] = 1.0;
            sense[i] = DAQP_ACTIVE | DAQP_IMMUTABLE;
        }
        for (int i = 0; i < nineq; ++i)
            A[(neq+i)*n+(i+neq)%n] = 1.0;
        qp = {n, m, 0, H.data(), f.data(), A.data(), bu.data(), bl.data(),
              sense.data(), nullptr, 0, 0};
    }
};

// Set up a workspace, as for a one-shot solve (such as daqp_quadprog) or as a
// workspace that is to be updated and solved repeatedly
template <class Problem>
bool setup_reduces(Problem& p, int policy, DAQPWorkspace& work,
                   bool one_shot = true) {
    work = {};
    allocate_daqp_settings(&work);
    work.settings->eq_reduction = policy;
    assert(setup_daqp_main(&p.qp, &work, nullptr,
                           one_shot ? DAQP_UPDATE_eliminate : 0) > 0);
    return DAQP_IS_REDUCED(&work);
}

void cleanup(DAQPWorkspace& work) {
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
}

} // namespace

int main() {
    DAQPWorkspace work{};

    // AUTO reduces a one-shot problem, but not a workspace that is to be
    // updated; ON reduces both
    TestProblem dense(false);
    assert(setup_reduces(dense, DAQP_EQ_REDUCTION_AUTO, work));
    cleanup(work);
    assert(!setup_reduces(dense, DAQP_EQ_REDUCTION_AUTO, work, false));
    cleanup(work);
    assert(setup_reduces(dense, DAQP_EQ_REDUCTION_ON, work, false));
    cleanup(work);

    // Updates of a workspace are only reduced with ON (or when an update
    // itself is marked as one-shot), and OFF gives up a reduction in place
    assert(setup_reduces(dense, DAQP_EQ_REDUCTION_AUTO, work));
    dense.bu.back() = 9.0;
    assert(daqp_update_ldp(DAQP_UPDATE_d, &work, &dense.qp) >= 0);
    assert(!DAQP_IS_REDUCED(&work));
    assert(daqp_update_ldp(DAQP_UPDATE_d | DAQP_UPDATE_eliminate, &work,
                           &dense.qp) >= 0);
    assert(DAQP_IS_REDUCED(&work));
    cleanup(work);

    assert(setup_reduces(dense, DAQP_EQ_REDUCTION_ON, work));
    assert(daqp_update_ldp(DAQP_UPDATE_d, &work, &dense.qp) >= 0);
    assert(DAQP_IS_REDUCED(&work));

    work.settings->eq_reduction = DAQP_EQ_REDUCTION_OFF;
    assert(daqp_update_ldp(DAQP_UPDATE_d, &work, &dense.qp) >= 0);
    assert(!DAQP_IS_REDUCED(&work));
    cleanup(work);

    TestProblem diagonal_with_inequalities(true);
    assert(setup_reduces(diagonal_with_inequalities,
                         DAQP_EQ_REDUCTION_AUTO, work));
    cleanup(work);

    SizedProblem tiny(10, 6, 10, false);
    assert(!setup_reduces(tiny, DAQP_EQ_REDUCTION_AUTO, work));
    cleanup(work);
    assert(setup_reduces(tiny, DAQP_EQ_REDUCTION_ON, work));
    cleanup(work);

    SizedProblem ratio_boundary(60, 6, 20, false);
    assert(!setup_reduces(ratio_boundary, DAQP_EQ_REDUCTION_AUTO, work));
    cleanup(work);
    assert(setup_reduces(ratio_boundary, DAQP_EQ_REDUCTION_ON, work));
    cleanup(work);

    SizedProblem diagonal_sparse(40, 6, 20, true);
    assert(!setup_reduces(diagonal_sparse, DAQP_EQ_REDUCTION_AUTO, work));
    cleanup(work);
    assert(setup_reduces(diagonal_sparse, DAQP_EQ_REDUCTION_ON, work));
    cleanup(work);

    TestProblem diagonal(true, false, false);
    assert(!setup_reduces(diagonal, DAQP_EQ_REDUCTION_AUTO, work));
    cleanup(work);
    assert(setup_reduces(diagonal, DAQP_EQ_REDUCTION_ON, work));
    work.settings->eq_reduction = DAQP_EQ_REDUCTION_OFF;
    std::vector<c_float> x(TestProblem::n), lam(TestProblem::m);
    DAQPResult result{};
    result.x = x.data();
    result.lam = lam.data();
    daqp_solve(&result, &work);
    assert(result.exitflag > 0);
    assert(!DAQP_IS_REDUCED(&work));
    cleanup(work);
    assert(!setup_reduces(dense, DAQP_EQ_REDUCTION_OFF, work));
    cleanup(work);

    // Solving a workspace repeatedly: whether the data or only the bounds are
    // updated, AUTO solves the full problem after the setup (also if the setup
    // was reduced) and ON reduces every update, with the same solutions
    const int structural = DAQP_UPDATE_Rinv | DAQP_UPDATE_M |
                           DAQP_UPDATE_v | DAQP_UPDATE_d |
                           DAQP_UPDATE_sense;
    for (int policy : {DAQP_EQ_REDUCTION_AUTO, DAQP_EQ_REDUCTION_ON}) {
        for (int mask : {structural, int(DAQP_UPDATE_d)}) {
            SizedProblem p(60, 30, 20, false);
            assert(setup_reduces(p, policy, work));
            std::vector<c_float> xs(p.n), lams(p.m);
            DAQPResult res{};
            res.x = xs.data();
            res.lam = lams.data();
            daqp_solve(&res, &work);
            assert(res.exitflag > 0);
            for (int i = 0; i < 3; ++i) {
                p.bu[p.m-1] = 9.0-i;
                assert(daqp_update_ldp(mask, &work, &p.qp) >= 0);
                assert(DAQP_IS_REDUCED(&work) ==
                       (policy == DAQP_EQ_REDUCTION_ON));
                daqp_solve(&res, &work);
                assert(res.exitflag > 0);
                for (int j = 0; j < p.neq; ++j) // The equalities hold
                    assert(std::abs(xs[j] - 1.0) < 1e-9);
            }
            cleanup(work);
        }
    }

    TestProblem singular(true, true);
    setup_reduces(singular, DAQP_EQ_REDUCTION_AUTO, work);
    assert(work.n_prox > 0);
    assert(work.eq != nullptr && work.eq->neq > 0);
    assert(daqp_eq_will_reduce(&work));
    cleanup(work);
}
