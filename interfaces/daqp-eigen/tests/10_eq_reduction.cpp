#include "api.h"
#include "utils.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
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

bool setup_reduces(TestProblem& p, int policy, DAQPWorkspace& work) {
    work = {};
    allocate_daqp_settings(&work);
    work.settings->eq_reduction = policy;
    assert(setup_daqp(&p.qp, &work, nullptr) > 0);
    return DAQP_IS_REDUCED(&work);
}

bool setup_reduces(SizedProblem& p, int policy, DAQPWorkspace& work) {
    work = {};
    allocate_daqp_settings(&work);
    work.settings->eq_reduction = policy;
    assert(setup_daqp(&p.qp, &work, nullptr) > 0);
    return DAQP_IS_REDUCED(&work);
}

void cleanup(DAQPWorkspace& work) {
    free_daqp_workspace(&work);
    free_daqp_ldp(&work);
}

} // namespace

int main() {
    DAQPWorkspace work{};

    TestProblem dense(false);
    assert(setup_reduces(dense, DAQP_EQ_REDUCTION_AUTO, work));
    assert(work.eq->rebuilds == 0);

    dense.bu.back() = 9.0;
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

    SizedProblem rebuild_candidate(60, 12, 20, false);
    assert(setup_reduces(rebuild_candidate, DAQP_EQ_REDUCTION_AUTO, work));
    const int structural = DAQP_UPDATE_Rinv | DAQP_UPDATE_M |
                           DAQP_UPDATE_v | DAQP_UPDATE_d |
                           DAQP_UPDATE_sense;
    for (int i = 1; i <= DAQP_EQ_MAX_REBUILDS; ++i) {
        assert(daqp_update_ldp(structural, &work,
                               &rebuild_candidate.qp) >= 0);
        assert(DAQP_IS_REDUCED(&work));
        assert(work.eq->rebuilds == i);
    }
    assert(daqp_update_ldp(structural, &work, &rebuild_candidate.qp) >= 0);
    assert(!DAQP_IS_REDUCED(&work));
    assert(work.eq->rebuilds == DAQP_EQ_MAX_REBUILDS);
    assert(daqp_update_ldp(DAQP_UPDATE_d, &work,
                           &rebuild_candidate.qp) >= 0);
    assert(DAQP_IS_REDUCED(&work));
    assert(work.eq->rebuilds == 1);
    assert(daqp_update_ldp(DAQP_UPDATE_d, &work,
                           &rebuild_candidate.qp) >= 0);
    assert(work.eq->rebuilds == 0);
    cleanup(work);

    // Many equalities: a dense Hessian keeps the reduction through rebuilds,
    // while a diagonal one makes the full constraints cheap enough to give up
    for (int diag = 0; diag < 2; ++diag) {
        SizedProblem dense_eq(60, 30, 20, diag == 1);
        assert(setup_reduces(dense_eq, DAQP_EQ_REDUCTION_AUTO, work));
        for (int i = 1; i <= DAQP_EQ_MAX_REBUILDS; ++i) {
            assert(daqp_update_ldp(structural, &work, &dense_eq.qp) >= 0);
            assert(DAQP_IS_REDUCED(&work));
        }
        assert(daqp_update_ldp(structural, &work, &dense_eq.qp) >= 0);
        assert(DAQP_IS_REDUCED(&work) == (diag == 0));
        cleanup(work);
    }

    TestProblem singular(true, true);
    setup_reduces(singular, DAQP_EQ_REDUCTION_AUTO, work);
    assert(work.n_prox > 0);
    assert(work.eq != nullptr && work.eq->neq > 0);
    assert(daqp_eq_will_reduce(&work));
    cleanup(work);
}
