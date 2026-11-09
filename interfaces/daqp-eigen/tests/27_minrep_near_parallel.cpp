#include "api.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <vector>

// Redundancy removal with nearly parallel constraints: daqp_minrep forces each
// constraint to equality (at its upper bound, the lower bound being -inf) and checks
// feasibility. The near-dependence triggers the refactorization in daqp_ldp, which
// must not move such a constraint to its (infinite) lower bound, since that gave
// false infeasibility and dropped constraints that are not redundant. The instances
// come from Fourier-Motzkin elimination of random polyhedra (PolyDAQP.jl); all three
// constraints in each are nonredundant (verified with HiGHS)

namespace {

void check_all_nonredundant(std::vector<c_float> A, std::vector<c_float> b, int n) {
    const int m = static_cast<int>(b.size());
    std::vector<int> is_redundant(m, -1);
    daqp_minrep(is_redundant.data(), A.data(), b.data(), n, m, 0);
    for (int i = 0; i < m; ++i) assert(is_redundant[i] == 0);
}

} // namespace

int main() {
    check_all_nonredundant({
        -0.1407115988511143, -0.5919173663814182, 0.3722899713576772, -0.16112264864186251, 0.6821097760961792,
        -0.14071425328600495, -0.5919232776583797, 0.37228572094793333, -0.16112627923310846, 0.6821055610439594,
        -0.14071610659701875, -0.5919174243534064, 0.37228550037321256, -0.16112040367474242, 0.6821117663715827},
        {0.029502198692654098, 0.029499387880622053, 0.029499723174929993}, 5);

    check_all_nonredundant({
        0.30319748635969795, 0.8433210005761512, 0.4395850058549491, 0.06038209072153924,
        0.3033251439684689, 0.8432983640228614, 0.4395324368032527, 0.06043974909501073,
        0.30333823563004414, 0.843304472922671, 0.43951214059384713, 0.0604364047911712},
        {0.01567729376025065, 0.015663658447281494, 0.015657969188361355}, 4);
    return 0;
}
