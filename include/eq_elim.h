#ifndef DAQP_EQ_ELIM_H
# define DAQP_EQ_ELIM_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"

// Whether the workspace currently holds a reduced LDP
#define DAQP_IS_REDUCED(work) ((work)->eq != NULL && (work)->eq->installed)

/*
 * Eliminate the equality constraints of the LDP in the workspace, and retrieve
 * the solution of the original problem from the reduced one. daqp_quadprog
 * applies these around an ordinary setup/solve; a workspace that is set up and
 * solved directly is reduced only if these are called.
 *
 * Returns the number of eliminated constraints (0 if the LDP was left intact)
 * or a negative exit flag if the equality constraints cannot be satisfied.
 */
int daqp_eq_eliminate(DAQPWorkspace* work);

// Reduce the LDP in the workspace by eliminating equality constraints.
// Returns 1 if a reduced LDP was installed, 0 if the LDP was left intact,
// and a negative exit flag if the equalities are infeasible.
int daqp_eq_reduce(DAQPWorkspace* work, const int mask);

// Whether the equality constraints will be eliminated
int daqp_eq_will_reduce(const DAQPWorkspace* work);

// Put the full LDP back into the workspace (no-op if it is not reduced)
void daqp_eq_restore(DAQPWorkspace* work);

// Expand the reduced iterate w in work->u into the full LDP iterate u,
// and extract the multipliers of the eliminated equality constraints.
void daqp_eq_expand(DAQPWorkspace* work);

void free_daqp_eq(DAQPWorkspace* work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_EQ_ELIM_H
