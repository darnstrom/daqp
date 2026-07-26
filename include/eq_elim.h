#ifndef DAQP_EQ_ELIM_H
# define DAQP_EQ_ELIM_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"

// Whether the workspace currently holds a reduced LDP
#define DAQP_IS_REDUCED(work) ((work)->eq != NULL && (work)->eq->installed)

/*
 * Whether the constraints of the full problem have not been formed. A reduction
 * only needs the equality rows of M, so the full M and d are skipped whenever
 * one is formed, and restoring only swaps the pointers back to them. Every path
 * that forms them again also clears neq, so it marks both.
 */
#define DAQP_EQ_IS_FULL_UNFORMED(work) \
    ((work)->eq != NULL && (work)->eq->neq != 0)

/*
 * Eliminate the equality constraints of the LDP in the workspace.
 * daqp_update_ldp applies the elimination when DAQP_UPDATE_eliminate is set;
 * daqp_extract_result retrieves the full solution after solving.
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

/*
 * Form the constraints of the full problem and put it in the workspace. An
 * elimination only forms the equality rows of M, so this is what makes the
 * workspace usable again once the reduction is dropped.
 */
int daqp_eq_form_full(DAQPWorkspace* work);

/*
 * Put back a reduction that a previous solve retrieved, so that solving the
 * same workspace twice is the same as solving it once. Returns 0 if there is
 * nothing to put back, and a negative exit flag if the equalities turn out to
 * be infeasible.
 */
int daqp_eq_reinstall(DAQPWorkspace* work);

// Expand the reduced iterate w in work->u into the full LDP iterate u,
// and extract the multipliers of the eliminated equality constraints.
void daqp_eq_expand(DAQPWorkspace* work);

void free_daqp_eq(DAQPWorkspace* work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_EQ_ELIM_H
