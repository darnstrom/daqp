#ifndef DAQP_FACTORIZATION_H
# define DAQP_FACTORIZATION_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"
#include "constants.h"

/* Dot product declared static inline so the compiler can inline it into every
 * call site and apply its own vectorisation.  No manual unrolling is needed:
 * with -O3 -fassociative-math the compiler auto-vectorises the simple loop
 * just as well as any hand-unrolled version.
 */
static inline c_float daqp_dot(const c_float* v1, const c_float* v2, const int n) {
    c_float sum = 0.0;
    int i;
    for (i = 0; i < n; i++) sum += v1[i] * v2[i];
    return sum;
}

void daqp_update_LDL_add(DAQPWorkspace *work, const int add_ind);
void daqp_update_LDL_remove(DAQPWorkspace *work, const int rm_ind);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_FACTORIZATION_H
