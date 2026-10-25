#ifndef DAQP_PRESOLVE_H
# define DAQP_PRESOLVE_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"
#include "constants.h"

DAQPProblem* daqp_presolve(DAQPWorkspace* work, DAQPProblem* qp);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_PRESOLVE_H
