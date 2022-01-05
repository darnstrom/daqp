#ifndef DAQP_UTILS_H
# define DAQP_UTILS_H
#include <time.h>
#include "daqp.h"
// Utils for transforming QP to LDP
int compute_Rinv_and_M(Workspace *work);
void update_v_and_d(c_float *f, Workspace *work) ;
//
// Utils for profiling
double time_diff(struct timespec tic, struct timespec toc);

#endif //ifndef DAQP_UTILS_H
