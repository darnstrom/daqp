#ifndef DAQP_UTILS_H
# define DAQP_UTILS_H
#include <time.h>
#include "daqp.h"
// Utils for transforming QP to LDP
int compute_Rinv_and_M(c_float *R, c_float *M, const c_float eps,const int n, const int m);
void update_v_and_d(c_float *f, c_float *bupper, c_float *blower, Workspace *work) ;
void pack_symmetric(c_float *S, c_float *Sp, int n);
//
// Utils for profiling
double time_diff(struct timespec tic, struct timespec toc);

#endif //ifndef DAQP_UTILS_H
