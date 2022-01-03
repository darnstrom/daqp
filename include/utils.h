#ifndef DAQP_UTILS_H
# define DAQP_UTILS_H
#include <time.h>
#include "daqp.h"
// Utils for transforming QP to LDP
int compute_Rinv_and_M(double *R, double *M, const double eps,const int n, const int m);
void update_v_and_d(double *f, double *bupper, double *blower, Workspace *work) ;
void pack_symmetric(double *S, double *Sp, int n);
//
// Utils for profiling
double time_diff(struct timespec tic, struct timespec toc);

#endif //ifndef DAQP_UTILS_H
