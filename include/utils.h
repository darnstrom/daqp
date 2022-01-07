#ifndef DAQP_UTILS_H
# define DAQP_UTILS_H
#include <time.h>
#include "daqp.h"
// Utils for transforming QP to LDP
int update_ldp(const int mask, Workspace *work);
int update_Rinv(Workspace *work);
void update_M(Workspace *work);
void update_v(c_float *f, Workspace *work);
void update_d(Workspace *work);
//
// Utils for profiling
double time_diff(struct timespec tic, struct timespec toc);

#endif //ifndef DAQP_UTILS_H
