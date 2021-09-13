#ifndef DAQP_API_H
# define DAQP_API_H
#include "daqp.h"
#include "daqp_prox.h"
void daqp_setup(Workspace** work, double* M, double* d,int n);

int daqp_feas(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound);
int daqp_feas_warmstart(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound, int* WS, const int n_active);

int minimal_Hrepresentation(Workspace *work);
int daqp_minHRep(Workspace* work, c_float* A, c_float*b, int *sense, const int m);

void daqp_lp_setup(ProxWorkspace** pws_ptr, double* f, double* A,double* b, int n, int m);
c_float daqp_lp_solve(ProxWorkspace* prox_work, c_float* f, c_float* A, c_float* b, const int m);
#endif //ifndef DAQP_API_H
