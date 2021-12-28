#ifndef DAQP_API_H
# define DAQP_API_H
#include "daqp.h"
#include "daqp_prox.h"
void daqp_setup(Workspace** work, double* M, double* d,int n);

int daqp_feas(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound);
int daqp_feas_warmstart(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound, int* WS, const int n_active);


int daqp_quadprog(double* x, double* H, double* f, double *A, double *b, int* sense, int n, int m, int packed);

void pack_H(double *H, double *R, int n);
int qp2ldp(double *R, double *v, double* M, double* d, int n, int m, double eps);
#endif //ifndef DAQP_API_H
