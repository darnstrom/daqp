#ifndef DAQP_API_H
# define DAQP_API_H
#include "daqp.h"
#include "daqp_prox.h"

typedef struct{
  c_float *x;
  c_float fval;
  int exitflag; 
  int iter;
  int outer_iter;
  double solve_time;
  double setup_time; 
}DAQPResult;

void daqp_setup(Workspace** work, double* M, double* d,int n);

int daqp_feas(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound);
int daqp_feas_warmstart(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound, int* WS, const int n_active);


void daqp_quadprog(DAQPResult* res, double* H, double* f, double *A, double *bupper, double *blower, int* sense, int n, int m, int packed,DAQPSettings* settings);

int qp2ldp(double *R, double *v, double* M, double* dupper, double* dlower, int n, int m, double eps);

#endif //ifndef DAQP_API_H
