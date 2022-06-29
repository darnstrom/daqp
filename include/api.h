#ifndef DAQP_API_H
# define DAQP_API_H
#include "daqp.h"
#include "daqp_prox.h"

typedef struct{
  c_float *x;
  c_float *lam;
  c_float fval;
  c_float soft_slack;  

  int exitflag; 
  int iter;
  c_float solve_time;
  c_float setup_time; 

}DAQPResult;

void daqp_solve(DAQPResult* res, Workspace *work);
void daqp_quadprog(DAQPResult* res, QP* qp,DAQPSettings* settings);

int setup_daqp(QP *qp, Workspace* work);
int setup_daqp_ldp(Workspace *work, QP* qp);
void allocate_daqp_settings(Workspace *work);
void allocate_daqp_workspace(Workspace *work, int n);

void free_daqp_ldp(Workspace *work);
void free_daqp_workspace(Workspace *work);

void daqp_extract_result(DAQPResult* res, Workspace* work);
void daqp_default_settings(DAQPSettings *settings);
#endif //ifndef DAQP_API_H
