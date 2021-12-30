#ifndef DAQP_PROX_H 
# define DAQP_PROX_H 

#include "types.h"
#include "constants.h"
#include "daqp.h"


typedef struct{
  c_float* f;
  c_float* b;
  
  c_float* xold;
  c_float* x;

  int inner_iterations;
  int outer_iterations;

  int m;
  int n;
  int cycle_counter;
  c_float epsilon;
  c_float fval;

  Workspace *work;
}ProxWorkspace;

int daqp_prox(ProxWorkspace *prox_work);
int gradient_step(ProxWorkspace *prox_work, Workspace* work);

void allocate_prox_workspace(ProxWorkspace *prox_work, int n, int m);
void free_prox_workspace(ProxWorkspace *prox_work);
void reset_prox_workspace(ProxWorkspace *prox_work);

#define ETA 1e-6
#define PROX_ITER_LIMIT 1000 
#endif //ifndef DAQP_PROX_H 
