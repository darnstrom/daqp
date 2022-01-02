#ifndef DAQP_H
# define DAQP_H

#include "factorization.h" 
#include "constants.h"
#include "auxiliary.h"

int daqp(Workspace *work);
void ldp2qp_solution(Workspace *work);

void warmstart_workspace(Workspace* work, int* WS, const int n_active); 
void allocate_daqp_workspace(Workspace *work, int n);
void free_daqp_workspace(Workspace *work);
void reset_daqp_workspace(Workspace *work);
void reset_daqp_workspace_warm(Workspace *work);

void allocate_daqp_ldp(Workspace *work, int n, int m);
void free_daqp_ldp(Workspace *work);
#endif //ifndef DAQP_H
