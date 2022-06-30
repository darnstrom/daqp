#ifndef DAQP_H
# define DAQP_H

#include "factorization.h" 
#include "constants.h"
#include "auxiliary.h"

int daqp_ldp(DAQPWorkspace *work);
void ldp2qp_solution(DAQPWorkspace *work);

void warmstart_workspace(DAQPWorkspace* work, int* WS, const int n_active); 
void reset_daqp_workspace(DAQPWorkspace *work);

#endif //ifndef DAQP_H
