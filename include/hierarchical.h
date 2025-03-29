#ifndef DAQP_HIERARCHICAL_H 
# define DAQP_HIERARCHICAL_H

#include "types.h"
#include "constants.h"
#include "daqp.h"

int daqp_hiqp(DAQPWorkspace *work);

typedef struct{
    int index;
    c_float* A;
    c_float* bu;
    c_float* bl;
    int m;
}DAQPTask;

void daqp_sort_tasks(DAQPTask* tasks, int nt);
DAQPProblem daqp_setup_hqp(DAQPTask* tasks, int nt, int n);


#endif //ifndef DAQP_HIERARCHICAL_H
