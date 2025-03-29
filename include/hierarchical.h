#ifndef DAQP_HIERARCHICAL_H 
# define DAQP_HIERARCHICAL_H

#include "types.h"
#include "constants.h"
#include "daqp.h"

int daqp_hiqp(DAQPWorkspace *work);

typedef struct{
    int level;
    c_float* A;
    c_float* bu;
    c_float* bl;
    int m;
}DAQPTask;

void daqp_sort_tasks(DAQPTask* tasks, int n_tasks);
DAQPProblem daqp_setup_hqp(DAQPTask* tasks, int n_tasks, int n);


#endif //ifndef DAQP_HIERARCHICAL_H
