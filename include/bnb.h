#ifndef DAQP_BNB_H 
# define DAQP_BNB_H 

#include "types.h"
#include "constants.h"
#include "daqp.h"


int daqp_bnb(DAQPWorkspace* work);
int get_branch_id(DAQPWorkspace* work);
void spawn_children(DAQPNode *node, int branch_id, DAQPWorkspace* work);
int process_node(DAQPNode *node, DAQPWorkspace* work);

int setup_daqp_bnb(DAQPWorkspace* work, int* bin_inds, int nb);
void free_daqp_bnb(DAQPWorkspace* work);

#endif //ifndef DAQP_BNB_H 
