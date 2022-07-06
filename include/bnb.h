#ifndef DAQP_BNB_H 
# define DAQP_BNB_H 

#include "types.h"
#include "constants.h"
#include "daqp.h"


int daqp_bnb(DAQPWorkspace* work);
int get_branch_id(DAQPWorkspace* work);
void spawn_children(DAQPNode *node, int branch_id, DAQPWorkspace* work);
int process_node(DAQPNode *node, DAQPWorkspace* work);

void setup_daqp_bnb(DAQPWorkspace* work, int* bin_inds, int nb);
void free_daqp_bnb(DAQPWorkspace* work);

//#define UNMASK_LOWER(x) (x>>15)
//#define UNMASK_BIN_ID(x) (x - (1<<16))
//#define MASK_LOWER(x) (x+(1<<16))


#endif //ifndef DAQP_BNB_H 
