#ifndef DAQP_BNB_H 
# define DAQP_BNB_H 

#include "types.h"
#include "constants.h"
#include "daqp.h"

typedef struct{
  int new_bin; 
  int depth;
}DAQPNode;

typedef struct{
  int* bin_ids;
  int nb;

  DAQPNode* tree;
  int n_tree;

  c_float* ustar;

  DAQPWorkspace* dwork;
}DAQPBnBWorkspace;

int daqp_bnb(DAQPBnBWorkspace* work);
int get_branch_id(DAQPBnBWorkspace* work);
void spawn_children(DAQPNode *node, int branch_id, DAQPBnBWorkspace* work);
int process_node(DAQPNode *node, DAQPBnBWorkspace* work);

void setup_daqp_bnb(DAQPWorkspace* work, DAQPBnBWorkspace* bnb_work);
void free_daqp_bnb(DAQPBnBWorkspace* bnb_work);

#define UNMASK_LOWER(x) (x>>30)
#define UNMASK_BIN_ID(x) (x & ~(1 << 31))
#define MASK_LOWER(x) (x+(1<<31))


#endif //ifndef DAQP_BNB_H 
