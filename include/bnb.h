#ifndef DAQP_BNB_H 
# define DAQP_BNB_H 

#include "types.h"
#include "constants.h"
#include "daqp.h"

int daqp_bnb(DAQPWorkspace *work);

typedef struct{
  int* bin_states; 
}DAQPNode;

typedef struct{
  int* bin_ids;
  int nb;

  DAQPNode* tree;
  int n_tree;
  int start_id;

  c_float  Jbar;
  c_float* ubar;

  DAQPWorkspace* dwork;
}DAQPBnBWorkspace;

int daqp_bnb(DAQPBnBWorkspace* work);
int get_branch_id(DAQPNode *node, DAQPBnBWorkspace* work);
int spawn_children(DAQPNode *node, int branch_id, DAQPBnBWorkspace* work);

#endif //ifndef DAQP_BNB_H 
