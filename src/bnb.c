#include <stdio.h>
#include <stdlib.h>
#include "bnb.h"
#include "daqp.h"

int daqp_bnb(DAQPBnBWorkspace* work){
  DAQPNode* current_node;
  int branch_ind;
  while(work->n_tree > 0){
	current_node = work->tree[work->start_id]; 
	work->n_tree--;
	work->start_id++; // TODO: do modulo...
	work->start_id = work->start_id > N ? 0 : work->start_id;

	process_node(current_node,work); // Solve relaxation

	// Cut conditions
	(work->dwork->fval >= work->Jbar) && continue; // Dominance cut 

	branch_ind = get_branch_id(current_node, work); 
	if(branch_id==-1){// Nothing to branch over => integer feasible
	  work->Jbar = fval;
	  work->ubar = work->dwork->u; // TODO: do pointer swap
	}
	else{
	  spawn_children(current_node, branch_id, work)
	}
  }
}

int get_branch_id(DAQPNode* node, DAQPBnBWorkspace* work){
  for(int i=0; i<work->nb; i++){
	if(node->bin_states[i]==0){
	  if(IS_INACTIVE(bin_ids[i])){ // TODO: do correctly...  
		return i 
	  }
	}
  }
  return -1;
}

int spawn_children(DAQPNode *node, int branch_id, DAQPBnBWorkspace* work){
  // TODO lifo vs fifo  
  DAQPNode* child1 = work->tree[start_id+n_tree] //TODO might be a bug when n_tree=0 (overwriting_node); 
  DAQPNode* child2 = work->tree[start_id+n_tree+1]
  for(int i =0;i<work->nb;i++){ // Copy
	child1->bin_ids[i] = node->bin_ids[i]; // TODO if lifo we can copy... 
	child2->bin_ids[i] = node->bin_ids[i];
  }
  child1->bin_ids[branch_id] = UPPER_BIN;
  child2->bin_ids[branch_id] = LOWER_BIN;

  work->n_tree+=2;
}
