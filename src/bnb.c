#include <stdio.h>
#include <stdlib.h>
#include "bnb.h"
#include "daqp.h"

int daqp_bnb(DAQPBnBWorkspace* work){
  DAQPNode current_node;
  int branch_id,exitflag;
  c_float *swp_ptr;

  // Setup root node
  work->tree[0].depth=-1;
  work->n_tree=1;
  
  // Start tree exploration
  while(work->n_tree > 0){
	current_node = work->tree[work->n_tree--]; 

	exitflag = process_node(&current_node,work); // Solve relaxation

	// Cut conditions
	if(exitflag==EXIT_INFEASIBLE) continue; // Dominance cut
	if(exitflag<0) return exitflag; // Inner solver failed => abort

	// Find index to branch over 
	branch_id = get_branch_id(work); 
	if(branch_id==-1){// Nothing to branch over => integer feasible
	  work->dwork->settings->fval_bound = work->dwork->fval;
	  // Set ubar = u
	  swp_ptr=work->ustar; work->ustar= work->dwork->u; work->dwork->u=swp_ptr;
	}
	else{
	  spawn_children(&current_node, branch_id, work);
	}
  }
  return 1;
}

int get_branch_id(DAQPBnBWorkspace* work){
  // TODO: pick the most fractional? 
  for(int i=0; i<work->nb; i++){
	// Branch on first inactive constraint 
	if((work->dwork->sense[work->bin_ids[i]]&1) != 1){
	  return i;
	}
  }
  return -1;
}

void spawn_children(DAQPNode *node, int branch_id, DAQPBnBWorkspace* work){
  // Since lifo current node can be reused as a new child
  DAQPNode* new_child = work->tree+(work->n_tree+1);

  // Update child1 (reuse node) 
  node->new_bin = MASK_LOWER(branch_id); // LOWER 
  node->depth++;

  // Update child2
  new_child->new_bin= branch_id; // UPPER
  new_child->depth=node->depth;
  
  work->n_tree+=2;
}

// TODO need to forbid pivoting
int process_node(DAQPNode *node, DAQPBnBWorkspace* bnb_work){
  int i,exitflag;
  DAQPWorkspace* work = bnb_work->dwork;

  if(node->depth>=0){
	// Cleanup sense 
	for(i=node->depth; i<work->n_active; i++)
	  work->sense[i]&=~(ACTIVE+IMMUTABLE);
	// Reset workspace, but keep previous binaries in WS 
	work->sing_ind=EMPTY_IND;
	work->n_active=node->depth;
	work->reuse_ind=node->depth; 

	// Setup new binary constraint
	i = UNMASK_BIN_ID(node->new_bin);
	add_constraint(work,i,1.0);
	if(work->sing_ind!=EMPTY_IND) return EXIT_INFEASIBLE; // Disregard singular node 
	work->sense[i] |= (IMMUTABLE+UNMASK_LOWER(node->new_bin));  

	// Possibly activate more constraints (Warm-start)...
	// TODO
  }

  // Solve relaxation
  exitflag = daqp_ldp(work);
  return exitflag;
}

void setup_daqp_bnb(DAQPWorkspace* work, DAQPBnBWorkspace* bnb_work){
  int i,nb;

  // Append workspace 
  bnb_work->dwork = work;
  // Setup binary ids 
  bnb_work->bin_ids= malloc(work->n*sizeof(int)); // nb>n => overdetermined system 
  nb = 0;
  for(i = 0; i<work->m;i++){
	if(work->sense[i]&BINARY){
	  bnb_work->bin_ids[nb++]=i; 
	  if(nb>work->n) return; // Cannot have more than n binaries 
	}

	// Allocate memory 
	bnb_work->ustar= malloc(work->n*sizeof(int));
	bnb_work->tree= malloc(bnb_work->nb*sizeof(DAQPNode));
	bnb_work->n_tree = 0; 
  }
}

void free_daqp_bnb(DAQPBnBWorkspace* bnb_work){
	free(bnb_work->ustar);
	free(bnb_work->tree);
	free(bnb_work->bin_ids);
	//free_daqp_workspace(bnb_work->dwork);
}
