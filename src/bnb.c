#include <stdio.h>
#include <stdlib.h>
#include "bnb.h"

int daqp_bnb(DAQPWorkspace* work){
  int branch_id, node_id, exitflag;
  c_float *swp_ptr;

  printf("REMOVE_LOWER_FLAG(18):%d\n",REMOVE_LOWER_FLAG(18));
  // Setup root node
  work->bnb->tree_depths[0]=-1;
  work->bnb->n_nodes=1;
  
  // Start tree exploration
  while(work->bnb->n_nodes > 0){

	node_id = --work->bnb->n_nodes;
	exitflag = process_node(node_id,work); // Solve relaxation

	// Cut conditions
	if(exitflag==EXIT_INFEASIBLE) continue; // Dominance cut
	if(exitflag<0) return exitflag; // Inner solver failed => abort

	// Find index to branch over 
	branch_id = get_branch_id(work); 
	if(branch_id==-1){// Nothing to branch over => integer feasible
	  work->settings->fval_bound = work->fval;
	  // Set ubar = u
	  swp_ptr=work->xold; work->xold= work->u; work->u=swp_ptr;
	}
	else{
	  spawn_children(node_id, branch_id, work);
	}
  }
  // Let work->u point to the best feasible solution 
  swp_ptr=work->u; work->u= work->xold; work->xold=swp_ptr;
  return 1;
}

int get_branch_id(DAQPWorkspace* work){
  // TODO: pick the most fractional? 
  for(int i=0; i < work->bnb->nb; i++){
	// Branch on first inactive constraint 
	if(IS_ACTIVE(work->bnb->bin_ids[i])) continue;
	return work->bnb->bin_ids[i];
  }
  return -1;
}

void spawn_children(const int node_id, const int branch_id, DAQPWorkspace* work){

  // Update child1 (reuse current node) 
  work->bnb->tree_bin_ids[node_id] = ADD_LOWER_FLAG(branch_id);
  work->bnb->tree_depths[node_id] +=1;

  // Update child2
  work->bnb->tree_bin_ids[node_id+1] = branch_id;
  work->bnb->tree_depths[node_id+1] = work->bnb->tree_depths[node_id];
  
  work->bnb->n_nodes+=2;
}

int process_node(const int node_id, DAQPWorkspace* work){
  int i,exitflag;
  const int depth = work->bnb->tree_depths[node_id]; 
  const int bin_id = work->bnb->tree_bin_ids[node_id]; 
  if(depth >=0){
	// Cleanup sense 
	for(i=depth; i<work->n_active; i++)
	  work->sense[work->WS[i]]&=~(ACTIVE+IMMUTABLE);
	// Reset workspace, but keep previous binaries in WS 
	work->sing_ind=EMPTY_IND;
	work->n_active=depth;
	work->reuse_ind=depth; 

	i = REMOVE_LOWER_FLAG(bin_id);
	add_constraint(work,i,1.0);
	// Setup new binary constraint
	if(work->sing_ind!=EMPTY_IND) return EXIT_INFEASIBLE; // Disregard singular node 
	work->sense[i] |= (IMMUTABLE+EXTRACT_LOWER_FLAG(bin_id));  

	// Possibly activate more constraints (Warm-start)...
	// TODO
  }

  // Solve relaxation
  exitflag = daqp_ldp(work);
  return exitflag;
}

int setup_daqp_bnb(DAQPWorkspace* work, int* bin_ids, int nb){
  if(nb > work->n) return EXIT_OVERDETERMINED_INITIAL;
  if((work->bnb == NULL) && (nb >0)){
	work->bnb= malloc(sizeof(DAQPBnB));

	work->bnb->nb = nb;
	work->bnb->bin_ids = bin_ids;

	// Setup tree
	work->bnb->tree_bin_ids= malloc((work->bnb->nb+1)*sizeof(int));
	work->bnb->tree_depths = malloc((work->bnb->nb+1)*sizeof(int));
	work->bnb->n_nodes = 0; 
  }
  return 1;
}

void free_daqp_bnb(DAQPWorkspace* work){
  if(work->bnb != NULL){
	free(work->bnb->tree_bin_ids);
	free(work->bnb->tree_depths);
	free(work->bnb);
	work->bnb = NULL;
  }
}
