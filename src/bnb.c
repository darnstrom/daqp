#include <stdio.h>
#include <stdlib.h>
#include "bnb.h"

int daqp_bnb(DAQPWorkspace* work){
  DAQPNode* current_node;
  int branch_id,exitflag;
  c_float *swp_ptr;

  // Setup root node
  work->bnb->tree[0].depth=-1;
  work->bnb->n_tree=1;
  
  // Start tree exploration
  int iter = 0;
  while((work->bnb->n_tree > 0) && (++iter < 6)){
	current_node = work->bnb->tree+(--work->bnb->n_tree); 

	exitflag = process_node(current_node,work); // Solve relaxation

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
	  spawn_children(current_node, branch_id, work);
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
	return i;
  }
  return -1;
}

void spawn_children(DAQPNode *node, int branch_id, DAQPWorkspace* work){
  // Since lifo current node can be reused as a new child
  DAQPNode* new_child = work->bnb->tree+(work->bnb->n_tree+1);

  // Update child1 (reuse node) 
  node->new_bin = branch_id; // LOWER 
  node->depth++;
  node->is_lower = 2;

  // Update child2
  new_child->new_bin= branch_id; // UPPER
  new_child->depth=node->depth;
  new_child->is_lower = 0;
  
  work->bnb->n_tree+=2;
}

// TODO need to forbid pivoting
int process_node(DAQPNode *node, DAQPWorkspace* work){
  int i,exitflag;
  if(node->depth >=0){
	// Cleanup sense 
	for(i=node->depth; i<work->n_active; i++)
	  work->sense[i]&=~(ACTIVE+IMMUTABLE);
	// Reset workspace, but keep previous binaries in WS 
	work->sing_ind=EMPTY_IND;
	work->n_active=node->depth;
	work->reuse_ind=node->depth; 

	add_constraint(work,node->new_bin,1.0);
	// Setup new binary constraint
	if(work->sing_ind!=EMPTY_IND) return EXIT_INFEASIBLE; // Disregard singular node 
	work->sense[node->new_bin] |= (IMMUTABLE+node->is_lower);  

	// Possibly activate more constraints (Warm-start)...
	// TODO
  }

  // Solve relaxation
  exitflag = daqp_ldp(work);
  return exitflag;
}

void setup_daqp_bnb(DAQPWorkspace* work, int* bin_ids, int nb){
  if((work->bnb == NULL) && (nb >0)){
	work->bnb= malloc(sizeof(DAQPBnB));

	work->bnb->nb = nb;
	work->bnb->bin_ids = bin_ids;

	// Setup tree
	work->bnb->tree= malloc((work->bnb->nb+1)*sizeof(DAQPNode));
	work->bnb->n_tree = 0; 
  }
}

void free_daqp_bnb(DAQPWorkspace* work){
  if(work->bnb != NULL){
	free(work->bnb->tree);
	free(work->bnb);
	work->bnb = NULL;
  }
}
