#include <stdio.h>
#include <stdlib.h>
#include "bnb.h"

int daqp_bnb(DAQPWorkspace* work){
  int branch_id, exitflag;
  DAQPNode* node;
  c_float *swp_ptr = NULL;

  work->bnb->itercount=0;
  work->bnb->nodecount=0;
  // Setup root node
  work->bnb->tree[0].depth=-1;
  work->bnb->tree[0].WS_start=0;
  work->bnb->tree[0].WS_end=0;
  work->bnb->n_nodes=1;
  work->bnb->tree[0].bin_id=0;
  
  // Start tree exploration

  while(work->bnb->n_nodes > 0){

	node = work->bnb->tree+(--work->bnb->n_nodes);
	exitflag = process_node(node,work); // Solve relaxation

	// Cut conditions
	if(exitflag==EXIT_INFEASIBLE) continue; // Dominance cut
	if(exitflag<0) return exitflag; // Inner solver failed => abort

	// Find index to branch over 
	save_warmstart(node,work);
	branch_id = get_branch_id(work); 
	if(branch_id==-1){// Nothing to branch over => integer feasible
	  work->settings->fval_bound = work->fval;
	  swp_ptr=work->xold; work->xold= work->u; work->u=swp_ptr; // Store feasible sol
	}
	else{
	  spawn_children(node,branch_id, work);
	}
  }

  // Exploration completed 
  work->iterations = work->bnb->itercount;
  if(swp_ptr==NULL)
	return EXIT_INFEASIBLE;
  else{
	// Let work->u point to the best feasible solution 
	swp_ptr=work->u; work->u= work->xold; work->xold=swp_ptr;
	return EXIT_OPTIMAL;
  }
}

int process_node(DAQPNode* node, DAQPWorkspace* work){
  int exitflag;
  work->bnb->nodecount+=1;
  if(node->depth >=0){
	if(work->bnb->n_nodes==0 || (node-1)->depth!=node->depth){ 
	  // Sibling has been processed => need to fix workspace state
	  node_cleanup_workspace(work);
	  warmstart_node(node->depth, node->WS_start,node->WS_end,work);
	}
	// Add binary constraint 
	add_upper_lower(node->bin_id,work);
	work->sense[REMOVE_LOWER_FLAG(node->bin_id)] |= IMMUTABLE; // Make equality
  }
  // Solve relaxation
  exitflag = daqp_ldp(work);
  work->bnb->itercount += work->iterations;

  return exitflag;
}

int get_branch_id(DAQPWorkspace* work){
  int i,disp;
  int branch_id = EMPTY_IND;
  for(i=0; i < work->bnb->nb; i++){
	// Branch on first inactive constraint 
	if(IS_ACTIVE(work->bnb->bin_ids[i])) continue;
	branch_id = work->bnb->bin_ids[i];
	break;
  }

  if(branch_id == EMPTY_IND) return EMPTY_IND; // Nothing to branch over (=>integer feasible)

  // Determine if upper or lower child should be processed first 
  // by computing whether the upper or lower bound is closer to be activated
  c_float diff = 0.5*(work->dupper[branch_id]+work->dlower[branch_id]);
  if(IS_SIMPLE(branch_id)){//Simple bound
	if(work->Rinv==NULL) diff-=work->u[branch_id]; //Hessian is identify 
	else{
	  for(i=branch_id,disp=branch_id+R_OFFSET(branch_id,NX);i<NX;i++) // 
		diff-=work->Rinv[disp++]*work->u[i];
	}
  }
  else{//General bound
	for(i=0,disp=NX*(branch_id-N_SIMPLE);i<NX;i++) 
	  diff-=work->M[disp++]*work->u[i];
  }
  branch_id = diff<0 ? branch_id : ADD_LOWER_FLAG(branch_id);
  return branch_id;
}

void spawn_children(DAQPNode* node, const int branch_id, DAQPWorkspace* work){

  // Update child1 (reuse current node) 
  node->bin_id = TOGGLE_LOWER_FLAG(branch_id);
  node->depth +=1;

  // Update child2
  (node+1)->bin_id = branch_id;
  (node+1)->depth = node->depth; 
  (node+1)->WS_start = node->WS_start;
  (node+1)->WS_end= node->WS_end;
  
  work->bnb->n_nodes+=2;
}

void node_cleanup_workspace(DAQPWorkspace* work){
  // Cleanup sense 
  for(int i=0; i<work->n_active; i++)
	work->sense[work->WS[i]]&= IS_BINARY(work->WS[i]) ? ~(ACTIVE+IMMUTABLE): ~ACTIVE;
  // Reset workspace
  work->sing_ind=EMPTY_IND;
  work->n_active=0;
  work->reuse_ind=0; 
}


void warmstart_node(const int depth, const int start_id, const int end_id, DAQPWorkspace* work){
  int i;
  for(i=start_id;i<end_id;i++){
	if(work->sing_ind != EMPTY_IND) break; // Abort warm start if singular basis 
	// Add the constraint
	add_upper_lower(work->bnb->tree_WS[i],work);
  }
  for(i=0;i<depth;i++) // Make binaries immutable 
	work->sense[REMOVE_LOWER_FLAG(work->WS[i])] |= IMMUTABLE; 
}

void save_warmstart(DAQPNode* node, DAQPWorkspace* work){
  // Save warmstart 
  if(work->bnb->n_nodes==0 || (node-1)->depth!=node->depth){ // Moving up the tree
    work->bnb->nWS = node->WS_start;
  }
  node->WS_start = work->bnb->nWS;

  int nb_added = 0;
  int id_to_add;
  work->bnb->nWS+=node->depth+1;
  for(int i =0; i<work->n_active;i++){
	id_to_add = (work->WS[i]+(IS_LOWER(work->WS[i]) << (LOWER_BIT-1)));
	if((work->sense[work->WS[i]]&(IMMUTABLE+BINARY))==IMMUTABLE+BINARY){
	  work->bnb->tree_WS[node->WS_start+(nb_added++)]= id_to_add; 
	}
	else{
	  work->bnb->tree_WS[work->bnb->nWS++]= id_to_add;
	}
  }
  node->WS_end = work->bnb->nWS;
}

int add_upper_lower(const int add_id, DAQPWorkspace* work){
  int true_add_id = REMOVE_LOWER_FLAG(add_id);
  // Setup new binary constraint
  if(EXTRACT_LOWER_FLAG(add_id)){
	SET_LOWER(true_add_id);
	add_constraint(work,true_add_id,-1.0);
  }
  else{
	SET_UPPER(true_add_id);
	add_constraint(work,true_add_id,1.0);
  }
  return 1;
}

int setup_daqp_bnb(DAQPWorkspace* work, int* bin_ids, int nb){
  if(nb > work->n) return EXIT_OVERDETERMINED_INITIAL;
  if((work->bnb == NULL) && (nb >0)){
	work->bnb= malloc(sizeof(DAQPBnB));

	work->bnb->nb = nb;
	work->bnb->bin_ids = bin_ids;

	// Setup tree
	work->bnb->tree= malloc((work->bnb->nb+1)*sizeof(DAQPNode));
	work->bnb->tree_WS= malloc((work->n+1)*(work->bnb->nb+1)*sizeof(int));
	work->bnb->n_nodes = 0; 
	work->bnb->nWS= 0; 
  }
  return 1;
}

void free_daqp_bnb(DAQPWorkspace* work){
  if(work->bnb != NULL){
	free(work->bnb->tree);
	free(work->bnb->tree_WS);
	free(work->bnb);
	work->bnb = NULL;
  }
}
