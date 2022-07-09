#include <stdio.h>
#include <stdlib.h>
#include "bnb.h"

int daqp_bnb(DAQPWorkspace* work){
  int branch_id, exitflag;
  DAQPNode* node;
  c_float *swp_ptr;

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
	  // Set ubar = u
	  swp_ptr=work->xold; work->xold= work->u; work->u=swp_ptr;
	}
	else{
	  spawn_children(node,branch_id, work);
	}
  }
  // Let work->u point to the best feasible solution 
  swp_ptr=work->u; work->u= work->xold; work->xold=swp_ptr;
  work->iterations = work->bnb->itercount;
  return 1;
}

int get_branch_id(DAQPWorkspace* work){
  int i,disp;
  int branch_id = -1;
  for(i=0; i < work->bnb->nb; i++){
	// Branch on first inactive constraint 
	if(IS_ACTIVE(work->bnb->bin_ids[i])) continue;
	branch_id = work->bnb->bin_ids[i];
	break;
  }

  if(branch_id == -1) return -1; // No index to branch over (=>integer feasible)

  // Determine if upper or lower child should be processed first 
  // by computing whether the upper or lower bound is closer to be activated
  c_float diff = 0.5*(work->dupper[branch_id]+work->dlower[branch_id]);
  if(branch_id < N_SIMPLE){//Simple bound
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


void save_warmstart(DAQPNode* node, DAQPWorkspace* work){
  // Save warmstart 
  if(work->bnb->n_nodes==0 || (node-1)->depth!=node->depth){ // Moving up the tree
    work->bnb->nWS = node->WS_start;
  }
  node->WS_start = work->bnb->nWS;

  for(int i =node->depth+1; i<work->n_active;i++){
    work->bnb->tree_WS[work->bnb->nWS++]=(work->WS[i]+(IS_LOWER(work->WS[i]) << (LOWER_BIT-1)));
  }
  node->WS_end = work->bnb->nWS;
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

int process_node(DAQPNode* node, DAQPWorkspace* work){
  int exitflag;
  work->bnb->nodecount+=1;
  if(node->depth >=0){
	// Setup workspace before processing node 
	add_new_binary(node->depth,node->bin_id,work);
	// Warm start 
	warmstart_node(node->WS_start,node->WS_end,work);
  }
  // Solve relaxation
  exitflag = daqp_ldp(work);
  work->bnb->itercount += work->iterations;

  return exitflag;
}

int add_new_binary(const int depth, const int bin_id, DAQPWorkspace* work){
  int i;
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
  work->sense[i] |= IMMUTABLE;
  if(EXTRACT_LOWER_FLAG(bin_id))
	SET_LOWER(i);
  else
	SET_UPPER(i);
  return 1;
}

void warmstart_node(const int start_id, const int end_id, DAQPWorkspace* work){
  int i,add_id;
  for(i=start_id;i<end_id;i++){
	if(work->sing_ind != EMPTY_IND) break; // Abort warm start if singular basis 
	// Add the constraint
	add_id = REMOVE_LOWER_FLAG(work->bnb->tree_WS[i]);
	if(EXTRACT_LOWER_FLAG(work->bnb->tree_WS[i])){
	  SET_LOWER(add_id);
	  add_constraint(work,add_id,-1.0);
	}
	else{
	  SET_UPPER(add_id);
	  add_constraint(work,add_id,1.0);
	}
  }
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
