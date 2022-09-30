#include <stdio.h>
#include <stdlib.h>
#include "bnb.h"

int daqp_bnb(DAQPWorkspace* work){
    int branch_id, exitflag;
    DAQPNode* node;
    c_float *swp_ptr = NULL;
    c_float fval_bound0 = work->settings->fval_bound;
    work->bnb->neq = work->n_active; 

    work->bnb->itercount=0;
    work->bnb->nodecount=0;
    // Setup root node
    work->bnb->tree[0].depth=-1;
    work->bnb->tree[0].WS_start=0;
    work->bnb->tree[0].WS_end=0;
    work->bnb->tree[0].bin_id=0;
    work->bnb->n_nodes=1;
    work->bnb->n_clean=work->bnb->neq;

    // Start tree exploration

    while( work->bnb->n_nodes > 0 ){

        node = work->bnb->tree+(--work->bnb->n_nodes);
        exitflag = process_node(node,work); // Solve relaxation

        // Cut conditions
        if(exitflag==EXIT_INFEASIBLE) continue; // Dominance cut
        if(exitflag<0) return exitflag; // Inner solver failed => abort

        // Find index to branch over 
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
    // Correct fval
    work->fval = work->settings->fval_bound;
    work->settings->fval_bound = fval_bound0;
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
            work->bnb->n_clean += (node->depth-(node+1)->depth);
            node_cleanup_workspace(work->bnb->n_clean,work);
            warmstart_node(node,work);
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

    save_warmstart(node,work);

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

void node_cleanup_workspace(int n_clean, DAQPWorkspace* work){
    // Cleanup sense 
    for(int i=n_clean; i<work->n_active; i++)
        work->sense[work->WS[i]]&= IS_BINARY(work->WS[i]) ? ~(ACTIVE+IMMUTABLE): ~ACTIVE;
    // Reset workspace
    work->sing_ind=EMPTY_IND;
    work->n_active=n_clean;
    work->reuse_ind=n_clean; 
}


void warmstart_node(DAQPNode* node, DAQPWorkspace* work){
    int i,n_clean_old;
    n_clean_old = work->bnb->n_clean;
    work->bnb->n_clean = work->bnb->neq+node->depth; 
    for(i=node->WS_start + n_clean_old-work->bnb->neq; i< node->WS_end ;i++){
        if(work->sing_ind != EMPTY_IND) break; // Abort warm start if singular basis 
                                               // Add the constraint
        add_upper_lower(work->bnb->tree_WS[i],work);
    }
    for(i=n_clean_old;i<work->bnb->n_clean;i++) // Make binaries immutable 
        work->sense[work->WS[i]] |= IMMUTABLE; 
    work->bnb->nWS = node->WS_start; // always move up tree after warmstart 
}

void save_warmstart(DAQPNode* node, DAQPWorkspace* work){
    // Save warmstart 
    node->WS_start = work->bnb->nWS;

    int nb_added = 0;
    int id_to_add;
    work->bnb->nWS+=(node->depth+1);
    for(int i =work->bnb->neq; i<work->n_active;i++){
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
