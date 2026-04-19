#include "bnb.h"

int daqp_bnb(DAQPWorkspace* work){
    int branch_id, exitflag;
    DAQPNode* node;
    c_float *swp_ptr = NULL;

    // Modify upper bound based on absolute/relative suboptimality tolerance
    c_float fval_bound0 = work->settings->fval_bound;
    c_float eps_r = 1/(1+work->settings->rel_subopt);
    work->settings->fval_bound = (fval_bound0 - work->settings->abs_subopt)*eps_r;

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

    exitflag = DAQP_EXIT_INFEASIBLE;
    // Start tree exploration
    while( work->bnb->n_nodes > 0 ){

        node = work->bnb->tree+(--work->bnb->n_nodes);
        exitflag = daqp_process_node(node,work); // Solve relaxation
        // Cut conditions
        if(exitflag==DAQP_EXIT_INFEASIBLE) continue; // Dominance cut
        if(exitflag<0) break; // Inner solver failed => exit loop

        // Find index to branch over
        branch_id = daqp_get_branch_id(work);
        if(branch_id==DAQP_EMPTY_IND){// Nothing to branch over => integer feasible
            work->settings->fval_bound = (0.5*work->fval - work->settings->abs_subopt)*eps_r;
            swp_ptr=work->xold; work->xold= work->u; work->u=swp_ptr; // Store feasible sol
        }
        else{
            daqp_spawn_children(node,branch_id, work);
        }
    }

    // Exploration completed
    work->iterations = work->bnb->itercount;
    if(swp_ptr==NULL){
        work->settings->fval_bound = fval_bound0;
        return exitflag < 0 ? exitflag : DAQP_EXIT_INFEASIBLE;
    }
    else{
        // Invert fval_bound = (0.5*fval_best - abs_subopt)*eps_r to recover fval_best
        work->fval = 2*work->settings->fval_bound/eps_r + 2*work->settings->abs_subopt;
        work->settings->fval_bound = fval_bound0;
        // Let work->u point to the best feasible solution
        swp_ptr=work->u; work->u= work->xold; work->xold=swp_ptr;
        return DAQP_EXIT_OPTIMAL;
    }
}

int daqp_process_node(DAQPNode* node, DAQPWorkspace* work){
    int exitflag;
    work->bnb->nodecount+=1;
    if(node->depth >=0){
        // Fix a binary constraints
        work->bnb->fixed_ids[node->depth] = node->bin_id;
        // Setup relaxation
        if(work->bnb->n_nodes==0 || (node-1)->depth!=node->depth){
            // Sibling has been processed => need to fix workspace state
            work->bnb->n_clean += (node->depth-(node+1)->depth);
            daqp_node_cleanup_workspace(work->bnb->n_clean,work);
            daqp_warmstart_node(node,work);
        }
        else{
            daqp_add_upper_lower(node->bin_id,work);
            work->sense[DAQP_REMOVE_LOWER_FLAG(node->bin_id)] |= DAQP_IMMUTABLE; //Equality
            if(work->sing_ind != DAQP_EMPTY_IND) // Need to cold start to not miss integer feasible
                daqp_setup_cold_bnb(node,work);

        }
        // Add binary constraint
    }
    // Solve relaxation
    exitflag = daqp_ldp(work);
    work->bnb->itercount += work->iterations;

    if(exitflag == DAQP_EXIT_CYCLE){// Try to repair (cold start)
        daqp_setup_cold_bnb(node,work);
        exitflag = daqp_ldp(work);
        work->bnb->itercount += work->iterations;
    }

    return exitflag;
}

int daqp_get_branch_id(DAQPWorkspace* work){
    int i, j, disp;
    int branch_id = DAQP_EMPTY_IND;
    c_float diff;

    for(i=0; i < work->bnb->nb; i++){
        branch_id = work->bnb->bin_ids[i];
        // Skip fixed binary constraints
        if(DAQP_IS_ACTIVE(branch_id)) continue;

        // Compute signed distance from midpoint between bounds
        diff = 0.5*(work->dupper[branch_id]+work->dlower[branch_id]);
        if(branch_id < work->ms){//Simple bound
            if(work->Rinv==NULL) diff -= work->u[branch_id]; //Hessian is identity
            else{
                for(j=branch_id,disp=branch_id+DAQP_R_OFFSET(branch_id,work->n);j<work->n;j++)
                    diff -= work->Rinv[disp++]*work->u[j];
            }
        }
        else{//General bound
            for(j=0,disp=work->n*(branch_id-work->ms);j<work->n;j++)
                diff -= work->M[disp++]*work->u[j];
        }

        // Branch on the first infeasible binary variable
        return diff < 0 ? branch_id : DAQP_ADD_LOWER_FLAG(branch_id);
    }

    return DAQP_EMPTY_IND;
}

void daqp_spawn_children(DAQPNode* node, const int branch_id, DAQPWorkspace* work){

    daqp_save_warmstart(node,work);

    // Update child1 (reuse current node)
    node->bin_id = DAQP_TOGGLE_LOWER_FLAG(branch_id);
    node->depth +=1;

    // Update child2
    (node+1)->bin_id = branch_id;
    (node+1)->depth = node->depth;
    (node+1)->WS_start = node->WS_start;
    (node+1)->WS_end= node->WS_end;

    work->bnb->n_nodes+=2;
}

void daqp_node_cleanup_workspace(int n_clean, DAQPWorkspace* work){
    int i;
    // Cleanup sense
    for(i=n_clean; i<work->n_active; i++)
        work->sense[work->WS[i]]&= DAQP_IS_BINARY(work->WS[i]) ?
            ~(DAQP_ACTIVE+DAQP_IMMUTABLE): ~DAQP_ACTIVE;
    // Reset workspace
    work->sing_ind=DAQP_EMPTY_IND;
    work->n_active=n_clean;
    work->reuse_ind=n_clean;
}


void daqp_warmstart_node(DAQPNode* node, DAQPWorkspace* work){
    int i;
    // Add fixed constraints
    for(i=work->bnb->n_clean - work->bnb->neq; i< node->depth+1;i++){
        daqp_add_upper_lower(work->bnb->fixed_ids[i],work);
        DAQP_SET_IMMUTABLE(DAQP_REMOVE_LOWER_FLAG(work->bnb->fixed_ids[i]));
    }
    work->bnb->n_clean = work->bnb->neq+node->depth;
    // Add free constraints
    for(i=node->WS_start; i < node->WS_end; i++){
        daqp_add_upper_lower(work->bnb->tree_WS[i],work);
        if(work->sing_ind != DAQP_EMPTY_IND) {
            work->n_active--;
            DAQP_SET_INACTIVE(work->WS[work->n_active]);
            work->sing_ind = DAQP_EMPTY_IND;
            break; // Abort warm start if singular basis
        }
    }
    work->bnb->nWS = node->WS_start; // always move up tree after warmstart
}

void daqp_save_warmstart(DAQPNode* node, DAQPWorkspace* work){
    int id_to_add, i;
    // Save warmstart
    node->WS_start = work->bnb->nWS;

    for(i =work->bnb->neq; i<work->n_active;i++){
        id_to_add = (work->WS[i]+(DAQP_IS_LOWER(work->WS[i]) << (DAQP_LOWER_BIT-1)));
        if((work->sense[work->WS[i]]&(DAQP_IMMUTABLE+DAQP_BINARY))!=DAQP_IMMUTABLE+DAQP_BINARY)
            work->bnb->tree_WS[work->bnb->nWS++]= id_to_add;
    }
    node->WS_end = work->bnb->nWS;
}

int daqp_add_upper_lower(const int add_id, DAQPWorkspace* work){
    int true_add_id = DAQP_REMOVE_LOWER_FLAG(add_id);
    // Setup new binary constraint
    if(DAQP_EXTRACT_LOWER_FLAG(add_id)){
        DAQP_SET_LOWER(true_add_id);
        daqp_add_constraint(work,true_add_id,-1.0);
    }
    else{
        DAQP_SET_UPPER(true_add_id);
        daqp_add_constraint(work,true_add_id,1.0);
    }
    return 1;
}

void daqp_setup_cold_bnb(DAQPNode* node,DAQPWorkspace *work){
    int i;
    daqp_node_cleanup_workspace(work->bnb->n_clean,work);
    for(i=work->bnb->n_clean - work->bnb->neq; i< node->depth+1;i++){
        daqp_add_upper_lower(work->bnb->fixed_ids[i],work);
        DAQP_SET_IMMUTABLE(DAQP_REMOVE_LOWER_FLAG(work->bnb->fixed_ids[i]));
    }
    work->bnb->n_clean = work->bnb->neq+node->depth;
}
