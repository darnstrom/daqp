#include "bnb.h"
#ifdef PROFILING
#include "utils.h"
#endif

static c_float daqp_binary_diff(const int id, DAQPWorkspace* work){
    int j, disp;
    c_float diff = 0.5*(work->dupper[id]+work->dlower[id]);

    if(id < work->ms){//Simple bound
        if(work->Rinv==NULL) diff -= work->u[id]; //Hessian is identity
        else{
            for(j=id,disp=id+DAQP_R_OFFSET(id,work->n);j<work->n;j++)
                diff -= work->Rinv[disp++]*work->u[j];
        }
    }
    else{//General bound; daqp_add_infeasible already computed M*u
        diff -= work->Mu[id-work->ms];
    }
    return diff;
}

// Store the free part of the working set in ids. Returns the number stored.
static int daqp_bnb_store_ws(int* ids, DAQPWorkspace* work){
    int i, n_ids = 0;
    for(i=work->bnb->neq; i<work->n_active;i++){
        if((work->sense[work->WS[i]]&(DAQP_IMMUTABLE+DAQP_BINARY))!=DAQP_IMMUTABLE+DAQP_BINARY)
            ids[n_ids++] = work->WS[i]+(DAQP_IS_LOWER(work->WS[i]) << (DAQP_LOWER_BIT-1));
    }
    return n_ids;
}

// Add the constraints in ids to the working set (aborted if the basis gets singular)
static void daqp_bnb_load_ws(const int* ids, const int n_ids, DAQPWorkspace* work){
    int i;
    for(i=0; i < n_ids; i++){
        // E.g., a binary that the new bounds fix is already active (as an equality)
        if(DAQP_IS_ACTIVE(DAQP_REMOVE_LOWER_FLAG(ids[i]))) continue;
        daqp_add_upper_lower(ids[i],work);
        if(work->sing_ind != DAQP_EMPTY_IND) {
            work->n_active--;
            DAQP_SET_INACTIVE(work->WS[work->n_active]);
            work->sing_ind = DAQP_EMPTY_IND;
            break;
        }
    }
}

// Add fixed binary id; a dependent one is left out: returns 1 if added, 0 if
// consistent, and -1 if inconsistent (the node is infeasible)
static int daqp_bnb_add_fixed(const int id, DAQPWorkspace* work){
    int j;
    c_float residual = 0, term;
    // Same tolerance as for the binary in daqp_get_branch_id
    const c_float tol = work->settings->primal_tol*
        (work->scaling != NULL ? work->scaling[DAQP_REMOVE_LOWER_FLAG(id)] : 1);
    daqp_add_upper_lower(id,work);
    if(work->sing_ind == DAQP_EMPTY_IND){
        DAQP_SET_IMMUTABLE(DAQP_REMOVE_LOWER_FLAG(id));
        return 1;
    }
    daqp_compute_singular_direction(work);
    for(j = 0; j < work->n_active; j++){
        const int k = work->WS[j];
        term = work->lam_star[j]*(DAQP_IS_LOWER(k) ? work->dlower[k] : work->dupper[k]);
        residual += term;
    }
    daqp_drop_singular_last(work);
    return (residual > tol || residual < -tol) ? -1 : 0;
}

// Add the fixed binaries of node (n_clean falls back to neq if one is left
// out). Returns -1 if the node is infeasible
static int daqp_bnb_add_fixed_ids(DAQPNode* node, DAQPWorkspace* work){
    int i, flag, dropped = 0;
    for(i=work->bnb->n_clean - work->bnb->neq; i< node->depth+1;i++){
        flag = daqp_bnb_add_fixed(work->bnb->fixed_ids[i],work);
        if(flag < 0){
            work->bnb->n_clean = work->bnb->neq;
            return -1;
        }
        dropped |= flag == 0;
    }
    work->bnb->n_clean = dropped ? work->bnb->neq : work->bnb->neq+node->depth;
    return 0;
}

// The immutable constraints (e.g., equalities) are kept fixed as a prefix of
// the working set throughout the tree. Returns the length of that prefix.
static int daqp_bnb_setup_root(DAQPWorkspace* work){
    int i, j, error_flag;
    for(i = 0; i < work->n_active; i++)
        if(!DAQP_IS_IMMUTABLE(work->WS[i])) break;
    for(j = i; j < work->n_active; j++)
        if(DAQP_IS_IMMUTABLE(work->WS[j])) break;
    if(j == work->n_active) return i; // Mutable constraints are a warm start

    // Immutable after mutable => only activate immutable constraints
    for(i = 0; i < work->n_active; i++)
        if(!DAQP_IS_IMMUTABLE(work->WS[i])) DAQP_SET_INACTIVE(work->WS[i]);
    reset_daqp_workspace(work);
    error_flag = daqp_activate_constraints(work);
    return error_flag < 0 ? error_flag : work->n_active;
}

// Use the candidate in work->x as incumbent. Returns its objective
// (internal scale) with u = R*x+v stored in work->xold or -1 if infeasible.
static c_float daqp_bnb_incumbent(DAQPWorkspace* work){
    int i, j, disp;
    const int n = work->n;
    c_float val, fval = 0, tol = work->settings->primal_tol;
    c_float *x = work->x, *u = work->xold;
    DAQPProblem* qp = work->qp;
    if(qp == NULL) return -1;

    // Check feasibility (soft constraints are treated as hard)
    for(i=0, disp=0; i < work->m; i++){
        if(i < work->ms) val = x[i];
        else for(j=0, val=0; j < n; j++) val += qp->A[disp++]*x[j];
        if(DAQP_IS_IMMUTABLE(i) && !DAQP_IS_ACTIVE(i) && !DAQP_IS_BINARY(i)) continue; // Ignored
        if(val > qp->bupper[i]+tol || val < qp->blower[i]-tol) return -1;
        if(DAQP_IS_BINARY(i) && val > qp->blower[i]+tol && val < qp->bupper[i]-tol) return -1;
    }

    // Invert ldp2qp_solution: u = R*x + v
    for(i=0; i < n; i++) u[i] = x[i];
    if(work->Rinv != NULL){
        if(work->scaling != NULL)
            for(i=0; i < work->ms; i++) u[i] *= work->scaling[i];
        for(i=n-1; i >= 0; i--){ // Back substitution with upper triangular Rinv
            disp = i+DAQP_R_OFFSET(i,n);
            for(j=i+1; j < n; j++) u[i] -= work->Rinv[disp+j-i]*u[j];
            u[i] /= work->Rinv[disp];
        }
    }
    else if(work->RinvD != NULL)
        for(i=0; i < n; i++) u[i] /= work->RinvD[i];
    if(work->v != NULL)
        for(i=0; i < n; i++) u[i] += work->v[i];

    for(i=0; i < n; i++) fval += u[i]*u[i];
    return fval;
}

// The regularization sum_j reg_j/2*x_j^2 of daqp_bnb_prox at u, which is not part of
// the objective J (only needed if rel_subopt > 0 or for the lower bound lb)
static c_float daqp_bnb_reg(const c_float* u, const c_float lb, DAQPWorkspace* work){
    int i;
    c_float val = 0, *x = work->bnb->xtmp;
    if(work->n_prox == 0 || x == NULL || (work->settings->rel_subopt <= 0 && lb <= -DAQP_INF))
        return 0;
    for(i = 0; i < work->n; i++) x[i] = work->v != NULL ? u[i]-work->v[i] : u[i];
    daqp_apply_Rinv(work,x);
    for(i = 0; i < work->n; i++) val += 0.5*work->bnb->reg[i]*x[i]*x[i];
    return val;
}

// Nodes must improve on the objective J by more than abs_subopt + rel_subopt*|J|.
// The bound is in the scale 0.5*fval = J + c used by daqp_ldp. If J is within this of
// the lower bound lb, no node can, and all are pruned
static c_float daqp_bnb_bound(c_float fval, c_float c, const c_float lb, DAQPSettings* settings){
    const c_float J = 0.5*fval - c;
    const c_float tol = settings->abs_subopt + settings->rel_subopt*(J < 0 ? -J : J);
    return (lb > -DAQP_INF && J - lb <= tol) ? -DAQP_INF : 0.5*fval - tol;
}

// Fixing binaries that dominate a row of M (e.g., big-M or small curvature on the binaries)
// leaves the fraction r of its squared norm that is on the other variables. Pivots of such
// rows are hence about r, and the singularity tolerance is made relative to the smallest r
static c_float daqp_bnb_sing_tol(DAQPWorkspace* work){
    int i, k;
    c_float rmin = 1, sb, tol;
    if(work->n_prox > 0) return work->settings->sing_tol; // Not for daqp_bnb_prox
    for(i = 0; i < work->m-work->ms; i++){ // The rows of M have unit norm => r = 1-sb
        const c_float* Mi = work->M+(size_t)i*work->n;
        for(k = 0, sb = 0; k < work->bnb->nb; k++)
            if(work->bnb->bin_ids[k] < work->ms) sb += Mi[work->bnb->bin_ids[k]]*Mi[work->bnb->bin_ids[k]];
        if(sb > 0 && 1-sb > 0 && 1-sb < rmin) rmin = 1-sb;
    }
    tol = rmin*work->settings->sing_tol;
    if(tol < DAQP_BNB_MIN_SING_TOL) // Limit by the accuracy of the pivots
        tol = work->settings->sing_tol < DAQP_BNB_MIN_SING_TOL ? work->settings->sing_tol : DAQP_BNB_MIN_SING_TOL;
    return tol;
}

int daqp_bnb(DAQPWorkspace* work, const c_float lb){
    int branch_id, exitflag;
    DAQPNode* node;
    c_float *swp_ptr = NULL;

    exitflag = daqp_bnb_setup_root(work);
    if(exitflag < 0) return exitflag;
    work->bnb->neq = exitflag;
    const c_float sing_tol0 = work->settings->sing_tol;
    work->settings->sing_tol = daqp_bnb_sing_tol(work);

    // Warm start the root with the root working set of the previous solve
    // (unless a warm start has been provided)
    if(work->n_active == work->bnb->neq)
        daqp_bnb_load_ws(work->bnb->root_WS,work->bnb->n_root_WS,work);

    // Objective J = 0.5*fval - c, with c = 0.5*||v||^2 - fp, where fp is the objective
    // removed by equality elimination (only needed for rel_subopt, fval_bound and lb)
    c_float fval_bound0 = work->settings->fval_bound, fval_best = 0, c = 0;
    if(work->settings->rel_subopt > 0 || fval_bound0 < DAQP_INF || lb > -DAQP_INF){
        if(work->v != NULL) for(int i=0; i < work->n; i++) c += 0.5*work->v[i]*work->v[i];
        if(work->eq != NULL && work->eq->installed) c -= work->eq->fp;
    }
    work->settings->fval_bound = fval_bound0 + c;

    // Start from a user-provided integer-feasible solution
    if(work->state & DAQP_STATE_INCUMBENT){
        work->state &= ~DAQP_STATE_INCUMBENT;
        fval_best = daqp_bnb_incumbent(work);
        if(fval_best >= 0 && 0.5*fval_best < work->settings->fval_bound){
            work->settings->fval_bound = daqp_bnb_bound(fval_best,c+daqp_bnb_reg(work->xold,lb,work),
                    lb,work->settings);
            swp_ptr = work->xold; // Marks that a feasible solution is stored in xold
        }
    }

    work->bnb->itercount=0;
    work->bnb->nodecount=0;
    // Setup root node
    work->bnb->tree[0].depth=-1;
    work->bnb->tree[0].WS_start=0;
    work->bnb->tree[0].WS_end=0;
    work->bnb->tree[0].bin_id=0;
    work->bnb->tree[0].lb=-DAQP_INF; // Always process the root
    work->bnb->n_nodes=1;
    work->bnb->n_clean=work->bnb->neq;
    work->bnb->nWS=0;

    exitflag = DAQP_EXIT_INFEASIBLE;
    // Start tree exploration
    while( work->bnb->n_nodes > 0 ){

        node = work->bnb->tree+(--work->bnb->n_nodes);
        if(node->lb > 2*work->settings->fval_bound) continue;
        exitflag = daqp_process_node(node,work); // Solve relaxation
        if(node->depth < 0 && exitflag > 0 && work->bnb->root_WS != NULL)
            work->bnb->n_root_WS = daqp_bnb_store_ws(work->bnb->root_WS,work);
#ifdef PROFILING
        // Individual relaxations are often too short to reach the timer check
        // in daqp_ldp, so also enforce the limit across the BnB tree.
        if(work->timer != NULL && (work->bnb->nodecount&31)==0){
            toc((DAQPtimer*)work->timer);
            if(get_time((DAQPtimer*)work->timer) > work->settings->time_limit){
                exitflag = DAQP_EXIT_TIMELIMIT;
                break;
            }
        }
#endif
        // Cut conditions
        if(exitflag==DAQP_EXIT_INFEASIBLE) continue; // Dominance cut
        if(exitflag<0) break; // Inner solver failed => exit loop

        // Find index to branch over
        branch_id = daqp_get_branch_id(work);
        if(branch_id==DAQP_EMPTY_IND){// Nothing to branch over => integer feasible
            fval_best = work->fval;
            work->settings->fval_bound = daqp_bnb_bound(fval_best,c+daqp_bnb_reg(work->u,lb,work),
                    lb,work->settings);
            swp_ptr=work->xold; work->xold= work->u; work->u=swp_ptr; // Store feasible sol
        }
        else{
            daqp_spawn_children(node,branch_id, work);
        }
    }

    // Exploration completed
    work->iterations = work->bnb->itercount;
    // Restore the root state (unfix binaries etc.) so that the workspace can
    // be reused for subsequent solves
    daqp_node_cleanup_workspace(work->bnb->neq,work);
    work->bnb->n_clean = work->bnb->neq;
    work->settings->fval_bound = fval_bound0;
    work->settings->sing_tol = sing_tol0;
    if(swp_ptr==NULL){
        return exitflag < 0 ? exitflag : DAQP_EXIT_INFEASIBLE;
    }
    else{
        work->fval = fval_best;
        // Let work->u point to the best feasible solution
        swp_ptr=work->u; work->u= work->xold; work->xold=swp_ptr;
        // Exploration ended early (time/iteration limit, cycling) => not proven optimal
        return exitflag < DAQP_EXIT_INFEASIBLE ? DAQP_EXIT_FEASIBLE : DAQP_EXIT_OPTIMAL;
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
            // Retain only fixed constraints of common ancestors (robust to pruned nodes)
            if(work->bnb->n_clean > work->bnb->neq+node->depth)
                work->bnb->n_clean = work->bnb->neq+node->depth;
            daqp_node_cleanup_workspace(work->bnb->n_clean,work);
            if(daqp_warmstart_node(node,work) < 0) return DAQP_EXIT_INFEASIBLE;
        }
        else{
            daqp_add_upper_lower(node->bin_id,work);
            work->sense[DAQP_REMOVE_LOWER_FLAG(node->bin_id)] |= DAQP_IMMUTABLE; //Equality
            if(work->sing_ind != DAQP_EMPTY_IND){ // Need to cold start to not miss integer feasible
                if(daqp_setup_cold_bnb(node,work) < 0) return DAQP_EXIT_INFEASIBLE;
            }
        }
        // Add binary constraint
    }
    // Solve relaxation
    exitflag = daqp_ldp(work);
    work->bnb->itercount += work->iterations;

    if(exitflag == DAQP_EXIT_CYCLE){// Try to repair (cold start)
        // A cycle can be caused by stale cached forward-substitution data.
        // Force the retained fixed prefix to be recomputed during repair.
        work->reuse_ind=0;
        if(daqp_setup_cold_bnb(node,work) < 0) return DAQP_EXIT_INFEASIBLE;
        exitflag = daqp_ldp(work);
        work->bnb->itercount += work->iterations;
    }

    return exitflag;
}

int daqp_get_branch_id(DAQPWorkspace* work){
    int i, j;
    int id = DAQP_EMPTY_IND;
    c_float diff, dist, tol, dist_far;
    DAQPNode* node = work->bnb->tree+work->bnb->n_nodes; // Node being processed

    for(i=0; i < work->bnb->nb; i++){
        id = work->bnb->bin_ids[i];
        // Skip fixed binary constraints (also those dropped from the working set since
        // their bounds are implied, which may still be off by rounding)
        if(DAQP_IS_ACTIVE(id)) continue;
        for(j = 0; j <= node->depth; j++)
            if(DAQP_REMOVE_LOWER_FLAG(work->bnb->fixed_ids[j]) == id) break;
        if(j <= node->depth) continue;

        // Compute signed distance from midpoint between bounds
        diff = daqp_binary_diff(id,work);

        // A zero-dual binary constraint can lie at an endpoint without being
        // active. It is already integer feasible and does not need branching.
        dist = 0.5*(work->dupper[id]-work->dlower[id])
            -(diff < 0 ? -diff : diff);
        tol = work->settings->primal_tol;
        if(work->scaling != NULL) tol *= work->scaling[id];
        if(dist <= tol) continue;

        // Bounds for the children (spawned in node/node+1) from KKT at u*:
        // ||u||^2 >= ||u*||^2 + (m_i'(u-u*))^2 since ||m_i|| = 1 (unless soft)
        dist_far = DAQP_IS_SOFT(id) ? 0 : work->dupper[id]-work->dlower[id]-dist-tol;
        dist = DAQP_IS_SOFT(id) ? 0 : dist-tol;
        (node+1)->lb = work->fval + dist*dist; // Nearest endpoint
        node->lb = work->fval + dist_far*dist_far;
        // Explore the endpoint nearest to the relaxation first.
        return diff < 0 ? id : DAQP_ADD_LOWER_FLAG(id);
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
    // Only the retained prefix can still have a valid cached substitution.
    if(work->reuse_ind > n_clean)
        work->reuse_ind=n_clean;
}


int daqp_warmstart_node(DAQPNode* node, DAQPWorkspace* work){
    // Add fixed constraints
    if(daqp_bnb_add_fixed_ids(node,work) < 0) return -1;
    // Add free constraints
    daqp_bnb_load_ws(work->bnb->tree_WS+node->WS_start,node->WS_end-node->WS_start,work);
    work->bnb->nWS = node->WS_start; // always move up tree after warmstart
    return 0;
}

void daqp_save_warmstart(DAQPNode* node, DAQPWorkspace* work){
    node->WS_start = work->bnb->nWS;
    work->bnb->nWS += daqp_bnb_store_ws(work->bnb->tree_WS+work->bnb->nWS,work);
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

int daqp_setup_cold_bnb(DAQPNode* node,DAQPWorkspace *work){
    daqp_node_cleanup_workspace(work->bnb->n_clean,work);
    return daqp_bnb_add_fixed_ids(node,work);
}
