#include "auxiliary.h"
#include "factorization.h"
void daqp_remove_constraint(DAQPWorkspace* work, const int rm_ind){
    int i;
    // Update data structures
    DAQP_SET_INACTIVE(work->WS[rm_ind]); 
    daqp_update_LDL_remove(work,rm_ind);
    (work->n_active)--;

    for(i=rm_ind;i<work->n_active;i++){
        work->WS[i] = work->WS[i+1]; 
        work->lam[i] = work->lam[i+1]; 
    }
    // Can only reuse work less than the ind that was removed 
    if(rm_ind < work->reuse_ind)
        work->reuse_ind = rm_ind;

    // Check if the removal lead to singularity (can happen due to numerics)
    if(work->n_active > 0 && work->D[work->n_active-1] < work->settings->sing_tol){
        work->sing_ind = work->n_active-1;
        work->D[work->n_active-1] = 0;
    }
    else{ // Pivot for improved numerics
        daqp_pivot_last(work);
    }
}
void daqp_add_constraint(DAQPWorkspace *work, const int add_ind, c_float lam){
    // Update data structures  
    DAQP_SET_ACTIVE(add_ind);
#ifdef SOFT_WEIGHTS
    if((DAQP_IS_LOWER(add_ind) && lam <= -work->d_ls[add_ind])||
            (DAQP_IS_LOWER(add_ind)==0 && lam >=  work->d_us[add_ind]))
        DAQP_SET_SLACK_FREE(add_ind);
    else
        DAQP_SET_SLACK_FIXED(add_ind);
#endif
    daqp_update_LDL_add(work, add_ind);
    work->WS[work->n_active] = add_ind;
    work->lam[work->n_active] = lam;
    work->n_active++;

    // Pivot for improved numerics
    daqp_pivot_last(work);
}

void daqp_compute_primal_and_fval(DAQPWorkspace *work){
    int i,j,disp,id;
    c_float fval=0;
    // Reset u & soft slack
    for(j=0;j<work->n;j++)
        work->u[j]=0;
    work->soft_slack = 0;
    //u[m] <-- Mk'*lam_star (zero if empty set)
    for(i=0;i<work->n_active;i++){
        id = work->WS[i];
        c_float lam_i = work->lam_star[i];
        if(id < work->ms){
            // Simple constraint 
            if(work->Rinv!=NULL){ // Hessian is not identity
                for(j=id, disp=DAQP_R_OFFSET(id,work->n);j<work->n;++j)
                    work->u[j]-=work->Rinv[disp+j]*lam_i;
            }
            else work->u[id]-=lam_i; // Hessian is identity
        }
        else{ // General constraint
            for(j=0,disp=work->n*(id-work->ms);j<work->n;j++)
                work->u[j]-=work->M[disp++]*lam_i;
        }
        if(DAQP_IS_SOFT(id)){
#ifdef SOFT_WEIGHTS
            if(DAQP_IS_LOWER(id))
                fval+= lam_i*lam_i*work->rho_ls[id];
            else
                fval+= lam_i*lam_i*work->rho_us[id];
#else
            fval+= lam_i*lam_i;
#endif
        }
    }
    // Check for progress 
#ifndef SOFT_WEIGHTS
    fval=fval*work->settings->rho_soft;
#endif
    work->soft_slack=fval;// XXX: keep this for now to return SOFT_OPTIMAL
    fval += daqp_dot(work->u, work->u, work->n);
    work->fval = fval;
}
int daqp_add_infeasible(DAQPWorkspace *work){
    int j,disp;
    c_float ep = -work->settings->primal_tol;
    c_float min_val = ep;
    c_float Mu,min_cand;
    int isupper=0, add_ind=DAQP_EMPTY_IND;
    // Simple bounds 
    for(j=0, disp=0;j<work->ms;j++){
        // Never activate immutable or already active constraints 
        if(work->sense[j]&(DAQP_ACTIVE+DAQP_IMMUTABLE)){ 
            disp+=work->n-j;
            continue;
        }
        if(work->Rinv==NULL){// Hessian is identify
            Mu=work->u[j]; 
        }
        else{
            Mu = daqp_dot(work->Rinv+disp,work->u+j,work->n-j);
        }
        disp+=work->n-j;
        min_cand = work->dupper[j]-Mu;
        if(min_cand < min_val && (work->scaling == NULL || min_cand < ep*work->scaling[j])){
            add_ind = j; isupper = 1;
            min_val = min_cand;
        }
        else{
            min_cand = Mu - work->dlower[j];
            if(min_cand < min_val && (work->scaling == NULL || min_cand < ep*work->scaling[j])){
                add_ind = j; isupper = 0;
                min_val = min_cand;
            }
        }
    }
    /* General two-sided constraints */
    for(j=work->ms, disp=0;j<work->m;j++){
        // Never activate immutable or already active constraints 
        if(work->sense[j]&(DAQP_ACTIVE+DAQP_IMMUTABLE)){ 
            disp+=work->n;// Skip ahead in M
            continue;
        }
        Mu = daqp_dot(work->M+disp,work->u,work->n);
        disp+=work->n;

        min_cand = work->dupper[j]-Mu;
        if(min_cand < min_val && (work->scaling == NULL || min_cand < ep*work->scaling[j])){
            add_ind = j; isupper = 1;
            min_val = min_cand;
        }
        else{
            min_cand = Mu - work->dlower[j];
            if(min_cand < min_val && (work->scaling == NULL || min_cand < ep*work->scaling[j])){
                add_ind = j; isupper = 0;
                min_val = min_cand;
            }
        }
    }
    // No constraint is infeasible => return
    if(add_ind == DAQP_EMPTY_IND) return 0;
    // Otherwise add infeasible constraint to working set 
    if(isupper)
        DAQP_SET_UPPER(add_ind);
    else
        DAQP_SET_LOWER(add_ind);
    // Set lam = lam_star
    c_float *swp_ptr;
    swp_ptr=work->lam; work->lam = work->lam_star; work->lam_star=swp_ptr;
    // Add the constraint
    if(isupper)
        daqp_add_constraint(work,add_ind,1);
    else
        daqp_add_constraint(work,add_ind,-1);
    return 1;
}
#ifdef SOFT_WEIGHTS
int daqp_remove_blocking(DAQPWorkspace *work){
    int i, ind, rm_ind = DAQP_EMPTY_IND;
    c_float alpha=DAQP_INF;
    c_float alpha_cand, lam_slack;
    const c_float dual_tol = work->settings->dual_tol;
    c_float p;
    for(i=0;i<work->n_active;i++){
        ind = work->WS[i];
        if(DAQP_IS_IMMUTABLE(ind)) continue;
        lam_slack = work->lam[i];
        p = (work->sing_ind == DAQP_EMPTY_IND) ? work->lam_star[i]-work->lam[i] : work->lam_star[i];
        if(DAQP_IS_LOWER(ind)){
            if(DAQP_IS_SLACK_FREE(ind)){ // lam <= -d_ls
                lam_slack += work->d_ls[ind];
                // lam* <= -d_ls for lower -> nothing to do
                if(p < dual_tol || work->lam_star[i] <= -work->d_ls[ind]+dual_tol) continue;
            }
            else{ // slack bound active (implying that -d_ls <= lam <= 0)
                // Remaining within the bound -d_ls <=lam* <= 0 -> nothing to do
                if(work->lam_star[i] <= dual_tol && (work->lam_star[i]+dual_tol >= -work->d_ls[ind]) &&
                        work->sing_ind == DAQP_EMPTY_IND) continue;
                if(p < 0){ // lam* < -d_ls
                    lam_slack += work->d_ls[ind];
                }
            }
        }
        else{ // IS_UPPER
            if(DAQP_IS_SLACK_FREE(ind)){ // lam >= d_us
                lam_slack -= work->d_us[ind];
                //lam* >= d_us for upper -> nothing to do
                if(p > -dual_tol || work->lam_star[i] >= work->d_us[ind]) continue;
            }
            else{ // slack bound active (implying that 0 <= lam <=d_us)
                // Remaining within the bound 0 <=lam* <=d_us -> nothing to do
                if(work->lam_star[i] >= -dual_tol && (work->lam_star[i] <= dual_tol+work->d_us[ind])
                        && work->sing_ind == DAQP_EMPTY_IND) continue;
                if(p > 0) // lam* > d_us
                    lam_slack -= work->d_us[ind];
            }
        }

        alpha_cand = -lam_slack/p;
        if(alpha_cand < alpha){
            alpha = alpha_cand;
            rm_ind = i;
        }
    }
    if(rm_ind == DAQP_EMPTY_IND) return 0; // Either dual feasible or primal infeasible
    // If blocking constraint -> update lambda
    alpha *= 1.001;
    if(work->sing_ind == DAQP_EMPTY_IND)
        for(i=0;i<work->n_active;i++){
            work->lam[i]+=alpha*(work->lam_star[i]-work->lam[i]);
        }
    else
        for(i=0;i<work->n_active;i++){
            work->lam[i]+=alpha*work->lam_star[i];
        }



    // Remove the constraint from the working set and update LDL
    work->sing_ind=DAQP_EMPTY_IND;
    int abs_rm_id = work->WS[rm_ind];
    lam_slack = work->lam[rm_ind];

    daqp_remove_constraint(work,rm_ind);
    if(DAQP_IS_SOFT(abs_rm_id)==0 || work->sing_ind != DAQP_EMPTY_IND) return 1;

    // Reactivate
    if((DAQP_IS_LOWER(abs_rm_id))){
        if(lam_slack  > 0) return 1;
    }
    else{//IS_UPPER
        if(lam_slack  < 0) return 1;
    }
    // Reactive and toggle sense to correctly update factorization
    daqp_add_constraint(work,abs_rm_id,lam_slack);

    return 1;
}
#else // not SOFT_WEIGHTS
int daqp_remove_blocking(DAQPWorkspace *work){
    int i,rm_ind = DAQP_EMPTY_IND; 
    c_float alpha=DAQP_INF;
    c_float alpha_cand;
    const c_float dual_tol = work->settings->dual_tol;
    // sing_ind is not mutated inside the search loop below, so it is safe to cache.
    const int is_normal = (work->sing_ind == DAQP_EMPTY_IND);
    for(i=0;i<work->n_active;i++){
        if(DAQP_IS_IMMUTABLE(work->WS[i])) continue;
        if(DAQP_IS_LOWER(work->WS[i])){
            if(work->lam_star[i]<dual_tol) continue; //lam <= 0 for lower -> dual feasible
        }
        else if(work->lam_star[i]>-dual_tol) continue; //lam* >= 0 for upper-> dual feasible

        alpha_cand = is_normal ? -work->lam[i]/(work->lam_star[i]-work->lam[i])
                               : -work->lam[i]/work->lam_star[i];
        if(alpha_cand < alpha){
            alpha = alpha_cand; 
            rm_ind = i;
        }
    }
    if(rm_ind == DAQP_EMPTY_IND) return 0; // Either dual feasible or primal infeasible
    // If blocking constraint -> update lambda.
    // Rewrite as lam = (1-alpha)*lam + alpha*lam_star to enable vectorization
    // (two independent scale-and-add ops instead of a data-dependent subtract-then-add).
    if(is_normal){
        c_float c = 1.0 - alpha;
        for(i=0;i<work->n_active;i++)
            work->lam[i] = c*work->lam[i] + alpha*work->lam_star[i];
    }
    else
        for(i=0;i<work->n_active;i++)
            work->lam[i]+=alpha*work->lam_star[i];

    // Remove the constraint from the working set and update LDL
    work->sing_ind=DAQP_EMPTY_IND;
    daqp_remove_constraint(work,rm_ind);
    return 1;
}
#endif // SOFT_WEIGHTS

void daqp_compute_CSP(DAQPWorkspace *work){
    int i,j,disp;
    c_float sum;
    // Forward substitution (xi <-- L\d) fused with D-scaling (zi = xi/di)
    for(i=work->reuse_ind,disp=DAQP_ARSUM(work->reuse_ind); i<work->n_active; i++){
        // Setup RHS
        if(DAQP_IS_LOWER(work->WS[i])){
            sum = -work->dlower[work->WS[i]];
#ifdef SOFT_WEIGHTS
            if(DAQP_IS_SOFT(work->WS[i]) && DAQP_IS_SLACK_FREE(work->WS[i])) 
                sum-= work->d_ls[work->WS[i]]*work->rho_ls[work->WS[i]]; 
#endif
        }
        else{
            sum = -work->dupper[work->WS[i]];
#ifdef SOFT_WEIGHTS
            if(DAQP_IS_SOFT(work->WS[i]) && DAQP_IS_SLACK_FREE(work->WS[i])) 
                sum+= work->d_us[work->WS[i]]*work->rho_us[work->WS[i]]; 
#endif
        }
        sum -= daqp_dot(work->L + disp, work->xldl, i);
        disp += i + 1; // advance past i sub-diagonal elements + 1 diagonal skip
        work->xldl[i] = sum;
        work->zldl[i] = sum / work->D[i]; // fuse: scale while D[i] is in cache
    }
    // Backward substitution (lam_star <-- L'\z): column-sweep for sequential L access.
    // Each outer iteration i finalizes lam_star[i] then scatters into lam_star[0..i-1],
    // reading L[ARSUM(i)..ARSUM(i)+i-1] sequentially.
    for(i=0; i<work->n_active; i++)
        work->lam_star[i] = work->zldl[i];
    for(i=work->n_active-1; i>0; i--){
        c_float li = work->lam_star[i];
        disp = DAQP_ARSUM(i);
        for(j=0; j<i; j++)
            work->lam_star[j] -= li * work->L[disp+j];
    }
    work->reuse_ind = work->n_active; // Save forward substitution information 
}

//TODO this could probably be directly calculated in L
void daqp_compute_singular_direction(DAQPWorkspace *work){
    // Step direction is stored in lam_star
    int i,j,disp;
    int offset_L = DAQP_ARSUM(work->sing_ind);

    // Initialize from the sing_ind column of L, then set the sing_ind entry to 1
    for(i=0; i<work->sing_ind; i++)
        work->lam_star[i] = -work->L[offset_L+i];
    work->lam_star[work->sing_ind] = 1;

    // Backwards substitution (p_tidle <-- L'\(-l)): column-sweep for sequential L access.
    for(i=work->sing_ind-1; i>0; i--){
        c_float li = work->lam_star[i];
        disp = DAQP_ARSUM(i);
        for(j=0; j<i; j++)
            work->lam_star[j] -= li * work->L[disp+j];
    }

    if(DAQP_IS_LOWER(work->WS[work->sing_ind])) //Flip to ensure descent direction 
        for(i=0;i<=work->sing_ind;i++)
            work->lam_star[i] =-work->lam_star[i];
}


void daqp_pivot_last(DAQPWorkspace *work){
    const int rm_ind = work->n_active-2; 
    if(work->n_active > 1 && 
            work->D[rm_ind] < work->settings->pivot_tol && // element in D small enough
            work->D[rm_ind] < work->D[work->n_active-1]){ // element in D smallar than neighbor
        const int ind_old = work->WS[rm_ind];
        // Ensure that binaries never swap order (since this order is exploited) 
        if(DAQP_IS_BINARY(ind_old) && DAQP_IS_BINARY(work->WS[work->n_active-1])) return; 
        if(work->bnb != NULL && rm_ind < work->bnb->n_clean) return;

        c_float lam_old = work->lam[rm_ind];
        daqp_remove_constraint(work,rm_ind); // pivot_last might be recursively called here 

        if(work->sing_ind!=DAQP_EMPTY_IND) return; // Abort if D becomes singular

        daqp_add_constraint(work,ind_old,lam_old);
    }	
}

// Activate constrainte that are marked active in sense
int daqp_activate_constraints(DAQPWorkspace *work){
    //TODO prioritize inequalities?
    int i;
    for(i =0;i<work->m;i++){
        if(DAQP_IS_ACTIVE(i)){
#ifdef SOFT_WEIGHTS
            if(DAQP_IS_LOWER(i)){
                if(DAQP_IS_SLACK_FREE(i))
                    daqp_add_constraint(work,i, -(work->d_ls[i]+1));
                else
                    daqp_add_constraint(work,i, -0.9*work->d_ls[i]);
            }
            else{ //IS Upper
                if(DAQP_IS_SLACK_FREE(i))
                    daqp_add_constraint(work,i, work->d_us[i]+1);
                else
                    daqp_add_constraint(work,i, 0.9*work->d_us[i]);
            }
#else
            if(DAQP_IS_LOWER(i))
                daqp_add_constraint(work,i, -1.0);
            else
                daqp_add_constraint(work,i, 1.0);
#endif
        }
        if(work->sing_ind != DAQP_EMPTY_IND){
            int exitflag = 1;
            for(;i<work->m;i++){ 
                // 1. Check if there are equalities that couldn't be activated
                // 2. Make sure that sense is clean for unactivated constraints 
                if(DAQP_IS_ACTIVE(i)){
                    if(DAQP_IS_IMMUTABLE(i)) 
                        exitflag = DAQP_EXIT_OVERDETERMINED_INITIAL;
                    else
                        DAQP_SET_INACTIVE(i);
                }
            }
            // Remove the last constraint that lead to singularity
            work->n_active--;
            work->sing_ind = DAQP_EMPTY_IND;
            return exitflag;
        }
    }
    return 1;
}

// Deactivate all active constraints that are mutable (i.e., not equality constraints)
void daqp_deactivate_constraints(DAQPWorkspace *work){
    int i;
    for(i =0;i<work->n_active;i++){
        if(DAQP_IS_IMMUTABLE(work->WS[i])) continue; 
        DAQP_SET_INACTIVE(work->WS[i]); 
    }
}
