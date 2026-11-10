#include "avi.h"
#include "daqp.h"
#include "utils.h"
#include "auxiliary.h"

// Max pivots to repair a warm-started working set before the DR iterations
#ifndef DAQP_AVI_MAX_REPAIR
#define DAQP_AVI_MAX_REPAIR 4
#endif

// Solve AVI
int daqp_solve_avi(DAQPWorkspace *work) {
    DAQPAVI* avi = work->avi;
    int n = work->n;
    int i,j,k,disp;
    int exitflag = -10;
    c_float val,sum,sum2;
    int tot_iter = 0;
    int counter = 0;
    int terminate_limit = 2;
    int retry_requested = 0;
    c_float minimum_newton_residual = DAQP_INF;

    work->outer_iter = 0;
    avi->Hinv_f_valid = 0;

    // Shortcut: return directly if the (warm-started) working set is optimal
    // (not if it has a dependent constraint, which daqp_ldp drops first)
    if(work->sing_ind == DAQP_EMPTY_IND && work->n_active <= n){
        for(i = 0; i < work->n_active; i++) avi->lam_buffer[i] = work->lam_star[i];
        // Otherwise, try to repair a warm-started working set by pivoting on violated constraints
        int optimal, pivots = 0, max_pivots = work->n_active > 0 ? DAQP_AVI_MAX_REPAIR : 0;
        daqp_solve_avi_kkt(work);
        while(!(optimal = daqp_check_optimal_avi(work)) && pivots++ < max_pivots){
            if(avi->viol < 0) daqp_remove_constraint(work, -avi->viol-1);
            else{
                if(avi->viol & 1) DAQP_SET_UPPER(avi->viol >> 1); else DAQP_SET_LOWER(avi->viol >> 1);
                daqp_add_constraint(work, avi->viol >> 1, (avi->viol & 1) ? 1.0 : -1.0);
            }
            if(work->sing_ind != DAQP_EMPTY_IND) break; // (daqp_ldp drops the dependent constraint)
            daqp_solve_avi_kkt(work);
        }
        if(optimal){
            for(i = 0; i < n; i++) work->x[i] = avi->x[i];
            work->iterations = 1 + pivots;
            return 1;
        }
        for(i = 0; i < work->n_active; i++) work->lam_star[i] = avi->lam_buffer[i];
    }

    // Initial avi iterate
    for(i=0; i < work->n; i++) work->avi->x[i] = work->x[i];

    // Start the iterations
    // TODO iter_limit should be the for tot_iter...
    for (k = 0; k < work->settings->iter_limit; k++) {
        work->outer_iter++;
        // Compute xtemp = H*x + f - (Hsym + I)x
        for(i=0, disp=0; i < work->n; i++){
            sum = sum2 = 0.0;
            for(j=0; j < work->n; j++) {
                sum += work->qp->H[disp]*avi->x[j];
                sum2 += avi->Hs_rho[disp++]*avi->x[j];
            }
            avi->Hx[i] = sum;
            avi->xtemp[i] = sum + work->qp->f[i] - sum2;
        }

        // Update linear term
        daqp_update_v(avi->xtemp,work);
        daqp_update_d(work, work->qp->bupper,work->qp->blower);

        exitflag = daqp_ldp(work);

        if (exitflag < 0) break;
        ldp2qp_solution(work);
        tot_iter += work->iterations;

        if (counter == terminate_limit){ // Check if Newton step made progress
            // Compute ||natural residual||^2
            sum = 0.0;
            for(i=0; i < n; i++){
                val = work->avi->x[i]-work->x[i];
                sum += val*val;
            }
            // If no decrease since last Newton iterate -> revert Newton step
            if(sum > minimum_newton_residual){
                for(i = 0; i < n; i++) work->avi->x[i] = work->xold[i];
                if(terminate_limit == 30 && avi->retry_rho_needed){
                    retry_requested = 1;
                    break;
                }
                terminate_limit += 5; // Increase terminate limit to give DR more time to converge
                if(terminate_limit > 30) terminate_limit = 30;
            }
            else{
                minimum_newton_residual = sum;
                for(i=0; i < n; i++) avi->y[i]=work->x[i];
            }
        }
        else{ // Update y iterate
            for(i=0; i < n; i++) avi->y[i]=work->x[i];
        }

        // AS has not changed -> Check KKT conditions
        if (work->iterations == 1) {
            if (++counter == terminate_limit) {
                for(i=0; i < n; i++) work->xold[i] = work->avi->x[i]; // Store xold in case Newton step fails
                // Find KKT point
                daqp_solve_avi_kkt(work);
                if (daqp_check_optimal_avi(work)) {
                    for(i=0; i < n; i++) work->x[i] = work->avi->x[i]; // TODO no need for local x in avi
                    exitflag = 1;
                    break;
                }
                continue;
            }
        }
        else {
            counter = 0;
        }

        for (i = 0; i < work->n; i++){
            avi->xtemp[i] = avi->rho*avi->y[i]+avi->Hx[i];
            avi->y[i] -= avi->x[i];
        }
        for (i = 0; i < n; i++) {
            const c_float* Hsym_i = avi->Hsym + i*n;
            for (j = 0, sum = 0.0; j < n; j++) sum += Hsym_i[j]*avi->y[j];
            avi->xtemp[i] += 0.5*sum;
        }
        daqp_lu_solve(avi->H_rho, avi->P_H2, avi->xtemp, avi->x, n);
    }
    if(retry_requested){
        int original_limit = work->settings->iter_limit;
        int previous_outer_iterations = work->outer_iter;
        int retry_flag = daqp_retry_avi_with_reduced_rho(work);
        if(retry_flag < 0) return retry_flag;
        if(retry_flag > 0 && k+1 < original_limit){
            work->settings->iter_limit = original_limit-(k+1);
            retry_flag = daqp_solve_avi(work);
            work->settings->iter_limit = original_limit;
            work->iterations += tot_iter;
            work->outer_iter += previous_outer_iterations;
            return retry_flag;
        }
        work->iterations = tot_iter;
        return DAQP_EXIT_ITERLIMIT;
    }
    if(k==work->settings->iter_limit) exitflag = -4;
    work->iterations = tot_iter;
    return exitflag;
}

// H^{-1}a_i (a_i = e_i for simple bounds), computed when first needed
static const c_float* daqp_avi_hinv_row(DAQPWorkspace* work, const int row){
    DAQPAVI* avi = work->avi;
    const int n = work->n;
    c_float* Yi = avi->HinvAT + (size_t)row*n;
    if(!avi->HinvAT_valid[row]){
        if(row < work->ms){ // Simple bound
            int j;
            for(j = 0; j < n; j++) avi->unit[j] = 0;
            avi->unit[row] = 1.0;
            daqp_lu_solve(avi->LU_H, avi->P_H, avi->unit, Yi, n);
        }
        else
            daqp_lu_solve(avi->LU_H, avi->P_H, work->qp->A + (size_t)(row-work->ms)*n, Yi, n);
        avi->HinvAT_valid[row] = 1;
    }
    return Yi;
}

void daqp_solve_avi_kkt(DAQPWorkspace* work) {
    DAQPAVI* avi = work->avi;
    int i, j, k, disp,row_idx;
    c_float sum;
    int nAS = work->n_active;
    const int n = work->n;
    c_float* lambda = work->lam_star;
    c_float* x = avi->x;
    c_float* S = avi->kkt_buffer;
    c_float* rhs = avi->kkt_buffer+nAS*nAS;
    c_float* temp = avi->Hinv_f;
    int* WS = work->WS;

    // Form and factorize S = A_WS * H^-1 * A_WS^T (unless cached for this WS)
    int reuse = (avi->kkt_nAS == nAS);
    for (i = 0; reuse && i < nAS; i++) reuse = (avi->kkt_WS[i] == WS[i]);
    if(!reuse){
        int has_soft = 0;
        // Reuse entries of S for constraints that were in the previous working set
        int* pos = avi->kkt_pos;
        for (i = 0; i < nAS; i++)
            for (pos[i] = -1, j = 0; j < avi->kkt_nAS; j++)
                if(avi->kkt_WS[j] == WS[i]){ pos[i] = j; break; }
        for (i = 0; i < nAS; i++) {
            const c_float* Yi = NULL;
            for (j = 0; j < nAS; j++) {
                if(pos[i] >= 0 && pos[j] >= 0)
                    sum = avi->kkt_raw_old[pos[j]*n + pos[i]];
                else{
                    if(Yi == NULL) Yi = daqp_avi_hinv_row(work, WS[i]);
                    row_idx = WS[j];
                    if(row_idx < work->ms) // Simple bound
                        sum = Yi[row_idx];
                    else{ // General constraint
                        const c_float* a_j = work->qp->A + (size_t)(row_idx-work->ms)*n;
                        for (k = 0, sum = 0.0; k < n; k++) sum += a_j[k]*Yi[k];
                    }
                }
                avi->kkt_raw[j*n + i] = sum;
                S[j * nAS + i] = sum;
            }
            // Soft constraints -> diagonal of S gets regularized
            if(DAQP_IS_SOFT(WS[i])){
                S[i*(nAS+1)] += work->settings->rho_soft/(work->scaling[WS[i]] * work->scaling[WS[i]]);
                has_soft = 1;
            }
        }
        c_float* swp = avi->kkt_raw; avi->kkt_raw = avi->kkt_raw_old; avi->kkt_raw_old = swp;
        daqp_lu(S, avi->P_S, nAS);
        for (i = 0; i < nAS; i++) avi->kkt_WS[i] = WS[i];
        avi->kkt_nAS = has_soft ? -1 : nAS; // (soft weights may change)
    }

    // Form rhs for Schur system: rhs = -A_WS * H^-1 * f - b_WS
    if(!avi->Hinv_f_valid){
        daqp_lu_solve(avi->LU_H, avi->P_H, work->qp->f, temp, n);
        avi->Hinv_f_valid = 1;
    }

    for (i = 0; i < nAS; i++) {
        row_idx = WS[i];
        sum = DAQP_IS_LOWER(row_idx) ?
            work->qp->blower[row_idx] : work->qp->bupper[row_idx];
        if(row_idx < work->ms) // Simple bound
            sum += temp[row_idx];
        else{ // General constraint
            disp = (row_idx-work->ms)*n;
            for (k = 0; k < n; k++) {
                sum += work->qp->A[disp++] * temp[k];
            }
        }
        rhs[i] = -sum;
    }

    // Get lambda by solving S * lambda = rhs
    daqp_lu_solve(S, avi->P_S, rhs, lambda, nAS);

    // x = -H^-1 f - sum_j lambda_j H^-1 a_j
    for (i = 0; i < n; i++) x[i] = -temp[i];
    for (j = 0; j < nAS; j++) {
        const c_float* Yj = daqp_avi_hinv_row(work, WS[j]);
        const c_float lj = lambda[j];
        for (i = 0; i < n; i++) x[i] -= lj*Yj[i];
    }
}


int daqp_check_optimal_avi(DAQPWorkspace* work){
    int i,j,disp;
    c_float dual_tol = work->settings->dual_tol;
    c_float primal_tol = work->settings->primal_tol;
    // First check dual variables
    for(i=0; i < work->n_active; i++){
        if(DAQP_IS_IMMUTABLE(work->WS[i])) continue;
        if(DAQP_IS_LOWER(work->WS[i])){
            if(work->lam_star[i] > dual_tol){ work->avi->viol = -1-i; return 0; }
        }
        else{
            if(work->lam_star[i] < -dual_tol){ work->avi->viol = -1-i; return 0; }
        }
    }
    // Simple constraints
    for(i=0; i < work->ms; i++){
        if(DAQP_IS_ACTIVE(i)) continue;
        if(work->avi->x[i] > work->qp->bupper[i] + primal_tol){ work->avi->viol = 2*i+1; return 0; }
        if(work->avi->x[i] < work->qp->blower[i] - primal_tol){ work->avi->viol = 2*i; return 0; }
    }

    // General constraints
    for(disp=0; i < work->m; i++){
        if(DAQP_IS_ACTIVE(i)){
            disp += work->n;
            continue;
        }
        c_float Ax = 0.0;
        for(j=0; j < work->n; j++) Ax += work->qp->A[disp++]*work->avi->x[j];
        if(Ax > work->qp->bupper[i] + primal_tol){ work->avi->viol = 2*i+1; return 0; }
        if(Ax < work->qp->blower[i] - primal_tol){ work->avi->viol = 2*i; return 0; }
    }
    // All checks passed -> optimal KKT point found
    return 1;
}
