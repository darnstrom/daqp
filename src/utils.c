#include "daqp.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

int daqp_update_ldp(const int mask, DAQPWorkspace *work, DAQPProblem* qp){
    // TODO: copy dimensions from work->qp? 
    int error_flag, i;
    int do_activate = 0;

    // Add original qp to workspace
    work->qp = qp;

    // Update dimensions of problem
    work->n = qp->n;
    work->m = qp->m;
    work->ms = qp->ms;

    // Reset unconstrained_optimal flag (re-evaluated below whenever relevant data changes)
    work->unconstrained_optimal = 0;

    /** Update constraint sense **/
    if(mask&DAQP_UPDATE_sense){
        if(work->qp->sense == NULL) // Assume all constraints are "normal" inequality constraints
            for(i=0;i<work->m;i++) work->sense[i] = 0;
        else{
            for(i=0;i<work->m;i++) work->sense[i] = qp->sense[i];
            do_activate = 1;
        }
    }

    /** Update Rinv **/
    if(mask&DAQP_UPDATE_Rinv){
        if(work->avi == NULL)
            error_flag = daqp_update_Rinv(work, qp->H, qp->problem_type==2 ? 1 : 0);
        else{
            daqp_update_avi(work->avi,qp);
            error_flag = daqp_update_Rinv(work, work->avi->Hs_rho,0);
        }
        if(error_flag<0)
            return error_flag;
    }

    /** Update v (moved before M to enable early-exit check below) **/
    if(mask&DAQP_UPDATE_Rinv||mask&DAQP_UPDATE_v){
        daqp_update_v(qp->f,work,mask);
    }

    /** Check if unconstrained optimum is primal feasible.
     *  Compute x_unc = Rinv * (-v) (using the raw, un-normalized Rinv) and
     *  verify dupper_unnorm = bupper - A*x_unc >= 0 and
     *         dlower_unnorm = blower - A*x_unc <= 0 for every constraint.
     *  If so, x_unc is optimal and we can skip the expensive M = A*Rinv step.
     *  Only applicable for standard QPs (no BnB, no hierarchy, no AVI, no prox). **/
    if((mask&DAQP_UPDATE_Rinv||mask&DAQP_UPDATE_M||
        mask&DAQP_UPDATE_v||mask&DAQP_UPDATE_d) &&
       work->bnb == NULL && work->nh <= 1 && work->avi == NULL){
        int j, disp;
        c_float sum;
        int feasible = 1;
        const int n = work->n;
        const c_float primal_tol = work->settings->primal_tol;

        // Compute x_unc = Rinv_unnorm * (-v) stored temporarily in work->u.
        // work->u is calloc'd to 0; after this block it is always reset to 0.
        if(work->v != NULL){
            if(work->Rinv != NULL){
                // Upper-triangular back-substitution: u[i] = sum_{j>=i} Rinv[i,j]*(-v[j])
                for(i = 0, disp = 0; i < n; i++){
                    for(j = i, sum = 0; j < n; j++)
                        sum += work->Rinv[disp++] * work->v[j];
                    work->u[i] = -sum;
                }
            } else if(work->RinvD != NULL){
                for(i = 0; i < n; i++)
                    work->u[i] = -work->RinvD[i] * work->v[i];
            } else {
                for(i = 0; i < n; i++)
                    work->u[i] = -work->v[i];
            }
        }
        // If v == NULL (f == 0), x_unc = 0 and work->u is already 0.

        // Check simple bounds: blower[i] <= x_unc[i] <= bupper[i]
        for(i = 0; i < work->ms && feasible; i++){
            if(work->u[i] > qp->bupper[i] + primal_tol ||
               work->u[i] < qp->blower[i] - primal_tol)
                feasible = 0;
        }

        // Check general constraints: blower[i] <= A[i,:]*x_unc <= bupper[i]
        if(feasible){
            for(i = work->ms, disp = 0; i < work->m && feasible; i++){
                for(j = 0, sum = 0; j < n; j++)
                    sum += qp->A[disp++] * work->u[j];
                if(sum > qp->bupper[i] + primal_tol ||
                   sum < qp->blower[i] - primal_tol)
                    feasible = 0;
            }
        }

        // Reset u back to 0 (ldp2qp_solution needs u = 0 for the unconstrained solution)
        for(i = 0; i < n; i++) work->u[i] = 0;

        if(feasible){
            work->unconstrained_optimal = 1;
            // Normalize Rinv so that ldp2qp_solution recovers the correct x = Rinv*(-v)
            if(mask&DAQP_UPDATE_Rinv)
                daqp_normalize_Rinv(work);
            // Skip M and d computation and return early
            return do_activate;
        }
    }

    /** Update M **/
    if(mask&DAQP_UPDATE_Rinv||mask&DAQP_UPDATE_M){
        error_flag = daqp_update_M(work,qp->A,mask);
        if(error_flag<0)
            return error_flag;
    }

    // Normalize Rinv
    if(mask&DAQP_UPDATE_Rinv){
        daqp_normalize_Rinv(work);
    }

    /** Update d **/
    if(mask&DAQP_UPDATE_Rinv||mask&DAQP_UPDATE_M||mask&DAQP_UPDATE_v||mask&DAQP_UPDATE_d){
        error_flag = daqp_update_d(work,qp->bupper,qp->blower);
        if(error_flag<0) return error_flag;
        if(error_flag==1) do_activate = 1;
    }

#ifdef SOFT_WEIGHTS
    // TODO: Use mask or something to avoid scaling something more times... 
    if(work->d_ls != NULL && work->scaling !=NULL){
        for(i=0;i<work->m; i++){
            work->d_ls[i]*=work->scaling[i];
            work->d_us[i]*=work->scaling[i];
            work->rho_ls[i]/=work->scaling[i]*work->scaling[i];
            work->rho_us[i]/=work->scaling[i]*work->scaling[i];
        }
    }
#endif

    /** Update hierarchy **/
    if(mask&DAQP_UPDATE_hierarchy){
        work->nh = qp->nh;
        work->break_points = qp->break_points;
    }

    // Make sure activate constraints are activated
    if(do_activate == 1){
        reset_daqp_workspace(work);
        if(work->nh < 2)
            error_flag = daqp_activate_constraints(work);
        else{// Activate the first level (since those constraints are hard)
            int m_tmp = work->m;
            work->m = work->break_points[0];
            error_flag = daqp_activate_constraints(work);
            work->m = m_tmp;
        }
        if(error_flag<0)
            return error_flag;
    }

    return 0;
}

int daqp_update_Rinv(DAQPWorkspace *work, c_float* H, int is_factored){
    int i, j, k, disp, disp2, disp3;
    const int n = work->n; 
    c_float eps = work->settings->eps_prox;
    c_float zero_tol = work->settings->zero_tol;

    // Reset the semi-proximal mask for this factorization
    if(work->prox_mask != NULL){
        for(i = 0; i < n; i++) work->prox_mask[i] = 0;
    }
    work->n_prox = 0;

    // Check if Diagonal
    int is_diagonal = 1;
    for (i=0, disp=1; i<n; i++, disp += (is_factored ? 1 : (i+1))){
        for (j=1; j<n-i; j++, disp++) {
            if(H[disp] > 1e-12 || H[disp] < -1e-12){
                is_diagonal = 0; break;
            }
        }
        if(!is_diagonal) break;
    }

    // Diagonal Case
    if(is_diagonal){
        if(work->Rinv != NULL){ work->RinvD = work->Rinv; work->Rinv = NULL; }

        disp = 0;
        for(i=0; i<n; i++){
            c_float Hi = H[disp];
            // If factored, skip eps and sqrt, else apply them
            if(!is_factored){
                // Only add eps if this direction would be singular without it
                if(Hi <= zero_tol){
                    if(work->prox_mask != NULL) work->prox_mask[i] = 1;
                    work->n_prox++;
                    Hi += eps;
                }
                if (Hi <= zero_tol) return DAQP_EXIT_NONCONVEX;
                Hi = sqrt(Hi);
                disp += n+1;
            } else {
                if(Hi <= zero_tol){
                    if(work->prox_mask != NULL) work->prox_mask[i] = 1;
                    work->n_prox++;
                    Hi = sqrt(Hi*Hi + eps); // Regularization for factors
                }
                disp += n-i;
            }

            work->RinvD[i] = 1/Hi;
            if(work->scaling != NULL && i < work->ms) work->scaling[i] = Hi;
        }
        return 1;
    }

    // Prepare Rinv for dense data
    if(work->RinvD != NULL){ work->Rinv = work->RinvD; work->RinvD = NULL; }

    // Cholesky 
    if (is_factored) {
        for(i=0, disp=0; i<n; i++){
            work->Rinv[disp] = 1/H[disp]; // Store 1/rii
            for (j=1, disp++; j<n-i; j++, disp++) 
                work->Rinv[disp] = H[disp];
        }
    } else {
        // Standard Cholesky (H -> R), adding eps only where the diagonal
        // of the factor would otherwise be non-positive (semi-proximal).
        for (i=0, disp=0, disp3=0; i<n; disp+=n-i, i++, disp3+=i) {
            // Accumulate diagonal contribution without eps first
            c_float diag_i = H[disp3++];
            for (k=0, disp2=i; k<i; k++, disp2+=n-k)
                diag_i -= work->Rinv[disp2] * work->Rinv[disp2];

            // Only regularize if this direction is singular
            if(diag_i <= zero_tol){
                if(work->prox_mask != NULL) work->prox_mask[i] = 1;
                work->n_prox++;
                diag_i += eps;
            }

            if (diag_i <= zero_tol) return DAQP_EXIT_NONCONVEX;
            work->Rinv[disp] = sqrt(diag_i);

            for (j=1; j<n-i; j++) {
                work->Rinv[disp+j] = H[disp3++];
                for (k=0, disp2=i; k<i; k++, disp2+=n-k)
                    work->Rinv[disp+j] -= work->Rinv[disp2] * work->Rinv[disp2+j];
                work->Rinv[disp+j] /= work->Rinv[disp];
            }
            work->Rinv[disp] = 1/work->Rinv[disp];
        }
    }

    // R -> Rinv
    for(k=0, disp=0; k<n; k++){
        disp2 = disp;
        work->Rinv[disp] = work->Rinv[disp2++];
        for(j=k+1; j<n; j++) work->Rinv[disp2++] *= -work->Rinv[disp];
        for(i=k+1, disp++; i<n; i++, disp++){
            work->Rinv[disp] *= work->Rinv[disp2++];
            for(j=1; j<n-i; j++) work->Rinv[disp+j] -= work->Rinv[disp2++] * work->Rinv[disp];
        }
    }
    return 1;
}

int daqp_update_M(DAQPWorkspace *work, c_float *A, const int mask){
    int i,j,k,disp,disp2;
    const int n = work->n;
    const int mA = work->m-work->ms;
    int stop_id =  (mask & DAQP_UPDATE_Rinv) ? n : n-work->ms;
    if(work->Rinv != NULL){
        for(k = 0,disp2=n*mA-1;k<mA;k++,disp2-=n){
            disp=DAQP_ARSUM(n);
            for(j = 0; j< stop_id ; ++j){
                for(i=0;i<j;++i)
                    work->M[disp2-i] += work->Rinv[--disp]*A[disp2-j];
                work->M[disp2-j]=work->Rinv[--disp]*A[disp2-j];
            }
            for(; j<n; ++j){// Take into account scaling in Rinv 
                for(i=0;i<j;++i)
                    work->M[disp2-i] += (work->Rinv[--disp]/work->scaling[n-j-1])*A[disp2-j];
                work->M[disp2-j]=(work->Rinv[--disp]/work->scaling[n-j-1])*A[disp2-j];
            }
        }
    }
    else{
        if(work->RinvD == NULL){ // Copy A to M 
            for(k = 0,disp=0;k<mA;k++){
                for(i=0;i<n;i++,disp++)
                    work->M[disp] = A[disp];
            }
        }
        else{
            for(k = 0,disp=0;k<mA;k++){
                for(i=0;i<n;i++,disp++)
                    work->M[disp] = A[disp]*work->RinvD[i];
            }
        }
    }

    reset_daqp_workspace(work); // Internal factorizations need to be redone!
    return daqp_normalize_M(work);
}

void daqp_update_v(c_float *f, DAQPWorkspace *work, const int mask){
    int i,j,disp;
    const int n = work->n;
    if(work->v == NULL || f == NULL) return;
    if(work->Rinv == NULL){// Rinv = I => v = R'\v = f
        if(work->RinvD != NULL)
            for(i=0;i<n;++i) work->v[i] = f[i]*work->RinvD[i];
        else
            for(i=0;i<n;++i) work->v[i] = f[i];
        return;
    }
    int stop_id =  (mask & DAQP_UPDATE_Rinv) ? 0 : work->ms;
    for(j=n-1,disp=DAQP_ARSUM(n);j>=stop_id;j--){
        for(i=n-1;i>j;i--)
            work->v[i] +=work->Rinv[--disp]*f[j];
        work->v[j]=work->Rinv[--disp]*f[j];
    }
    for(;j>=0;j--){// Take into accoutn scaling in Rinv
        for(i=n-1;i>j;i--)
            work->v[i] +=(work->Rinv[--disp]/work->scaling[j])*f[j];
        work->v[j]=(work->Rinv[--disp]/work->scaling[j])*f[j];
    }
}

int daqp_update_d(DAQPWorkspace *work, c_float *bupper, c_float *blower){
    /* Compute d  = b+M*v */
    int i,j,disp;
    int do_activate = 0;
    c_float sum;

#ifndef DAQP_ASSUME_VALID
    c_float diff;
    for(i =0;i<work->m;i++){
        if(DAQP_IS_IMMUTABLE(i)) continue;
        diff = bupper[i] - blower[i];
        // Check for trivial infeasibility
        if ( diff < -work->settings->primal_tol ){
            return DAQP_EXIT_INFEASIBLE;
        }
        // Check for unmarked equality constraint (blower == bupper)
        else if ( diff < work->settings->zero_tol ){
            work->sense[i] |= DAQP_ACTIVE + DAQP_IMMUTABLE;
            do_activate = 1;
        }
        // TODO: Make innactive here
    }
#endif

    const int n = work->n;
    work->reuse_ind = 0; // RHS of KKT system changed => cannot reuse intermediate results
    // Take into scaling of constraints
    if(work->scaling != NULL){
        for(i = 0;i<work->m;i++){
            work->dupper[i] = bupper[i]*work->scaling[i];
            work->dlower[i] = blower[i]*work->scaling[i];
        }
    }
    else{
        for(i = 0;i<work->m;i++){
            work->dupper[i] = bupper[i];
            work->dlower[i] = blower[i];
        }
    }

    if(work->v == NULL) return do_activate;
    // Simple bounds 
    if(work->Rinv !=NULL){
        for(i = 0,disp=0;i<work->ms;i++){
            for(j=i, sum=0;j<n;j++)
                sum+=work->Rinv[disp++]*work->v[j];
            work->dupper[i]+=sum;
            work->dlower[i]+=sum;
        }
    }else{
        for(i = 0,disp=0;i<work->ms;i++){
            work->dupper[i]+=work->v[i];
            work->dlower[i]+=work->v[i];
        }
    }
    //General bounds
    for(i = work->ms, disp=0;i<work->m;i++){
        for(j=0, sum=0;j<n;j++)
            sum+=work->M[disp++]*work->v[j];
        work->dupper[i]+=sum;
        work->dlower[i]+=sum;
    }
    return do_activate;
}

void daqp_normalize_Rinv(DAQPWorkspace* work){
    int i,j,disp;
    c_float scaling_i;
    // Normalize simple constraints
    if(work->Rinv !=NULL){
        for(i=0, disp=0; i < work->ms;i++){
            scaling_i = 0;
            for(j=i; j < work->n; j++,disp++){
                scaling_i+=work->Rinv[disp]*work->Rinv[disp];
            }
            scaling_i = 1/sqrt(scaling_i);
            work->scaling[i] = scaling_i; // Need to save to correctly retrieve solution
            for(j=i,disp-=(work->n-i); j < work->n; j++,disp++)
                work->Rinv[disp]*= scaling_i;
        }
    }
}
int daqp_normalize_M(DAQPWorkspace* work){
    int i,j,disp;
    c_float scaling_i;
    c_float zero_tol = work->settings->zero_tol;
    // Normalize general constraints 
    for(i=work->ms, disp=0;i<work->m;i++){
        scaling_i = 0;
        for(j=0;j<work->n;disp++,j++)
            scaling_i+=work->M[disp]*work->M[disp];
        if(scaling_i < zero_tol){
#ifndef DAQP_ASSUME_VALID
            if(work->qp->bupper[i] < -zero_tol || work->qp->blower[i] > zero_tol)
                if(work->sense[i] != DAQP_IMMUTABLE)
                    return DAQP_EXIT_INFEASIBLE;
#endif
            work->sense[i] = DAQP_IMMUTABLE; // ignore zero-row constraint
            continue; // TODO: mark infeasibility if dupper & dlower are nonzero
        }
        scaling_i = 1/sqrt(scaling_i);
        work->scaling[i]=scaling_i;
        for(j=0, disp-=work->n;j<work->n;j++,disp++)
            work->M[disp]*=scaling_i;
    }
    return 0;
}

int daqp_update_avi(DAQPAVI* avi, DAQPProblem* p){
    const int n = p->n;
    // Setup matrices Hsym, Hs_rho, and H_rho, LU_H
    int i,j,disp;
    c_float val;
    avi->rho = 0.0;
    for (i = 0, disp=0; i < n; i++) {
        for (j = 0; j < n; j++, disp++) {
            val = (p->H[disp] + p->H[j * n + i]) * 0.5;
            avi->Hsym[disp] = val;
            avi->Hs_rho[disp] = val;
            avi->H_rho[disp] = p->H[disp];
            avi->LU_H[disp] = p->H[disp];
            avi->rho += p->H[disp] * p->H[disp];
        }
    }
    // Regularization 
    avi->rho = sqrt(avi->rho)/2;
    for(i=0,disp=0; i<n;i++){
        avi->Hs_rho[disp] += avi->rho;
        avi->H_rho[disp] += avi->rho;
        disp += n+1;
    }

    // Factorize H and H_rho 
    daqp_lu(avi->LU_H, avi->P_H, n);
    daqp_lu(avi->H_rho, avi->P_H2, n);
    return 1;
}

int daqp_lu(c_float* A, int* P, int n) {
    c_float max_val,pA;
    for (int i = 0; i < n; i++) P[i] = i; // Initialize permutation vector
    for (int i = 0; i < n; i++) {
        // Pivot
        max_val = 0.0;
        int pivot = i;
        for (int j = i; j < n; j++) {
            pA = A[j*n+i];
            pA = (pA < 0) ? -pA : pA; // |pA| 
                if (pA > max_val) {
                max_val = pA;
                pivot = j;
            }
        }

        // Check for singularity
        if (max_val < 1e-12) return -1;

        // Swap rows in A
        for (int k = 0; k < n; k++) {
            c_float temp = A[i * n + k];
            A[i * n + k] = A[pivot * n + k];
            A[pivot * n + k] = temp;
        }
        // Swap elements in permutation vector
        int tempP = P[i];
        P[i] = P[pivot];
        P[pivot] = tempP;

        // Elimination
        for (int j = i + 1; j < n; j++) {
            A[j * n + i] /= A[i * n + i]; // Store multiplier in L part
            for (int k = i + 1; k < n; k++) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }
    return 0;
}

void daqp_lu_solve(c_float* LU, int* P, c_float* b, c_float* x, int n) {
    // Solve Ly = Pb
    for (int i = 0; i < n; i++) {
        x[i] = b[P[i]]; // Apply permutation to b
        for (int j = 0; j < i; j++) {
            x[i] -= LU[i * n + j] * x[j];
        }
    }
    // Solve Ux = y
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            x[i] -= LU[i * n + j] * x[j];
        }
        x[i] /= LU[i * n + i];
    }
}

/* Remove Minrep */
void daqp_minrep_work(int* is_redundant, DAQPWorkspace* work){
    int i,j,exitflag;

    for(i=0; i < work->m; i++)
        is_redundant[i] = -1;

    for(i=0; i < work->m; i++){
        if(is_redundant[i] != -1 || DAQP_IS_IMMUTABLE(i)) continue;
        reset_daqp_workspace(work);
        work->sense[i] = 5;
        daqp_add_constraint(work,i,1.0);
        //work->dupper[i] += tol_weak; TODO support weaky infeasible constraints
        exitflag = daqp_ldp(work);
        if(exitflag== DAQP_EXIT_INFEASIBLE){
            is_redundant[i] = 1;
            work->sense[i] &=~DAQP_ACTIVE; // deactive (remains immutable -> ignored)
        }
        else{
            is_redundant[i] = 0;
            work->sense[i] &=~DAQP_IMMUTABLE;
            if(exitflag==DAQP_EXIT_OPTIMAL)
                for(j=0; j < work->n_active; j++) // all active constraint must also be nonredundant
                    is_redundant[work->WS[j]] = 0;
        }
        // work->dupper[i] -= tol_weak; // TODO support weakly infeasible constraints
        daqp_deactivate_constraints(work);
    }
}

/* Profiling */
#ifdef PROFILING
#ifdef _WIN32
void tic(DAQPtimer *timer){
    QueryPerformanceCounter(&(timer->start));
}
void toc(DAQPtimer *timer){
    QueryPerformanceCounter(&(timer->stop));
}
double get_time(DAQPtimer *timer){
    LARGE_INTEGER f;
    QueryPerformanceFrequency(&f);
    return (double)(timer->stop.QuadPart - timer->start.QuadPart)/f.QuadPart;
}
#else // not _WIN32 (assume that time.h works) 

void tic(DAQPtimer *timer){
    clock_gettime(CLOCK_MONOTONIC, &(timer->start));
}
void toc(DAQPtimer *timer){
    clock_gettime(CLOCK_MONOTONIC, &(timer->stop));
}

double get_time(DAQPtimer *timer){
    struct timespec diff;
    if ((timer->stop.tv_nsec - timer->start.tv_nsec) < 0) {
        diff.tv_sec  = timer->stop.tv_sec - timer->start.tv_sec - 1;
        diff.tv_nsec = 1e9 + timer->stop.tv_nsec - timer->start.tv_nsec;
    } else {
        diff.tv_sec  = timer->stop.tv_sec - timer->start.tv_sec;
        diff.tv_nsec = timer->stop.tv_nsec - timer->start.tv_nsec;
    }
    return (double)diff.tv_sec + (double )diff.tv_nsec / 1e9;
}
#endif // _WIN32
#endif // PROFILING 
