#include "daqp.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

int update_ldp(const int mask, DAQPWorkspace *work, DAQPProblem* qp){
    // TODO: copy dimensions from work->qp? 
    int error_flag, i;
    int do_activate = 0;

    // Add original qp to workspace
    work->qp = qp;

    // Update dimensions of problem
    work->n = qp->n;
    work->m = qp->m;
    work->ms = qp->ms;

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
            error_flag = update_Rinv(work, qp->H);
        else{
            update_avi(work->avi,qp);
            error_flag = update_Rinv(work, work->avi->Hs_rho);
        }
        if(error_flag<0)
            return error_flag;
    }
    /** Update M **/
    if(mask&DAQP_UPDATE_Rinv||mask&DAQP_UPDATE_M){
        error_flag = update_M(work,qp->A,mask);
        if(error_flag<0)
            return error_flag;
    }

    /** Update v **/
    if(mask&DAQP_UPDATE_Rinv||mask&DAQP_UPDATE_v){
        update_v(qp->f,work,mask);
    }

    // Normalize Rinv
    if(mask&DAQP_UPDATE_Rinv){
        normalize_Rinv(work);
    }

    /** Update d **/
    if(mask&DAQP_UPDATE_Rinv||mask&DAQP_UPDATE_M||mask&DAQP_UPDATE_v||mask&DAQP_UPDATE_d){
        error_flag = update_d(work,qp->bupper,qp->blower);
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
            error_flag = activate_constraints(work);
        else{// Activate the first level (since those constraints are hard)
            int m_tmp = work->m;
            work->m = work->break_points[0];
            error_flag = activate_constraints(work);
            work->m = m_tmp;
        }
        if(error_flag<0)
            return error_flag;
    }

    return 0;
}

int update_Rinv(DAQPWorkspace *work, c_float *H){
    int i,j,k,disp,disp2,disp3;
    const int n = work->n; 
        // Check if diagonal
    int is_diagonal = 1;
    for (i=0,disp=1; i<n; i++, disp+=i+1){
        for (j=1; j<n-i; j++,disp++) {
            if(H[disp] > 1e-12 || H[disp] < -1e-12){
                is_diagonal=0;
                break;
            }
        }
        if(is_diagonal == 0) break;
    }

    // If diagonal, just keep track of variable scaling and use Rinv = I
    if(is_diagonal==1){
        if(work->Rinv != NULL){
            work->RinvD = work->Rinv;
            work->Rinv = NULL;
        }
        c_float Hi;
        i=0; disp=0;
        if(work->scaling != NULL){
            for(;i<work->ms;i++,disp+=n){ // Combine with settings scaling
                Hi = H[disp++]+work->settings->eps_prox;
                if (Hi <= 0) return DAQP_EXIT_NONCONVEX;
                Hi = sqrt(Hi);
                work->RinvD[i] = 1/Hi;
                work->scaling[i] = Hi;
            }
        }
        for(;i<n;i++,disp+=n){
            Hi = H[disp++] + work->settings->eps_prox;
            if (Hi <= 0) return DAQP_EXIT_NONCONVEX;
            Hi = sqrt(Hi);
            work->RinvD[i] = 1/Hi;
        }
        return 1;
    }
    // Make sure Rinv can be assinged if not diagonal
    //(necessary if H change from diagonal to non-diagonal)
    if(work->RinvD != NULL && work->Rinv ==NULL){
        work->Rinv= work->RinvD;
        work->RinvD = NULL;
    }


    // Cholesky
    for (i=0,disp=0,disp3=0; i<n; disp+=n-i,i++,disp3+=i) {
        // Diagonal element
        work->Rinv[disp] = H[disp3++]+work->settings->eps_prox;// Add regularization
        for (k=0,disp2=i; k<i; k++,disp2+=n-k) 
            work->Rinv[disp] -= work->Rinv[disp2]*work->Rinv[disp2];
        if (work->Rinv[disp] <= 0) return DAQP_EXIT_NONCONVEX; // Not positive definite 

        work->Rinv[disp] = sqrt(work->Rinv[disp]);

        // Off-diagonal elements
        for (j=1; j<n-i; j++) {
            work->Rinv[disp+j]=H[disp3++];
            for (k=0,disp2=i; k<i; k++,disp2+=n-k)
                work->Rinv[disp+j] -= work->Rinv[disp2]*work->Rinv[disp2+j];
            work->Rinv[disp+j] /= work->Rinv[disp];
        }
        // Store 1/r_ii instead of r_ii 
        // to get multiplication instead division when forward/backward substituting
        work->Rinv[disp] = 1/work->Rinv[disp]; 
    }

    // Compute Rinv (store in R) by Rinv = R\I 
    for(k=0,disp=0;k<n;k++){
        disp2=disp;
        work->Rinv[disp]=work->Rinv[disp2++]; // Break out first iteration to get rhs
        for(j=k+1;j<n;j++)
            work->Rinv[disp2++]*=-work->Rinv[disp];
        disp++;
        for(i=k+1;i<n;i++,disp++){
            work->Rinv[disp]*=work->Rinv[disp2++];
            for(j=1;j<n-i;j++)
                work->Rinv[disp+j]-=work->Rinv[disp2++]*work->Rinv[disp];
        }
    }
    return 1;
}

int update_M(DAQPWorkspace *work, c_float *A, const int mask){
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
    return normalize_M(work);
}

void update_v(c_float *f, DAQPWorkspace *work, const int mask){
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

int update_d(DAQPWorkspace *work, c_float *bupper, c_float *blower){
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

void normalize_Rinv(DAQPWorkspace* work){
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
int normalize_M(DAQPWorkspace* work){
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

int update_avi(DAQPAVI* avi, DAQPProblem* p){
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

    // Set x0 to zero 
    for(i=0;i<n;i++) avi->x[i] = 0; // TODO separate from setup...
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
        add_constraint(work,i,1.0);
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
        deactivate_constraints(work);
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
