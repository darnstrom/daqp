#include "avi.h"
#include "daqp.h"
#include "utils.h"

// Solve AVI
int _daqp_avi(DAQPAVI *avi) {
    DAQPWorkspace* work = avi->work;
    int n = work->n;
    int i,j,k,disp;
    c_float val,sum,sum2;
    int tot_iter = 0;
    int counter = 0;
    int terminate_limit = 1;
    int exitflag = -4;

    // Set starting iterate
    for(i=0;i<n;i++) avi->x[i] = 0;

    // Start the iterations
    // TODO iter_limit should be the for tot_iter...
    for (k = 0; k < work->settings->iter_limit; k++) {
        // Compute xtemp = H*x + f - (Hsym + I)x
        for(i=0, disp=0; i < work->n; i++){
            sum = sum2 = 0.0;
            for(j=0; j < work->n; j++) {
                sum += avi->problem.H[disp]*avi->x[j];
                sum2 += avi->H1pI[disp++]*avi->x[j];
            }
            avi->Hx[i] = sum;
            avi->xtemp[i] = sum + avi->problem.f[i] - sum2;
        } 

        // Update linear term 
        update_v(avi->xtemp,work,0);
        update_d(work, avi->problem.bupper,avi->problem.blower);

        exitflag = daqp_ldp(work);

        if (exitflag < 0) break;
        ldp2qp_solution(work);
        tot_iter += work->iterations;
        for(i=0; i < n; i++) avi->y[i]=work->x[i];

        // AS has not changed -> Check KKT conditions
        if (work->iterations == 1) {
            if (++counter == terminate_limit) {
                daqp_solve_avi_kkt(avi);
                if (daqp_check_optimal_avi(avi)) {
                    exitflag = 1;
                    break;
                } 
                terminate_limit *= 2;
                if(terminate_limit > 10) terminate_limit = 10;
                continue;
            }
        }
        else {
            counter = 0;
        }

        for (i = 0; i < work->n; i++){
            avi->xtemp[i] = 0.5*avi->y[i]+avi->Hx[i];
            avi->y[i] -= avi->x[i];
        }
        for (i = 0; i < n; i++) {
            avi->xtemp[i] += 0.5*avi->Hsym[i * n + i] * avi->y[i]; // diagonal
            for (j = i + 1; j < n; j++) {
                val = 0.5*avi->Hsym[i * n + j];
                avi->xtemp[i] += val * avi->y[j];
                avi->xtemp[j] += val * avi->y[i];
            }
        }
        daqp_lu_solve(avi->H2pI, avi->P_H2, avi->xtemp, avi->x, n);
    }
    work->iterations = tot_iter;
    return exitflag;
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
void daqp_solve_avi_kkt(DAQPAVI* avi) {
    // TODO handle simple bounds...
    int i, j, k, disp;
    int nAS = avi->work->n_active;
    int n = avi->work->n;
    c_float* lambda = avi->work->lam_star;
    c_float* x = avi->x;
    c_float* S = avi->kkt_buffer;
    c_float* rhs_S = avi->kkt_buffer+nAS*nAS;
    c_float* temp = avi->xtemp; // TODO make sure this does not collide with other
    int* WS = avi->work->WS;

    // Compute S = A_WS * H^-1 * A_WS^T
    for (i = 0; i < nAS; i++) {
        // temp = H^-1 * A_row_WS[i]^T
        daqp_lu_solve(avi->LU_H, avi->P_H, avi->problem.A + WS[i] * n, temp, n);

        for (j = 0; j < nAS; j++) {
            double sum = 0;
            disp = WS[j] * n;
            for (k = 0; k < n; k++) {
                sum += avi->problem.A[disp + k] * temp[k];
            }
            S[j * nAS + i] = sum;
        }
    }

    // Form rhs for Schur system: rhs_S = -A_WS * H^-1 * f - b_WS
    daqp_lu_solve(avi->LU_H, avi->P_H, avi->problem.f, temp, n);

    for (i = 0; i < nAS; i++) {
        int row_idx = WS[i];
        double sum = 0;
        disp = row_idx * n;
        for (j = 0; j < n; j++) {
            sum -= avi->problem.A[disp++] * temp[j];
        }
        sum -= avi->work->sense[row_idx]&2 ? avi->problem.blower[row_idx] : avi->problem.bupper[row_idx];
        rhs_S[i] = sum;
    }

    // Get lambda by solving S * lambda = rhs_S
    daqp_lu(S, avi->P_S, nAS); 
    daqp_lu_solve(S, avi->P_S, rhs_S, lambda, nAS);

    // Get x by solving for: H * x = -f - A_WS^T * lambda
    for (i = 0; i < n; i++) temp[i] = -avi->problem.f[i];
    for (j = 0; j < nAS; j++) {
        c_float lj = lambda[j];
        disp = WS[j] * n;
        for (i = 0; i < n; i++) temp[i] -= avi->problem.A[disp++] * lj;
    }

    // Final back-substitution to get x
    daqp_lu_solve(avi->LU_H, avi->P_H, temp, x, n);
}


int daqp_check_optimal_avi(DAQPAVI* avi){
    DAQPWorkspace* work = avi->work;
    int i,j,disp;
    c_float dual_tol = work->settings->dual_tol;
    c_float primal_tol = work->settings->primal_tol;
    // First check dual variables
    for(i=0; i < work->n_active; i++){
        if(IS_LOWER(work->WS[i])){
            if(work->lam_star[i] > dual_tol) return 0;
        }
        else{
            if(work->lam_star[i] < -dual_tol) return 0;
        }
    }
    // Simple constraints
    for(i=0; i < work->ms; i++){
        if(IS_ACTIVE(i)) continue;
        if(avi->x[i] > avi->problem.bupper[i] + primal_tol) return 0; 
        if(avi->x[i] < avi->problem.blower[i] -primal_tol) return 0;
    }

    // General constraints
    for(disp=0; i < work->m; i++){
        if(IS_ACTIVE(i)){
            disp += work->n;
            continue;
        }
        c_float Ax = 0.0;
        for(j=0; j < work->n; j++) Ax += avi->problem.A[disp++]*avi->x[j];
        if(Ax > avi->problem.bupper[i]+primal_tol) return 0; 
        if(Ax < avi->problem.blower[i] - primal_tol) return 0;
    }
    // All checks passed -> optimal KKT point found
    return 1;
}
