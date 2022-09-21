
/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
#include <stdlib.h>

#include "daqp/daqp.h"
#include "daqp/api.h"
#include "daqp/types.h"
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 2
#define y_width 1

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
#define N_VARS 2
#define N_CONSTRS 3 
#define N_BOUNDS 0
#define N_CONTROLS 2 
#define WARM_START 0
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Start function
 *
 */
void daqp_sfun_Start_wrapper(void **pW)
{
/* %%%-SFUNWIZ_wrapper_Start_Changes_BEGIN --- EDIT HERE TO _END */
int n = N_VARS;
    int m = N_CONSTRS;
    int ms = N_BOUNDS;

    // Allocate work struct 
    DAQPWorkspace* work = calloc(1,sizeof(DAQPWorkspace));
    // Allocate settings
    allocate_daqp_settings(work);
    daqp_default_settings(work->settings);

    allocate_daqp_workspace(work,n,0);

    // Allocate qp 
    DAQPProblem *qp = calloc(1,sizeof(DAQPProblem));
    work->qp = qp;
    work->qp->n = n;
    work->qp->m = m;
    work->qp->ms = ms;

    work->qp->H = malloc(n*n*sizeof(c_float));
    work->qp->A = malloc(n*(m-ms)*sizeof(c_float));
    work->qp->f = malloc(n*sizeof(c_float));
    work->qp->bupper = malloc(m*sizeof(c_float));
    work->qp->blower = malloc(m*sizeof(c_float));
    work->qp->sense = NULL; 

    // Setup daqp workspace 
    work->Rinv = malloc(((n+1)*n/2)*sizeof(c_float));
    work->M = malloc(n*(m-ms)*sizeof(c_float));
    work->scaling= malloc(m*sizeof(c_float));
    work->dupper = malloc(m*sizeof(c_float));
    work->dlower = malloc(m*sizeof(c_float));
    work->v = malloc(n*sizeof(c_float));
    work->sense = calloc(1,m*sizeof(int));

    // Add to parameter 
    pW[0] = work;
/* %%%-SFUNWIZ_wrapper_Start_Changes_END --- EDIT HERE TO _BEGIN */
}
/*
 * Output function
 *
 */
void daqp_sfun_Outputs_wrapper(const real_T *H,
			const real_T *f,
			const real_T *A,
			const real_T *bupper,
			const real_T *blower,
			real_T *u,
			void **pW)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
DAQPWorkspace* work = pW[0];

    int ii,jj;
    int n = N_VARS;
    int m = N_CONSTRS;
    int mA = m-N_BOUNDS;

    int update_mask=0;

    // Update H and A
    if(WARM_START < 2 || work->iterations == 0){
        for(ii=0; ii < n*n; ii++) work->qp->H[ii] = H[ii];
        for (ii=0; ii<mA; ++ii )
            for ( jj=0; jj<n; ++jj )
                work->qp->A[ii*n + jj] = A[jj*mA + ii];
        update_mask += UPDATE_Rinv+UPDATE_M;
    }

    // Update f & bupper/blower
    update_mask += UPDATE_v + UPDATE_d;
    for(ii=0; ii < n; ii++) work->qp->f[ii] = f[ii];
    for(ii=0; ii < m; ii++){ 
        work->qp->bupper[ii] = bupper[ii];
        work->qp->blower[ii] = blower[ii];
    }
    // Update internal problem 
    update_ldp(update_mask,work);

    // Correct working set
    if(WARM_START == 0){
        deactivate_constraints(work);
        reset_daqp_workspace(work);
    }
    else if(WARM_START == 1)
        activate_constraints(work);

    // Solve problem
    daqp_ldp(work);
    ldp2qp_solution(work);

    // Extract output
    for(ii = 0; ii < N_CONTROLS; ii++)
        u[ii] = work->x[ii];
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
 * Terminate function
 *
 */
void daqp_sfun_Terminate_wrapper(void **pW)
{
/* %%%-SFUNWIZ_wrapper_Terminate_Changes_BEGIN --- EDIT HERE TO _END */
DAQPWorkspace* work = pW[0]; 
    free_daqp_workspace(work);
    free_daqp_ldp(work);

    free(work->qp->H);
    free(work->qp->f);
    free(work->qp->A);
    free(work->qp->bupper);
    free(work->qp->blower);
    free(work->qp);

    free(work);
/* %%%-SFUNWIZ_wrapper_Terminate_Changes_END --- EDIT HERE TO _BEGIN */
}

