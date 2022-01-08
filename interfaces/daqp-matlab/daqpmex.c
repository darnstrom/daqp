#include "mex.h"
#include "api.h"
#include "utils.h"
#include <string.h>

const char* INFO_FIELDS[] = {
  "setup_time",           
  "solve_time",           
  "iter",           
  "outer_iter",
  "soft_slack"}; 

const char* SETTINGS_FIELDS[] = {
  "primal_tol",           
  "dual_tol",           
  "zero_tol",           
  "pivot_tol",
  "progress_tol",
  "cycle_tol",
  "iter_limit",
  "eps_prox",
  "eta_prox",
  "prox_iter_limit",
  "rho_soft"}; 


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{

  // RHS 

  // Extract command
  char cmd[64];
  mxGetString(prhs[0], cmd, sizeof(cmd));
  
  // Extract workspace pointer 
  Workspace *work;
  long long *work_i;
  union{long long i; void *ptr;} work_ptr; // Used for int64 & pointer juggling..
  if(nrhs>1){// Pointer always second argument (stored as int64)
	work_i = (long long *)mxGetData(prhs[1]);
	work_ptr.i  = *work_i; 
	work = work_ptr.ptr;
  }

	if (!strcmp("new", cmd)) {
	  // Allocate new workspace and return pointer to it
	  work_ptr.ptr = calloc(1,sizeof(Workspace));
	  // Return pointer as int64
	  plhs[0] = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	  work_i = (long long *) mxGetData(plhs[0]);
	  *work_i = work_ptr.i;
	  return;
	}
	else if (!strcmp("delete", cmd)) {
	  // Free workspace 
	  free_daqp_workspace(work);
	  free_daqp_ldp(work);
	  //if(work->qp && work->qp->f) free(work->qp->f);
	  //if(work->qp && work->qp->bupper) free(work->qp->bupper);
	  //if(work->qp && work->qp->blower) free(work->qp->blower);
	  if(work->qp) free(work->qp);
	  free(work->settings);
	  free(work);
	  return;
	}
	else if (!strcmp("setup", cmd)) {

	  if(work->qp) mexErrMsgTxt("Setup already completed.");
	  QP *qp = calloc(1,sizeof(QP));
	  // Extract data
	  double *H,*f,*A,*bupper,*blower;
	  int *sense,error_flag;
	  
	  // Get dimensions 
	  int n = mxGetM(prhs[4]);
	  int m = mxGetM(prhs[5]);
	  int ms = m-mxGetN(prhs[4]);
	  
	  // Setup QP struct
	  qp->n = n;
	  qp->m = m;
	  qp->ms = ms;
	  qp->H= mxGetPr(prhs[2]);
	  qp->f= mxGetPr(prhs[3]);
	  qp->A= mxGetPr(prhs[4]);
	  qp->bupper= mxGetPr(prhs[5]);
	  qp->blower= mxGetPr(prhs[6]);
	  qp->sense= (int *)mxGetPr(prhs[7]);
	  // Copy f and b since these might be used after mex ends
	  //qp->f = calloc(n,sizeof(c_float));
	  //qp->bupper = calloc(m,sizeof(c_float));
	  //qp->blower = calloc(m,sizeof(c_float));
	  //for(int i=0;i<n;i++) qp->f[i] = f[i];
	  //for(int i=0;i<m;i++) qp->bupper[i]=bupper[i];	
	  //for(int i=0;i<m;i++) qp->blower[i]=blower[i];	

	  // Setup QP struct
	  //QP qp={n,m,ms,H,f,A,bupper,blower,sense};
	  // Call C API 
	  error_flag = setup_daqp(qp,work->settings,work);
	  if(error_flag < 0){
		printf("Setup failed (%d)\n",error_flag);
		free(work->qp);
		work->qp = NULL;
	  }

	  return;
	}
	else if (!strcmp("solve", cmd)) {
	  int error_flag;
	  double *fval;
	  int *exitflag;
	  DAQPResult result;
	  if(work->qp == NULL) mexErrMsgTxt("No problem to solve");
	  // Update QP pointers 
	  work->qp->H= mxGetPr(prhs[2]);
	  work->qp->f= mxGetPr(prhs[3]);
	  work->qp->A= mxGetPr(prhs[4]);
	  work->qp->bupper= mxGetPr(prhs[5]);
	  work->qp->blower= mxGetPr(prhs[6]);
	  work->qp->sense= (int *)mxGetPr(prhs[7]);
	  // Setup output 
	  plhs[0] = mxCreateDoubleMatrix((mwSize)work->n,1,mxREAL); // x_star
	  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); // fval
	  plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL); //Exit flag
	  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); //CPU time
	  plhs[0] = mxCreateDoubleMatrix((mwSize)work->n,1,mxREAL); // x_star
	  result.x = mxGetPr(plhs[0]);
	  fval= mxGetPr(plhs[1]);
	  exitflag = (int *)mxGetPr(plhs[2]);
	  
	  // Solve problem
	  daqp_solve(&result,work); 
	  // Package solution information
	  exitflag[0] = result.exitflag; 
	  fval[0] = result.fval;
	  //
	  int n_info = sizeof(INFO_FIELDS)/sizeof(INFO_FIELDS[0]);
	  mxArray* info_struct = mxCreateStructMatrix(1,1,n_info,INFO_FIELDS);
	  mxSetField(info_struct, 0, "solve_time", mxCreateDoubleScalar(result.solve_time));
	  mxSetField(info_struct, 0, "setup_time", mxCreateDoubleScalar(result.setup_time));
	  mxSetField(info_struct, 0, "iter", mxCreateDoubleScalar(result.iter));
	  mxSetField(info_struct, 0, "outer_iter", mxCreateDoubleScalar(result.outer_iter));
	  mxSetField(info_struct, 0, "soft_slack", mxCreateDoubleScalar(result.soft_slack));
	  plhs[3] = info_struct;
	}
	else if (!strcmp("set_default_settings", cmd)){
	  if(work->settings == 0) work->settings = malloc(sizeof(DAQPSettings));
	  daqp_default_settings(work->settings);
	}
	else if (!strcmp("get_settings", cmd)) {
	  if(work->settings != NULL){
		int n_settings = sizeof(SETTINGS_FIELDS)/sizeof(SETTINGS_FIELDS[0]);
		mxArray* s = mxCreateStructMatrix(1,1,n_settings,SETTINGS_FIELDS);
		mxSetField(s, 0, "primal_tol", mxCreateDoubleScalar(work->settings->primal_tol));
		mxSetField(s, 0, "dual_tol", mxCreateDoubleScalar(work->settings->dual_tol));
		mxSetField(s, 0, "zero_tol", mxCreateDoubleScalar(work->settings->zero_tol));
		mxSetField(s, 0, "pivot_tol", mxCreateDoubleScalar(work->settings->pivot_tol));
		mxSetField(s, 0, "progress_tol", mxCreateDoubleScalar(work->settings->progress_tol));
		mxSetField(s, 0, "cycle_tol", mxCreateDoubleScalar(work->settings->cycle_tol));
		mxSetField(s, 0, "iter_limit", mxCreateDoubleScalar(work->settings->iter_limit));
		mxSetField(s, 0, "eps_prox", mxCreateDoubleScalar(work->settings->eps_prox));
		mxSetField(s, 0, "eta_prox", mxCreateDoubleScalar(work->settings->eta_prox));
		mxSetField(s, 0, "prox_iter_limit", mxCreateDoubleScalar(work->settings->prox_iter_limit));
		mxSetField(s, 0, "rho_soft", mxCreateDoubleScalar(work->settings->rho_soft));
		plhs[0] = s;
	  }
	}
	else if (!strcmp("set_settings", cmd)) {
	  const mxArray* s = prhs[2];
	  work->settings->primal_tol = (c_float)mxGetScalar(mxGetField(s, 0, "primal_tol"));
	  work->settings->dual_tol =  (c_float)mxGetScalar(mxGetField(s, 0, "dual_tol"));
	  work->settings->zero_tol = (c_float)mxGetScalar(mxGetField(s, 0, "zero_tol"));
	  work->settings->pivot_tol = (c_float)mxGetScalar(mxGetField(s, 0, "pivot_tol"));
	  work->settings->progress_tol = (c_float)mxGetScalar(mxGetField(s, 0, "progress_tol"));
	  work->settings->cycle_tol = (int)mxGetScalar(mxGetField(s, 0, "cycle_tol"));
	  work->settings->iter_limit= (int)mxGetScalar(mxGetField(s, 0, "iter_limit"));
	  work->settings->eps_prox = (c_float)mxGetScalar(mxGetField(s, 0, "eps_prox"));
	  work->settings->eta_prox= (c_float)mxGetScalar(mxGetField(s, 0, "eta_prox"));
	  work->settings->prox_iter_limit= (int)mxGetScalar(mxGetField(s, 0, "prox_iter_limit"));
	  work->settings->rho_soft= (c_float)mxGetScalar(mxGetField(s, 0, "rho_soft"));
	}
	else if (!strcmp("update", cmd)) {
	  if(work->qp == NULL) mexErrMsgTxt("No problem to update");
	  // Update QP pointers 
	  work->qp->H= mxGetPr(prhs[2]);
	  work->qp->f= mxGetPr(prhs[3]);
	  work->qp->A= mxGetPr(prhs[4]);
	  work->qp->bupper= mxGetPr(prhs[5]);
	  work->qp->blower= mxGetPr(prhs[6]);
	  work->qp->sense= (int *)mxGetPr(prhs[7]);
	  // Update LDP with new QP data
	  const int update_mask = (int)mxGetScalar(prhs[8]);
	  update_ldp(update_mask,work);
	}

  // RHS
}
