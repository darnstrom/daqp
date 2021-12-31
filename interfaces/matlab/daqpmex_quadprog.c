#include "mex.h"
#include "api.h"

const char* INFO_FIELDS[] = {
  "setup_time",           
  "solve_time",           
  "iter",           
  "outer_iter"}; 

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
  double *x_star, *fval;
  double *H,*f,*A,*b;
  int *exitflag, *sense;
  DAQPResult res;
  /* check for proper number of arguments */
  if(nrhs!=6) {
	mexErrMsgIdAndTxt("DAQP:nrhs","6 inputs required.");
  }
  if(nlhs!=4) {
	mexErrMsgIdAndTxt("DAQP:nlhs","4 output required.");
  }

  int n = mxGetM(prhs[2]);
  int m = mxGetN(prhs[2]);

  // RHS
  if(mxIsEmpty(prhs[0])) H=NULL; else H= mxGetPr(prhs[0]);
  f= mxGetPr(mxDuplicateArray(prhs[1]));
  A= mxGetPr(prhs[2]);
  b= mxGetPr(mxDuplicateArray(prhs[3]));
  sense= (int *)mxGetPr(mxDuplicateArray(prhs[4]));
  
  // LHS
  plhs[0] = mxCreateDoubleMatrix((mwSize)n,1,mxREAL); // x_star
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); // fval
  plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL); //Exit flag
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); //CPU time

  res.x = mxGetPr(plhs[0]);
  fval= mxGetPr(plhs[1]);
  exitflag = (int *)mxGetPr(plhs[2]);

  // Setup settings
  DAQPSettings settings;
  if(mxIsEmpty(prhs[5]))
	daqp_default_settings(&settings);
  else{
	settings.primal_tol = (c_float)mxGetScalar(mxGetField(prhs[5], 0, "primal_tol"));
	settings.dual_tol =  (c_float)mxGetScalar(mxGetField(prhs[5], 0, "dual_tol"));
	settings.zero_tol = (c_float)mxGetScalar(mxGetField(prhs[5], 0, "zero_tol"));
	settings.pivot_tol = (c_float)mxGetScalar(mxGetField(prhs[5], 0, "pivot_tol"));
	settings.progress_tol = (c_float)mxGetScalar(mxGetField(prhs[5], 0, "progress_tol"));
	settings.cycle_tol = (int)mxGetScalar(mxGetField(prhs[5], 0, "cycle_tol"));
	settings.iter_limit= (int)mxGetScalar(mxGetField(prhs[5], 0, "iter_limit"));
	settings.eps_prox = (c_float)mxGetScalar(mxGetField(prhs[5], 0, "eps_prox"));
	settings.eta_prox= (c_float)mxGetScalar(mxGetField(prhs[5], 0, "eta_prox"));
	settings.prox_iter_limit= (int)mxGetScalar(mxGetField(prhs[5], 0, "prox_iter_limit"));
  }

  // Solve problem 
  daqp_quadprog(&res,H,f,A,b,sense,n,m,0,&settings);
  
  // Extract info
  exitflag[0] = res.exitflag; 
  fval[0] = res.fval;

  // Create info struct
  int n_info = sizeof(INFO_FIELDS)/sizeof(INFO_FIELDS[0]);
  mxArray* info_struct = mxCreateStructMatrix(1,1,n_info,INFO_FIELDS);
  mxSetField(info_struct, 0, "solve_time", mxCreateDoubleScalar(res.solve_time));
  mxSetField(info_struct, 0, "setup_time", mxCreateDoubleScalar(res.setup_time));
  mxSetField(info_struct, 0, "iter", mxCreateDoubleScalar(res.iter));
  mxSetField(info_struct, 0, "outer_iter", mxCreateDoubleScalar(res.outer_iter));
  plhs[3] = info_struct;
}
