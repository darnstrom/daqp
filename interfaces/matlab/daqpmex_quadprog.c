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
  double *H,*f,*A,*bupper,*blower,*ub,*lb;
  int *exitflag, *sense;
  DAQPResult res;
  /* check for proper number of arguments */
  if(nrhs!=9) {
	mexErrMsgIdAndTxt("DAQP:nrhs","9 inputs required.");
  }
  if(nlhs!=4) {
	mexErrMsgIdAndTxt("DAQP:nlhs","4 output required.");
  }

  int n = mxGetM(prhs[2]);
  int m = mxGetM(prhs[3]);
  int ms = m-mxGetN(prhs[2]);

  // RHS
  if(mxIsEmpty(prhs[0])) H=NULL; else H= mxGetPr(prhs[0]);
  f= mxGetPr(mxDuplicateArray(prhs[1]));
  A= mxGetPr(prhs[2]);
  bupper= mxGetPr(mxDuplicateArray(prhs[3]));
  blower= mxGetPr(mxDuplicateArray(prhs[4]));
  ub= mxGetPr(mxDuplicateArray(prhs[5]));
  lb= mxGetPr(mxDuplicateArray(prhs[6]));
  sense= (int *)mxGetPr(mxDuplicateArray(prhs[7]));
  
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
  const mxArray* mex_settings = prhs[8];
  if(mxIsEmpty(mex_settings))
	daqp_default_settings(&settings);
  else{
	settings.primal_tol = (c_float)mxGetScalar(mxGetField(mex_settings, 0, "primal_tol"));
	settings.dual_tol =  (c_float)mxGetScalar(mxGetField(mex_settings, 0, "dual_tol"));
	settings.zero_tol = (c_float)mxGetScalar(mxGetField(mex_settings, 0, "zero_tol"));
	settings.pivot_tol = (c_float)mxGetScalar(mxGetField(mex_settings, 0, "pivot_tol"));
	settings.progress_tol = (c_float)mxGetScalar(mxGetField(mex_settings, 0, "progress_tol"));
	settings.cycle_tol = (int)mxGetScalar(mxGetField(mex_settings, 0, "cycle_tol"));
	settings.iter_limit= (int)mxGetScalar(mxGetField(mex_settings, 0, "iter_limit"));
	settings.eps_prox = (c_float)mxGetScalar(mxGetField(mex_settings, 0, "eps_prox"));
	settings.eta_prox= (c_float)mxGetScalar(mxGetField(mex_settings, 0, "eta_prox"));
	settings.prox_iter_limit= (int)mxGetScalar(mxGetField(mex_settings, 0, "prox_iter_limit"));
  }

  // Solve problem 
  daqp_quadprog(&res,H,f,A,bupper,blower,sense,n,m,ms,0,&settings);
  
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
