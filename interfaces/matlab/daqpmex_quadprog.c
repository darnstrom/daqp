#include "mex.h"
#include "api.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
  double *x_star, *cpuTime, *fval;
  double *H,*f,*A,*b, *eps;
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
  eps = mxGetPr(mxDuplicateArray(prhs[5]));
  
  // LHS
  plhs[0] = mxCreateDoubleMatrix((mwSize)n,1,mxREAL); // x_star
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); // fval
  plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL); //Exit flag
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); //CPU time

  res.x = mxGetPr(plhs[0]);
  fval= mxGetPr(plhs[1]);
  exitflag = (int *)mxGetPr(plhs[2]);
  cpuTime = mxGetPr(plhs[3]);
  
  
  daqp_quadprog(&res,H,f,A,b,sense,n,m,0,eps[0]);
  
  exitflag[0] = res.exitflag; 
  fval[0] = res.fval;
  cpuTime[0] = res.solve_time; 
}
