#include "mex.h"
#include "api.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
  double *x_star, *cpuTime, *fval;
  double *H,*f,*A,*b;
  int *exitflag, *sense;
  DAQPResult res;
  /* check for proper number of arguments */
  if(nrhs!=5) {
	mexErrMsgIdAndTxt("DAQP:nrhs","5 inputs required.");
  }
  if(nlhs!=4) {
	mexErrMsgIdAndTxt("DAQP:nlhs","4 output required.");
  }

  int n = mxGetM(prhs[2]);
  int m = mxGetN(prhs[2]);
  /* create a pointer to the real data in the input matrix  */
#if MX_HAS_INTERLEAVED_COMPLEX
  H = mxGetDoubles(prhs[0]);
  f= mxGetDoubles(mxDuplicateArray(prhs[1]));
  A= mxGetDoubles(prhs[2]);
  b= mxGetDoubles(mxDuplicateArray(prhs[3]));
  sense = =mxGetInt32s(mxDuplicateArray(prhs[4]));
#else
  H= (c_float *)mxGetPr(prhs[0]);
  f= (c_float *)mxGetPr(mxDuplicateArray(prhs[1]));
  A= (c_float *)mxGetPr(prhs[2]);
  b= (c_float *)mxGetPr(mxDuplicateArray(prhs[3]));
  sense= (int *)mxGetPr(mxDuplicateArray(prhs[4]));
#endif
  // LHS
  plhs[0] = mxCreateDoubleMatrix((mwSize)n,1,mxREAL); // x_star
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); // fval
  plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL); //Exit flag
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); //CPU time
#if MX_HAS_INTERLEAVED_COMPLEX
  res.x = mxGetDoubles(plhs[0]);
#else
  res.x = mxGetPr(plhs[0]);
#endif
  fval= mxGetPr(plhs[1]);
  exitflag = (int *)mxGetPr(plhs[2]);
  cpuTime = mxGetPr(plhs[3]);
  
  
  daqp_quadprog(&res,H,f,A,b,sense,n,m,0);
  
  exitflag[0] = res.exitflag; 
  fval[0] = res.fval;
  cpuTime[0] = res.solve_time; 
}
