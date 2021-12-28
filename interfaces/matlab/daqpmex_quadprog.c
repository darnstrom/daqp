#include "mex.h"
#include "api.h"
#include <time.h>

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
  double *x_star, *cpuTime, *fval;
  double *H,*f,*A,*b;
  int *exitflag, *sense;
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
  x_star = mxGetDoubles(plhs[0]);
#else
  x_star = mxGetPr(plhs[0]);
#endif
  fval= mxGetPr(plhs[1]);
  exitflag = (int *)mxGetPr(plhs[2]);
  cpuTime = mxGetPr(plhs[3]);
  
  //TIC
  struct timespec tic,toc;
  clock_gettime(CLOCK_MONOTONIC, &tic);
  
  exitflag[0] = daqp_quadprog(x_star,H,f,A,b,sense,n,m,0);
  
  //TOC
  clock_gettime(CLOCK_MONOTONIC, &toc);

  struct timespec diff;
  if ((toc.tv_nsec - tic.tv_nsec) < 0) {
    diff.tv_sec  = toc.tv_sec - tic.tv_sec - 1;
    diff.tv_nsec = 1e9 + toc.tv_nsec - tic.tv_nsec;
  } else {
    diff.tv_sec  = toc.tv_sec - tic.tv_sec;
    diff.tv_nsec = toc.tv_nsec - tic.tv_nsec;
  }
  double solve_time = (double)diff.tv_sec + (double )diff.tv_nsec / 1e9;
  cpuTime[0] = solve_time;

  // TODO compute fval
  fval[0] = 0;
}
