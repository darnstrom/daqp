#include "mex.h"
#include "daqp.h"
#include <time.h>

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
  double *x_star, *cpuTime, *fval;
  int *exitflag;
  /* check for proper number of arguments */
  if(nrhs!=3) {
	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","3 inputs required.");
  }
  if(nlhs!=4) {
	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","3 output required.");
  }

  int n = mxGetM(prhs[0]);
  int m = mxGetN(prhs[0]);
  Workspace work;
  work.n=n;
  work.m=m;
  /* create a pointer to the real data in the input matrix  */
#if MX_HAS_INTERLEAVED_COMPLEX
  work.M = mxGetDoubles(prhs[0]);
  work.d= mxGetDoubles(prhs[1]);
  work.sense=mxGetInt32s(prhs[2])
#else
  work.M= (c_float *)mxGetPr(prhs[0]);
  work.d= (c_float *)mxGetPr(prhs[1]);
  work.sense= (int *)mxGetPr(prhs[2]);
#endif
  // LHS
  plhs[0] = mxCreateDoubleMatrix((mwSize)n,1,mxREAL); // x_star
#if MX_HAS_INTERLEAVED_COMPLEX
  x_star = mxGetDoubles(plhs[0]);
#else
  x_star = mxGetPr(plhs[0]);
#endif

  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); // fval
  plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL); //Exit flag
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); //CPU time
  exitflag = (int *)mxGetPr(plhs[2]);
  cpuTime = mxGetPr(plhs[3]);
  fval= mxGetPr(plhs[1]);
  
  // Allocate daqp
  allocate_daqp_workspace(&work, n);

  //TIC
  struct timespec tic,toc;
  clock_gettime(CLOCK_MONOTONIC, &tic);
  // DAQP
  add_equality_constraints(&work);
  exitflag[0] = daqp(&work);
  
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

  double func_val;
  fval[0] = work.fval;

  // Allocate final value
  for(int i =0;i<n;i++)
	x_star[i] = work.u[i];
  
  // Free memory
  //freeWorkspaceIters(&work);
}
