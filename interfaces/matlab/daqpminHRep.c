#include "mex.h"
#include "daqp.h"
#include "api.h"
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
  if(nlhs!=2) {
	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","2 output required.");
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
  plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL); //Exit flag
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); //CPU time
  exitflag = (int *)mxGetPr(plhs[0]);
  cpuTime = mxGetPr(plhs[1]);
  
  // Allocate daqp
  allocate_daqp_workspace(&work, n);

  //TIC
  struct timespec tic,toc;
  clock_gettime(CLOCK_MONOTONIC, &tic);
  // DAQP
  exitflag[0] = minimal_Hrepresentation(&work); //TODO, change exitflag...
  
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

  // Free memory
  //freeWorkspaceIters(&work);
}
