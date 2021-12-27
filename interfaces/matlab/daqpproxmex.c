#include "mex.h"
#include "daqp.h"
#include "daqp_prox.h"
#include <time.h>

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
  double *x_star, *cpuTime, *fval;
  int i,*exitflag, *iters;
  /* check for proper number of arguments */
  if(nrhs!=6) {
	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","6 inputs required.");
  }
  if(nlhs!=5) {
	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","5 output required.");
  }

  int n = mxGetM(prhs[0]);
  int m = mxGetN(prhs[0]);
  int n_R= mxGetN(prhs[2]);
  ProxWorkspace prox_work;
  /* create a pointer to the real data in the input matrix  */
#if MX_HAS_INTERLEAVED_COMPLEX
  prox_work.M=mxGetDoubles(prhs[0]);
  prox_work.b = mxGetDoubles(prhs[1]);
  if(n_R>0)
	prox_work.R= mxGetDoubles(prhs[2]);
  else
	prox_work.R= NULL; 
  prox_work.f= mxGetDoubles(prhs[3]);
  prox_work.epsilon = mxGetDoubles(prhs[4])[0];
#else
  prox_work.M= (c_float *)mxGetPr(prhs[0]);
  prox_work.b= (c_float *)mxGetPr(prhs[1]);
  if(n_R>0)
	prox_work.R= (c_float *)mxGetPr(prhs[2]);
  else
	prox_work.R= NULL;
  prox_work.f= (c_float *)mxGetPr(prhs[3]);
  prox_work.epsilon = mxGetPr(prhs[4])[0];
#endif
  // LHS
  plhs[0] = mxCreateDoubleMatrix((mwSize)n,1,mxREAL); // x_star
#if MX_HAS_INTERLEAVED_COMPLEX
  x_star = mxGetDoubles(plhs[0]);
#else
  x_star = mxGetPr(plhs[0]);
#endif


  // Setup output
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); // fval
  plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL); //Exit flag
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); //CPU time
  plhs[4] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL); //Exit flag
  exitflag = (int *)mxGetPr(plhs[2]);
  iters = (int *)mxGetPr(plhs[4]);
  cpuTime = mxGetPr(plhs[3]);
  fval= mxGetPr(plhs[1]);
  
  // Allocate workspaces 
  allocate_prox_workspace(&prox_work,n,m);
  prox_work.work->M = prox_work.M;
#if MX_HAS_INTERLEAVED_COMPLEX
  prox_work.work->sense=mxGetInt32s(prhs[5])
#else
  prox_work.work->sense= (int *)mxGetPr(prhs[5]);
#endif
  
  //TIC
  struct timespec tic,toc;
  clock_gettime(CLOCK_MONOTONIC, &tic);

  add_equality_constraints(prox_work.work);
  // DAQP
  exitflag[0] = daqp_prox(&prox_work);
  iters[0] = prox_work.inner_iterations;
  //printf("daqp_prox leave\n");
  
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


  // Allocate final value
  for(int i =0;i<n;i++)
	x_star[i] = prox_work.x[i];
  
  double func_val=0;
  for(int i = 0;i<n;i++)
	func_val += prox_work.f[i]*x_star[i];
  fval[0] = func_val;
  
  // Free memory
  //free_daqp_workspace(&work);
  //free_prox_workspace(&prox_work);
  //printf("Leaving daqpproxmex\n");
}
