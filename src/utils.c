#include "daqp.h"
#include "utils.h"
#include <math.h>
// Compute upper cholesky factor for H+eps*I and M = A*Rinv
// H is stored in R (packed form) 
// A is stored in M 
// Works inplace, so H and A are overwritten
int compute_Rinv_and_M(c_float *R, c_float *M, const c_float eps,const int n, const int m){
  int i,j,k,disp,disp2;
  // Cholesky
  for (i=0,disp=0; i<n; disp+=n-i,i++) {
	// Diagonal element
	R[disp] += eps;// Add regularization
	for (k=0,disp2=i; k<i; k++,disp2+=n-k) 
	  R[disp] -= R[disp2]*R[disp2];
	if (R[disp] <= 0)
	  return -1; // Not positive definite
	R[disp] = sqrt(R[disp]);

	// Off-diagonal elements
	for (j=1; j<n-i; j++) {
	  for (k=0,disp2=i; k<i; k++,disp2+=n-k)
		R[disp+j] -= R[disp2]*R[disp2+j];
	  R[disp+j] /= R[disp];
	}
	// Store 1/r_ii instead of r_ii 
	// to get multiplication instead division when forward/backward substituting
	R[disp] = 1/R[disp]; 
  }
  // Compute M = A\R  
  for(k = 0,disp=0;k<m;k++)
	for(i = 0,disp2=0; i<n;i++,disp++){
	  M[disp]*=R[disp2++];// Final divide
	  for(j=1;j<n-i;j++)
		M[disp+j] -= R[disp2++]*M[disp];
	}


  // Compute Rinv (store in R) by Rinv = R\I 
  for(k=0,disp=0;k<n;k++){
    disp2=disp;
    R[disp]=R[disp2++]; // Break out first iteration to get rhs
    for(j=k+1;j<n;j++)
      R[disp2++]*=-R[disp];
    disp++;
    for(i=k+1;i<n;i++,disp++){
      R[disp]*=R[disp2++];
      for(j=1;j<n-i;j++)
    	R[disp+j]-=R[disp2++]*R[disp];
    }
  }
  return 0;
}

void update_v_and_d(c_float *f, c_float *bupper, c_float *blower, Workspace *work) 
{
  int i,j,disp;
  c_float sum;

  if(f!=NULL && work->Rinv!=NULL){
	// Compute v = R'\f  
	for(j=work->n-1,disp=ARSUM(work->n);j>=0;j--){
	  for(i=work->n-1;i>j;i--)
		work->v[i] +=work->Rinv[--disp]*f[j];
	  work->v[j]=work->Rinv[--disp]*f[j];
	}
  }

  /* Compute d  = b+M*v */
  // Simple bounds 
  for(i = 0,disp=0;i<N_SIMPLE;i++){
	for(j=i, sum=0;j<work->n;j++)
	  sum+=work->Rinv[disp++]*work->v[j];
	work->dupper[i]=bupper[i]+sum;
	work->dlower[i]=blower[i]+sum;
  }
  //General bounds
  for(i = N_SIMPLE, disp=0;i<N_CONSTR;i++){
	for(j=0, sum=0;j<work->n;j++)
	  sum+=work->M[disp++]*work->v[j];
	work->dupper[i]=bupper[i]+sum;
	work->dlower[i]=blower[i]+sum;
  }
}

// Convert symmetric matrix to packed form
// Example = [1 2 3; 4 5 6; 7 8 9] -> [1 2 3 5 6 9] 
void pack_symmetric(c_float *S, c_float *Sp, const int n){
  if(S==NULL) return;
  int i,j,disp,disp2;
  for(i=0,disp=0,disp2=0;i<n;i++,disp2+=i)
	for(j=i;j<n;j++){
	  Sp[disp++] = S[disp2++];
	}
}

/* Profiling */
double time_diff(struct timespec tic, struct timespec toc){
  struct timespec diff;
  if ((toc.tv_nsec - tic.tv_nsec) < 0) {
    diff.tv_sec  = toc.tv_sec - tic.tv_sec - 1;
    diff.tv_nsec = 1e9 + toc.tv_nsec - tic.tv_nsec;
  } else {
    diff.tv_sec  = toc.tv_sec - tic.tv_sec;
    diff.tv_nsec = toc.tv_nsec - tic.tv_nsec;
  }
  return (double)diff.tv_sec + (double )diff.tv_nsec / 1e9;
}
