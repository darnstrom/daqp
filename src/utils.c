#include "daqp.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
// Compute upper cholesky factor for H+eps*I and M = A*Rinv
// H is stored in R (packed form) 
// A is stored in M 
int compute_Rinv_and_M(Workspace *work){
  // TODO: Allow M to be computed separately
  // maybe check work->qp->H == NULL instead ? 
  if(work->Rinv == NULL) return 1; //Assume identity Hessian 
  int i,j,k,disp,disp2,disp3;
  const int n = work->n; 
  // Cholesky
  for (i=0,disp=0,disp3=0; i<n; disp+=n-i,i++,disp3+=i) {
	// Diagonal element
	work->Rinv[disp] = work->qp->H[disp3++]+work->settings->eps_prox;// Add regularization
	for (k=0,disp2=i; k<i; k++,disp2+=n-k) 
	  work->Rinv[disp] -= work->Rinv[disp2]*work->Rinv[disp2];
	if (work->Rinv[disp] <= 0) return -1; // Not positive definite
	work->Rinv[disp] = sqrt(work->Rinv[disp]);

	// Off-diagonal elements
	for (j=1; j<n-i; j++) {
	  work->Rinv[disp+j]=work->qp->H[disp3++];
	  for (k=0,disp2=i; k<i; k++,disp2+=n-k)
		work->Rinv[disp+j] -= work->Rinv[disp2]*work->Rinv[disp2+j];
	  work->Rinv[disp+j] /= work->Rinv[disp];
	}
	// Store 1/r_ii instead of r_ii 
	// to get multiplication instead division when forward/backward substituting
	work->Rinv[disp] = 1/work->Rinv[disp]; 
  }

  // Compute Rinv (store in R) by Rinv = R\I 
  for(k=0,disp=0;k<n;k++){
    disp2=disp;
    work->Rinv[disp]=work->Rinv[disp2++]; // Break out first iteration to get rhs
    for(j=k+1;j<n;j++)
      work->Rinv[disp2++]*=-work->Rinv[disp];
    disp++;
    for(i=k+1;i<n;i++,disp++){
      work->Rinv[disp]*=work->Rinv[disp2++];
      for(j=1;j<n-i;j++)
    	work->Rinv[disp+j]-=work->Rinv[disp2++]*work->Rinv[disp];
    }
  }
  
  // Compute M = A*Rinv (by using M' = Rinv'*A')
  const int mA = work->m-work->ms;
  for(k = 0,disp2=n*mA-1;k<mA;k++,disp2-=n){
	for(j = 0, disp=ARSUM(work->n); j<n; ++j){
	  for(i=0;i<j;++i)
		work->M[disp2-i] += work->Rinv[--disp]*work->qp->A[disp2-j];
	  work->M[disp2-j]=work->Rinv[--disp]*work->qp->A[disp2-j];
	}
  }
  return 0;
}

void update_v_and_d(c_float *f, Workspace *work) 
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
	work->dupper[i]=work->qp->bupper[i]+sum;
	work->dlower[i]=work->qp->blower[i]+sum;
  }
  //General bounds
  for(i = N_SIMPLE, disp=0;i<N_CONSTR;i++){
	for(j=0, sum=0;j<work->n;j++)
	  sum+=work->M[disp++]*work->v[j];
	work->dupper[i]=work->qp->bupper[i]+sum;
	work->dlower[i]=work->qp->blower[i]+sum;
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
