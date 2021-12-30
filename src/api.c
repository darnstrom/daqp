#include "api.h" 
#include <stdlib.h>
//#include <stdio.h>
#include <math.h>
// Check feasibility of Ax <=b 
// (ret = 1 if feasible) 
int daqp_feas(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound){
  int ret;
  reset_daqp_workspace(work); // Reset workspace
  work->m = m;
  work->M = A;
  work->d = b;
  work->sense = sense;
  work->fval_bound = fval_bound;
  ret = daqp(work);
  // Cleanup sense for reuse...
  for(int j=0;j<work->n_active;j++)
	work->sense[work->WS[j]]=INACTIVE_INEQUALITY; // TODO handle equality constraints...
  return ret;
}

// Check feasibility of Ax <= b with the sovler started with working set WS 
// (ret = 1 if feasible)
int daqp_feas_warmstart(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound, int* WS, const int n_active){
  int ret;
  reset_daqp_workspace(work); // Reset workspace
  work->m = m;
  work->M = A;
  work->d = b;
  work->sense = sense;
  warmstart_workspace(work,WS,n_active);
  work->fval_bound = fval_bound;

  ret = daqp(work);

  // Return working set in WS
  for(int i=0;i<work->n_active;i++) {
	WS[i] = work->WS[i];
  }
  WS[work->n] = work->n_active; // Number of active in last element
  
  for(int j=0;j<work->n_active;j++)
	work->sense[work->WS[j]]=INACTIVE_INEQUALITY; // TODO handle equality constraints...
  return ret;
}

// Create workspace, including allocated iterates  
// ws_ptr is pointing to the created workspace
void daqp_setup(Workspace** ws_ptr, double* M, double* d, int n){
  Workspace* work=malloc(sizeof(Workspace));
  work->n = n;
  work->M = M;
  work->d = d;
  allocate_daqp_workspace(work, n);
  *ws_ptr = work;
}

int daqp_quadprog(double* x, double* H, double* f, double *A, double *b, int* sense, int n, int m, int packed){
  int ret; 
  Workspace work;

  // Transform QP to ldp (R,v,M,d is stored inplace of H,f,A,b)
  if(!packed) pack_symmetric(H,H,n); 
  qp2ldp(H,f,A,b,n,m,0);

  // Setup daqp workspace 
  allocate_daqp_workspace(&work,n);
  work.n = n; work.m = m;
  work.R = H; work.v = f;
  work.M = A; work.d = b;
  work.x = x;
  work.sense = sense;
  add_equality_constraints(&work);
  
  // Solve LDP
  ret = daqp(&work);

  // Copy solution and cleanup
  free_daqp_workspace(&work);
  return ret;
}

int qp2ldp(double *R, double *v, double* M, double* d, int n, int m, double eps)
{
  // Assumption: 
  //  H is stored in R (in packed form)
  //  A is stored in M
  //  f is stored in v 
  //  b is stored in d
  int i, j, k;
  int disp,disp2;

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
  // Compute v = R'\f  
  for(i = 0,disp=0;i<n;i++){
	v[i]*=R[disp++];
	for(j=i+1;j<n;j++)
	  v[j] -= R[disp++]*v[i];
  }

  // Compute d  = b+M*v
  for(i = 0, disp=0;i<m;i++)
	for(j=0;j<n;j++)
	  d[i]+=M[disp++]*v[j];

  return 0;
}

void pack_symmetric(double *S, double *Sp, int n){
  int i,j,disp,disp2;
  for(i=0,disp=0,disp2=0;i<n;i++,disp2+=i)
	for(j=i;j<n;j++){
	  Sp[disp++] = S[disp2++];
	}
}
