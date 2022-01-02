#include "api.h" 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// Check feasibility of Ax <=b 
// (ret = 1 if feasible) 
int daqp_feas(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound){
  int ret;
  reset_daqp_workspace(work); // Reset workspace
  work->m = m;
  work->M = A;
  work->dupper = b;
  //TODO: dlower...
  work->sense = sense;
  work->fval_bound = fval_bound;
  ret = daqp(work);
  // Cleanup sense for reuse...
  for(int j=0;j<work->n_active;j++){
	if(IS_IMMUTABLE(work->WS[j])) continue;
	SET_INACTIVE(work->WS[j]); 
  }
  return ret;
}

// Check feasibility of Ax <= b with the sovler started with working set WS 
// (ret = 1 if feasible)
int daqp_feas_warmstart(Workspace* work, c_float* A, c_float*b, int *sense, const int m, c_float fval_bound, int* WS, const int n_active){
  int ret;
  reset_daqp_workspace(work); // Reset workspace
  work->m = m;
  work->M = A;
  work->dupper = b;
  //TODO: set dlower...
  work->sense = sense;
  warmstart_workspace(work,WS,n_active);
  work->fval_bound = fval_bound;

  ret = daqp(work);

  // Return working set in WS
  for(int i=0;i<work->n_active;i++) {
	WS[i] = work->WS[i];
  }
  WS[work->n] = work->n_active; // Number of active in last element
  
  // Cleanup sense for reuse...
  for(int j=0;j<work->n_active;j++){
	if(IS_IMMUTABLE(work->WS[j])) continue;
	SET_INACTIVE(work->WS[j]); 
  }
  return ret;
}

// Create workspace, including allocated iterates  
// ws_ptr is pointing to the created workspace
void daqp_setup(Workspace** ws_ptr, double* M, double* d, int n){
  Workspace* work=malloc(sizeof(Workspace));
  work->n = n;
  work->M = M;
  work->dupper = d;
  // TODO: setup dlower...
  allocate_daqp_workspace(work, n);
  *ws_ptr = work;
}

void daqp_quadprog(DAQPResult *res, double* H, double* f, double *A, double *bupper, double* blower, int* sense, int n, int m, int packed, DAQPSettings *settings){
  struct timespec tstart,tsetup,tsol;
  //DAQPSettings settings;
  //daqp_default_settings(&settings);
  int i;
  // TIC start
  clock_gettime(CLOCK_MONOTONIC, &tstart);

  // Transform QP to ldp (R,v,M,d is stored inplace of H,f,A,b)
  if(!packed) pack_symmetric(H,H,n); 
  qp2ldp(H,f,A,bupper,blower,n,m,settings->eps_prox);

  if(settings->eps_prox==0){
	// Setup daqp workspace 
	Workspace work;
	allocate_daqp_workspace(&work,n);
	work.n = n; work.m = m;
	work.ms = 0;
	work.Rinv = H; work.v = f;
	work.M = A; work.dupper = bupper; work.dlower = blower;
	work.x = res->x;
	work.sense = sense;
	work.settings = settings;
	add_equality_constraints(&work);
	clock_gettime(CLOCK_MONOTONIC, &tsetup);


	// Solve LDP
	res->exitflag = daqp(&work);
	clock_gettime(CLOCK_MONOTONIC, &tsol);
	// Correct offset in fval 
	for(i=0;i<n;i++)
	  work.fval-=work.v[i]*work.v[i];// 
	work.fval *=0.5;
	res->fval = work.fval;
	res->iter = work.iterations;
	res->outer_iter = 0;
	
	free_daqp_workspace(&work);
  }
  else{
	// Setup workspace
	ProxWorkspace prox_work;
	allocate_prox_workspace(&prox_work,n,m);
	prox_work.work->ms = 0; // TODO: move this into allocate
	prox_work.work->Rinv=H; prox_work.work->M=A; 
	prox_work.bupper = bupper; prox_work.blower=blower; prox_work.f = f;
	prox_work.work->sense = sense;
	prox_work.epsilon = settings->eps_prox;
	prox_work.work->settings= settings;
	add_equality_constraints(prox_work.work);
	clock_gettime(CLOCK_MONOTONIC, &tsetup);

	// Solve problem
	res->exitflag = daqp_prox(&prox_work);
	clock_gettime(CLOCK_MONOTONIC, &tsol);

	for(i=0;i<n;i++) // Extract solution
	  res->x[i]=prox_work.x[i];
	res->fval = 0; // TODO: correct this
	res->iter = prox_work.inner_iterations;
	res->outer_iter = prox_work.outer_iterations;

	free_prox_workspace(&prox_work);
  }

  // Append time 
  res->solve_time = time_diff(tsetup,tsol);
  res->setup_time = time_diff(tstart,tsetup);

}

int qp2ldp(double *R, double *v, double* M, double* dupper, double*  dlower, int n, int m, double eps)
{
  // Assumptions: 
  //  H stored in R (in packed form); f is stored in v 
  //  A is stored in M; b is stored in d
  int i, j, k;
  int disp,disp2;
  c_float sum;
  
  // If LP -> no cholesky, just scale f with eps
  if(R==NULL){ 
	for(i=0;i<n;i++)
	  v[i]/=eps;
	return -1; 
  }

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


  // Compute Rinv (store in R...)
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

  if(eps != 0) return  0; // No need to compute v & d for prox
  
  // Compute v = R'\f  
  for(j=n-1,disp=ARSUM(n);j>=0;j--){
	for(i=n-1;i>j;i--)
	  v[i] +=R[--disp]*v[j];
	v[j]*=R[--disp];
  }

  // Compute d  = b+M*v
  for(i = 0, disp=0;i<m;i++){
	sum = 0;
	for(j=0;j<n;j++)
	  sum+=M[disp++]*v[j];
	dupper[i]+=sum;
	dlower[i]+=sum;
  }

  return 0;
}

void pack_symmetric(double *S, double *Sp, int n){
  if(S==NULL) return;
  int i,j,disp,disp2;
  for(i=0,disp=0,disp2=0;i<n;i++,disp2+=i)
	for(j=i;j<n;j++){
	  Sp[disp++] = S[disp2++];
	}
}

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
