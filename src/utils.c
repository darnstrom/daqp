#include "daqp.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

int update_ldp(const int mask, Workspace *work){
  // TODO: copy dimensions from work->qp? 
  int error_flag;
  /** Update Rinv **/
  if(mask&UPDATE_Rinv){
	error_flag = update_Rinv(work);
	if(error_flag<0)
	  return error_flag;
  }
  /** Update M **/
  if(mask&UPDATE_Rinv||mask&UPDATE_M){
	update_M(work);
  }

  /** Update v **/
  if(mask&UPDATE_Rinv||mask&UPDATE_v){
	update_v(work->qp->f,work);
  }
  
  /** Update d **/
  if(mask&UPDATE_Rinv||mask&UPDATE_M||mask&UPDATE_v||mask&UPDATE_d){
	update_d(work);
  }

  /** Update constraint sense **/
  if(mask&UPDATE_sense){
	if(work->qp->sense == NULL) // Assume all constraints are "normal" inequality constraints
	  for(int i=0;i<work->m;i++) work->sense[i] = 0;
	else
	  for(int i=0;i<work->m;i++) work->sense[i] = work->qp->sense[i];
  }
  
  return 0;
}

int update_Rinv(Workspace *work){
  int i,j,k,disp,disp2,disp3;
  const int n = work->n; 
  // Cholesky
  for (i=0,disp=0,disp3=0; i<n; disp+=n-i,i++,disp3+=i) {
	// Diagonal element
	work->Rinv[disp] = work->qp->H[disp3++]+work->settings->eps_prox;// Add regularization
	for (k=0,disp2=i; k<i; k++,disp2+=n-k) 
	  work->Rinv[disp] -= work->Rinv[disp2]*work->Rinv[disp2];
	if (work->Rinv[disp] <= 0) return EXIT_NONCONVEX; // Not positive definite 
	//TODO: handle singular case by regularization
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
  return 1;
}

void update_M(Workspace *work){
  int i,j,k,disp,disp2;
  const int n = work->n;
  const int mA = work->m-work->ms;
  for(k = 0,disp2=n*mA-1;k<mA;k++,disp2-=n){
	for(j = 0, disp=ARSUM(n); j<n; ++j){
	  for(i=0;i<j;++i)
		work->M[disp2-i] += work->Rinv[--disp]*work->qp->A[disp2-j];
	  work->M[disp2-j]=work->Rinv[--disp]*work->qp->A[disp2-j];
	}
  }
  reset_daqp_workspace(work); // Internal factorizations need to be redone!
}

void update_v(c_float *f, Workspace *work){
  int i,j,disp;
  const int n = work->n;
  if(work->Rinv == 0){// Rinv = I => v = R'\v = f
	for(i=0;i<n;++i) work->v[i] = f[i];
	return;
  }
  for(j=n-1,disp=ARSUM(n);j>=0;j--){
	for(i=n-1;i>j;i--)
	  work->v[i] +=work->Rinv[--disp]*f[j];
	work->v[j]=work->Rinv[--disp]*f[j];
  }
}

void update_d(Workspace *work){
  /* Compute d  = b+M*v */
  // Simple bounds 
  int i,j,disp;
  c_float sum;
  const int n = work->n;
  if(work->Rinv !=NULL){
	for(i = 0,disp=0;i<N_SIMPLE;i++){
	  for(j=i, sum=0;j<n;j++)
		sum+=work->Rinv[disp++]*work->v[j];
	  work->dupper[i]=work->qp->bupper[i]+sum;
	  work->dlower[i]=work->qp->blower[i]+sum;
	}
  }else{
	for(i = 0,disp=0;i<N_SIMPLE;i++){
	  work->dupper[i]=work->qp->bupper[i]+work->v[i];
	  work->dlower[i]=work->qp->blower[i]+work->v[i];
	}
  }
  //General bounds
  for(i = N_SIMPLE, disp=0;i<N_CONSTR;i++){
	for(j=0, sum=0;j<n;j++)
	  sum+=work->M[disp++]*work->v[j];
	work->dupper[i]=work->qp->bupper[i]+sum;
	work->dlower[i]=work->qp->blower[i]+sum;
  }

  work->reuse_ind = 0; // RHS of KKT system changed => cannot reuse intermediate results 
}

/* Profiling */
#ifdef _WIN32
void tic(DAQPtimer *timer){
  QueryPerformanceCounter(&(timer->start));
}
void toc(DAQPtimer *timer){
  QueryPerformanceCounter(&(timer->stop));
}
double get_time(DAQPtimer *timer){
  LARGE_INTEGER f;
  QueryPerformanceFrequency(&f);
  return (double)(timer->stop.QuadPart - timer->start.QuadPart)/f.QuadPart;
}
#else // not _WIN32 (assume that time.h works) 

void tic(DAQPtimer *timer){
  clock_gettime(CLOCK_MONOTONIC, &(timer->start));
}
void toc(DAQPtimer *timer){
  clock_gettime(CLOCK_MONOTONIC, &(timer->stop));
}

double get_time(DAQPtimer *timer){
  struct timespec diff;
  if ((timer->stop.tv_nsec - timer->start.tv_nsec) < 0) {
    diff.tv_sec  = timer->stop.tv_sec - timer->start.tv_sec - 1;
    diff.tv_nsec = 1e9 + timer->stop.tv_nsec - timer->start.tv_nsec;
  } else {
    diff.tv_sec  = timer->stop.tv_sec - timer->start.tv_sec;
    diff.tv_nsec = timer->stop.tv_nsec - timer->start.tv_nsec;
  }
  return (double)diff.tv_sec + (double )diff.tv_nsec / 1e9;
}
#endif // _WIN32
