#include "api.h" 
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void daqp_quadprog(DAQPResult *res, QP* qp, DAQPSettings *settings){
  struct timespec tstart,tsetup,tsol;
  // TIC start
  clock_gettime(CLOCK_MONOTONIC, &tstart);
  
  Workspace work;
  setup_daqp(qp,settings,&work);

  // Solve problem
  if(settings->eps_prox==0){
	clock_gettime(CLOCK_MONOTONIC, &tsetup);
	res->exitflag = daqp(&work);
	clock_gettime(CLOCK_MONOTONIC, &tsol);
	work.inner_iter = work.iterations;
  }
  else{//Prox
	clock_gettime(CLOCK_MONOTONIC, &tsetup);
	res->exitflag = daqp_prox(&work);
	clock_gettime(CLOCK_MONOTONIC, &tsol);
  }
  
  // Extract results
  daqp_extract_result(res,&work);

  // Add time to result
  res->solve_time = time_diff(tsetup,tsol);
  res->setup_time = time_diff(tstart,tsetup);
  
  // Free memory
  free_daqp_workspace(&work);
  free_daqp_ldp(&work);
}

void setup_daqp(QP* qp, DAQPSettings *settings, Workspace *work){
  // Check if QP is well-posed
  //validate_QP(qp);
  
  // Setup workspace
  work->settings = settings;
  allocate_daqp_workspace(work,qp->n);
  setup_daqp_ldp(work,qp);
  add_equality_constraints(work);
}

//  Setup LDP from QP  
void setup_daqp_ldp(Workspace *work, QP *qp){
  int i;

  work->n = qp->n;
  work->m = qp->m;
  work->ms = qp->ms;
  work->qp = qp;
 // Extract data for Rinv and M 
  if(qp->H!=NULL){ 
	work->Rinv = malloc(((qp->n+1)*qp->n/2)*sizeof(c_float));
	work->M = malloc(qp->n*(qp->m-qp->ms)*sizeof(c_float));
  }
  else{// H = I =>  no need to transform H->Rinv and M->A 
	work->Rinv = qp->H;
	work->M = qp->A;
  }

  // Extract data for v and d
  if(qp->f!=NULL || work->settings->eps_prox != 0){
	work->dupper = malloc(qp->m*sizeof(c_float));
	work->dlower = malloc(qp->m*sizeof(c_float));
	work->v = malloc(qp->n*sizeof(c_float));
  }
  else{ // f = 0 => no need to transform f->v and b->d 
	work->v=qp->f;
	work->dupper = qp->bupper; 
	work->dlower = qp->blower; 
  }
  
  // Transform QP to LDP
  if(qp->H != 0){
	// Copy H->Rinv and A->M
	pack_symmetric(qp->H,work->Rinv,qp->n);
	const int nA = qp->n*(qp->m-qp->ms);
	for(i=0; i<nA;i++) work->M[i] = qp->A[i];
	
	// Compute Rinv and M
	compute_Rinv_and_M(work->Rinv,work->M,work->settings->eps_prox,qp->n,qp->m-qp->ms);
  }
  // Compute v and M
  update_v_and_d(qp->f,qp->bupper,qp->blower,work);

  // Setup up constraint states
  if(qp->sense == NULL) // Assume all constraint are "normal" inequality constraints
	work->sense = calloc(qp->m,sizeof(int));
  else{
	work->sense = malloc(qp->m*sizeof(int));
	for(int i=0;i<qp->m;i++) work->sense[i] = qp->sense[i];
  }
}

// Free data for LDP 
void free_daqp_ldp(Workspace *work){
  free(work->sense);
  if(work->Rinv != NULL){
	free(work->Rinv);
	free(work->M);
  }
  if(work->v != NULL){
	free(work->v);
	free(work->dupper);
	free(work->dlower);
  }

}

void daqp_extract_result(DAQPResult* res, Workspace* work){
  // Extract optimal solution and correct fval offset
  res->fval = work->fval;
  for(int i=0;i<work->n;i++){
	res->x[i] = work->x[i];
	res->fval-=work->v[i]*work->v[i]; 
  }
  res->fval *=0.5;
  res->soft_slack = work->soft_slack;
  res->iter = work->inner_iter;
  res->outer_iter = work->outer_iter; // TODO add iter...
}
