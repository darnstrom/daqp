#include "api.h" 
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// DAQP + timing
void daqp_solve(DAQPResult *res, Workspace *work){
  struct timespec tstart,tsol;
  clock_gettime(CLOCK_MONOTONIC, &tstart); //TIC
  // Select algorithm
  if(work->settings->eps_prox==0){
	res->exitflag = daqp(work);
	work->inner_iter = work->iterations;
  }
  else{//Prox
	res->exitflag = daqp_prox(work);
  }

  clock_gettime(CLOCK_MONOTONIC, &tsol); // TOC
  
  // Package result
  daqp_extract_result(res,work);
  // Add time to result
  res->solve_time = time_diff(tstart,tsol);
  res->setup_time = 0; 
}

void daqp_quadprog(DAQPResult *res, QP* qp, DAQPSettings *settings){
  struct timespec tstart,tsetup;
  int setup_flag;
  
  clock_gettime(CLOCK_MONOTONIC, &tstart); //TIC
  Workspace work;
  setup_flag = setup_daqp(qp,settings,&work);
  clock_gettime(CLOCK_MONOTONIC, &tsetup); //TOC

  if(setup_flag >= 0)
	daqp_solve(res,&work);
  else
	res->exitflag = setup_flag;
  
  // Add setup time to result 
  res->setup_time = time_diff(tstart,tsetup);
  // Free memory
  free_daqp_workspace(&work);
  free_daqp_ldp(&work);
}

int setup_daqp(QP* qp, DAQPSettings *settings, Workspace *work){
  int errorflag;
  // Check if QP is well-posed
  //validate_QP(qp);
  
  // Setup workspace
  work->settings = settings;
  allocate_daqp_workspace(work,qp->n);
  errorflag = setup_daqp_ldp(work,qp);
  if(errorflag < 0){
	free_daqp_workspace(work);
	return errorflag;
  }
  errorflag = add_equality_constraints(work);
  if(errorflag < 0){
	free_daqp_workspace(work);
	return errorflag;
  }
  return 1;
}

//  Setup LDP from QP  
int setup_daqp_ldp(Workspace *work, QP *qp){
  int error_flag,update_mask=0;

  work->n = qp->n;
  work->m = qp->m;
  work->ms = qp->ms;
  work->qp = qp;
 
  // Allocate memory for Rinv and M 
  if(qp->H!=NULL){ 
	work->Rinv = malloc(((qp->n+1)*qp->n/2)*sizeof(c_float));
	work->M = malloc(qp->n*(qp->m-qp->ms)*sizeof(c_float));
	update_mask += UPDATE_Rinv+UPDATE_M;
  }
  else{// H = I =>  no need to transform H->Rinv and M->A 
	work->Rinv = NULL;
	work->M = qp->A;
  }

  // Allocate memory for d and v 
  if(qp->f!=NULL || work->settings->eps_prox != 0){
	work->dupper = malloc(qp->m*sizeof(c_float));
	work->dlower = malloc(qp->m*sizeof(c_float));
	work->v = malloc(qp->n*sizeof(c_float));
	update_mask+=UPDATE_v+UPDATE_d;
  }
  else{ // f = 0 => no need to transform f->v and b->d 
	work->v= NULL;
	work->dupper = qp->bupper; 
	work->dlower = qp->blower; 
  }
  
  // Allocate memory for local constraint states
  work->sense = malloc(qp->m*sizeof(int));
  update_mask += UPDATE_sense;
	
  // Form LDP
  error_flag = update_ldp(update_mask, work);
  if(error_flag<0){
	free_daqp_ldp(work);
	return error_flag;
  }
  return 1;
}

// Free data for LDP 
void free_daqp_ldp(Workspace *work){
  if(work->sense==NULL) return; // Already freed
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
  work->sense = NULL;
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
