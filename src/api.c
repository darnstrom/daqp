#include "api.h" 
#include "utils.h"
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

void daqp_quadprog(DAQPResult *res, QP* qp, DAQPSettings *settings){
  struct timespec tstart,tsetup,tsol;
  // TIC start
  clock_gettime(CLOCK_MONOTONIC, &tstart);
  // Transform QP to ldp (R,v,M,d is stored inplace of H,f,A,b)

  Workspace work;
  printf("Allocate iters (n:%d)\n",qp->n);
  allocate_daqp_workspace(&work,qp->n);
  work.settings = settings;
  printf("Setup LDP\n");
  setup_daqp_ldp(&work,qp);
  work.settings = settings; // TODO unnecessary? 
  printf("Add equalities\n");
  add_equality_constraints(&work);

  clock_gettime(CLOCK_MONOTONIC, &tsetup);
  printf("Solving\n");
  if(settings->eps_prox==0){
	res->exitflag = daqp(&work);

  }
  else{//Prox
	clock_gettime(CLOCK_MONOTONIC, &tsetup);
	res->exitflag = daqp_prox(&work);
	clock_gettime(CLOCK_MONOTONIC, &tsol);
  }
  clock_gettime(CLOCK_MONOTONIC, &tsol);
  
  printf("Extract results\n");
  // Extract results and free workspace 
  daqp_extract_result(res,&work);
  free_daqp_workspace(&work);
  free_daqp_ldp(&work);

  // Add time to result
  res->solve_time = time_diff(tsetup,tsol);
  res->setup_time = time_diff(tstart,tsetup);
}

void setup_daqp(QP* qp, DAQPSettings *settings, Workspace *work){
  // Check if QP is well-posed
  //validate_QP(qp);
  
  // Setup workspace
  allocate_daqp_workspace(work,qp->n);
  setup_daqp_ldp(work,qp);
}

//  Setup LDP from QP  
void setup_daqp_ldp(Workspace *work, QP *qp){
  int i;

  work->n = qp->n;
  work->m = qp->m;
  work->ms = qp->ms;
  work->qp = qp;
  printf("Extract Rinv M\n");
 // Extract data for Rinv and M 
  if(qp->H!=NULL){ 
	work->Rinv = malloc(((qp->n+1)*qp->n/2)*sizeof(c_float));
	work->M = malloc(qp->n*(qp->m-qp->ms)*sizeof(c_float));
  }
  else{// H = I =>  no need to transform H->Rinv and M->A 
	work->Rinv = qp->H;
	work->M = qp->A;
  }

  printf("Extract v and d\n");
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
  printf("Transform QP\n");
  if(qp->H != 0){
	// Copy H->Rinv and A->M
	printf("Copy matrices\n");
	pack_symmetric(qp->H,work->Rinv,qp->n);
	const int nA = qp->n*(qp->m-qp->ms);
	for(i=0; i<nA;i++) work->M[i] = qp->A[i];
	
	// Compute Rinv and M
	printf("Compute Rinv and m\n");
	compute_Rinv_and_M(work->Rinv,work->M,work->settings->eps_prox,qp->n,qp->m-qp->ms);
  }
  // Compute v and M
  update_v_and_d(qp->f,qp->bupper,qp->blower,work);

  // Setup up constraint states
  printf("Setup sense\n");
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
