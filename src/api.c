#include "api.h" 
#include <stdlib.h>
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

// Compute minimal Hrepresentation of Mx<=d
// Redundant constraints are marked in sense work->sense 
int minimal_Hrepresentation(Workspace *work){
  int i,j,exitflag,N_nonred= work->m;

  for(i=0; i<N_CONSTR;i++){
	// Add equality constraint
	work->add_ind = i;
	update_LDL_add(work);
	work->WS[work->n_active] = i;
	work->n_active++;
	work->sense[i] = EQUALITY;

	// Solve feasibility problem
	exitflag = daqp(work);
	//Clear active cosntraint
	for(j=0;j<work->n_active;j++)
	  work->sense[work->WS[j]]=INACTIVE_INEQUALITY; // TODO handle equality constraints...
	// constraint i redundant if LDP was infeasible
	if(exitflag==EXIT_INFEASIBLE){
	  N_nonred--;
	  work->sense[i] = FREE_CONSTRAINT; 
	}
	// Reset working set
	reset_daqp_workspace(work);
  }
  return N_nonred; 
}

// API to call minHRep... with A,b
int daqp_minHRep(Workspace* work, c_float* A, c_float*b, int *sense, const int m){
  reset_daqp_workspace(work); // Reset workspace
  work->m = m;
  work->M = A;
  work->d = b;
  work->sense = sense;
  return minimal_Hrepresentation(work);
}

// LP API
void daqp_lp_setup(ProxWorkspace** pws_ptr, double* f, double* A,double* b, int n, int m){
  ProxWorkspace* pws=malloc(sizeof(ProxWorkspace));
  pws->Rinv = NULL;
  pws->epsilon=1;
  pws->f = f;
  pws->M = A;
  pws->b = b;
  allocate_prox_workspace(pws, n,m);
  *pws_ptr = pws;
}

c_float daqp_lp_solve(ProxWorkspace* prox_work, c_float* f, c_float* A, c_float* b, const int m){
  // Solve min_x f' x s.t. A x <= b 

  reset_daqp_workspace(prox_work->work); // Reset daqp workspace
  reset_prox_workspace(prox_work); // Reset prox workspace
  
  prox_work->m = m;
  prox_work->work->m = m;
  prox_work->f = f;
  prox_work->M = A;
  prox_work->b = b;
 
  daqp_prox(prox_work);
  
  // Return function value 
  c_float fval=0;
  for(int i = 0;i<prox_work->n;i++)
	fval += f[i]*prox_work->x[i];
  return fval;
}
