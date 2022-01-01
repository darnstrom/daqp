#include <stdio.h>
#include <stdlib.h>
#include "daqp.h" 

int daqp(Workspace *work){
  while(1){
	if(work->iterations++>work->settings->iter_limit)
	  return EXIT_ITERLIMIT;
	if(work->sing_ind==EMPTY_IND){ 
	  compute_CSP(work);
	  if(work->fval > work->fval_bound)
		return EXIT_INFEASIBLE;
	  if(work->cycle_counter > work->settings->cycle_tol){
		if(work->tried_repair==1)
		  return EXIT_CYCLE;
		else{
		  // Cycling -> Try to reorder and refactorize LDL
		  reorder_LDL(work);
		  warmstart_workspace(work, work->WS,work->n_active);
		  work->tried_repair=1;
		  continue;
		}

	  }
	  // Check dual feasibility of CSP
	  find_blocking_constraints(work);
	  if(work->n_blocking==0){ //lam_star >= 0 (i.e., dual feasible)
		find_constraint_to_add(work);
		if(work->add_ind == EMPTY_IND){ //mu >= (i.e., primal feasible)
		  // All KKT-conditions satisfied -> optimum found 
		  // LDP-> QP solution (x* stored in u)
		  ldp2qp_solution(work->x,work->R,work->u,work->v,work->n); 
		  return EXIT_OPTIMAL; 
		}
		else{
		  // Set lam = lam_star
		  work->swp_pointer=work->lam;
		  work->lam = work->lam_star;
		  work->lam_star=work->swp_pointer;

		  add_constraint(work);
		}
	  }
	  else{// Blocking constraints -> remove constraint from working set
		for(int i=0; i<work->n_active;i++)// p = lam^*-lam (stored in lam^*)
		  work->lam_star[i]-=work->lam[i];
		compute_alpha_and_rm_blocking(work);
		remove_constraint(work);
	  }
	}
	else{// Singular case
	  compute_singular_direction(work);
	  find_blocking_constraints(work);
	  if(work->n_blocking==0){
		// Infeasible problem
		return EXIT_INFEASIBLE;
	  }
	  else{
		compute_alpha_and_rm_blocking(work);
		work->sing_ind=EMPTY_IND;
		remove_constraint(work);
	  }
	}
  }
}

// Compute x = -R\(u+v)
void ldp2qp_solution(double *x, double *R, double *u, double *v, int nx){
  int i,j,disp;
  for(i=0;i<nx;i++)
	x[i] = -(u[i]+v[i]);
  if(R != NULL) // Backwards substitution (Skip if LP since R = I)
	for(i=nx-1,disp=(nx+1)*nx/2-1;i>=0;i--){
	  for(j=nx-1;j>i;j--)
		x[i]-=R[disp--]*x[j];
	  x[i]*=R[disp--];
	}
}

void warmstart_workspace(Workspace* work, int* WS, const int n_active){
  // TODO, will probably be error with equality constraints here... (Make sure reorder always adds inequality constraints...)
  reset_daqp_workspace(work); // Reset workspace
  for(int i = 0; i<n_active; i++){
	if(work->sing_ind!=EMPTY_IND){ //If  
	  work->add_ind = WS[i];
	  add_constraint(work);
	  work->lam[i] = 1;
	}else{ //Make sure that the unadded constraints are inactive in sense
	  SET_INACTIVE(work->WS[i]);
	}
  }
}

// Allocate memory for iterates  
void allocate_daqp_workspace(Workspace *work, int n){
  work->n = n;
  work->R = NULL;
  work->v = NULL;

  work->lam = malloc((n+1)*sizeof(c_float));
  work->lam_star = malloc((n+1)*sizeof(c_float));
  
  work->WS= malloc((n+1)*sizeof(int));
  work->BS= malloc((n+1)*sizeof(int));
  
  work->L= malloc(((n+1)*(n+2)/2)*sizeof(c_float));
  work->D= malloc((n+1)*sizeof(c_float));

  work->xldl= malloc((n+1)*sizeof(c_float));
  work->zldl= malloc((n+1)*sizeof(c_float));
  
  work->u= malloc(n*sizeof(c_float));
  work->x = work->u; 

  reset_daqp_workspace(work);
}


// Free memory for iterates
void free_daqp_workspace(Workspace *work){
  free(work->lam);
  free(work->lam_star);

  free(work->WS);
  free(work->BS);

  free(work->L);
  free(work->D);

  free(work->xldl);
  free(work->zldl);

  free(work->u);
}

// Reset workspace to default values
void reset_daqp_workspace(Workspace *work){
  work->iterations=0;
  work->sing_ind=EMPTY_IND;
  work->add_ind=EMPTY_IND;
  work->rm_ind=EMPTY_IND;
  work->n_active =0;
  work->n_blocking=0;
  work->reuse_ind=0;
  work->cycle_counter=0;
  work->tried_repair=0;
  work->fval= -1;
  work->fval_bound= INF;
}

// Reset workspace in a state that allows for warmstarting
void reset_daqp_workspace_warm(Workspace *work){
  work->iterations=0;
  work->reuse_ind=0;
  work->fval= -1;
}

// Allocate memory for problem data
void allocate_daqp_ldp(Workspace *work,int n, int m){
  work->R = malloc(((n+1)*n/2)*sizeof(c_float));
  work->M = malloc(n*m*sizeof(c_float));
  work->dupper = malloc(m*sizeof(c_float));
  work->dlower = malloc(m*sizeof(c_float));
  work->v = malloc(n*sizeof(c_float));
}

// Free data for problem data
void free_daqp_ldp(Workspace *work){
  free(work->R);
  free(work->M);
  free(work->dupper);
  free(work->dlower);
  free(work->v);
}
