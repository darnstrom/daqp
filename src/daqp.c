#include <stdio.h>
#include <stdlib.h>
#include "daqp.h" 

int daqp(Workspace *work){
  while(1){
	if(work->iterations++>MAX_ITER)
	  return EXIT_ITERLIMIT;
	if(work->sing_ind==EMPTY_IND){ 
	  compute_CSP(work);
	  if(work->fval > work->fval_bound)
		return EXIT_INFEASIBLE;
	  if(work->cycle_counter > CYCLE_TOL){
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
		  return EXIT_OPTIMAL; 
		}
		else{
		  // Set lam = lam_star
		  work->swp_pointer=work->lam;
		  work->lam = work->lam_star;
		  work->lam_star=work->swp_pointer;
		  work->swp_pointer = NULL;

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

void warmstart_workspace(Workspace* work, int* WS, const int n_active){
  // TODO, will probably be error with equality constraints here... (Make sure reorder always adds inequality constraints...)
  reset_daqp_workspace(work); // Reset workspace
  for(int i = 0; i<n_active; i++){
	if(work->sing_ind!=EMPTY_IND){ //If  
	  work->add_ind = WS[i];
	  add_constraint(work);
	  work->lam[i] = 1;
	}else{ //Make sure that the unadded constraints are inactive in sense
	  work->sense[work->WS[i]]=INACTIVE_INEQUALITY;
	}
  }
}

// Allocate memory for iterates  
void allocate_daqp_workspace(Workspace *work, int n){
  work->n = n;
  work->lam = malloc((n+1)*sizeof(c_float));
  work->lam_star = malloc((n+1)*sizeof(c_float));
  
  work->WS= malloc((n+1)*sizeof(int));
  work->BS= malloc((n+1)*sizeof(int));
  
  work->L= malloc(((n+1)*(n+2)/2)*sizeof(c_float));
  work->D= malloc((n+1)*sizeof(c_float));

  work->xldl= malloc((n+1)*sizeof(c_float));
  work->zldl= malloc((n+1)*sizeof(c_float));
  
  work->u= malloc(n*sizeof(c_float));
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

  // Important that M and d has been cleared before running this!
  // or that pointers to M and d are retained outside the work-struct.
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
