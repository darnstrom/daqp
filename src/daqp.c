#include "daqp.h" 


int daqp(Workspace *work){
  c_float *swp_ptr;
  int tried_repair=0, cycle_counter=0;
  c_float best_fval = -1;
  while(1){
	/* CHECK TERMINAL CONDITIONS */
	if(++work->iterations>work->settings->iter_limit) return EXIT_ITERLIMIT;
	if(work->sing_ind==EMPTY_IND){ 
	  compute_CSP(work);
	  // Check dual feasibility of CSP
	  if(!remove_blocking(work)){ //lam_star >= 0 (i.e., dual feasible)
		compute_primal_and_fval(work);
		if(!add_infeasible(work)){ //mu >= (i.e., primal feasible)
		  // All KKT-conditions satisfied -> optimum found 
		  ldp2qp_solution(work); 
		  if(work->soft_slack > work->settings->primal_tol) 
			return EXIT_SOFT_OPTIMAL; 
		  else
			return EXIT_OPTIMAL; 
		}
		else{
		  // Set lam = lam_star
		  swp_ptr=work->lam; work->lam = work->lam_star; work->lam_star=swp_ptr;
		}

		/* Check fval terminal conditions */
		if(best_fval > work->fval+work->settings->progress_tol){ 
		  if(cycle_counter++ > work->settings->cycle_tol && tried_repair++ == 1) 
			return EXIT_CYCLE;
		  else{// Cycling -> Try to reorder and refactorize LDL
			reorder_LDL(work);
			warmstart_workspace(work, work->WS,work->n_active);
		  }
		}
		else{ // Progress was made
		  best_fval = work->fval;
		  cycle_counter = 0;
		  if(best_fval > work->fval_bound) return EXIT_INFEASIBLE;
		}
	  }
	}
	else{// Singular case
	  compute_singular_direction(work);
	  if(!remove_blocking(work)) return EXIT_INFEASIBLE;
	}
  }
}

// Compute x = -R\(u+v)
void ldp2qp_solution(Workspace *work){
  int i,j,disp;
  // x* = Rinv*(u-v)
  for(i=0;i<work->n;i++)
	work->x[i]=work->u[i]-work->v[i];
  if(work->Rinv != NULL) // (Skip if LP since R = I)
	for(i=0,disp=0;i<work->n;i++){
	  work->x[i]*=work->Rinv[disp++];
	  for(j=i+1;j<work->n;j++)
		work->x[i]+=work->Rinv[disp++]*work->x[j];
	}
}

void warmstart_workspace(Workspace* work, int* WS, const int n_active){
  // TODO, will probably be error with equality constraints here... (Make sure reorder always adds inequality constraints...)
  reset_daqp_workspace(work); // Reset workspace
  for(int i = 0; i<n_active; i++){
	if(work->sing_ind!=EMPTY_IND){ //If  
	  add_constraint(work,WS[i],0); // TODO: replace 0 with upper/lower...
	  work->lam[i] = 1;
	}else{ //Make sure that the unadded constraints are inactive in sense
	  SET_INACTIVE(work->WS[i]);
	}
  }
}


// Reset workspace to default values
void reset_daqp_workspace(Workspace *work){
  work->iterations=0;
  work->inner_iter=0;
  work->outer_iter=0;
  work->sing_ind=EMPTY_IND;
  work->n_active =0;
  work->reuse_ind=0;
  work->fval= -1;
  work->fval_bound= INF;
}

// Reset workspace in a state that allows for warmstarting
void reset_daqp_workspace_warm(Workspace *work){
  work->iterations=0;
  work->reuse_ind=0;
  work->fval= -1;
}
