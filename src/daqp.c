#include "daqp.h" 

int daqp_ldp(DAQPWorkspace *work){
  int exitflag=EXIT_ITERLIMIT,iter;
  int tried_repair=0, cycle_counter=0;
  c_float best_fval = -1;
  for(iter=1; iter < work->settings->iter_limit; ++iter){
	if(work->sing_ind==EMPTY_IND){ 
	  compute_CSP(work);
	  // Check dual feasibility of CSP
	  if(!remove_blocking(work)){ //lam_star >= 0 (i.e., dual feasible)
		compute_primal_and_fval(work);
		if(!add_infeasible(work)){ //mu >= (i.e., primal feasible)
		  // All KKT-conditions satisfied -> optimum found 
		  if(work->soft_slack > work->settings->primal_tol)
			exitflag = EXIT_SOFT_OPTIMAL; 
		  else
			exitflag = EXIT_OPTIMAL;
		  break;
		}

		/* Check fval terminal conditions */
		if(best_fval > work->fval-work->settings->progress_tol){ 
		  if(cycle_counter++ > work->settings->cycle_tol){
			if(tried_repair == 1){
			  exitflag = EXIT_CYCLE;
			  break;
			}
			else{// Cycling -> Try to reorder and refactorize LDL
			  tried_repair =1;
			  reset_daqp_workspace(work);
			  activate_constraints(work);
			  //reorder_LDL(work);
			  //warmstart_workspace(work, work->WS,work->n_active);
			  cycle_counter=0;
			  best_fval = -1;
			}
		  }
		}
		else{ // Progress was made
		  best_fval = work->fval;
		  cycle_counter = 0;
		  if(best_fval > work->settings->fval_bound){ 
			exitflag = EXIT_INFEASIBLE;
			break;
		  }
		}
	  }
	}
	else{// Singular case
	  compute_singular_direction(work);
	  if(!remove_blocking(work)){ 
		exitflag = EXIT_INFEASIBLE;
		break;
	  }
	}
  }
  // Finalize result before returning
  work->iterations = iter;
  return exitflag;
}

// Compute x = -R\(u+v)
void ldp2qp_solution(DAQPWorkspace *work){
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

void warmstart_workspace(DAQPWorkspace* work, int* WS, const int n_active){
  // TODO, will probably be error with equality constraints here... (Make sure reorder always adds inequality constraints...)
  reset_daqp_workspace(work); // Reset workspace
  for(int i = 0; i<n_active; i++){
	if(work->sing_ind!=EMPTY_IND){
	  add_constraint(work,WS[i],1.0); 
	}else{ //Make sure that the unadded constraints are inactive in sense
	  SET_INACTIVE(work->WS[i]);
	}
  }
}


// Reset workspace to default values
void reset_daqp_workspace(DAQPWorkspace *work){
  work->sing_ind=EMPTY_IND;
  work->n_active =0;
  work->reuse_ind=0;
}
