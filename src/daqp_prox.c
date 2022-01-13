#include <stdio.h>
#include <stdlib.h>
#include "daqp_prox.h"
#include "utils.h"

int daqp_prox(Workspace *work){
  int i,total_iter=0;
  const int nx=work->n;
  int exitflag;
  c_float *swp_ptr;
  c_float diff,eps=work->settings->eps_prox;

  while(total_iter  <  work->settings->iter_limit){
	// xold <-- x
	swp_ptr = work->xold; work->xold = work->x; work->x = swp_ptr;
	
	// ** Solve least-distance problem **
	work->u = work->x;
	exitflag = daqp_ldp(work);
	
	total_iter += work->iterations;
	if(exitflag<0) 
	  return exitflag; // Could not solve LDP -> return
	else 
	 ldp2qp_solution(work); // Get qp solution 

	if(eps==0) break; // No regularization -> no outer iterations
	
	// ** Check convergence **
	if(work->iterations==1){ // No changes to the working set 
	  for(i=0, diff= 0;i<nx;i++){ // ||x_old - x|| > eta  ?
		diff= work->x[i] - work->xold[i];
		if((diff> work->settings->eta_prox) || (diff< -work->settings->eta_prox)) break;
	  }
	  if(i==nx){
		exitflag = EXIT_OPTIMAL; // Fix point reached
		break;
	  }
	}
	
	// ** Perturb problem **
	// Compute v = R'\(f-eps*x) (FWS Skipped if LP since R = I) 
	if(work->Rinv== NULL){ 
	  eps*= (work->iterations==1) ? 10 : 0.9; // Adapt epsilon TODO: add to settings 
	  for(i = 0; i<nx;i++) 
		work->v[i] = work->qp->f[i]*eps-work->x[i];
	}
	else{
	  for(i = 0; i<nx;i++) 
		work->v[i] = work->qp->f[i]-eps*work->x[i];
	  update_v(work->v,work); // 
	}
	// Perturb RHS of constraints 
	update_d(work);
  }
  // Finalize results
  if(total_iter >= work->settings->iter_limit) exitflag = EXIT_ITERLIMIT; 
  work->iterations = total_iter;
  return exitflag;
}
