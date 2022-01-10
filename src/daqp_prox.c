#include <stdio.h>
#include <stdlib.h>
#include "daqp_prox.h"
#include "utils.h"

int daqp_prox(Workspace *work){
  int i;
  const int nx=work->n;
  int exitflag;
  c_float *swp_ptr;
  c_float diff,eps=work->settings->eps_prox;

  while(work->outer_iter++  <  work->settings->prox_iter_limit){
	// xold <-- x
	swp_ptr = work->xold; work->xold = work->x; work->x = swp_ptr;
	
	// ** Solve least-distance problem **
	work->u = work->x;
	reset_daqp_workspace_warm(work);
	exitflag = daqp(work);
	
	work->inner_iter+=work->iterations;
	if(exitflag!=EXIT_OPTIMAL) return exitflag;
	if(eps==0){
	  if(work->soft_slack > work->settings->primal_tol) return EXIT_SOFT_OPTIMAL;
	  return EXIT_OPTIMAL; 
	}
	
	// ** Check convergence **
	if(work->iterations==1){ // No changes to the working set 
	  for(i=0, diff= 0;i<nx;i++){ // ||x_old - x|| > eta  ?
		diff= work->x[i] - work->xold[i];
		if((diff> work->settings->eta_prox) || (diff< -work->settings->eta_prox)) break;
	  }
	  if(i==nx) return EXIT_OPTIMAL; // Fix point reached
	  // Take gradient step if LP (and we are not constrained to a vertex) 
	  if((work->Rinv == NULL)&&(work->n_active != work->n)){ 
		//if(gradient_step(work)==EMPTY_IND) return EXIT_UNBOUNDED;
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
  return EXIT_ITERLIMIT;
}
