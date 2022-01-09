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

// Gradient step
// TODO: could probably reuse code from daqp 
int gradient_step(Workspace* work){
  int j,k,disp,add_ind=EMPTY_IND;
  int add_isupper;
  const int nx=work->n;
  const int m=work->m;
  const int ms=work->ms;
  c_float Ax,delta_s, min_alpha=INF;
  // Find constraint j to add: j =  argmin_j s_j 
  // Simple bounds 
  for(j=0, disp=0;j<ms;j++){
	if(IS_ACTIVE(work->sense[j])) continue;
	delta_s = work->x[j]-work->xold[j];
	if(delta_s>0 && //Feasible descent direction
	   work->qp->bupper[j]<INF && // Not single-sided
	   work->qp->bupper[j]-work->x[j]<min_alpha*delta_s){
	  add_ind = j;
	  add_isupper = 1; 
	  min_alpha = (work->qp->bupper[j]-work->x[j])/delta_s;
	}
	else if(delta_s < 0 && //Feasible descent direction
			work->qp->blower[j]>-INF && // Not single-sided
			work->qp->blower[j]-work->x[j]>min_alpha*delta_s){
	  add_ind = j;
	  add_isupper=0;
	  min_alpha = (work->qp->blower[j]-work->x[j])/delta_s;
	}
  }
  //General bounds
  for(j=ms, disp=0;j<m;j++){
	if(IS_ACTIVE(work->sense[j])){
	  disp+=nx;// Skip ahead in A 
	  continue;
	}
	//delta_s[j] = A[j,:]*delta_x
	for(k=0,delta_s=0,Ax=0;k<nx;k++){ // compute s = A(x-xold) and Ax
	  Ax += work->M[disp]*work->x[k];
	  delta_s-=work->M[disp++]*work->xold[k]; 
	}
	delta_s +=Ax;

	if(delta_s>0 && // Feasible descent direction
	   work->qp->bupper[j]<INF && // Not single-sided
	   work->qp->bupper[j]-Ax < delta_s*min_alpha){
	  add_ind = j;
	  add_isupper=1;
	  min_alpha=(work->qp->bupper[j]-Ax)/delta_s;
	}
	else if(delta_s<0 && // Feasible descent direction
			work->qp->blower[j]>-INF && // Not single-sided
			work->qp->blower[j]-Ax > delta_s*min_alpha){
	  add_ind = j;
	  add_isupper=0;
	  min_alpha=(work->qp->blower[j]-Ax)/delta_s;
	}
  }
  if(add_ind == EMPTY_IND) return EMPTY_IND;
  // TODO Maybe don't add to AS?
  work->add_ind = add_ind;
  work->add_isupper= add_isupper;
  add_constraint(work); // Update working set and LDL'
  // TODO: Won't this be handled inside daqp anyways? 
  if(work->sing_ind!=EMPTY_IND){// Remove constraint if basis becomes singular 
	remove_constraint(work,work->sing_ind);
	work->sing_ind = EMPTY_IND;
  }
  else{
	for(k=0;k<nx;k++) // x <-- x+alpha deltax 
	  work->x[k]+=min_alpha*(work->x[k]-work->xold[k]);
  }
  return add_ind;
}
