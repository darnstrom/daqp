#include <stdio.h>
#include <stdlib.h>
#include "daqp_prox.h"
#include "utils.h"

int daqp_prox(Workspace *work){
  int i,cycle_counter=0;
  const int nx=work->n;
  int exitflag;
  c_float diff,fval=INF, *swp_ptr;

  while(work->outer_iter++  <  work->settings->prox_iter_limit){

	// xold <-- x
	swp_ptr = work->xold; work->xold = work->x; work->x = swp_ptr;
	
	// ** Solve least-distance problem **
	work->u = work->x;
	reset_daqp_workspace_warm(work);
	exitflag = daqp(work);
	
	work->inner_iter+=work->iterations;
	if(exitflag!=EXIT_OPTIMAL) return exitflag;
	if(work->settings->eps_prox == 0){
	  if(work->soft_slack > work->settings->primal_tol) return EXIT_SOFT_OPTIMAL;
	  return EXIT_OPTIMAL; 
	}
	
	// ** Check convergence **
	if(work->iterations==1 && work->outer_iter&1){ // No changes to the working set 
	  for(i=0, diff= 0;i<nx;i++){ // ||x_old - x|| > eta  ?
		diff= work->x[i] - work->xold[i];
		if((diff> work->settings->eta_prox) || (diff< -work->settings->eta_prox)) break;
	  }
	  if(i==nx) return EXIT_OPTIMAL; // Fix point reached
	  // Take gradient step if LP (and we are not constrained to a vertex) 
	  if((work->Rinv == NULL)&&(work->n_active != work->n)){ 
		if(gradient_step(work)==EMPTY_IND) return EXIT_UNBOUNDED;
	  }
	}
	// Compute objective function value to detect progress
	// (TODO: this is currently only valid for LPs...)
	if(work->Rinv == NULL){
	  for(i=0, diff=fval;i<nx;i++)
		diff-=work->qp->f[i]*work->x[i];
	  if(diff<0){ 
		if(cycle_counter++ > work->settings->cycle_tol) return EXIT_OPTIMAL;
	  }
	  else{ // Progress -> update objective function value
		fval -=diff ;
		cycle_counter=0;
	  }
	}
	
	// ** Perturb problem **
	// Compute v = R'\(f-eps*x) (FWS Skipped if LP since R = I) 
	if(work->Rinv== NULL) 
	  for(i = 0; i<nx;i++) 
		work->v[i] = work->qp->f[i]-work->x[i];
	else{
	  for(i = 0; i<nx;i++) 
		work->v[i] = work->qp->f[i]-work->settings->eps_prox*work->x[i];
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
  const int nx=work->n;
  const int m=work->m;
  c_float delta_s,alpha, min_alpha=INF;
  // Find constraint j to add: j =  argmin_j s_j 
  for(j=0, disp=0;j<m;j++){
	if(IS_ACTIVE(work->sense[j])){
	  disp+=nx;// Skip ahead in A 
	  continue;
	}
	//delta_s[j] = A[j,:]*delta_x
	delta_s= 0; 
	for(k=0;k<nx;k++) // 
	  delta_s+=work->M[disp++]*(work->x[k]-work->xold[k]); 
	if(delta_s>0){
	  // Compute alphaj = (b[j]-A[j,:]*x]/delta_s
	  alpha = work->qp->bupper[j];
	  // TODO: ALSO include blower here...
	  // (Note that bounded constraints->always bounded)
	  for(k=0, disp-=nx;k<nx;k++) // 
		alpha-=work->M[disp++]*work->x[k];
	  alpha/=delta_s;
	  if(alpha<min_alpha&&alpha>0){
		add_ind = j;
		min_alpha= alpha;
	  }
	}
  }
  if(add_ind == EMPTY_IND) return EMPTY_IND;
  // Maybe don't add to AS?
  work->add_ind = add_ind;
  add_constraint(work); // Update working set and LDL'
  if(work->sing_ind!=EMPTY_IND){// Remove constraint if basis becomes singular 
	work->rm_ind = work->sing_ind;
	remove_constraint(work);
	work->sing_ind = EMPTY_IND;
  }
  else{
	for(k=0;k<nx;k++) // x <-- x+alpha deltax 
	  work->x[k]+=min_alpha*work->xold[k];
  }
  return add_ind;
}
