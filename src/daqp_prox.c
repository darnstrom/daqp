#include <stdio.h>
#include <stdlib.h>
#include "daqp_prox.h"

int daqp_prox(ProxWorkspace *prox_work){
  int i,j,disp;
  const int nx=prox_work->n;
  int exitflag,fixpoint;
  c_float sum, *swp_pointer;
  Workspace *work = prox_work->work;

  while(prox_work->outer_iterations++<PROX_ITER_LIMIT){

	// ** Perturb problem **

	// Compute v = R'\(f-eps*x) (FWS Skipped if LP since R = I) 
	if(work->R== NULL) 
	  for(i = 0; i<nx;i++) 
		work->v[i] = prox_work->f[i]-prox_work->x[i];
	else{
	  for(i = 0; i<nx;i++) 
		work->v[i] = prox_work->f[i]-prox_work->epsilon*prox_work->x[i];
	  for(i = 0,disp=0; i<nx;i++){
		work->v[i]*=work->R[disp++];
		for(j=i+1;j<nx;j++)
		  work->v[j] -= work->R[disp++]*work->v[i];
	  }
	} 

	// Perturb d (d = b+M*v)
	for(i = 0,disp=0; i<prox_work->m;i++){
	  sum = prox_work->b[i];
	  for(j=0; j<nx;j++)
		sum += work->M[disp++]*work->v[j]; 
	  work->d[i] = sum;
	}

	// xold <-- x
	swp_pointer = prox_work->xold;
	prox_work->xold = prox_work->x;
	prox_work->x = swp_pointer;
	
	
	// ** Solve least-distance problem **
	reset_daqp_workspace_warm(work);
	work->x = prox_work->x;
	exitflag = daqp(work);
	
	prox_work->inner_iterations+=work->iterations;
	if(exitflag!=EXIT_OPTIMAL) return exitflag;
	if(prox_work->epsilon == 0) return EXIT_OPTIMAL; // No regularization -> optimal solution
	
	//Check convergence
	if(work->iterations==1 && prox_work->outer_iterations%2){ // No changes to the working set 
	  fixpoint = 1;
	  for(i=0;i<nx;i++){
		prox_work->xold[i]= prox_work->x[i] - prox_work->xold[i];
		if((prox_work->xold[i]> ETA) || (prox_work->xold[i]< -ETA)) // ||x_old - x|| > eta 
		  fixpoint = 0;
	  }
	  if(fixpoint==1) return EXIT_OPTIMAL; // Fix point reached
	  // Take gradient step if LP (and we are not constrained to a vertex) 
	  if((work->R == NULL)&&(work->n_active != work->n)){ 
		if(gradient_step(prox_work,work)==EMPTY_IND) return EXIT_UNBOUNDED;
	  }
	}
	// Compute objective function value to detect progress
	sum = 0;
	for(i=0;i<nx;i++)
	  sum+=prox_work->f[i]*prox_work->x[i];
	if(sum>prox_work->fval){ 
	  if(prox_work->cycle_counter++ > 10)
		return EXIT_OPTIMAL; // No progress -> fix-point
	}
	else{ // Progress -> update objective function value
	  prox_work->fval = sum;
	  prox_work->cycle_counter=0;
	}
  }
  return EXIT_ITERLIMIT;
}

// Gradient step
// TODO: could probably reuse code from daqp 
int gradient_step(ProxWorkspace* prox_work, Workspace* work){
  int j,k,disp,add_ind=EMPTY_IND;
  const int nx=prox_work->n;
  c_float delta_s,alpha, min_alpha=INF;
  // Find constraint j to add: j =  argmin_j s_j 
  for(j=0, disp=0;j<prox_work->m;j++){
	if(work->sense[j]!=INACTIVE_INEQUALITY){
	  disp+=nx;// Skip ahead in A 
	  continue;
	}
	//delta_s[j] = A[j,:]*delta_x
	delta_s= 0; 
	for(k=0;k<nx;k++) // 
	  delta_s+=work->M[disp++]*prox_work->xold[k];
	if(delta_s>0){
	  // Compute alphaj = (b[j]-A[j,:]*x]/delta_s
	  alpha = prox_work->b[j];
	  for(k=0, disp-=nx;k<nx;k++) // 
		alpha-=work->M[disp++]*prox_work->x[k];
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
	prox_work->x[k]+=min_alpha*prox_work->xold[k];
  }
  return add_ind;
}
		
// Utils
void allocate_prox_workspace(ProxWorkspace *prox_work, int n, int m){
  
  prox_work->x= malloc(n*sizeof(c_float));
  prox_work->xold= malloc(n*sizeof(c_float));

  prox_work->n=n;
  prox_work->m=m;
  reset_prox_workspace(prox_work);

  // Setup daqp workspace
  prox_work->work=malloc(sizeof(Workspace));
  allocate_daqp_workspace(prox_work->work, n);
  prox_work->work->m=m;
  prox_work->work->d= malloc(m*sizeof(c_float));
  prox_work->work->v= malloc(n*sizeof(c_float));
}

void free_prox_workspace(ProxWorkspace *prox_work){
  free(prox_work->x);
  free(prox_work->xold);

  free_daqp_workspace(prox_work->work);
  free(prox_work->work->d);
  free(prox_work->work->v);
  free(prox_work->work);
} 

void reset_prox_workspace(ProxWorkspace *prox_work){
  for(int i=0;i<prox_work->n;i++) // x0 = 0
	prox_work->x[i] = 0; 
  prox_work->inner_iterations=0;
  prox_work->outer_iterations=0;
  prox_work->fval=INF;
  prox_work->cycle_counter=0;
}

