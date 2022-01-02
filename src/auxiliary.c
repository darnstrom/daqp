#include "auxiliary.h"
#include "factorization.h"
void remove_constraint(Workspace* work){
  int i;
  update_LDL_remove(work);
  
  // Update data structures
  SET_INACTIVE(work->WS[work->rm_ind]); 
  (work->n_active)--;
  for(i=work->rm_ind;i<work->n_active;i++){
	work->WS[i] = work->WS[i+1]; 
	work->lam[i] = work->lam[i+1]; 
  }
  // Can only reuse work less than the ind that was removed 
  if(work->rm_ind < work->reuse_ind)
	work->reuse_ind = work->rm_ind;
  
  // Pivot for improved numerics
  pivot_last(work);
}
// Maybe take add_ind as input instead?
void add_constraint(Workspace *work){
  update_LDL_add(work);
  // Update data structures  
  SET_ACTIVE(work->add_ind);
  if(work->add_isupper)
	SET_UPPER(work->add_ind);
  else
	SET_LOWER(work->add_ind);
  work->WS[work->n_active] = work->add_ind;
  work->lam[work->n_active] = 0;
  work->n_active++;

  // Pivot for improved numerics
  pivot_last(work);
}

void find_constraint_to_add(Workspace *work){
  int i,j,k,disp;
  c_float min_val = -work->settings->primal_tol;
  c_float Mu;
  int isupper=0, add_ind=EMPTY_IND;
  // Reset u
  for(j=0;j<NX;j++)
	work->u[j]=0;
  //u[m] <-- Mk'*lam_star (zero if empty set)
  for(i=0;i<work->n_active;i++){
	if(IS_SIMPLE(work->WS[i])){
	  // Simple constraint 
	  if(work->Rinv!=NULL){ // Hessian is not identity
	  for(j=work->WS[i], disp=R_OFFSET(work->WS[i],NX);j<NX;j++)
		work->u[j]-=work->Rinv[disp++]*work->lam_star[i];
	  }
	  else work->u[j]-=work->lam_star[work->WS[i]]; // Hessian is identity
	}
	else{ 
	  // General constraint
	for(j=0,disp=NX*(work->WS[i]-N_SIMPLE);j<NX;j++)
	  work->u[j]-=work->M[disp++]*work->lam_star[i];
	}
  }
  // Check for progress 
  c_float fval_diff=work->fval;
  for(j=0;j<NX;j++)
	fval_diff-=work->u[j]*work->u[j];
  if(fval_diff<-work->settings->progress_tol){
	work->fval -= fval_diff;
	work->cycle_counter=0;
  }else work->cycle_counter++;
  
  // If empty working mu = d
  if(work->n_active==0){  
	for(j=0;j<N_CONSTR;j++)
	  if(work->dupper[j]<min_val){//dupper
		add_ind = j; isupper = 1; 
		min_val = work->dupper[j]; 
	  }
	  else if(-work->dlower[j]<min_val){//dlower
		add_ind = j; isupper = 0;
		min_val = -work->dlower[j];
	  }
  }
  else{// Non-empty working set 
	for(j=0, disp=0;j<N_SIMPLE;j++){
	  if(IS_ACTIVE(j)){
		disp+=NX-j;
		continue;
	  }
	  if(work->Rinv==NULL){// Hessian is identify
		Mu=work->u[i]; 
		disp+=NX-j;
	  }
	  else{
		for(k=j,Mu=0;k<NX;k++) // 
		  Mu+=work->Rinv[disp++]*work->u[k];
	  }
	  if(work->dupper[j]-Mu<min_val){
		add_ind = j; isupper = 1;
		min_val = work->dupper[j]-Mu;
	  }
	  else if(-(work->dlower[j]-Mu)<min_val){
		add_ind = j; isupper = 0;
		min_val = -(work->dlower[j]-Mu);
	  }
	}
	/* General two-sided constraints */
	for(j=N_SIMPLE, disp=0;j<N_CONSTR;j++){
	  if(IS_ACTIVE(j)){
		disp+=NX;// Skip ahead in M
		continue;
	  }
	  for(k=0,Mu=0;k<NX;k++) 
		Mu+=work->M[disp++]*work->u[k];

	  if(work->dupper[j]-Mu<min_val){
		add_ind = j; isupper = 1;
		min_val = work->dupper[j]-Mu;
	  }
	  else if(-(work->dlower[j]-Mu)<min_val){
		add_ind = j; isupper = 0;
		min_val = -(work->dlower[j]-Mu);
	  }
	}
  }
  work->add_ind = add_ind;
  work->add_isupper = isupper;
}
void find_blocking_constraints(Workspace *work){
  c_float violation;
  work->n_blocking=0;
  for(int i=0;i<work->n_active;i++){
	if(IS_IMMUTABLE(work->WS[i])) continue;
	if(IS_LOWER(work->WS[i])) 
	  violation = work->lam_star[i]-work->settings->dual_tol;
	else
	  violation = -work->lam_star[i]-work->settings->dual_tol;
	if(violation>0)
	  work->BS[work->n_blocking++]=i;
  }
}
void compute_alpha_and_rm_blocking(Workspace *work){
  // p is stored in lam_star...
  int i;
  c_float alpha_cand; 
  c_float alpha=-work->lam[work->BS[0]]/work->lam_star[work->BS[0]]; 
  work->rm_ind = work->BS[0];
  for(i=1;i<work->n_blocking;i++){
	alpha_cand = -work->lam[work->BS[i]]/work->lam_star[work->BS[i]];
	if(alpha_cand < alpha){
	  alpha = alpha_cand;
	  work->rm_ind =work->BS[i]; 
	}
  }
  // lam = lam + alpha p
  for(i=0;i<work->n_active;i++)
	work->lam[i]+=alpha*work->lam_star[i];
}
void compute_CSP(Workspace *work){
  int i,j,disp,start_disp;
  double sum;
  // Forward substitution (xi <-- L\d)
  for(i=work->reuse_ind,disp=ARSUM(work->reuse_ind); i<work->n_active; i++){
	// Setup RHS
	if(IS_LOWER(work->WS[i])) 
	  sum = -work->dlower[work->WS[i]];
	else 
	  sum = -work->dupper[work->WS[i]];
	for(j=0; j<i; j++)
	  sum -= work->L[disp++]*work->xldl[j];
	disp++; //Skip 1 in L 
	work->xldl[i] = sum;
  }
  // Scale with D  (zi = xi/di)
  for(i=work->reuse_ind; i<work->n_active; i++)
	work->zldl[i] = work->xldl[i]/work->D[i];
  //Backward substitution  (lam_star <-- L'\z)
  start_disp = ARSUM(work->n_active)-1;
  for(i = work->n_active-1;i>=0;i--){
	sum=work->zldl[i];
	disp = start_disp--;
	for(j=work->n_active-1;j>i;j--){
	  sum-=work->lam_star[j]*work->L[disp];
	  disp-=j;
	} 
	work->lam_star[i] = sum;
  }
  
  work->reuse_ind = work->n_active; // Save forward substitution information 
}

//TODO this could probably be directly calculated in L
void compute_singular_direction(Workspace *work){
  // Step direction is stored in lam_star
  int i,j,disp,offset_L= ARSUM(work->sing_ind);
  int start_disp= offset_L-1;

  // Backwards substitution (p_tidle <-- L'\(-l))
  for(i = work->sing_ind-1;i>=0;i--){
	work->lam_star[i] = -work->L[offset_L+i];
	disp = start_disp--;
	for(j=work->sing_ind-1;j>i;j--){
	  work->lam_star[i]-=work->lam_star[j]*work->L[disp];
	  disp-=j;
	} 
  }
  work->lam_star[work->sing_ind]=1;

  if(IS_LOWER(work->WS[work->sing_ind])) //Flip to ensure descent direction 
	for(i=0;i<=work->sing_ind;i++)
	  work->lam_star[i] =-work->lam_star[i];

  for(i =work->sing_ind+1;i<work->n_active;i++) //Probably uneccesary...
	work->lam_star[i] = 0;
}


void reorder_LDL(Workspace *work){
  // Extract first column l1,: of  L 
  // and store l1,:^2 in the beginning of L (since L will be overwritten anyways...)  
  // (a large value of l^2 signify linear dependence with the first constraint)
  int i,j,disp;
  c_float swp;
  for(i = 1, disp = 1; i < work->n_active; i++){
	work->L[i] = work->L[disp]*work->L[disp];  
	disp+=i+1;
  }
  // Sort l1,:^2 elements (and reorder the working set accordingly)
  // Bubble sort (use disp for swapping int)
  for(i=work->n_active-1; i>0; i--){
	for(j=1; j<i; j++){

	  if(work->L[j] > work->L[j+1]){
		// Swap
		swp= work->L[j];
		disp = work->WS[j]; 
		work->L[j] = work->L[j+1];
		work->WS[j] = work->WS[j+1];
		work->L[j+1]= swp; 
		work->WS[j+1]= disp; 
	  }
	}
  }
}

void pivot_last(Workspace *work){
  if(work->n_active > 1 && work->D[work->n_active-2] < work->settings->pivot_tol && 
	 work->D[work->n_active-2] < work->D[work->n_active-1]){
	work->rm_ind = work->n_active-2; 
	int ind_old = work->WS[work->rm_ind];
	int old_sense =work->sense[ind_old]; // Make sure that equality constraints are retained... 

	c_float lam_old = work->lam[work->rm_ind];
	remove_constraint(work);

	if(work->sing_ind!=EMPTY_IND) return; // Abort if D becomes singular

	work->add_ind = ind_old;
	add_constraint(work);
	work->sense[ind_old] = old_sense;
	work->lam[work->n_active-1] = lam_old;
  }	
}

void add_equality_constraints(Workspace *work){
  for(int i =0;i<N_CONSTR;i++){
	//TODO prioritize inequalities?
	if(IS_ACTIVE(i)){
	  work->add_ind=i; 
	  update_LDL_add(work);
	  work->WS[work->n_active] = i;
	  work->n_active++;
	  //TODO check singularity...
	  pivot_last(work);
	}
  }
}

void daqp_default_settings(DAQPSettings* settings){
  settings->primal_tol = DEFAULT_PRIM_TOL;
  settings->dual_tol = DEFAULT_DUAL_TOL; 
  settings->zero_tol = DEFAULT_ZERO_TOL;
  settings->pivot_tol = DEFAULT_PIVOT_TOL;
  settings->progress_tol= DEFAULT_PROG_TOL;

  settings->cycle_tol = DEFAULT_CYCLE_TOL;
  settings->iter_limit = DEFAULT_ITER_LIMIT;

  settings->eps_prox = 0;
  settings->eta_prox = DEFAULT_ETA;
  settings->prox_iter_limit = DEFAULT_PROX_ITER_LIMIT; 
}
