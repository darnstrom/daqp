#include <stdio.h>
#include <stdlib.h>
#include "daqp.h" 
#include "debug.h"

// Some printing
void printmatrix(c_float* matrix, int row ,int col){
  int n = 0;
  for (int i=0;i<row;i++){
	for(int j=0;j<col;j++){
	  printf("%f ",matrix[n]);
	  n++;
	} 
	printf("\n");		
  }
  printf("\n");		
}
void printvector(c_float* vec, int row){
  for (int i=0;i<row;i++){
	printf("%f ",vec[i]);
	printf("\n");		
  }
  printf("\n");		
}
void printWS(int* WS, int row){
  for (int i=0;i<row;i++){
	printf("%d ",WS[i]);
  }
  printf("\n");		
  printf("\n");		
}
void print_L(c_float* L, int n){
  int disp = 0;
  for(int i = 0; i<=n;i++){
	for(int j = 0;j<i;j++)
	  printf("%15.15f ",L[disp++]);
	printf("\n");
  }
}


// daqp with print
int daqp_debug(Workspace *work){
  //Begin loop
  printf("M: \n");
  printmatrix(work->M,work->m,work->n);
  printf("d:\n");
  printvector(work->d,work->m);
  printf("\n ==== DAQP starting ====\n\n");
  while(1){
	if(work->iterations++>MAX_ITER)
	  return EXIT_ITERLIMIT;
	printf("|fval: %10.2f | ", work->fval);
	printf("AS: { ");
	for(int i=0;i<work->n_active;i++)
	  printf("%d ",work->WS[i]+1);
	printf("}|\n");
	printf("lambda\n: ");
	//for(int i=0;i<work->n_active;i++)
	//  printf("%f\n ",work->lam[i]);
	//printf("\n");
	//printf("\nL: ");
	//print_L(work->L,work->n_active);

	//printf("\nD:\n ");
	//for(int i =0;i<work->n_active;i++)
	//  printf("%15.15f\n",work->D[i]);
	//printf("\n");

	if(work->sing_ind==EMPTY_IND){ 
	  compute_CSP(work);
	  if(work->fval > work->fval_bound)
		return EXIT_INFEASIBLE;
	  if(work->cycle_counter > CYCLE_TOL){
		if(work->tried_repair==1)
		  return EXIT_CYCLE;
		else{
		  // Try to reorder and refactorize LDL
		  //printf("Cycling detected, trying to reorder/refactorize\n"); 
		  reorder_LDL(work);
		  warmstart_workspace(work, work->WS,work->n_active);
		  work->tried_repair=1;
		  continue;
		}

	  }
	  //printf("\n lam_star: ");
	  //for(int i =0;i<work->n_active;i++)
	  //  printf("%f ",work->lam_star[i]);
	  //printf("\n");
	  //printf("D: ");
	  //for(int i =0;i<work->n_active;i++)
	  //  printf("%f ",work->D[i]);
	  //printf("\n");
	  //printf("L:");
	  //print_L(work->L,work->n_active);
	  find_blocking_constraints(work);
	  if(work->n_blocking==0){ //lam_star >= 0
		find_constraint_to_add(work);
		if(work->add_ind == EMPTY_IND){
		  //printf("Optimal\n");
		  //computePrimalSolution(work);
		  //printf("x_star: ");
		  //  for(int i =0; i<work->n;i++)
		  //    printf(" %f ",work->x[i]);
		  //printf("\n");
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
	  else{// Blocking constraints
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
		//printf("Infeasible\n");
		return EXIT_INFEASIBLE;
		//if(validate_farkas(work)==1) return EXIT_INFEASIBLE;
		//if(work->tried_repair==1) return EXIT_INFEASIBLE;
		//printf("Farkas lemma did not hold, trying repair\n");
		//reorder_LDL(work);
		//warmstart_workspace(work, work->WS,work->n_active-1);
		//work->tried_repair=1;
		//printf("Trying again...\n");
	  }
	  else{
		compute_alpha_and_rm_blocking(work);
		work->sing_ind=EMPTY_IND;
		remove_constraint(work);
	  }
	}
  }
}
