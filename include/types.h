#ifndef DAQP_TYPES_H
# define DAQP_TYPES_H
#include "constants.h"
typedef struct{
 int n; // Number of primal variable  
 int m; // Number of constraints  
 c_float *M; // M' M is the Hessian of the dual objective function (dimensions: n x m)  
 c_float *d; // Lienar part of dual objective function (dimensions: m x 1) 
 c_float *u; // Unconstrained primal solution
 c_float *H_chol_inv; // H = H_chol * H_chol' (Diagonals inverse to avoid division)
}QPData; // TODO Use H,f,A,b instead...

typedef struct{
}SolverSettings;

typedef struct{
  // Problem data 
  int n; // Number of primal variable  
  int m; // Number of constraints  
  c_float *M; // M' M is the Hessian of the dual objective function (dimensions: n x m)  
  c_float *d; // Lienar part of dual objective function (dimensions: m x 1) 
  int *sense; // Denotes inequality or equality constraints

  // Iterates
  c_float* lam; // Dual iterate 
  c_float* lam_star; // Current constrained stationary point 
  c_float* u; // Stores Mk' lam_star
  c_float fval;
  c_float fval_bound;

  // LDL factors (Mk Mk' = L D L')
  c_float *L;
  c_float *D;

  int *WS; // Working set
  int *BS; // Set of blocking constraints 

  int iterations;
  int sing_ind; // Flag for denoting whether Mk Mk' is singular or not 
  int add_ind; // Index to add to the working set
  int rm_ind; // Index to remove from the working set
  int n_active; // Number of active contraints 
  int n_blocking; // Number of blocking constraints  
  int reuse_ind; // How much work that can be saved when solving Mk Mk' lam* = -dk
  int cycle_counter; // Number of iterations when no progress has been made
  int tried_repair; // Flag to mark if repair has taken place 

  c_float *swp_pointer;

  // Intermittent variables (LDL') lam_star = -dk
  c_float* xldl; // Solution to L xdldl = -dk
  c_float* zldl; // zldl_i = xldl_i/D_i
}Workspace;
#endif //ifndef DAQP_TYPES_H
