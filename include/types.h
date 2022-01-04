#ifndef DAQP_TYPES_H
# define DAQP_TYPES_H
#include "constants.h"

typedef struct{

// Data for the QP problem
//
// min		0.5 x'*H*x + f'x
// s.t 	blower <= A*x <= bupper 
// 			lb <=  x  <= ub
//
// n  - dimension of x 
// m  - # of rows in b 
// ms - # of simple bounds (assumed to be the first ms elements in b) 
// (The number of rows in A is hence m-ms)

// sense define the state of the constraints 
// (active, immutable, upper/lower, soft). 

int n;
int m;
int ms;

c_float* H;
c_float* f;

c_float* A;
c_float* bupper;
c_float* blower;

int* sense; 
}QP;
typedef struct{
  c_float primal_tol; 
  c_float dual_tol; 
  c_float zero_tol; 
  c_float pivot_tol;
  c_float progress_tol;
  
  int cycle_tol;
  int iter_limit;
 
  c_float eps_prox;
  c_float eta_prox;
  int prox_iter_limit;

  c_float rho_soft;
}DAQPSettings;

typedef struct{
  QP* qp;
  // LDP data 
  int n; // Number of primal variable  
  int m; // Number of constraints  
  int ms; // Number of simple bounds
  c_float *M; // M' M is the Hessian of the dual objective function (dimensions: n x m)  
  c_float *dupper; // Linear part of dual objective function (dimensions: m x 1) 
  c_float *dlower; // Linear part of dual objective function (dimensions: m x 1) 
  c_float *Rinv; // Inverse of upper cholesky factor of primal Hessian 
  c_float *v; // v = R'\f (used to transform QP to LDP 
  int *sense; // State of constraints  


  // Iterates
  c_float *x; // The final primal solution
  c_float *xold; // The latest primal solution (used for proximal-point iteratios)

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
  int outer_iter;
  int inner_iter;
  int sing_ind; // Flag for denoting whether Mk Mk' is singular or not 
  int add_ind; // Index to add to the working set
  int add_isupper; // Marks if index to add is upper or lower 
  int rm_ind; // Index to remove from the working set
  int n_active; // Number of active contraints 
  int n_blocking; // Number of blocking constraints  
  int reuse_ind; // How much work that can be saved when solving Mk Mk' lam* = -dk
  int cycle_counter; // Number of iterations when no progress has been made
  int tried_repair; // Flag to mark if repair has taken place 

  // Intermittent variables (LDL') lam_star = -dk
  c_float* xldl; // Solution to L xdldl = -dk
  c_float* zldl; // zldl_i = xldl_i/D_i

  // Minimum slack for soft constraints
  c_float soft_slack; 

  DAQPSettings* settings;
}Workspace;

#endif //ifndef DAQP_TYPES_H
