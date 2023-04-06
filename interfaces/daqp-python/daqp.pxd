cdef extern from "types.h":
    ctypedef struct DAQPProblem:
        int n
        int m
        int ms 
    
        double* H
        double* f
    
        double* A
        double* bupper
        double* blower 
    
        int* sense
    
        int* bin_ids
        int nb
    
    ctypedef struct DAQPSettings:
        double primal_tol;
        double dual_tol;
        double zero_tol;
        double pivot_tol;
        double progress_tol;
    
        int cycle_tol;
        int iter_limit;
        double fval_bound;
    
        double eps_prox;
        double eta_prox;
    
        double rho_soft;
    
        double rel_subopt;
        double abs_subopt;

cdef extern from "api.h":
    ctypedef struct DAQPResult:
        double *x;
        double *lam;
        double fval;
        double soft_slack;
    
        int exitflag;
        int iter;
        int nodes;
        double solve_time;
        double setup_time;
    
    cdef extern nogil:
        int daqp_quadprog(DAQPResult *res, DAQPProblem *prob, DAQPSettings *settings)
    cdef extern nogil:
        int daqp_default_settings(DAQPSettings *settings)
