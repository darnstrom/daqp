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
        int daqp_minrep(int *is_redundant, double *A, double*b, int n, int m, int ms)
    cdef extern nogil:
        int daqp_default_settings(DAQPSettings *settings)

cdef extern from "constants.h":
    cdef double DAQP_INF
    cdef double DEFAULT_PRIM_TOL
    cdef double DEFAULT_DUAL_TOL
    cdef double  DEFAULT_ZERO_TOL
    cdef double  DEFAULT_PROG_TOL
    cdef double  DEFAULT_PIVOT_TOL
    cdef int DEFAULT_CYCLE_TOL
    cdef double DEFAULT_ETA
    cdef int DEFAULT_ITER_LIMIT
    cdef double DEFAULT_RHO_SOFT
    cdef double  DEFAULT_REL_SUBOPT
    cdef double  DEFAULT_ABS_SUBOPT
