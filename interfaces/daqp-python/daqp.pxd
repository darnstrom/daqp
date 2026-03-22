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
    
        int* break_points
        int nh
        int problem_type

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

        double sing_tol;
        double refactor_tol;
        double time_limit;

    ctypedef struct DAQPWorkspace:
        int n
        int m
        int ms
        int nh
        DAQPSettings* settings

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
    cdef extern nogil:
        void daqp_primal_init_active(DAQPProblem *qp, double *x)
    cdef extern nogil:
        void daqp_dual_init_active(DAQPProblem *qp, double *lam)
    cdef extern nogil:
        void daqp_set_primal_start(DAQPWorkspace *work, double *x)
    cdef extern nogil:
        int setup_daqp(DAQPProblem *qp, DAQPWorkspace *work, double *setup_time)
    cdef extern nogil:
        void daqp_solve(DAQPResult *res, DAQPWorkspace *work)
    cdef extern nogil:
        void free_daqp_workspace(DAQPWorkspace *work)
    cdef extern nogil:
        void free_daqp_ldp(DAQPWorkspace *work)
    cdef extern nogil:
        void allocate_daqp_settings(DAQPWorkspace *work)

cdef extern from "utils.h":
    cdef extern nogil:
        int daqp_update_ldp(int mask, DAQPWorkspace *work, DAQPProblem *qp)

cdef extern from "constants.h":
    cdef double DAQP_INF
    cdef double DAQP_DEFAULT_PRIM_TOL
    cdef double DAQP_DEFAULT_DUAL_TOL
    cdef double  DAQP_DEFAULT_ZERO_TOL
    cdef double  DAQP_DEFAULT_PROG_TOL
    cdef double  DAQP_DEFAULT_PIVOT_TOL
    cdef int DAQP_DEFAULT_CYCLE_TOL
    cdef double DAQP_DEFAULT_ETA
    cdef int DAQP_DEFAULT_ITER_LIMIT
    cdef double DAQP_DEFAULT_RHO_SOFT
    cdef double  DAQP_DEFAULT_REL_SUBOPT
    cdef double  DAQP_DEFAULT_ABS_SUBOPT
    cdef double  DAQP_DEFAULT_SING_TOL
    cdef double  DAQP_DEFAULT_REFACTOR_TOL
    cdef double  DAQP_DEFAULT_EPS_PROX
