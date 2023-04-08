import numpy as np
cimport daqp

def solve(double[:, :] H, double[:] f, double[:, :] A,
        double[:] bupper, double[:] blower=None, int[:] sense =None,
        primal_tol = DEFAULT_PRIM_TOL, dual_tol = DEFAULT_DUAL_TOL, zero_tol = DEFAULT_ZERO_TOL,
        pivot_tol = DEFAULT_PIVOT_TOL, progress_tol = DEFAULT_PROG_TOL, 
        cycle_tol = DEFAULT_CYCLE_TOL, iter_limit =  DEFAULT_ITER_LIMIT, fval_bound = DAQP_INF,
        eps_prox= 0, eta_prox = DEFAULT_ETA, rho_soft = DEFAULT_RHO_SOFT,
        rel_subopt = DEFAULT_REL_SUBOPT, abs_subopt = DEFAULT_ABS_SUBOPT):
    # Setup problem
    cdef int n,m,ms,nb
    cdef int* bin_ids_cand 
    A = np.ascontiguousarray(A)
    mA, n = np.shape(A)
    m = np.size(bupper)

    if blower is None:
        blower = np.fill(m, -DAQP_INF)
    if sense is None:
        sense = np.zeros(m, dtype=int)

    cdef DAQPProblem problem = [n,m,m-mA, &H[0,0], &f[0], 
            &A[0,0], &bupper[0], &blower[0], &sense[0]]

    # Setup settings
    cdef DAQPSettings settings = [primal_tol, dual_tol, zero_tol, pivot_tol,
            progress_tol, cycle_tol, iter_limit, fval_bound,
            eps_prox, eta_prox, rho_soft, rel_subopt, abs_subopt]
        
    # Setup output
    cdef double[::1] x = np.zeros(n)
    cdef double[::1] lam = np.zeros(m)
    cdef DAQPResult res  = [&x[0],&lam[0],0,0,0,0,0,0,0] 

    # Solve 
    with nogil:
        daqp_quadprog(&res, &problem, &settings)
    info = {'solve_time':res.solve_time,
            'setup_time': res.setup_time,
            'iterations': res.iter,
            'nodes': res.nodes,
            'lam': np.asarray(lam)}
    return np.asarray(x), res.fval, res.exitflag, info
