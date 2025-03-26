import numpy as np
cimport daqp

def solve(double[:, :] H, double[:] f, double[:, :] A,
        double[:] bupper, double[:] blower=None, int[:] sense =None,
        primal_tol = DEFAULT_PRIM_TOL, dual_tol = DEFAULT_DUAL_TOL, zero_tol = DEFAULT_ZERO_TOL,
        pivot_tol = DEFAULT_PIVOT_TOL, progress_tol = DEFAULT_PROG_TOL, 
        cycle_tol = DEFAULT_CYCLE_TOL, iter_limit =  DEFAULT_ITER_LIMIT, fval_bound = DAQP_INF,
        eps_prox= 0, eta_prox = DEFAULT_ETA, rho_soft = DEFAULT_RHO_SOFT,
        rel_subopt = DEFAULT_REL_SUBOPT, abs_subopt = DEFAULT_ABS_SUBOPT):
    """
    Solve the quadratic program      minimize       0.5 x'*H*x + f' x
                                    subject to   blower <= A x <= bupper
    Example calls:
    x, fval, exitflag, info = daqp.solve(H, f, A, bupper)
    x, fval, exitflag, info = daqp.solve(H, f, A, bupper, blower)
    x, fval, exitflag, info = daqp.solve(H, f, A, bupper, blower, sense)

    Equality and binary constraints are supported (see information about sense below)

    If bupper and blower have more elements than rows of A, the first elements of
    bupper and bupper are interpreted as simple bounds.
    For example
        A = [[7.0 9.0]]
        bupper = [2.0 3.0]
        blower = [-4.0, -5.0]
    is interpreted as
        -4.0 <= x1 <= 2.0
        -5.0 <= 7.0 x1 + 9.0 x2 <= 3.0.

    Parameters
    ---
    H :
        Cost matrix (should be positive definite).
        For a singular cost matrix, e.g. for LPs, use a positive value for the setting eps_prox
        (see information of how to change settings below)
    f :
        Cost vector.
    A :
        Linear constraint matrix.
    bupper :
        Upper bound on linear constraints.
    blower:
        Lower bound on linear constraints (default = -inf).
    sense:
        Integer array that determines constraint types. For example:
          * 0  ->  Inequality constraint (default)
          * 1  ->  Active inequality constraint (used as warm start)
          * 5  ->  Equality constraint
          * 8  ->  Soft constraint (allowed to be violated if necessary) 
          * 16 ->  Binary constraint (the upper or lower bound should hold with equality)
    See <https://darnstrom.github.io/daqp/parameters>`_ for more information

    Keyword arguments are used for settings. For example:
       daqp.solve(H, f, A, bupper, blower, sense, primal_tol=1e-6, iter_limit=1000)
    Available settings include:
       * iter_limit : Maximum numer of allowed iterations
       * primal_tol : Primal feasibility tolerance
       * dual_tol   : Dual feasibility tolerance
       * eps_pro x  : Regularization used for proximal-point iterations
    See <https://darnstrom.github.io/daqp/parameters>`_ for all available settings.

    Returns
    -------
    x :
        Optimal primal solution.
    fval :
        Optimal cost.
    exitflag :
        Termination status; 1 signifies that an optimal solution was found.
        See <https://darnstrom.github.io/daqp/parameters>`_ for a complete list of exitflags
    info :
        Contains additional information from the solver:
           * setup_time : Time for settings up the problem
           * solve_time : Time for solving the problem
           * iterations : Number of performed iterations
           * nodes      : Explored nodes in branch-and-bound tree
           * lam        : Optimal dual solution
    """

    # Setup problem
    cdef int n,m,ms,nb
    cdef int* bin_ids_cand 
    cdef double *A_ptr, *bu_ptr, *bl_ptr
    cdef int* sense_ptr
    A = np.ascontiguousarray(A)
    mA, n = np.shape(A)
    m = np.size(bupper)

    # By default, set lower bounds to -inf and interpret constraints as inequalities 
    if blower is None:
        blower = np.full(m, -DAQP_INF)
    if sense is None:
        sense = np.zeros(m, dtype=np.intc)

    H_ptr = NULL if H is None else &H[0,0]
    f_ptr = NULL if f is None else &f[0]
    A_ptr = NULL if mA == 0 else &A[0,0]

    if m == 0:
        bu_ptr, bl_ptr, sense_ptr = NULL, NULL, NULL
    else:
        bu_ptr, bl_ptr, sense_ptr  = &bupper[0], &blower[0], &sense[0]


    cdef DAQPProblem problem = [n,m,m-mA, H_ptr, f_ptr, A_ptr, bu_ptr, bl_ptr, sense_ptr]

    # Setup settings
    cdef DAQPSettings settings = [primal_tol, dual_tol, zero_tol, pivot_tol,
            progress_tol, cycle_tol, iter_limit, fval_bound,
            eps_prox, eta_prox, rho_soft, rel_subopt, abs_subopt]
        
    # Setup output
    cdef double[::1] x = np.zeros(n)
    cdef double[::1] lam = np.zeros(m)
    cdef double *lam_ptr
    lam_ptr = NULL if m == 0 else &lam[0]
    cdef DAQPResult res  = [&x[0],lam_ptr,0,0,0,0,0,0,0]

    # Solve 
    with nogil:
        daqp_quadprog(&res, &problem, &settings)
    info = {'solve_time':res.solve_time,
            'setup_time': res.setup_time,
            'iterations': res.iter,
            'nodes': res.nodes,
            'lam': np.asarray(lam)}
    return np.asarray(x), res.fval, res.exitflag, info
def minrep(double[:,:] A, double[:] b):
    # Setup problem
    cdef int n,m,ms
    A = np.ascontiguousarray(A)
    mA, n = np.shape(A)
    m = np.size(b)
    ms = m-mA

    &A[0,0]

    # Setup output
    cdef int[::1] is_redundant = np.zeros(m, dtype=np.intc)

    # Solve 
    with nogil:
        daqp_minrep(&is_redundant[0], &A[0,0], &b[0],n,m,ms)
    return np.asarray(is_redundant)
