import numpy as np
cimport daqp
from libc.stdlib cimport calloc, free

def solve(double[:, :] H, double[:] f, double[:, :] A,
          double[:] bupper, double[:] blower=None,
          int[:] sense =None, int[:] break_points=np.zeros(0,dtype=np.intc),
          is_avi=False,
          primal_tol = DAQP_DEFAULT_PRIM_TOL, dual_tol = DAQP_DEFAULT_DUAL_TOL,
          zero_tol = DAQP_DEFAULT_ZERO_TOL, pivot_tol = DAQP_DEFAULT_PIVOT_TOL,
          progress_tol = DAQP_DEFAULT_PROG_TOL, cycle_tol = DAQP_DEFAULT_CYCLE_TOL,
          iter_limit =  DAQP_DEFAULT_ITER_LIMIT, fval_bound = DAQP_INF,
          eps_prox= 0, eta_prox = DAQP_DEFAULT_ETA, rho_soft = DAQP_DEFAULT_RHO_SOFT,
          rel_subopt = DAQP_DEFAULT_REL_SUBOPT, abs_subopt = DAQP_DEFAULT_ABS_SUBOPT,
          sing_tol = DAQP_DEFAULT_SING_TOL, refactor_tol = DAQP_DEFAULT_REFACTOR_TOL,
          primal_start=None, dual_start=None):
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
    break_points:
        Break points that define prioritized constraints
    is_avi:
        flag for interpreting the problem as an AVI (relevant for asymmetric H) 
    primal_start:
        Initial guess of primal iterate (used for warm starting). Sets the initial
        active set based on which constraints are nearly active at this iterate, and
        also sets the initial primal iterate in the workspace.
    dual_start:
        Initial guess of dual iterate / Lagrange multipliers (used for warm starting).
        Sets the initial active set based on the sign of the multipliers.

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
    cdef int n,m,ms,nb,problem_type
    cdef int* bin_ids_cand 
    cdef double *A_ptr, *bu_ptr, *bl_ptr
    cdef int* sense_ptr
    A = np.ascontiguousarray(A)
    mA, n = np.shape(A)
    m = np.size(bupper)
    nh = np.size(break_points)

    # By default, set lower bounds to -inf and interpret constraints as inequalities 
    if blower is None:
        blower = np.full(m, -DAQP_INF)
    if sense is None:
        sense = np.zeros(m, dtype=np.intc)
    problem_type = 1 if is_avi else 0

    H_ptr = NULL if H is None else &H[0,0]
    f_ptr = NULL if f is None else &f[0]
    A_ptr = NULL if mA == 0 else &A[0,0]
    bp_ptr = NULL if nh == 0 else &break_points[0]

    if m == 0:
        bu_ptr, bl_ptr, sense_ptr = NULL, NULL, NULL
    else:
        bu_ptr, bl_ptr, sense_ptr  = &bupper[0], &blower[0], &sense[0]


    cdef DAQPProblem problem = [n,m,m-mA, H_ptr, f_ptr, A_ptr, bu_ptr, bl_ptr, sense_ptr, bp_ptr, nh, problem_type]

    # Setup settings
    cdef DAQPSettings settings = [primal_tol, dual_tol, zero_tol, pivot_tol,
            progress_tol, cycle_tol, iter_limit, fval_bound,
            eps_prox, eta_prox, rho_soft, rel_subopt, abs_subopt, sing_tol, refactor_tol]
        
    # Setup output
    cdef double[::1] x = np.zeros(n)
    cdef double[::1] lam = np.zeros(m)
    cdef double *lam_ptr
    lam_ptr = NULL if m == 0 else &lam[0]
    cdef DAQPResult res  = [&x[0],lam_ptr,0,0,0,0,0,0,0]

    # Warm starting: initialize active set from primal or dual iterate
    cdef int[::1] sense_copy
    cdef double[::1] primal_start_c, dual_start_c
    cdef DAQPWorkspace* work
    cdef int setup_flag
    cdef double setup_time_c

    if dual_start is not None or primal_start is not None:
        # Copy sense to avoid modifying the user's array when init functions mark constraints active
        if m > 0:
            sense_copy = np.array(sense, dtype=np.intc, copy=True)
            problem.sense = &sense_copy[0]

        if dual_start is not None:
            dual_start_c = np.ascontiguousarray(dual_start, dtype=np.double)
            with nogil:
                daqp_dual_init_active(&problem, &dual_start_c[0])
        else:
            primal_start_c = np.ascontiguousarray(primal_start, dtype=np.double)
            with nogil:
                daqp_primal_init_active(&problem, &primal_start_c[0])

    if primal_start is not None:
        # Use separate setup/solve to allow setting the initial primal iterate
        work = <DAQPWorkspace*>calloc(1, sizeof(DAQPWorkspace))
        # Point workspace settings to stack-allocated settings struct
        work.settings = &settings
        setup_time_c = 0.0
        with nogil:
            setup_flag = setup_daqp(&problem, work, &setup_time_c)
        res.setup_time = setup_time_c
        res.exitflag = setup_flag
        if setup_flag >= 0:
            with nogil:
                daqp_set_primal_start(work, &primal_start_c[0])
                daqp_solve(&res, work)
        # Nullify settings pointer before freeing to prevent free of stack memory
        work.settings = NULL
        with nogil:
            free_daqp_workspace(work)
            free_daqp_ldp(work)
        free(work)
    else:
        # Solve (active set already initialized above if dual_start was given)
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
