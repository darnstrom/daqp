import numpy as np
cimport daqp
from libc.stdlib cimport calloc, free
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _solve_warm_start(DAQPProblem* problem, DAQPSettings* settings,
                            DAQPResult* res, int[:] sense, int m,
                            object primal_start, object dual_start) except *:
    """Handle warm-started solve; called only when a warm start is requested.

    Called from solve() when primal_start or dual_start is not None.
    All expensive Python work is isolated here so the common (cold-start)
    path of solve() stays as lean as possible.
    """
    cdef int[::1]    sense_copy
    cdef double[::1] primal_start_c, dual_start_c
    cdef DAQPWorkspace* work
    cdef int    setup_flag
    cdef double setup_time_c

    # Copy sense so the init helpers don't mutate the caller's array
    if m > 0:
        sense_copy = np.array(sense, dtype=np.intc, copy=True)
        problem.sense = &sense_copy[0]

    if primal_start is not None:
        primal_start_c = np.ascontiguousarray(primal_start, dtype=np.double)

    if dual_start is not None:
        dual_start_c = np.ascontiguousarray(dual_start, dtype=np.double)
        with nogil:
            daqp_dual_init_active(problem, &dual_start_c[0])
    else:
        # primal_start is not None (checked by caller)
        with nogil:
            daqp_primal_init_active(problem, &primal_start_c[0])

    if primal_start is not None:
        # Need a workspace to inject the initial primal iterate before solving
        work = <DAQPWorkspace*>calloc(1, sizeof(DAQPWorkspace))
        work.settings = settings
        setup_time_c = 0.0
        with nogil:
            setup_flag = setup_daqp_main(problem, work, &setup_time_c,1)
        res.setup_time = setup_time_c
        res.exitflag  = setup_flag
        if setup_flag >= 0:
            with nogil:
                daqp_set_primal_start(work, &primal_start_c[0])
                daqp_solve(res, work)
        # Nullify settings pointer before freeing to prevent free of stack memory
        work.settings = NULL
        with nogil:
            free_daqp_workspace(work)
            free_daqp_ldp(work)
        free(work)
    else:
        # dual_start only: active set already initialized above
        with nogil:
            daqp_quadprog(res, problem, settings)

@cython.boundscheck(False)
@cython.wraparound(False)
def solve(double[:, :] H, double[:] f, double[:, :] A,
          double[:] bupper, double[:] blower=None,
          int[:] sense =None, int[:] break_points=np.zeros(0,dtype=np.intc),
          is_avi=False,
          primal_tol = DAQP_DEFAULT_PRIM_TOL, dual_tol = DAQP_DEFAULT_DUAL_TOL,
          zero_tol = DAQP_DEFAULT_ZERO_TOL, pivot_tol = DAQP_DEFAULT_PIVOT_TOL,
          progress_tol = DAQP_DEFAULT_PROG_TOL, cycle_tol = DAQP_DEFAULT_CYCLE_TOL,
          iter_limit =  DAQP_DEFAULT_ITER_LIMIT, fval_bound = DAQP_INF,
          eps_prox= DAQP_DEFAULT_EPS_PROX, eta_prox = DAQP_DEFAULT_ETA, 
          rho_soft = DAQP_DEFAULT_RHO_SOFT,
          rel_subopt = DAQP_DEFAULT_REL_SUBOPT, abs_subopt = DAQP_DEFAULT_ABS_SUBOPT,
          sing_tol = DAQP_DEFAULT_SING_TOL, refactor_tol = DAQP_DEFAULT_REFACTOR_TOL,
          time_limit = 0,
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

    # Setup problem dimensions using fast memoryview shape attributes
    cdef int mA = A.shape[0]
    cdef int n = A.shape[1]
    cdef int m = bupper.shape[0]
    cdef int nh = break_points.shape[0]
    cdef int problem_type = 1 if is_avi else 0

    # Ensure A is C-contiguous; only copy if necessary (common case is free)
    if not A.is_c_contig():
        A = np.ascontiguousarray(A)

    # By default, set lower bounds to -inf and interpret constraints as inequalities
    if blower is None:
        blower = np.full(m, -DAQP_INF)
    if sense is None:
        sense = np.zeros(m, dtype=np.intc)

    # Setup output buffers (np.empty avoids zeroing; the solver fills all entries)
    cdef double[::1] x = np.empty(n)
    cdef double[::1] lam = np.empty(m) if m > 0 else np.empty(0)

    # Build C-level structs — all pointer arithmetic from this point on
    cdef double *H_ptr = NULL if H is None else &H[0,0]
    cdef double *f_ptr = NULL if f is None else &f[0]
    cdef double *A_ptr = NULL if mA == 0 else &A[0,0]
    cdef int   *bp_ptr = NULL if nh == 0 else &break_points[0]
    cdef double *bu_ptr
    cdef double *bl_ptr
    cdef int   *sense_ptr
    cdef double *lam_ptr
    if m == 0:
        bu_ptr = bl_ptr = NULL
        sense_ptr = NULL
        lam_ptr = NULL
    else:
        bu_ptr = &bupper[0]
        bl_ptr = &blower[0]
        sense_ptr = &sense[0]
        lam_ptr = &lam[0]

    cdef DAQPProblem problem = [n, m, m-mA, H_ptr, f_ptr, A_ptr, bu_ptr, bl_ptr, sense_ptr, bp_ptr, nh, problem_type]
    cdef DAQPSettings settings = [primal_tol, dual_tol, zero_tol, pivot_tol,
            progress_tol, cycle_tol, iter_limit, fval_bound,
            eps_prox, eta_prox, rho_soft, rel_subopt, abs_subopt, sing_tol, refactor_tol,
            time_limit]
    cdef DAQPResult res = [&x[0], lam_ptr, 0, 0, 0, 0, 0, 0, 0]

    if primal_start is None and dual_start is None:
        # Fast path: no warm start — enter C as quickly as possible
        with nogil:
            daqp_quadprog(&res, &problem, &settings)
    else:
        # Warm-start path: initialize active set, then solve
        _solve_warm_start(&problem, &settings, &res, sense, m,
                          primal_start, dual_start)

    return (np.asarray(x), res.fval, res.exitflag,
            {'solve_time': res.solve_time, 'setup_time': res.setup_time,
             'iterations': res.iter, 'nodes': res.nodes,
             'lam': np.asarray(lam)})
cdef class Model:
    """
    Wraps a DAQP workspace for efficient reuse across similar problems.

    Allocates memory once during ``setup`` and reuses it for all subsequent
    ``solve`` calls — no heap allocations occur during solving — making this
    the recommended interface for embedded or real-time applications.

    Methods
    -------
    setup(H, f, A, bupper, blower=None, sense=None, break_points=None,
          primal_start=None, dual_start=None, is_avi=False)
        Set up the QP and allocate workspace memory.
    solve()
        Solve the currently set-up problem.
    update(H=None, f=None, A=None, bupper=None, blower=None,
           sense=None, break_points=None)
        Update problem data without re-running the full setup.
    settings
        Property (dict) for getting and setting solver parameters.

    Example
    -------
    >>> d = daqp.Model()
    >>> d.setup(H, f, A, bupper, blower, sense)
    >>> x, fval, exitflag, info = d.solve()
    >>> d.update(f=f_new)
    >>> x, fval, exitflag, info = d.solve()
    """

    cdef DAQPWorkspace* _work
    cdef DAQPProblem    _qp

    # Typed memoryviews keep the underlying numpy buffers alive and provide
    # C-level pointers into them.
    cdef double[:, ::1] _H
    cdef double[::1]    _f
    cdef double[:, ::1] _A
    cdef double[::1]    _bupper
    cdef double[::1]    _blower
    cdef int[::1]       _sense
    cdef int[::1]       _break_points

    # Pre-allocated result buffers
    cdef double[::1]    _x
    cdef double[::1]    _lam

    cdef bint _has_model

    def __cinit__(self):
        self._work = <DAQPWorkspace*>calloc(1, sizeof(DAQPWorkspace))
        if self._work == NULL:
            raise MemoryError("Failed to allocate DAQPWorkspace")
        allocate_daqp_settings(self._work)
        self._has_model = False

    def __dealloc__(self):
        if self._work != NULL:
            free_daqp_workspace(self._work)
            free_daqp_ldp(self._work)
            free(self._work)
            self._work = NULL

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def setup(self, double[:, :] H, double[:] f, double[:, :] A,
              double[:] bupper, double[:] blower=None,
              int[:] sense=None, int[:] break_points=None,
              primal_start=None, dual_start=None, is_avi=False):
        """
        Set up the QP problem and allocate the solver workspace.

        Parameters follow the same convention as ``daqp.solve``.  Call this
        once before the first ``solve``; use ``update`` for subsequent data
        changes so that already-allocated memory is reused.

        Parameters
        ----------
        H : 2-D array_like
            Cost matrix (positive definite).  Pass ``None`` for LPs.
        f : 1-D array_like
            Cost vector.
        A : 2-D array_like
            Linear constraint matrix, shape ``(mA, n)``.
        bupper : 1-D array_like
            Upper bounds (length ``m >= mA``).  The first ``m - mA``
            elements are interpreted as simple variable bounds.
        blower : 1-D array_like, optional
            Lower bounds (default: ``-inf``).
        sense : 1-D int array_like, optional
            Constraint types (default: all 0 = inequality).
        break_points : 1-D int array_like, optional
            Break points for hierarchical QPs.
        primal_start : 1-D array_like, optional
            Initial primal iterate for warm starting.
        dual_start : 1-D array_like, optional
            Initial dual iterate for warm starting.

        Returns
        -------
        exitflag : int
            >= 0 on success, < 0 on failure.
        setup_time : float
            Time (seconds) spent in setup.
        """
        cdef int n, m, mA, ms, nh
        cdef double setup_time_c = 0.0
        cdef int setup_flag
        cdef DAQPSettings old_settings_val
        cdef bint restore_settings = False
        cdef double[::1] primal_start_c
        cdef double[::1] dual_start_c
        cdef int[::1]    sense_copy

        A = np.ascontiguousarray(A)
        mA, n = np.shape(A)
        m  = np.size(bupper)
        ms = m - mA

        if blower is None:
            blower = np.full(m, -DAQP_INF)
        if sense is None:
            sense = np.zeros(m, dtype=np.intc)
        if break_points is None:
            break_points = np.zeros(0, dtype=np.intc)
        nh = np.size(break_points)

        # ---- store persistent typed memoryviews ----
        self._H            = np.ascontiguousarray(np.asarray(H), dtype=np.double)
        self._f            = np.ascontiguousarray(np.asarray(f), dtype=np.double)
        self._A            = np.ascontiguousarray(np.asarray(A), dtype=np.double)
        self._bupper       = np.ascontiguousarray(np.asarray(bupper), dtype=np.double)
        self._blower       = np.ascontiguousarray(np.asarray(blower), dtype=np.double)
        self._sense        = np.ascontiguousarray(np.asarray(sense),  dtype=np.intc)
        self._break_points = np.ascontiguousarray(np.asarray(break_points), dtype=np.intc)

        # Pre-allocate result buffers
        self._x   = np.zeros(n, dtype=np.double)
        self._lam = np.zeros(m, dtype=np.double)

        # ---- save current settings before any teardown ----
        if self._work.settings != NULL:
            old_settings_val = self._work.settings[0]
            restore_settings = True

        # ---- free old workspace internals if re-setting up ----
        if self._has_model:
            free_daqp_workspace(self._work)   # also frees settings
            free_daqp_ldp(self._work)
        else:
            if self._work.settings != NULL:
                free(self._work.settings)
                self._work.settings = NULL

        # ---- pre-allocate settings with user values so setup_daqp sees them ----
        # This is critical: daqp_update_Rinv (called inside setup_daqp) reads
        # eps_prox from settings to decide which Cholesky diagonals need
        # regularisation.  If settings were NULL or held default (eps_prox=0)
        # when the Cholesky runs, a singular Hessian would be rejected as
        # non-convex even when the caller supplied a non-zero eps_prox.
        allocate_daqp_settings(self._work)   # fresh allocation, defaults
        if restore_settings:
            self._work.settings[0] = old_settings_val  # apply user values NOW

        # ---- build the problem struct ----
        cdef double* H_ptr  = NULL if H is None else &self._H[0, 0]
        cdef double* f_ptr  = NULL if f is None else &self._f[0]
        cdef double* A_ptr  = NULL if mA == 0   else &self._A[0, 0]
        cdef double* bu_ptr = NULL if m  == 0   else &self._bupper[0]
        cdef double* bl_ptr = NULL if m  == 0   else &self._blower[0]
        cdef int*    s_ptr  = NULL if m  == 0   else &self._sense[0]
        cdef int*    bp_ptr = NULL if nh == 0   else &self._break_points[0]

        self._qp = [n, m, ms, H_ptr, f_ptr, A_ptr, bu_ptr, bl_ptr,
                    s_ptr, bp_ptr, nh, 1 if is_avi else 0]

        # ---- warm-start active-set initialisation ----
        # Initialising the active set modifies sense, so work on a copy.
        if dual_start is not None:
            dual_start_c = np.ascontiguousarray(dual_start, dtype=np.double)
        if primal_start is not None:
            primal_start_c = np.ascontiguousarray(primal_start, dtype=np.double)

        if dual_start is not None or primal_start is not None:
            sense_copy = np.array(np.asarray(self._sense), dtype=np.intc, copy=True)
            self._sense = sense_copy
            s_ptr = &self._sense[0]
            self._qp.sense = s_ptr

            if dual_start is not None:
                with nogil:
                    daqp_dual_init_active(&self._qp, &dual_start_c[0])
            else:
                with nogil:
                    daqp_primal_init_active(&self._qp, &primal_start_c[0])

        # ---- call setup_daqp ----
        with nogil:
            setup_flag = setup_daqp_main(&self._qp, self._work, &setup_time_c,0)

        if setup_flag < 0:
            # setup_daqp freed settings on failure; re-allocate and restore.
            allocate_daqp_settings(self._work)
            if restore_settings:
                self._work.settings[0] = old_settings_val
            self._has_model = False
            return setup_flag, setup_time_c

        # ---- restore user-modified settings ----
        if restore_settings:
            self._work.settings[0] = old_settings_val

        self._has_model = True

        # ---- set primal iterate ----
        if primal_start is not None:
            with nogil:
                daqp_set_primal_start(self._work, &primal_start_c[0])

        return setup_flag, setup_time_c

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def solve(self):
        """
        Solve the currently set-up problem.

        Returns
        -------
        x : ndarray
            Optimal primal solution.
        fval : float
            Optimal cost.
        exitflag : int
            Termination status (> 0 success, < 0 failure).
        info : dict
            Additional solver information with keys ``solve_time``,
            ``setup_time``, ``iterations``, ``nodes``, and ``lam``.
        """
        if not self._has_model:
            raise RuntimeError("Model has not been set up. Call setup() first.")

        cdef int     m       = self._work.m
        cdef double* lam_ptr = NULL if m == 0 else &self._lam[0]
        cdef DAQPResult res  = [&self._x[0], lam_ptr, 0, 0, 0, 0, 0, 0, 0]

        with nogil:
            daqp_solve(&res, self._work)

        info = {
            'solve_time': res.solve_time,
            'setup_time': 0.0,
            'iterations': res.iter,
            'nodes':      res.nodes,
            'lam':        np.asarray(self._lam).copy(),
        }
        return np.asarray(self._x).copy(), res.fval, res.exitflag, info

    def update(self, H=None, f=None, A=None, bupper=None, blower=None,
               sense=None, break_points=None):
        """
        Update problem data and refresh the internal LDP representation.

        Only components that are explicitly provided (not ``None``) are
        updated.  The updated arrays must have the same dimensions as those
        passed to ``setup``.

        Parameters
        ----------
        H : 2-D array_like or None
            Updated cost matrix.
        f : 1-D array_like or None
            Updated cost vector.
        A : 2-D array_like or None
            Updated constraint matrix.
        bupper : 1-D array_like or None
            Updated upper bounds.
        blower : 1-D array_like or None
            Updated lower bounds.
        sense : 1-D int array_like or None
            Updated constraint sense flags.
        break_points : 1-D int array_like or None
            Updated break points.

        Returns
        -------
        exitflag : int
            0 on success, negative on failure.
        """
        if not self._has_model:
            raise RuntimeError("Model has not been set up. Call setup() first.")

        cdef int update_mask = 0
        cdef int n  = self._work.n
        cdef int m  = self._work.m
        cdef int ms = self._work.ms
        cdef int mA = m - ms
        cdef int nh = self._work.nh
        cdef int exitflag

        if H is not None:
            H_arr = np.ascontiguousarray(H, dtype=np.double)
            if H_arr.shape[0] == n and H_arr.shape[1] == n:
                self._H = H_arr
                self._qp.H = &self._H[0, 0]
                update_mask |= 1   # DAQP_UPDATE_Rinv

        if A is not None:
            A_arr = np.ascontiguousarray(A, dtype=np.double)
            if A_arr.shape[0] == mA and A_arr.shape[1] == n:
                self._A = A_arr
                self._qp.A = &self._A[0, 0]
                update_mask |= 2   # DAQP_UPDATE_M

        if f is not None:
            f_arr = np.ascontiguousarray(f, dtype=np.double)
            if f_arr.shape[0] == n:
                self._f = f_arr
                self._qp.f = &self._f[0]
                update_mask |= 4   # DAQP_UPDATE_v

        if bupper is not None:
            bu_arr = np.ascontiguousarray(bupper, dtype=np.double)
            if bu_arr.shape[0] == m:
                self._bupper = bu_arr
                self._qp.bupper = &self._bupper[0]
                update_mask |= 8   # DAQP_UPDATE_d

        if blower is not None:
            bl_arr = np.ascontiguousarray(blower, dtype=np.double)
            if bl_arr.shape[0] == m:
                self._blower = bl_arr
                self._qp.blower = &self._blower[0]
                update_mask |= 8   # DAQP_UPDATE_d

        if sense is not None:
            s_arr = np.ascontiguousarray(sense, dtype=np.intc)
            if s_arr.shape[0] == m:
                self._sense = s_arr
                self._qp.sense = &self._sense[0]
                update_mask |= 16  # DAQP_UPDATE_sense

        if break_points is not None:
            bp_arr = np.ascontiguousarray(break_points, dtype=np.intc)
            if bp_arr.shape[0] == nh:
                self._break_points = bp_arr
                self._qp.break_points = &self._break_points[0]
                update_mask |= 32  # DAQP_UPDATE_hierarchy

        with nogil:
            exitflag = daqp_update_ldp(update_mask, self._work, &self._qp)
        return exitflag

    @property
    def settings(self):
        """
        Current solver settings as a ``dict``.

        Assign a ``dict`` (with a subset of keys) to update individual
        settings without touching the others.

        Available keys: ``primal_tol``, ``dual_tol``, ``zero_tol``,
        ``pivot_tol``, ``progress_tol``, ``cycle_tol``, ``iter_limit``,
        ``fval_bound``, ``eps_prox``, ``eta_prox``, ``rho_soft``,
        ``rel_subopt``, ``abs_subopt``, ``sing_tol``, ``refactor_tol``,
        ``time_limit``.
        """
        if self._work.settings == NULL:
            return {}
        cdef DAQPSettings* s = self._work.settings
        return {
            'primal_tol':   s.primal_tol,
            'dual_tol':     s.dual_tol,
            'zero_tol':     s.zero_tol,
            'pivot_tol':    s.pivot_tol,
            'progress_tol': s.progress_tol,
            'cycle_tol':    s.cycle_tol,
            'iter_limit':   s.iter_limit,
            'fval_bound':   s.fval_bound,
            'eps_prox':     s.eps_prox,
            'eta_prox':     s.eta_prox,
            'rho_soft':     s.rho_soft,
            'rel_subopt':   s.rel_subopt,
            'abs_subopt':   s.abs_subopt,
            'sing_tol':     s.sing_tol,
            'refactor_tol': s.refactor_tol,
            'time_limit':   s.time_limit,
        }

    @settings.setter
    def settings(self, new_settings):
        """Update solver settings from a dict (only provided keys are changed)."""
        if self._work.settings == NULL:
            return
        cdef DAQPSettings* s = self._work.settings
        if 'primal_tol'   in new_settings: s.primal_tol   = new_settings['primal_tol']
        if 'dual_tol'     in new_settings: s.dual_tol     = new_settings['dual_tol']
        if 'zero_tol'     in new_settings: s.zero_tol     = new_settings['zero_tol']
        if 'pivot_tol'    in new_settings: s.pivot_tol    = new_settings['pivot_tol']
        if 'progress_tol' in new_settings: s.progress_tol = new_settings['progress_tol']
        if 'cycle_tol'    in new_settings: s.cycle_tol    = new_settings['cycle_tol']
        if 'iter_limit'   in new_settings: s.iter_limit   = new_settings['iter_limit']
        if 'fval_bound'   in new_settings: s.fval_bound   = new_settings['fval_bound']
        if 'eps_prox'     in new_settings: s.eps_prox     = new_settings['eps_prox']
        if 'eta_prox'     in new_settings: s.eta_prox     = new_settings['eta_prox']
        if 'rho_soft'     in new_settings: s.rho_soft     = new_settings['rho_soft']
        if 'rel_subopt'   in new_settings: s.rel_subopt   = new_settings['rel_subopt']
        if 'abs_subopt'   in new_settings: s.abs_subopt   = new_settings['abs_subopt']
        if 'sing_tol'     in new_settings: s.sing_tol     = new_settings['sing_tol']
        if 'refactor_tol' in new_settings: s.refactor_tol = new_settings['refactor_tol']
        if 'time_limit'   in new_settings: s.time_limit   = new_settings['time_limit']


@cython.boundscheck(False)
@cython.wraparound(False)
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
