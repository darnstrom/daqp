import numpy as np
cimport daqp

def solve(double[:, :] H, double[:] f, double[:, :] A, double[:] bupper, double[:] blower, int[:] sense):
    cdef int n,m,ms,nb
    cdef int* bin_ids_cand 
    A = np.ascontiguousarray(A)
    mA, n = np.shape(A)
    m = np.size(bupper)
    cdef DAQPProblem problem = [n,m,m-mA, &H[0,0], &f[0], 
            &A[0,0], &bupper[0], &blower[0], &sense[0], NULL,0]
    #cdef DAQPSettings settings
    #settings.iter_limit = 10000
    #with nogil:
    #    daqp_default_settings(&settings)
        
    cdef double[::1] x = np.zeros(n)
    cdef double[::1] lam = np.zeros(m)
    cdef DAQPResult res  = [&x[0],&lam[0],0,0,0,0,0,0,0] 

    with nogil:
        daqp_quadprog(&res, &problem, NULL)
    info = {'solve_time':res.solve_time,
            'setup_time': res.setup_time,
            'iterations': res.iter,
            'nodes': res.nodes,
            'lam': np.asarray(lam)}
    return np.asarray(x), res.fval, res.exitflag, info
