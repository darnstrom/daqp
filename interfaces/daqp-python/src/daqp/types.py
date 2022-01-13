from ctypes import * 
import ctypes.util

class QP(Structure):
    _fields_ = [('n', c_int),
            ('m', c_int),
            ('ms',c_int),
            ('H', POINTER(c_double)),
            ('f', POINTER(c_double)),
            ('A', POINTER(c_double)),
            ('bupper', POINTER(c_double)),
            ('blower', POINTER(c_double)),
            ('sense', POINTER(c_int))]

class DAQPSettings(Structure):
    def __init__(self, **kwargs):
        values = type(self)._defaults_.copy()
        for (key, val) in kwargs.items():
            values[key] = val
        super(DAQPSettings, self).__init__(**values)     
    
    _fields_= [('primal_tol', c_double),
            ('dual_tol', c_double),
            ('zero_tol', c_double),
            ('pivot_tol', c_double),
            ('progress_tol', c_double),
            ('cycle_tol', c_int),
            ('iter_limit', c_int),
            ('eps_prox', c_double),
            ('eta_prox', c_double),
            ('rho_soft', c_double)]

    _defaults_ = { "primal_tol" : 1e-6,
            "dual_tol" : 1e-12, 
            "zero_tol" : 1e-14,
            "pivot_tol" : 1e-2,
            "progress_tol" : 1e-6,
            "cycle_tol" : 10,
            "iter_limit" : 1000,
            "eps_prox" : 0,
            "eta_prox" : 1e-6,
            "rho_soft" :1e-3 
            }

class DAQPResult(Structure):
    _fields_ = [('x', POINTER(c_double)),
            ('fval',c_double),
            ('soft_slack',c_double),
            ('exitflag',c_int),
            ('iter',c_int),
            ('solve_time',c_double),
            ('setup_time',c_double)] 

