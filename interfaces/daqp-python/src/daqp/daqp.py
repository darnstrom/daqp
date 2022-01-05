import numpy as np
from ctypes import * 
import ctypes.util
import os.path
import daqp.types as types

class daqp:
    def __init__(self):

        # load library
        try:
            self._daqp=CDLL("libdaqp.so")
        except:
            print("Could not locate libdaqp.so. Make sure DAQP is installed correctly.")
    def solve(self):
        self._daqp.daqp_solve(self.work)

    def quadprog(self, H=None,f=None,
            A=None,bupper=None,blower=None,sense=None, **settings):
        (mA, n) = np.shape(A)
        m = np.size(bupper)
        ms = mA-m
        # Setup qp, settings and result
        qp = types.QP(n,m,ms,
                np.ascontiguousarray(H).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(f).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(A).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(bupper).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(blower).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(sense).ctypes.data_as(POINTER(c_int)))
        daqp_options = types.DAQPSettings()
        # Create struct to put result in
        result = types.DAQPResult()
        x = np.zeros([n,1]);
        result.x = np.ascontiguousarray(x).ctypes.data_as(POINTER(c_double));
        # Call C api
        self._daqp.daqp_quadprog(byref(result),byref(qp),byref(daqp_options))
        # Collect results 
        profiling = {'solve_time':result.solve_time,
                'setup_time': result.setup_time,
                'inner_iter': result.iter,
                'outer_iter': result.outer_iter}
        return x, result.fval, result.exitflag, profiling 

    def linprog(self, f=None,
            A=None,bupper=None,blower=None,sense=None, **settings):
        self.quadprog(None,f,A,bupper,blower,sense, settings)
