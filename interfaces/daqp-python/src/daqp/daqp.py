import numpy as np
from ctypes import * 
import ctypes.util
import os.path
import platform
import daqp.types as types

class daqp:
    def __init__(self):
        # load library
        ROOT_PATH = os.path.dirname(os.path.abspath(__file__))
        libdaqp = [os.path.join(ROOT_PATH,f) for f in os.listdir(ROOT_PATH)
                if f.startswith('libdaqp')]
        if not libdaqp:
            print("Could not find dynamic library")
        else:
            try: 
                self._daqp=CDLL(os.path.join(ROOT_PATH,libdaqp[0]))
            except:
                print("CDLL failed")


    def solve(self):
        self._daqp.daqp_solve(self.work)

    def default_settings(self):
        opts = types.DAQPSettings()
        self._daqp.daqp_default_settings(byref(opts))
        return opts

    def quadprog(self, H=None,f=None,
            A=None,bupper=None,blower=None,sense=None, settings={}):
        (mA, n) = np.shape(A)
        m = np.size(bupper)
        ms = m-mA
        bin_ids_cand = np.where(sense&16)[0]
        nb = np.size(bin_ids_cand) 
        if nb > 0:
            bin_ids = np.array(bin_ids_cand,dtype=c_int)
        else:
            bin_ids = None

        # Setup qp
        qp = types.QP(n,m,ms,
                np.ascontiguousarray(H).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(f).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(A).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(bupper).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(blower).ctypes.data_as(POINTER(c_double)),
                np.ascontiguousarray(sense).ctypes.data_as(POINTER(c_int)),
                np.ascontiguousarray(bin_ids).ctypes.data_as(POINTER(c_int)),
                nb)
        # Setup settings
        if settings:
            opts = self.default_settings()
            for key, val in settings.items():
                setattr(opts,key,val)
            opts_ptr = byref(opts)
        else:
            opts_ptr = None

        # Create struct to put result in
        result = types.DAQPResult()
        x = np.zeros([n,1]);
        lam = np.zeros([m,1]);
        result.x = np.ascontiguousarray(x).ctypes.data_as(POINTER(c_double));
        result.lam = np.ascontiguousarray(lam).ctypes.data_as(POINTER(c_double));
        # Call C api
        self._daqp.daqp_quadprog(byref(result),byref(qp),opts_ptr)
        # Collect results 
        info = {'solve_time':result.solve_time,
                'setup_time': result.setup_time,
                'iterations': result.iter,
                'nodes': result.nodes,
                'lam': lam}
        return x, result.fval, result.exitflag, info

    def linprog(self, f=None,
            A=None,bupper=None,blower=None,sense=None, settings={}):
        self.quadprog(None,f,A,bupper,blower,sense, settings)
