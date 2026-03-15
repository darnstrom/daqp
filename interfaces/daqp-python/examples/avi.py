import daqp
import numpy as np
from ctypes import * 
import ctypes.util

H = np.array([[1, 1.75], [0, 1]],dtype=c_double)
f = np.array([2, 2],dtype=c_double)
A = np.array([[1, 0], [0, 1]],dtype=c_double)
bupper = np.array([1,1],dtype=c_double)
blower= np.array([-1,-1],dtype=c_double)

x,fval,exitflag,info = daqp.solve(H,f,A,bupper,blower,is_avi=True)
print("Optimal solution:")
print(x)
print("Exit flag:",exitflag)
print("Info:",info)


