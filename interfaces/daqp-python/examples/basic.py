## Test quadprog
import daqp
import numpy as np
from ctypes import * 
import ctypes.util

H = np.array([[1, 0], [0, 1]],dtype=c_double)
f = np.array([2, 2],dtype=c_double)
A = np.array([[1, 0], [0, 1]],dtype=c_double)
bupper = np.array([1,1],dtype=c_double)
blower= np.array([-1,-1],dtype=c_double)
sense = np.array([0,0],dtype=c_int)

x,fval,exitflag,info = daqp.solve(H,f,A,bupper,blower,sense)
print("Optimal solution:")
print(x)
print("Exit flag:",exitflag)
print("Info:",info)

## Solve persistent using Model (allocates workspace once, reuses across solves)
d = daqp.Model()
exitflag, setup_time = d.setup(H, f, A, bupper, blower, sense)
print("\n--- Model interface ---")
print("Setup exit flag:", exitflag)

x, fval, exitflag, info = d.solve()
print("Optimal solution:", x)
print("Exit flag:", exitflag)
print("Info:", info)

# Update cost vector and re-solve (no reallocation)
f = -1 * f
d.update(f=f)
x, fval, exitflag, info = d.solve()
print("\nAfter updating f:")
print("Optimal solution:", x)
print("Exit flag:", exitflag)

# Change solver settings
d.settings = {'iter_limit': 100, 'primal_tol': 1e-8}
print("\nUpdated settings:", d.settings)
