---
layout: page
title: Python 
permalink: /start/python
nav_order: 4
nav_icon: python
parent: Interfaces 
math: mathjax3
---


## Setting up the problem
In Python we define the problem as 
```python
# Import relevant modules 
import daqp
import numpy as np
from ctypes import * 
import ctypes.util

# Define the problem
H = np.array([[1, 0], [0, 1]],dtype=c_double)
f = np.array([1, 1],dtype=c_double)
A = np.array([[1, 2], [1, -1]],dtype=c_double)
bupper = np.array([1,2,3,4],dtype=c_double)
blower = np.array([-1,-2,-3,-4],dtype=c_double)
sense = np.array([0,0,0,0],dtype=c_int)

```
`sense` determines the type of the constraints (more details are given [here](/daqp/parameters/#constraint-classification)).

Note: When $$b_u$$ and $$b_l$$ has more elements than the number of rows in $$A$$, the first elements in $$b_u$$ and $$b_l$$ are interpreted as simple bounds. 

## Calling DAQP
DAQP can be called as: 
```python
(xstar,fval,exitflag,info) = daqp.solve(H,f,A,bupper,blower,sense)
```

## Using a Workspace
The `daqp.solve` function above allocates memory on every call. For embedded or real-time
applications where the same problem structure is solved repeatedly (e.g., an MPC loop),
`daqp.Model` offers a workspace-based interface that allocates memory once during `setup`
and reuses it across all subsequent `solve` calls — no heap allocations occur during solving.

```python
d = daqp.Model()
d.setup(H, f, A, bupper, blower, sense)
xstar, fval, exitflag, info = d.solve()
```

If the problem data changes between solves, use `update` to push new data into the existing
workspace without re-running the full factorisation:
```python
# Update cost vector and bounds, then re-solve
d.update(f=f_new, bupper=bupper_new, blower=blower_new)
xstar, fval, exitflag, info = d.solve()
```
