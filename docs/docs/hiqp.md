---
layout: page
title: Hierarchical QP
permalink: /start/advanced/hiqp
nav_order: 2
parent: Advanced Problem Types
math: mathjax3
---

In a hierarchical QP, constraints are organized into $$n_h$$ priority levels. DAQP first minimizes
infeasibility at the highest-priority level, then satisfies lower-priority constraints to the extent
permitted by the higher-priority solution.

The hierarchy is specified through a `break_points` vector of length $$n_h + 1$$, where
`break_points[i]` is the (0-indexed) row of `A` at which priority level $$i$$ begins.

### Julia
```julia
using DAQP

# Two priority levels:
#   Level 1 (high): -1 <= x1 + x2 <= 1
#   Level 2 (low):  -3 <= x1 - x2 <= 3
A = [1.0 1.0; 1.0 -1.0];
bupper = [1.0; 3.0];
blower = [-1.0; -3.0];
sense  = zeros(Cint, 2);
break_points = Cint[0; 1; 2];   # level 1: row 0, level 2: row 1

d = DAQP.Model();
DAQP.setup(d, zeros(0,0), zeros(0), A, bupper, blower, sense;
           break_points = break_points);
x, fval, exitflag, info = DAQP.solve(d);
```

### MATLAB
```matlab
% hidaqp takes cell arrays of A, bu, bl for each priority level
As  = {[1 1]; [1 -1]};
bus = {1;  3};
bls = {-1; -3};

[x, es, exitflag, info] = daqp.hidaqp(As, bus, bls);
% es{i} contains the constraint violation at priority level i
```

### Python
```python
import daqp, numpy as np
from ctypes import c_double, c_int

A       = np.array([[1, 1], [1, -1]], dtype=c_double)
bupper  = np.array([1, 3], dtype=c_double)
blower  = np.array([-1, -3], dtype=c_double)
sense   = np.zeros(2, dtype=c_int)
bp      = np.array([0, 1, 2], dtype=np.intc)   # break points

x, fval, exitflag, info = daqp.solve(
    None, np.zeros(2), A, bupper, blower, sense, break_points=bp)
```
