---
layout: page
title: Mixed-Integer QP
permalink: /start/advanced/miqp
nav_order: 1
parent: Advanced Problem Types
math: mathjax3
---

DAQP can solve mixed-integer QPs with binary constraints of the form

$$Ax \in \lbrace b_l,\, b_u \rbrace,$$

with $$x_i \in \lbrace 0, 1 \rbrace$$ as a special case. Binary variables are declared by setting the
corresponding entry in the `sense` vector to `16` (`DAQP_BINARY`). DAQP then applies
branch-and-bound to find the optimal integer solution.

### Julia
```julia
using DAQP

# A small MIQP: minimize 0.5 x'x + f'x  s.t. x ∈ {0,1}^2
H = [1.0 0; 0 1];
f = [-1.5; -2.5];
A = zeros(0, 2);           # no general constraints
bupper = [1.0; 1.0];      # simple upper bounds
blower = [0.0; 0.0];      # simple lower bounds
sense  = Cint[16; 16];    # both variables are binary

x, fval, exitflag, info = DAQP.quadprog(H, f, A, bupper, blower, sense)
println("Solution: ", x, "  (nodes explored: ", info.nodes, ")")
```

### MATLAB
```matlab
H = eye(2);
f = [-1.5; -2.5];
A = zeros(0, 2);
bupper = [1; 1];
blower = [0; 0];
sense  = int32([16; 16]);   % binary

[x, fval, exitflag, info] = daqp.quadprog(H, f, A, bupper, blower, sense);
```

### Python
```python
import daqp, numpy as np
from ctypes import c_double, c_int

H       = np.eye(2, dtype=c_double)
f       = np.array([-1.5, -2.5], dtype=c_double)
A       = np.zeros((0, 2), dtype=c_double)
bupper  = np.array([1, 1], dtype=c_double)
blower  = np.array([0, 0], dtype=c_double)
sense   = np.array([16, 16], dtype=c_int)   # binary

x, fval, exitflag, info = daqp.solve(H, f, A, bupper, blower, sense)
```

The `info` struct contains a `nodes` field with the number of branch-and-bound nodes explored.
