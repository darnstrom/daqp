---
layout: page
title: Advanced Problem Types
permalink: /start/advanced
nav_order: 5
math: mathjax3
---

Beyond standard convex QPs, DAQP supports three additional problem classes: mixed-integer QPs,
hierarchical QPs, and affine variational inequalities. All are accessed through the same interface
as regular QPs by supplying additional arguments.

## Mixed-Integer QP (MIQP)

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

---

## Hierarchical QP (HQP)

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

---

## Affine Variational Inequality (AVI)

DAQP can find a solution to the **affine variational inequality** (AVI):

$$\text{find } x^\star \in C \text{ such that } \langle Hx^\star + f,\, y - x^\star\rangle \geq 0 \quad \forall\, y \in C,$$

where $$C = \lbrace x \mid b_l \leq Ax \leq b_u \rbrace$$. Unlike a QP, $$H$$ need not be symmetric
or positive definite; AVI generalizes both QPs and complementarity problems.

### Julia
```julia
using DAQP

H = [1.0 1.75; 0.0 1.0];   # asymmetric
f = [2.0; 2.0];
A = [1.0 0.0; 0.0 1.0];
bupper = [1.0; 1.0];
blower = [-1.0; -1.0];

x, _, exitflag, info = DAQP.avi(H, f, A, bupper, blower)
```

### MATLAB
```matlab
H = [1 1.75; 0 1];
f = [2; 2];
A = eye(2);
bupper = [1; 1];
blower = [-1; -1];
sense  = int32(zeros(2, 1));

[x, fval, exitflag, info] = daqp.avi(H, f, A, bupper, blower, sense);
```

### Python
```python
import daqp, numpy as np
from ctypes import c_double

H      = np.array([[1, 1.75], [0, 1]], dtype=c_double)
f      = np.array([2, 2], dtype=c_double)
A      = np.eye(2, dtype=c_double)
bupper = np.array([1, 1], dtype=c_double)
blower = np.array([-1, -1], dtype=c_double)

x, fval, exitflag, info = daqp.solve(H, f, A, bupper, blower, is_avi=True)
```
