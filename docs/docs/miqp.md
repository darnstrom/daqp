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

### <img src="{{ '/assets/icons/c.svg' | relative_url }}" class="nav-icon" alt="C"> C
```c
#include "api.h"

int n = 2, m = 2, ms = 2;          // ms == m: all bounds are simple
double H[4]     = {1, 0, 0, 1};
double f[2]     = {-1.5, -2.5};
double bupper[2]= {1.0, 1.0};
double blower[2]= {0.0, 0.0};
int    sense[2] = {16, 16};        // DAQP_BINARY = 16

DAQPProblem qp = {n, m, ms, H, f, NULL, bupper, blower, sense};

double x[2], lam[2];
DAQPResult result = {x, lam};
daqp_quadprog(&result, &qp, NULL);
```

### <img src="{{ '/assets/icons/julia.svg' | relative_url }}" class="nav-icon" alt="Julia"> Julia
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

### <img src="{{ '/assets/icons/matlab.svg' | relative_url }}" class="nav-icon" alt="MATLAB"> MATLAB
```matlab
H = eye(2);
f = [-1.5; -2.5];
A = zeros(0, 2);
bupper = [1; 1];
blower = [0; 0];
sense  = int32([16; 16]);   % binary

[x, fval, exitflag, info] = daqp.quadprog(H, f, A, bupper, blower, sense);
```

### <img src="{{ '/assets/icons/python.svg' | relative_url }}" class="nav-icon" alt="Python"> Python
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

### Singular Hessians
Branch and bound requires a strictly convex objective. If the Hessian is singular (for example, a
mixed-integer LP, or binary variables without a quadratic cost), the singular directions are
regularized by a small quadratic term, scaled by the bounds of the variables. The solution is then
optimal up to a small bias from this term, and the exit flag is `4` (optimal, inexact) unless the
bias is within the allowed suboptimality `abs_subopt + rel_subopt*abs(J)` (see
[settings]({{ '/parameters/' | relative_url }})). If eliminating the equality constraints makes
the Hessian positive definite, the problem is solved exactly (exit flag `1`).

### Presolve
Mixed-integer problems with big-M constraints, such as $$y \le M b$$ with a binary $$b$$, often have
weak relaxations that make branch and bound explore many nodes. Setting `presolve = 1` tightens
such problems before branch and bound:

- bound propagation tightens the simple bounds (and fixes binaries) using the constraints;
- coefficient tightening reduces the coefficients of binaries in one-sided constraints (e.g., a
  loose $$M$$) as far as the bounds of the other variables allow.

Neither changes the integer-feasible solutions, so the solution is the same, but the number of
nodes can decrease substantially. The problem data is not changed: the presolve is applied to the
internal representation of the problem when it is set up and whenever it is updated. Soft
constraints are not relied on, and immutable constraints are not changed. Equality constraints
are eliminated (see `eq_reduction`) from the presolved problem.

Presolve is enabled by default. Besides reducing the number of nodes, it makes branch and bound
more robust for loose big-M constants, which otherwise make the problem badly scaled. Its cost
is about that of setting up the problem, and it is done whenever the problem is updated (an update
of only the bounds then also forms the constraint matrix of the LDP again). For problems that are
solved repeatedly with small trees and tight big-M constants (e.g., hybrid MPC), this cost can
exceed what it saves, and `presolve = 0` can be faster. Generated code (see code generation) is
formed from the problem itself, since the presolve depends on the data at the time of generation.
