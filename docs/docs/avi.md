---
layout: page
title: Affine Variational Inequalities
permalink: /start/advanced/avi
nav_order: 3
parent: Advanced Problem Types
math: mathjax3
---

DAQP can find a solution to the **affine variational inequality** (AVI):

$$\text{find } x^\star \in C \text{ such that } \langle Hx^\star + f,\, y - x^\star\rangle \geq 0 \quad \forall\, y \in C,$$

where $$C = \lbrace x \mid b_l \leq Ax \leq b_u \rbrace$$. Unlike a QP, $$H$$ need not be symmetric; AVI generalizes both QPs and complementarity problems.

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
