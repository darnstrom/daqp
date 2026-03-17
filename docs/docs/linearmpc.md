---
layout: page
title: LinearMPC
permalink: /start/linearmpc
nav_order: 7
nav_icon: linearmpc
parent: Interfaces
math: mathjax3
---

DAQP can be used as the underlying QP solver in code-generation frameworks for **linear Model
Predictive Control (MPC)**. These tools let you design an MPC controller symbolically and then
automatically generate a self-contained, library-free C implementation that can be deployed on
embedded hardware.

## LinearMPC.jl (Julia)

[LinearMPC.jl](https://github.com/darnstrom/LinearMPC.jl) is a Julia package for designing and
deploying linear MPC controllers. It uses DAQP to solve the underlying condensed QP at each
sampling instant and can export the controller — including the pre-factored workspace — as
standalone C code.

### Installation
```julia
] add LinearMPC
```

### Basic usage

```julia
using LinearMPC

# Discrete-time double integrator
A = [1.0 1.0; 0.0 1.0];
B = [0.5; 1.0];
C = [1.0 0.0];

# MPC parameters
Np = 10;   # prediction horizon
Q  = [1.0 0.0; 0.0 0.0];   # state cost
R  = [0.1;;];               # input cost

ctrl = LinearMPC.LQRController(A, B, C, Q, R, Np);

# Set input and state constraints
LinearMPC.set_constraints!(ctrl; u_min=[-1.0], u_max=[1.0])

# Compute optimal input for current state x0
x0  = [1.0; 0.0];
u   = LinearMPC.control(ctrl, x0)
```

### Code generation

Once the controller is set up, export it as C code with:

```julia
LinearMPC.codegen(ctrl; dir="mpc_code")
```

This writes a self-contained C implementation (using DAQP) to the `mpc_code/` directory, ready
to be compiled and deployed on embedded targets.

---

## lmpc (Python)

[lmpc](https://github.com/darnstrom/lmpc) is the Python counterpart to LinearMPC.jl. It provides
the same workflow — design an MPC controller in Python, then generate embedded C code powered by
DAQP.

### Installation
```bash
pip install lmpc
```

### Basic usage

```python
import numpy as np
from lmpc import LQRController

# Discrete-time double integrator
A = np.array([[1.0, 1.0], [0.0, 1.0]])
B = np.array([[0.5], [1.0]])

# MPC parameters
Np = 10
Q  = np.diag([1.0, 0.0])
R  = np.array([[0.1]])

ctrl = LQRController(A, B, Q, R, Np)

# Set constraints
ctrl.set_constraints(u_min=[-1.0], u_max=[1.0])

# Compute optimal input for current state
x0 = np.array([1.0, 0.0])
u  = ctrl.control(x0)
```

### Code generation

```python
ctrl.codegen(dir="mpc_code")
```

As with LinearMPC.jl, this produces a standalone C implementation using DAQP that can be compiled
and flashed to embedded hardware without any external dependencies.
