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

## <img src="/daqp/assets/icons/julia.svg" height="18" alt="Julia"> LinearMPC.jl (Julia)

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
# Continuous time system dx = A x + B u, y = C x
A = [0 1 0 0; 0 -10 9.81 0; 0 0 0 1; 0 -20 39.24 0]; 
B = 100*[0;1.0;0;2.0;;];
C = [1.0 0 0 0; 0 0 1.0 0];

# create an MPC control with sample time 0.01, prediction/control horizon 50/5
Ts = 0.01
mpc = LinearMPC.MPC(A,B,Ts;C,Np=50,Nc=5);

# set the objective functions weights
set_objective!(mpc;Q=[1.2^2,1], R=[0.0], Rr=[1.0])

# set actuator limits
set_bounds!(mpc; umin=[-2.0],umax=[2.0])
# additional functions for adding constraints: set_output_bounds!, add_constraint!

```

### Code generation

Once the controller is set up, generate embeddable C code with:

```julia
LinearMPC.codegen(mpc;dir="codgen_dir")
```

This writes a self-contained C implementation (using DAQP) to the `mpc_code/` directory, ready
to be compiled and deployed on embedded targets.

---

## <img src="/daqp/assets/icons/python.svg" height="18" alt="Python"> lmpc (Python)

[lmpc](https://github.com/darnstrom/lmpc) is the Python counterpart to LinearMPC.jl. It provides
the same workflow — design an MPC controller in Python, then generate embedded C code powered by
DAQP.

### Installation
```bash
pip install lmpc
```

### Basic usage

```python
import numpy
from lmpc import MPC,ExplicitMPC

# Continuous time system dx = A x + B u
A = numpy.array([[0, 1, 0, 0], [0, -10, 9.81, 0], [0, 0, 0, 1], [0, -20, 39.24, 0]])
B = 100*numpy.array([0,1.0,0,2.0])
C = numpy.array([[1.0, 0, 0, 0], [0, 0, 1.0, 0]])


# create an MPC control with sample time 0.01, prediction horizon 10 and control horizon 5 
Np,Nc = 10,5
Ts = 0.01
mpc = MPC(A,B,Ts,C=C,Nc=Nc,Np=Np);

# set the objective functions weights
mpc.set_objective(Q=[1.44,1],R=[0.0],Rr=[1.0])

# set actuator limits
mpc.set_bounds(umin=[-2.0],umax=[2.0])
```

### Code generation

```python
mpc.codegen(dir="codgen_dir")
```

As with LinearMPC.jl, this produces a standalone C implementation using DAQP that can be compiled and flashed to embedded hardware without any external dependencies.
