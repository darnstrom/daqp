---
layout: page
title: Julia 
permalink: /start/julia
nav_order: 2
nav_icon: julia
parent: Interfaces 
math: mathjax3
---


## Setting up the problem
In Julia we define the problem as 
```julia
# Define the problem
H = [1.0 0; 0 1];
f = [1.0 ;1];
A = [1.0 2; 1 -1];
bupper = [1.0 ; 2; 3; 4];
blower = [-1.0; -2; -3; -4];
sense = zeros(Cint,4);
```
`sense` determines the type of the constraints (more details are given [here](/daqp/parameters/#constraint-classification)).

Note: When $$b_u$$ and $$b_l$$ has more elements than the number of rows in $$A$$, the first elements in $$b_u$$ and $$b_l$$ are interpreted as simple bounds. 

## Calling DAQP
There are two ways of calling DAQP in Julia. The first way is through a quadprog call: 
```julia
x,fval,exitflag,info = DAQP.quadprog(H,f,A,bupper,blower,sense);
```
This will solve the problem with default settings. A more flexible interface is also offered, where we first setup the problem and then solve it 
```julia
d = DAQP.Model();
DAQP.setup(d,H,f,A,bupper,blower,sense);
x,fval,exitflag,info = DAQP.solve(d);
```
Using `DAQP.Model` is the recommended approach for embedded or real-time applications because it
allocates memory only once during `setup` and reuses it across all subsequent `solve` calls — no
heap allocations occur during solving.

If the problem data changes between solves (e.g., updated cost vector or bounds in an MPC loop),
use `update` to push the new data into the existing workspace without re-running the full setup:
```julia
# Update bounds and cost vector, then re-solve
DAQP.update(d, nothing, f_new, nothing, bupper_new, blower_new, nothing)
x, fval, exitflag, info = DAQP.solve(d)
```

## Changing settings
If we, for example, want to change the maximum number of iterations to 2000 we can do so by
```julia
DAQP.settings(d, Dict(:iter_limit =>2000))
```

A full list of available settings is provided [here](/daqp/parameters/#settings)

## Using a locally built library

`DAQPBase` normally loads the shared library supplied by `DAQP_jll`. To use C changes from a DAQP checkout, build the library and its matching Julia wrapper together. Run these commands from the checkout root.

```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DJULIA=ON
cmake --build build --config Release
```

On Windows, put MinGW-w64 and Ninja on your `PATH` and add `-G Ninja -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++` to the first command.

CMake stages the wrapper and shared library in `build/interfaces/daqp-julia`. Start Julia from the checkout root with your application's project active, then develop the staged `DAQPBase` package in that project. The active project must allow the `DAQPBase` version listed in the staged `Project.toml`.

```julia
using Pkg
Pkg.develop(path=abspath("build", "interfaces", "daqp-julia"))
```

After switching or rebuilding the library, start a fresh Julia session with the same project active. Verify that Julia selects the library in `build/interfaces/daqp-julia`.

```julia
using DAQPBase
println(realpath(DAQPBase.libdaqp))
```

## Using DAQP in JuMP
DAQP can also be interfaced to [JuMP](https://jump.dev/). The following code sets up and solves the problem considered above

```julia
using DAQP
using JuMP

## Setup problem
model = Model(DAQP.Optimizer)
@variable(model, -1<= x1 <=1)
@variable(model, -2<= x2 <=2)
@objective(model, Min, 0.5*(x1^2+x2^2)+x1+x2)
@constraint(model, c1, -3 <= x1 + 2*x2 <= 3)
@constraint(model, c2, -4 <= x1 - x2 <= 4)

## Solve problem
optimize!(model)
```
