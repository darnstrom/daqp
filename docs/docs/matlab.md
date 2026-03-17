---
layout: page
title: MATLAB 
permalink: /start/matlab
nav_order: 3
nav_icon: matlab
parent: Interfaces 
math: mathjax3
---


## Setting up the problem
In MATLAB we define the problem as 
```matlab
% Define the problem
H = eye(2);
f = [1;1]; 
A = [1 2; 1 -1];
bupper = [1; 2; 3; 4];
blower = [-1; -2; -3; -4];
sense = zeros(4,1,'int32');
```
`sense` determines the type of the constraints (more details are given [here](/daqp/parameters/#constraint-classification)).

Note: When $$b_u$$ and $$b_l$$ have more elements than the number of rows in $$A$$, the first elements in $$b_u$$ and $$b_l$$ are interpreted as simple bounds. 

## Calling DAQP
There are two ways of calling DAQP in MATLAB. The first way is through a quadprog call: 
```matlab
[x,fval,exitflag,info] = daqp.quadprog(H,f,A,bupper,blower,sense);
```
This will solve the problem with default settings. A more flexible interface is also offered, where we first setup the problem and then solve it:
```matlab
d = daqp();
d.setup(H,f,A,bupper,blower,sense);
[x,fval,exitflag,info] = d.solve();
```
Using the `daqp` object is the recommended approach for embedded or real-time applications because
it allocates memory only once during `setup` and reuses it across all subsequent `solve` calls —
no heap allocations occur during solving.

If the problem data changes between solves (e.g., updated cost vector or bounds in an MPC loop),
call `setup` again on the same object to update the internal workspace:
```matlab
% Update bounds and re-solve (workspace is reused)
bupper_new = [1; 2; 2; 4];
blower_new = -bupper_new;
d.setup(H, f, A, bupper_new, blower_new, sense);
[x, fval, exitflag, info] = d.solve();
```

## Changing settings
If we, for example, want to change the maximum number of iterations to 2000 we can do so by
```matlab
d.settings('iter_limit',2000)
```

A full list of available settings is provided [here](/daqp/parameters/#settings).

## Using DAQP in YALMIP
DAQP can also be interfaced to [YALMIP](https://yalmip.github.io/). The following code sets up and solves the problem considered above

```matlab
% Setup problem
x = sdpvar(2,1);
cons = [-1 <= x(1) <=1,...
        -2 <= x(2) <=2,...
        -3 <= x(1)+2*x(2) <= 3,...
        -4 <= x(1)-x(2) <= 4];
obj = 0.5*x'*x+x(1)+x(2);
options = sdpsettings('solver','daqp');

% Solve problem
sol = optimize(cons,obj,options);
```
