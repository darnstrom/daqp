---
layout: page
title: C/C++
permalink: /start/c
nav_order: 1
nav_icon: c
parent: Interfaces 
math: mathjax3
---


## Setting up the problem
In C we define the problem as 
```c
// Define the problem
int n = 2; // Number of decision variables
int m = 4; // Number of constraints (general + simple)
int ms= 2; // Number of simple bounds
double H[4] = {1, 0, 0, 1};
double f[2] = {1,1}; 
double A[4] = {1, 2, 1, -1};
double bupper[4] = {1, 2, 3, 4};
double blower[4] = {-1, -2, -3, -4};
int sense[4] = {0,0,0,0}; // Only inequality constraints
DAQPProblem qp = {n,m,ms,H,f,A,bupper,blower,sense};
```
`sense` determines the type of the constraints (more details are given [here](/daqp/parameters/#constraint-classification)).

Note: When $$b_u$$ and $$b_l$$ have more elements than the number of rows in $$A$$, the first elements in $$b_u$$ and $$b_l$$ are interpreted as simple bounds. 

## Calling DAQP
A high-level function `daqp_quadprog` can be used to solve the problem.
```c
double x[2],lam[4];
DAQPResult result;
result.x = x; // primal variable
result.lam = lam; // dual variable
daqp_quadprog(&result,&qp,NULL);
```
The last argument is a pointer to a `DAQPSettings` struct, but passing a null-pointer will result in the default settings being used. The first argument is a pointer to a `DAQPResults` struct in which solution information will be populated. 

The optimal solution can be found in `result.x`, the optimal function value in `result.fval`, and the exit flag in `result.exitflag`. The struct `result` also contains some profiling information such as solve time and number of iterations.

## Using a Workspace
The `daqp_quadprog` function above allocates memory on every call, which may be undesirable in
real-time or embedded applications. To avoid this, DAQP provides a workspace-based interface where
memory is allocated once and reused across multiple solves.

```c
// Allocate and set up the workspace once
DAQPWorkspace work = {0};
setup_daqp(&qp, &work, NULL);   // allocates internal memory and factorizes H

// Solve
double x[2], lam[4];
DAQPResult result;
result.x   = x;
result.lam = lam;
daqp_solve(&result, &work);
```

Once the workspace is set up, problem data (such as the cost vector or bounds) can be modified and
the problem re-solved without any additional allocation:

```c
// Update the cost vector and re-solve
double f_new[2] = {-1, 0};
qp.f = f_new;
daqp_update_ldp(DAQP_UPDATE_v, &work, &qp);  // DAQP_UPDATE_v recomputes v = R'\f
daqp_solve(&result, &work);
```

When the workspace is no longer needed, free it with:
```c
free_daqp_workspace(&work);
free_daqp_ldp(&work);
```

## Changing settings
If we, for example, want to change the maximum number of iterations to 2000 we can do so by
```c
DAQPSettings settings;
daqp_default_settings(&settings); // Populate settings with default values
settings.iter_limit = 2000;
```

A full list of available settings is provided [here](/daqp/parameters/#settings).
## Using DAQP with Eigen
DAQP also has an interface to matrices/vectors defined with [Eigen](https://eigen.tuxfamily.org/). The following code sets up and solves the problem considered above

```cpp
// Setup problem
Eigen::MatrixXd H = Eigen::MatrixXd::Identity(2, 2);
Eigen::VectorXd f = Eigen::VectorXd::Ones(2);
Eigen::MatrixXd A = (Eigen::MatrixXd(2, 2) << 1, 2, 1, -1).finished();
Eigen::VectorXd bu = (Eigen::VectorXd(4) << 1,2,3,4).finished();
Eigen::VectorXd bl = (Eigen::VectorXd(4) << -1,-2,-3,-4).finished();

// Solve problem
EigenDAQPResult result = daqp_solve(H,f,A,bu,bl);
```
The primal solution can obtain with `result.get_primal()`, and the dual solution can be obtain with `result.get_dual()`.
