---
layout: page
title: C/C++
permalink: /start/c
nav: 3 
parent: Interfaces 
grand_parent: Getting started 
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
