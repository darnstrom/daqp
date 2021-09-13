---
layout: page
title: Getting started 
permalink: /start/
nav: 3 
---

## Example in C 

The following data is assumed to be available before calling daqp 
```c
// Assumed given and initalized:
int n // decision variables 
int m // number of constraints 
double* M // Stored in row major
double* d 
int* sense // 0 if inequality constraint, 1 if equality constraint 
```
This requires the QP to be transformed to an LDP (see [Background](/) for details).

If the above-mentioned problem data is available, the following code solves the corresponding LDP:

```c
// Setup workspace with problem data
Workspace work 
work.n = n;
work.m = m;
work.M = M;
work.d = d; 
work.sense = sense;

// Allocate memory for tempory variables (for example, iterates) 
allocate_daqp_workspace(&work, n);

// Activate equality constraints
add_equality_constraints(&work);

// Solve LDP
int exitflag = daqp(&work);

// optimal solution to LDP is available in work.u
```
Note that `sense` is modified during execution to contain the the optimal active set at termination. This allows for easy warm starts, but it also requires `sense` to be reset if the following problem should be cold-started. 

The memory allocated in `allocate_daqp_workspace` can be freed with the function `freeWorkspaceIters(&work)`. 
