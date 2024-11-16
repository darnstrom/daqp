---
layout: page
title: Simulink 
permalink: /start/simulink
nav: 3 
parent: Interfaces 
grand_parent: Getting started 
math: mathjax3
---


## Setting up the problem

In Simulink, an s-function block can be used to solve a quadratic program. The s-function needs to link to the compiled c-code called "daqp_sfunc". To compile the c-code, the following matlab-script can be used in the terminal:

```matlab
make_sfunc()
```

_Note_: Please make sure to have a C-compiler installed on your system and linked in matlab. This can be checked by running ``mex -setup`` in the terminal:

After the c-code has been compiled, a new file called ``daqp_sfunc.mexw64`` (windows) or ``daqp_sfunc.mexa64`` (linux, mac) should be found.

The corresponding S-function block has the following inputs:

| Inputs | Description | Size | 
| --- | --- | --- |
| $$H$$ | The Hessian matrix of the quadratic objective function | [ $$n$$ x $$n$$ ] |
| $$f$$ | The linear part of the objective function | [ $$n$$ x 1 ] |
| $$A$$ | The matrix of linear constraints | [ $$m_g$$ x $$n$$ ] |
| $$b_l$$ | The lower bound of the linear constraints | [ $$m$$ x 1 ] |
| $$b_u$$ | The upper bound of the linear constraints | [ $$m$$ x 1 ] |
| sense | The type of the constraints | [ $$m$$ x 1 ] |

_Note_: When $$b_u$$ and $$b_l$$ have more elements than the number of rows in $$A$$ the first elements in $$b_u$$ and $$b_l$$ are interpreted as simple bounds. 

The following parameters are used to set up the problem and need to be set in the mask of the s-function block:

| Parameters | Description | Size |
| --- | --- | --- |
| maxIter | Maximum number of iterations | 1 |

The block has the following outputs:

| Outputs | Description | Size |
| --- | --- | --- |
| $$x$$ | The optimal solution | [ $$n$$ x 1 ] |
| $$\lambda$$ | The optimal Lagrange multipliers | [ $$m$$ x 1 ] |
| $$f_{val}$$ | The optimal value of the objective function | 1 |
| exitflag | The exit flag of the solver | 1 |
| iter | The number of iterations used by the solver | 1 |
