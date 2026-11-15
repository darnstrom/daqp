---
layout: page
title: Simulink 
permalink: /start/simulink
nav_order: 5
nav_icon: simulink
parent: Interfaces 
math: mathjax3
---


## Setting up the problem

In Simulink, an S-Function block can be used to solve a quadratic program at every time step. The S-function is called `DAQP_sfunc` and is compiled by running the following in MATLAB (from `interfaces/daqp-simulink`):

```matlab
make_sfunc()
```

_Note_: Please make sure to have a C-compiler installed on your system and linked in MATLAB. This can be checked by running ``mex -setup``.

After the compilation, a file called ``DAQP_sfunc.mexw64`` (Windows), ``DAQP_sfunc.mexa64`` (Linux) or ``DAQP_sfunc.mexmaca64``/``DAQP_sfunc.mexmaci64`` (Mac) is created. Make sure that `interfaces/daqp-simulink` is on the MATLAB path, and set the S-function name of the S-Function block to `DAQP_sfunc`. For code generation, `rtwmakecfg.m` adds the DAQP sources to the build.

An example model, where a point tracks a moving reference subject to constraints, is built and simulated by running

```matlab
addpath examples
daqp_example()      % daqp_example(true) also saves it as daqp_example.slx
```

The S-function block solves

$$
\begin{aligned}
\min_{x} &\quad \frac{1}{2} x^T H x + f^T x \\
\text{s.t.} &\quad b_l \leq \begin{bmatrix} I & 0 \\ & A\end{bmatrix} x \leq b_u
\end{aligned}
$$

and has the following inputs:

| Inputs | Description | Size |
| --- | --- | --- |
| $$H$$ | The Hessian matrix of the quadratic objective function | [ $$n$$ x $$n$$ ] |
| $$f$$ | The linear part of the objective function | [ $$n$$ x 1 ] |
| $$A$$ | The matrix of general linear constraints (optional) | [ $$m_g$$ x $$n$$ ] |
| $$b_l$$ | The lower bounds of the constraints (optional, $$-\infty$$ if not connected) | [ $$m$$ x 1 ] |
| $$b_u$$ | The upper bounds of the constraints (optional, $$+\infty$$ if not connected) | [ $$m$$ x 1 ] |
| sense | The type of the constraints (optional, 0 if not connected), see [constraint classification](/daqp/parameters/#constraint-classification) | [ $$m$$ x 1 ] |

_Note_: When $$b_u$$ and $$b_l$$ have more elements than the number of rows in $$A$$, the first $$m-m_g$$ elements in $$b_u$$ and $$b_l$$ are interpreted as simple bounds on the first $$m-m_g$$ elements of $$x$$.

The following parameters are set in the S-function parameters field of the block, separated by commas (for example `1000, 1`):

| Parameters | Description | Default |
| --- | --- | --- |
| maxIter | Maximum number of iterations | (required) |
| warm_start | 0: every solve starts from the sense input (cold start). 1: every solve starts from the working set of the previous solve (warm start). 2: as 1, but $$H$$ and $$A$$ are assumed to be constant and are only read at the first time step.  | 0 |
| break_points | Break points for [hierarchical QPs](/daqp/start/advanced/hiqp): `break_points(i)` is the index of the last constraint in priority level `i` (the constraints in the first level are hard), and the last element must equal $$m$$. Leave empty (`[]`) for normal QPs. | [] |

A change in the sense input always restarts the solver from the new sense, also when warm starting.

The block has the following outputs:

| Outputs | Description | Size |
| --- | --- | --- |
| $$x$$ | The optimal solution | [ $$n$$ x 1 ] |
| $$\lambda$$ | The optimal Lagrange multipliers | [ $$m$$ x 1 ] |
| $$f_{val}$$ | The optimal value of the objective function | 1 |
| exitflag | The exit flag of the solver (see [exit flags](/daqp/parameters/)) | 1 |
| iter | The number of iterations used by the solver | 1 |

The solution is only meaningful when exitflag is positive.
