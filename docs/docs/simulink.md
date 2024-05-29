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
In Simulink, an s-function block can be used to solve a quadratic program. The s-function needs to link to the compiled c-code called "DAQP_sfunc". To compile the c-code, the following matlab-script can be used in the terminal:

```matlab
make_sfunc()
```

_Note_: Please make sure to have a C-compiler installed on your system and linked in matlab. This can be checked by running ``mex -setup`` in the terminal:

After the c-code has been compiled, a new file called ``DAQP_sfunc.mexw64`` (windows) or ``DAQP_sfunc.mexa64`` (linux,mac) should be found.

The s-function block has the following inputs:
| Inputs | Description | Original Size | Vectorized Size |
| --- | --- | --- | --- |
| $H$ | The Hessian matrix of the quadratic objective function | [$n$ x $n$] | [$n^2$ x 1] |
| $f$ | The linear part of the objective function | [$n$ x 1] | - |
| $A$ | The matrix of linear constraints | [$m_g$ x $n$] | [$m_g\cdot n$ x 1] |
| $b_l$ | The lower bound of the linear constraints | [$m$ x 1] | - |
| $b_u$ | The upper bound of the linear constraints | [$m$ x 1] | - |
| sense | The type of the constraints | [$m$ x 1] | - |

_Note_: When $b_u$ and $b_l$ have more elements than the number of rows in $A$, the first elements in $b_u$ and $b_l$ are interpreted as simple bounds. 

_Note_: The $A$ matrix (and $H$, however $H$ is typically symmetric) is **row-major-order**, meaning that the first $n$ elements are the first row of $A$, the next $n$ elements are the second row of $A$, and so on.

The following parameters are used to set up the problem and need to be set in the mask of the s-function block:

| Parameters | Description | Size |
| --- | --- | --- |
| $n$ | Number of decision variables | 1 |
| $m_g$ | Number of general constraints | 1 |
| maxIter | Maximum number of iterations | 1 |

The block has the following outputs:

| Outputs | Description | Size |
| --- | --- | --- |
| $x$ | The optimal solution | [$n$ x 1] |
| $\lambda$ | The optimal Lagrange multipliers | [$m_g$ x 1] |
| fval | The optimal value of the objective function | 1 |
| exitflag | The exit flag of the solver | 1 |
| iter | The number of iterations used by the solver | 1 |


## Example

An example can be found in the simulink file ``simulink_example.slx``. The example solves the following quadratic program:

$$\begin{aligned}
&\underset{x}{\text{minimize}}&& \frac{1}{2} \left(x_1^2 + x_2^2\right)\\
&\text{subject to} && 0.5 \leq \:x_1 \:\leq inf, \\
& && 2 \leq \:x_2 \:\leq inf, \\
& && 2.6 \leq \:x_1 + x_2 \: \leq inf, \\
\end{aligned}$$

The solution of the problem is $x = [0.6, 2]^T$.