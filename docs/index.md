---
layout: default
title: About 
nav_order: 1
description: "DAQP is a dual active-set solver for convex quadratic programs"
permalink: /
math: mathjax3
---

<h2>Background</h2>
DAQP is an active-set solver for quadratic programs (QPs) in the form 

\begin{equation} 
\begin{aligned}
&\underset{x}{\text{minimize}}&& \frac{1}{2} x^T H x + f^T x\\
&\text{subject to} && A x \leq b,
\end{aligned}
\end{equation}

where $$H\succ 0$$. The case when $$H\succeq 0$$ can, however, also be handled by extending DAQP with proximal-point iterations. 

<h3>Transforming a QP into a least-distance problem (LDP)</h3> 
DAQP is developed to solve QPs that arise in real-time linear MPC applications. Therefore, since $$ H $$ and $$ A $$ are often constant between problem instances in those applications, most of the code is written to solve the LDP in the form  

\begin{equation} 
\begin{aligned}
&\underset{u}{\text{minimize}}&& \frac{1}{2}\|u\|_2^2 \\
&\text{subject to} && M u \leq d.
\end{aligned}
\end{equation}

Hence, the QP in (1) is assumed to be transformed into an LDP in the form (2), which can be done by using a Cholesky factor $$ R $$ for $$ H $$ (that is, an upper-triangular matrix $$ R $$ such that $$H = R^T R $$).  We then get the LDP by defining $$ M \triangleq A R^{-1} $$ and $$ d \triangleq b+M v$$, with $$v \triangleq R^{-T} f$$. 
Finally, a solution $$ x^*$$ to (1) can be obtained by transforming a solution $$ u^*$$ to the LDP in (2) through  

$$
x^* = R^{-1}(u^* - v).
$$
