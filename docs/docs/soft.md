---
layout: page
title: Soft Constraints
permalink: /start/advanced/soft
nav_order: 4
parent: Advanced Problem Types
math: mathjax3
---

A constraint that is marked **soft** may be violated, at a cost. Softening a constraint guarantees
that the problem stays feasible, and the penalty decides how the unavoidable violation is shared
among the soft constraints.

## The penalty

With soft constraints being given by the index set $$\mathcal{S}$$, DAQP solves

$$
\begin{aligned}
\min_{x,\,s_l,\,s_u}\quad & \tfrac{1}{2}x^\top H x + f^\top x +
\sum_{i\in\mathcal{S}} \left( w_{l,i}\, s_{l,i} + \frac{s_{l,i}^2}{2\rho_{l,i}}
+ w_{u,i}\, s_{u,i} + \frac{s_{u,i}^2}{2\rho_{u,i}} \right) \\
\text{subject to}\quad & b_l - s_l \leq A x \leq b_u + s_u, \qquad s_l,\, s_u \geq 0,
\end{aligned}
$$

so each side of a soft constraint has two weights:

| | meaning | effect |
|:--|:--|:--|
| $$w$$ | linear (L1) weight | the multiplier has to exceed $$w$$ before the constraint is violated at all, so a large enough $$w$$ makes the penalty *exact*: the constraint is satisfied whenever that is possible |
| $$\rho$$ | **reciprocal** quadratic (L2) weight | a larger $$\rho$$ permits more violation; $$\rho \to 0$$ recovers a hard constraint |

The penalty is continuously differentiable, and the slack is zero while $$\lvert \lambda \rvert \leq w$$
and $$\rho(\lvert \lambda \rvert - w)$$ beyond that, where $$\lambda$$ is the multiplier of the constraint.

## Marking a constraint as soft

Set bit 8 (`DAQP_SOFT`) in `sense` for the constraints that may be violated.

### <img src="{{ '/assets/icons/c.svg' | relative_url }}" class="nav-icon" alt="C"> C
```c
int sense[3] = {0, DAQP_SOFT, DAQP_SOFT}; // rows 1 and 2 may be violated
```
### <img src="{{ '/assets/icons/julia.svg' | relative_url }}" class="nav-icon" alt="Julia"> Julia
```julia
sense = Cint[0, 8, 8]
```
### <img src="{{ '/assets/icons/matlab.svg' | relative_url }}" class="nav-icon" alt="MATLAB"> MATLAB
```matlab
d.soften_constraints([2 3]); % or set sense(i) = 8 directly
```
### <img src="{{ '/assets/icons/python.svg' | relative_url }}" class="nav-icon" alt="Python"> Python
```python
sense = np.array([0, 8, 8], dtype=np.intc)
```

## Uniform weights

By default every soft constraint uses the settings `rho_soft` and `w_soft`
(see [Settings]({{ '/start/settings' | relative_url }})). With the default `w_soft = 0` the penalty is
purely quadratic, and `rho_soft = 1e-6` keeps the violation small, i.e. the constraint is
*almost* hard.

```julia
DAQP.settings(d, Dict(:rho_soft => 1e-4, :w_soft => 1e3))
```

> These two are given in the **normalized** formulation that the solver works in, where the rows of
> the constraint matrix have unit norm. The individual weights below are instead given in the scale
> of the original problem, so the same number means the same thing only when the rows of $$A$$
> (after the Hessian factor) already have unit norm.

## Individual weights

To weight the constraints differently, set one weight per constraint and side. Only the entries that
are nonzero take effect; a zero entry falls back to `rho_soft`/`w_soft`.

### <img src="{{ '/assets/icons/c.svg' | relative_url }}" class="nav-icon" alt="C"> C
```c
#include "api.h"

DAQPWorkspace work = {0};
setup_daqp(&qp, &work, NULL);

double rho_l[3] = {0, 1e-2, 1e-4};  // reciprocal quadratic weights
double rho_u[3] = {0, 1e-2, 1e-4};
double w_l[3]   = {0, 0.5,  10.0};  // linear weights
double w_u[3]   = {0, 0.5,  10.0};
if(!daqp_set_soft_weights(&work, rho_l, rho_u, w_l, w_u))
    printf("built with DAQP_NO_SOFT_WEIGHTS\n");

DAQPResult result = {x, lam};
daqp_solve(&result, &work);
```
A `NULL` argument leaves that weight alone, so `daqp_set_soft_weights(&work, rho, rho, NULL, NULL)`
sets only the quadratic weights. The weights may also be written directly into `work.rho_ls`,
`work.rho_us`, `work.w_ls` and `work.w_us` after `daqp_allocate_soft_weights(&work)`.

### <img src="{{ '/assets/icons/julia.svg' | relative_url }}" class="nav-icon" alt="Julia"> Julia
```julia
d = DAQP.Model()
DAQP.setup(d, H, f, A, bupper, blower, sense)
soft_weights(d; rho_l = [0, 1e-2, 1e-4], rho_u = [0, 1e-2, 1e-4],
                w_l   = [0, 0.5,  10.0], w_u   = [0, 0.5,  10.0])
x, fval, exitflag, info = DAQP.solve(d)
```
### <img src="{{ '/assets/icons/matlab.svg' | relative_url }}" class="nav-icon" alt="MATLAB"> MATLAB
```matlab
d.soft_weights([0 1e-2 1e-4], [0 1e-2 1e-4], [0 0.5 10], [0 0.5 10]);
% [] leaves that weight at its default, e.g. d.soft_weights(rho_l, rho_u);
```
### <img src="{{ '/assets/icons/python.svg' | relative_url }}" class="nav-icon" alt="Python"> Python
```python
d = daqp.Model()
d.setup(H, f, A, bupper, blower, sense)
d.soft_weights(rho_l=[0, 1e-2, 1e-4], rho_u=[0, 1e-2, 1e-4],
               w_l=[0, 0.5, 10.0],    w_u=[0, 0.5, 10.0])
x, fval, exitflag, info = d.solve()
```
### C++ (Eigen)
```cpp
DAQP solver(n, m, m);
solver.update(H, f, A, bupper, blower, sense, break_points);
Eigen::VectorXd rho(3), w(3);
rho << 0, 1e-2, 1e-4;
w   << 0, 0.5,  10.0;
solver.set_soft_weights(rho, rho, w, w);  // false if built without support
solver.solve();
```

Individual weights are part of every default build. They can be compiled out with
`-DDAQP_NO_SOFT_WEIGHTS` (`cmake -DSOFT_WEIGHTS=OFF`), which makes every soft constraint use
`rho_soft`/`w_soft`; the functions above then report that the weights are unavailable instead of
silently ignoring them.

## The resulting violation

If a soft constraint ends up violated, the solver returns the exit flag `2` (`SOFT_OPTIMAL`) instead
of `1`, and `soft_slack` in the result holds the largest violation, in the units of the original
problem.

## Mapping from a nominal slack bound

A formulation that penalizes a slack with its own lower bound, as in
[acados]({{ '/start/acados' | relative_url }}),

$$\min\ \ldots + z\,s + \tfrac{1}{2}Z s^2 \quad \text{s.t.}\quad A x \leq b_u + s,\quad s \geq d,$$

is obtained by substituting $$s = d + s_u$$, which gives

$$b_u \mathrel{+}= d, \qquad w_u = \max(0,\, z + Z d), \qquad \rho_u = 1/Z.$$

The $$\max$$ only guards against a negative linear weight, which would make the slack unbounded.
