---
layout: page
title: Multistage Problems
permalink: /start/advanced/multistage
nav_order: 5
parent: Advanced Problem Types
math: mathjax3
---

DAQP has a structure-exploiting solver for multistage QPs, as in optimal control (MPC, moving horizon estimation). It runs the same dual active-set method as for dense QPs, but computes everything by Riccati recursions, so the cost grows linearly with the horizon. It is built with the CMake option `DAQP_MULTISTAGE` and requires [BLASFEO](https://github.com/giaf/blasfeo).

## Problem
The data is given per stage $$t$$ in the HPIPM/acados form, extended with constraints between consecutive stages:

$$
\begin{aligned}
\min_{x,u}\;& \sum_{t=0}^{N-1} \left(\frac{1}{2} \begin{bmatrix} u_t \\ x_t \end{bmatrix}^T \begin{bmatrix} R_t & S_t \\ S_t^T & Q_t \end{bmatrix} \begin{bmatrix} u_t \\ x_t \end{bmatrix} + \begin{bmatrix} r_t \\ q_t \end{bmatrix}^T \begin{bmatrix} u_t \\ x_t \end{bmatrix}\right) + \frac{1}{2} x_N^T Q_N x_N + q_N^T x_N\\
\text{s.t. }\;& \begin{aligned}[t]
& x_{t+1} = A_t x_t + B_t u_t + b_t, && t = 0,\dots,N-1\\
& lb_t \leq \begin{bmatrix} u_t \\ x_t \end{bmatrix}_{\text{idxb}_t} \leq ub_t, && t = 0,\dots,N\\
& lg_t \leq C_t x_t + D_t u_t \leq ug_t, && t = 0,\dots,N\\
& lc_t \leq E^x_t x_t + E^u_t u_t + F^x_t x_{t+1} + F^u_t u_{t+1} \leq uc_t, && t = 0,\dots,N-1
\end{aligned}
\end{aligned}
$$

There is no $$u_N$$, dimensions may vary over the horizon, and the problem must be convex.

## Data
The problem is passed in a `DAQPMultistage` struct (`include/multistage.h`). Per-stage arrays have length $$N+1$$.

| Field | Content |
|-------|---------|
| `N` | Horizon length |
| `nx`, `nu` | State and input dimensions per stage |
| `nb`, `ng`, `nc` | Number of bounds, general constraints, and coupling constraints per stage (`NULL`: none) |
| `A`, `B`, `b` | Dynamics (`NULL` `A`: zero, `NULL` `B`: identity, `NULL` `b`: zero) |
| `Q`, `S`, `R`, `q`, `r` | Cost (`NULL`: zero) |
| `idxb`, `lb`, `ub` | Bounds; `idxb[t]` indexes $$[u_t; x_t]$$ |
| `C`, `D`, `lg`, `ug` | General constraints |
| `Ex`, `Eu`, `Fx`, `Fu`, `lc`, `uc` | Coupling constraints between stage $$t$$ and $$t+1$$ |
| `sense` | Constraint types (`NULL`: all hard) |
| `BAbt`, `RSQrq`, `DCt`, `Et`, `Ft` | The matrices as BLASFEO matrices |

The solver reads the BLASFEO matrices, in HPIPM/acados format with stage variables ordered $$[u_t; x_t]$$: `BAbt[t]` $$= [B_t^T; A_t^T]$$, `RSQrq[t]` $$= [R_t\ S_t;\ S_t^T\ Q_t]$$ (lower triangle), `DCt[t]` $$= [D_t^T; C_t^T]$$, `Et[t]` $$= [E^{u\,T}_t; E^{x\,T}_t]$$, `Ft[t]` $$= [F^{u\,T}_t; F^{x\,T}_t]$$. Only the first $$n_u+n_x$$ rows are read; $$b$$, $$q$$, $$r$$ and all other vectors are plain arrays. The matrices are referenced, not copied.

For dense row-major data, `daqp_multistage_pack` forms the BLASFEO matrices from `A`, `B`, `Q`, `S`, `R`, `C`, `D`, `Ex`, `Eu`, `Fx`, `Fu` (only the lower triangles of `Q` and `R` are read), and `daqp_multistage_pack_free` frees them. Stages with the same pointers (time-invariant problems) share one BLASFEO matrix.

### Constraint order and `sense`
Constraints are numbered stage by stage, and within a stage as [bounds; general; coupling]. The total is `daqp_multistage_m(&prob)`. `sense` uses this order:

* `DAQP_SOFT`: [soft constraint]({{ '/start/advanced/soft' | relative_url }}) (weights via `daqp_set_soft_weights`).
* `DAQP_BINARY`: binary constraint (see [Mixed-integer problems](#mixed-integer-problems)).

### Initial state and equalities
If hard bounds at stage 0 fix all of $$x_0$$ at setup (as in MPC), $$x_0$$ is treated as given. Otherwise (as in MHE), $$x_0$$ is optimized, and the cost must be strictly convex in its free components.

Bounds and general constraints with equal finite bounds at setup are eliminated in the Riccati recursion. Other equalities are kept in the working set.

## Usage
```c
#include "multistage.h"

DAQPMultistage prob = {0};
// ... N, nx, nu, nb, idxb, lb, ub, A, B, Q, R, q, r, ...
daqp_multistage_pack(&prob);                  // BLASFEO matrices from the row-major data

DAQPWorkspace work; work.settings = NULL;
setup_daqp_multistage(&prob, &work);

double z[NZ], lam[M], pi[NPI];        // NZ = daqp_multistage_nz(&prob), M = daqp_multistage_m(&prob), NPI = nx_1+...+nx_N
DAQPResult res; res.x = z; res.lam = lam;
daqp_multistage_solve(&res, &work);          // z = [x_0; u_0; x_1; u_1; ...; x_N]
daqp_multistage_duals(&work, lam, pi);       // (if needed) complete multipliers and costates

// MPC loop: new initial state
lb[0][...] = ub[0][...] = x0_new;
daqp_multistage_update(&work, 0);
daqp_multistage_solve(&res, &work);

free_daqp_multistage(&work);
daqp_multistage_pack_free(&prob);
```

`setup_daqp_multistage` returns 1 on success, `DAQP_EXIT_NONCONVEX` for a nonconvex problem, and `DAQP_EXIT_OVERDETERMINED_INITIAL` for inconsistent equalities in the working set; the workspace is then freed.

### Updates and warm starts
`daqp_multistage_update(&work, mask)` refreshes the workspace after data changes and warm starts from the previous working set. Changed vectors (`b`, `q`, `r`, bounds) are detected automatically; `mask` flags what is not:

| Flag | Meaning |
|------|---------|
| `DAQP_MULTISTAGE_UPDATE_MATRICES` | The matrices changed (including after `daqp_multistage_pack`) |

As in MPC, the working set is shifted one stage along with the problem if this fit the previous solve better than the unshifted one (no flag is needed). If the working sets alternate between solves (chattering), the one from two solves back is used instead. If $$x_0$$ was fixed at setup but no longer is, the update returns `DAQP_EXIT_UNSUPPORTED` and a new setup is needed.

### Dual variables
`res.lam` holds the multipliers of the working set, except those of eliminated equalities and of bounds fixing $$x_0$$, which are zero. `daqp_multistage_duals(&work, lam, pi)` (either pointer may be `NULL`) computes the complete `lam` and the costates $$\pi = [\pi_0; \dots; \pi_{N-1}]$$, with $$\pi_t \in \mathbb{R}^{n_{x,t+1}}$$ the multiplier of the dynamics $$x_{t+1} = \dots$$. They are valid until the next update or solve.

If an elimination of equalities turns out inaccurate, the equalities are instead kept in the working set and the problem is re-solved; `daqp_multistage_duals` then returns 2.

## Semidefinite problems
If the cost is not strictly convex in the inputs, proximal-point iterations are used, controlled by the [settings]({{ '/parameters/' | relative_url }}#settings) `eps_prox` (negative, the default: when needed; positive: always; zero: never, returning `DAQP_EXIT_NONCONVEX`) and `eta_prox`.

## Mixed-integer problems
Binary constraints in `sense` are handled by DAQP's [branch and bound]({{ '/start/advanced/miqp' | relative_url }}), with the following differences:
* `res.lam` is zero and `daqp_multistage_duals` returns 0.
* Proximal-point iterations are not supported (`DAQP_EXIT_NONCONVEX`).
* Bounds on a fixed $$x_0$$ cannot be binary, and there can be at most as many binary constraints as inputs in total.
* `binary_horizon` $$> 0$$ keeps only the binary constraints of stages $$t <$$ `binary_horizon` and relaxes the rest.
* The binary assignment of the previous solution (at which bound each binary constraint was), as is and shifted by one stage, is fixed and solved, and the better feasible objective prunes the tree. The assignment does not depend on the coordinates of the problem (as the steps of an SQP iteration).

## How it works
DAQP solves QPs as least-distance problems $$\min \tfrac{1}{2}\|v\|^2$$ s.t. $$d_l \leq Mv \leq d_u$$. For a multistage problem, $$v$$ parametrizes perturbations from the unconstrained (LQR) optimum through the Riccati recursion:

$$\Delta u_t = L_{uu,t}^{-T}(v_t - L_{xu,t}^T \Delta x_t),\qquad \Delta x_{t+1} = A_t \Delta x_t + B_t \Delta u_t,$$

with $$L_t$$ the Cholesky factor of the stage's cost-to-go Hessian. $$M$$ is never formed (`work->M` stays `NULL`):
* products $$Mv$$ come from a forward recursion;
* the row of a constraint at stage $$t$$ comes from a backward recursion from $$t$$ to 0, and is only computed when the constraint enters the working set.

Equalities are eliminated stage by stage through a pivoted QR factorization, a free $$x_0$$ is handled by an extra leading stage, and branch and bound runs unchanged on the implicit $$M$$.

The approach follows [DAOCP](https://arxiv.org/abs/2609.33451) (Zaupa et al., 2026). The coupling constraints follow the multistage formulation of [PIQP](https://github.com/PREDICT-EPFL/piqp).
