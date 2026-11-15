---
layout: page
title: Multistage Problems
permalink: /start/advanced/multistage
nav_order: 5
parent: Advanced Problem Types
math: mathjax3
---

DAQP has a structure-exploiting solver for multistage QPs, in which the variables come in stages that are only coupled to the next stage, as in optimal control (for example MPC and moving horizon estimation). It runs the same dual active-set method as for dense QPs, but never forms the condensed problem: everything it needs is computed by Riccati recursions over the horizon, whose cost grows linearly with the horizon length. The solver is built with the CMake option `DAQP_MULTISTAGE` and requires [BLASFEO](https://github.com/giaf/blasfeo).

## Problem
The data is given per stage $$t$$, in the form of HPIPM and acados, with constraints between consecutive stages as an extension:

$$
\begin{aligned}
\min_{x,u}\;& \sum_{t=0}^{N-1} \left(\tfrac{1}{2} x_t^T Q_t x_t + u_t^T S_t x_t + \tfrac{1}{2} u_t^T R_t u_t + q_t^T x_t + r_t^T u_t\right) + \tfrac{1}{2} x_N^T Q_N x_N + q_N^T x_N\\
\text{s.t. }\;& x_{t+1} = A_t x_t + B_t u_t + b_t, && t = 0,\dots,N-1\\
& lb_t \leq [u_t; x_t]_{\text{idxb}_t} \leq ub_t, && t = 0,\dots,N\\
& lg_t \leq C_t x_t + D_t u_t \leq ug_t, && t = 0,\dots,N\\
& lc_t \leq E^x_t x_t + E^u_t u_t + F^x_t x_{t+1} + F^u_t u_{t+1} \leq uc_t, && t = 0,\dots,N-1
\end{aligned}
$$

There is no $$u_N$$. The dimensions can vary over the horizon. The problem has to be convex.

## Data
The problem is passed in a `DAQPMultistage` struct (`include/multistage.h`). All arrays of per-stage data have length $$N+1$$.

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
| `sense` | Constraint types (`NULL`: all hard), see below |
| `BAbt`, `RSQrq`, `DCt`, `Et`, `Ft` | The matrices as BLASFEO matrices |

The solver reads the matrices as BLASFEO matrices in the HPIPM/acados format, with the variables of a stage ordered as $$[u_t; x_t]$$: `BAbt[t]` $$= [B_t^T; A_t^T]$$, `RSQrq[t]` $$= [R_t\ S_t;\ S_t^T\ Q_t]$$ (lower triangle), `DCt[t]` $$= [D_t^T; C_t^T]$$, `Et[t]` $$= [E^{u\,T}_t; E^{x\,T}_t]$$ and `Ft[t]` $$= [F^{u\,T}_t; F^{x\,T}_t]$$. Only their first $$n_u+n_x$$ rows are read (HPIPM keeps $$b$$, $$r$$ and $$q$$ in an extra row; here they are always given as arrays). They are referenced, not copied, and have to stay valid. If you have dense row-major matrices instead, `daqp_multistage_pack` forms the BLASFEO matrices from `A`, `B`, `Q`, `S`, `R`, `C`, `D`, `Ex`, `Eu`, `Fx` and `Fu` (allocated at the first call, repacked in place afterwards; `daqp_multistage_pack_free` frees them). Only the lower triangles of `Q` and `R` are read. Consecutive stages given the same row-major matrices (the same pointers, as in a time-invariant problem) share one BLASFEO matrix: it is packed once, takes memory once, and the recursions read it from cache (about 15% faster products $$Mv$$ at the size of the Atlas benchmark). If that sharing changes (a stage gets its own matrix, or the pointers of stages become equal), the next pack allocates the matrices again. Vectors are always given as plain arrays.

### Constraints and `sense`
The constraints are numbered stage by stage, and within stage $$t$$ as [bounds; general constraints; coupling constraints]. Their total number is `daqp_multistage_m(&prob)`. The optional `sense` array uses this order:

* `DAQP_SOFT` marks a [soft constraint]({{ '/start/advanced/soft' | relative_url }}) (weights through `daqp_set_soft_weights`).
* `DAQP_BINARY` marks a binary constraint: its value is at its lower or its upper bound (see [Mixed-integer problems](#mixed-integer-problems)).

### Initial state
If the hard bounds at stage 0 fix all of $$x_0$$ at setup (equal lower and upper bounds), as for the initial state in MPC, $$x_0$$ is treated as given, and the other bounds on $$x_0$$ alone are only checked. Otherwise, as in moving horizon estimation, $$x_0$$ is optimized, which requires the cost to be strictly convex in its free components.

### Equality constraints
Bounds and general constraints with equal (finite) lower and upper bounds at setup are hard equalities that are eliminated in the Riccati recursion. Other equalities (coupling constraints, or bounds that only become equal after setup) are kept in the working set.

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

`daqp_multistage_update` refreshes the workspace after the data has changed, keeping the working set as a warm start. The vectors (`b`, `q`, `r` and the bounds, including the initial state) are always read again, and only the stages where they changed are recomputed; the mask (0, or a combination of the flags) gives what the update cannot detect:

| Flag | Meaning |
|------|---------|
| `DAQP_MULTISTAGE_UPDATE_MATRICES` | The matrices have changed (refactorizes) |
| `DAQP_MULTISTAGE_UPDATE_SHIFT` | The problem has moved by one stage, as in MPC: the warm start may be shifted with it |

A change of which bounds are equal is handled as a change of the matrices, and a new value of an eliminated equality is detected. If $$x_0$$ was fixed at setup but no longer is, the update returns `DAQP_EXIT_UNSUPPORTED` and a new setup is needed. If you repack with `daqp_multistage_pack`, pass `DAQP_MULTISTAGE_UPDATE_MATRICES`.

The warm start is the working set of the latest solve, or that of two solves back if the working sets alternate: in a closed loop that chatters (for example at saturated inputs), the latest working set is then the worst guess. They are taken to alternate if the working set of two solves back fit the latest two solutions clearly better than the previous one did (differing in at most half as many constraints, from at least two).

In MPC, the next problem is often the latest one moved by one stage: its initial state is (close to) the next state of the latest solution. Pass `DAQP_MULTISTAGE_UPDATE_SHIFT` in the mask of such updates (in SQP, only at the first QP of a time step). The shifted working set gives a constraint of stage $$t$$ the state in the working set (and the multiplier) of the corresponding constraint of stage $$t+1$$ of the latest solution: a bound on the same entry of $$[u; x]$$, and the general and coupling constraints with the same index, if the two stages have as many; the others, and those of the last stage, keep theirs. Whether this is the better guess depends on the problem (whether the active constraints move along the horizon with time), so the shift is only made if the shifted working set of the previous solve fit the latest solution better than the unshifted one did (and the working sets do not alternate). Where it applies, it can save most of the iterations, in particular with soft constraints: for a chain of masses with soft state bounds (12 states, 3 inputs, $$N = 30$$), about 18 instead of 148 iterations per solve with a quadratic penalty, and 90 instead of 365 with a large linear one. A warm start that is activated (after a factorization, a shift, or a change of the working set) starts from the multipliers of the latest solve, where they fit the state of each constraint, instead of from the bounds of their ranges.

After a factorization, the working set of the warm start is factorized at once, as the Cholesky factorization of $$M_A M_A^T$$ from the rows of its constraints (by blocked BLASFEO routines, and using that a row of a constraint at stage $$t$$ ends with $$v_t$$), which is several times faster than adding the constraints one by one for large working sets (for example in SQP iterations with many active constraints). This takes no extra memory: the Cholesky factor is formed in the memory of the LDL factorization (for up to about $$0.7 n_v$$ constraints), else in the row cache with the rows moved to the memory of the LDL factorization; if neither fits, the longest leading part of the working set that does is factorized at once and the rest is added one by one. A small pivot (where the additions one by one would reorder or drop a constraint) makes all constraints be added one by one. Soft constraints are included (a free slack adds its weight to the diagonal, so that the working set can hold more than $$n_v$$ constraints); small working sets and branch and bound use the one by one additions.

DAQP never forms $$M$$: it only computes the rows of $$M$$ that the working set needs, and products $$Mv$$, by Riccati recursions (`work->M` stays `NULL`).

`setup_daqp_multistage` returns 1 on success, `DAQP_EXIT_NONCONVEX` if the problem is not convex, and `DAQP_EXIT_OVERDETERMINED_INITIAL` for inconsistent equalities in the working set (inconsistent eliminated equalities make the solve return `DAQP_EXIT_INFEASIBLE`). On failure, the workspace is freed.

`res.lam` contains the multipliers of the working set, in the order of the constraints. The multipliers of eliminated equalities and of the bounds that fix $$x_0$$ are zero there: they follow from one backward recursion, which the solve leaves out so that a controller that only uses the primal solution does not pay for it. `daqp_multistage_duals(&work, lam, pi)` (either pointer can be `NULL`) runs this recursion at its first call after a solve and gives the complete multipliers `lam` and the multipliers of the dynamics $$\pi = [\pi_0; \dots; \pi_{N-1}]$$, with $$\pi_t$$ of length $$n_{x,t+1}$$ and $$\nabla_{x_{t+1}} L = \pi_t$$ for the Lagrangian $$L$$ (objective plus $$\lambda^T$$constraints plus $$\pi^T$$dynamics). They are valid until the next update or solve. An SQP method, which needs the complete dual solution, calls it after every solve. The same pass checks the eliminated equalities: the stationarity in the inputs they leave free shows whether the elimination was accurate (it is not for a large gain from the state to the inputs that the equalities determine, amplified along the horizon). The first solve after a factorization with eliminated equalities runs this pass itself (a later `daqp_multistage_duals` reuses it), and if the elimination was inaccurate, it keeps the equalities in the working set (until they change) and solves again. `daqp_multistage_duals` then returns 2; if it finds a later solve inaccurate, the equalities are kept in the working set from the next update on.

## Semidefinite problems
If the problem is not strictly convex in the inputs (the Riccati recursion breaks down), it is solved by proximal-point iterations: $$\epsilon I$$ is added to the Hessian, and the problem is re-solved around the latest solution until it no longer changes. This is controlled by the [settings]({{ '/parameters/' | relative_url }}#settings) `eps_prox` (negative: only when needed, the default; positive: always; zero: never, the solve then returns `DAQP_EXIT_NONCONVEX`) and `eta_prox`. The solution is then stationary to about `1e-6`.

## Mixed-integer problems
With binary constraints in `sense`, the problem is solved by DAQP's [branch and bound]({{ '/start/advanced/miqp' | relative_url }}). Then:
* `res.lam` is zero and `daqp_multistage_duals` returns 0 (no multipliers or costates).
* Problems that need proximal-point iterations are not supported (`DAQP_EXIT_NONCONVEX`).
* Bounds on a fixed $$x_0$$ cannot be binary, and the number of binary constraints can be at most the total number of inputs.
* `binary_horizon` $$> 0$$ keeps the binary constraints of the stages $$t <$$ `binary_horizon` only, and relaxes those of the later stages to their ranges (set at setup). In MPC, where only the first input is applied, this reduces the size of the tree with the horizon length (a horizon-tail relaxation).
* The inputs of the latest solution (one stage earlier with `DAQP_MULTISTAGE_UPDATE_SHIFT`, the last one repeated), applied from the new $$x_0$$, are a candidate: if its binary constraints are at a bound and the other constraints hold, its objective bounds the tree, which prunes the nodes that cannot be better (without eliminated equalities and a free $$x_0$$). For a closed loop of a double integrator with thrusters that have a minimum thrust ($$N = 30$$, 60 binaries), this saved about 25-70% of the time per solve.

## How it works
**The LDP.** DAQP solves QPs as least-distance problems in variables $$v$$ in which the objective is $$\tfrac{1}{2}\|v\|^2$$ and the constraints are $$d_l \leq Mv \leq d_u$$. For a multistage problem, $$v$$ comes from the Riccati recursion of the problem without inequality constraints: a perturbation of the unconstrained optimum $$z_{\text{lqr}}$$ that satisfies the dynamics (with $$\Delta x_0 = 0$$) is

$$\Delta u_t = L_{uu,t}^{-T}(v_t - L_{xu,t}^T \Delta x_t),\qquad \Delta x_{t+1} = A_t \Delta x_t + B_t \Delta u_t,$$

where $$L_t$$ is the Cholesky factor of $$R_t + B_t^T P_{t+1} B_t$$ (and the coupling to the state), and the objective increases by $$\tfrac{1}{2}\|v\|^2$$. A constraint $$a^T z$$ is then $$a^T z_{\text{lqr}} + M_a v$$. $$M$$ is never formed:
* the products $$Mv$$ for all constraints follow from a forward recursion (simulating $$v$$ through the horizon);
* the row $$M_a$$ of a constraint at stage $$t$$ follows from a backward recursion from $$t$$ to 0 (it is only nonzero in $$v_0,\dots,v_t$$). Rows are only computed for constraints that enter the working set, and are cached.

The Riccati recursion uses the square-root form of HPIPM, and the classical form when the cost-to-go is not positive definite (decided again at each factorization, as in HPIPM). The products P_{t+1} b_t are kept until b_t or the factorization changes (as HPIPM's use_Pb). After a change of the vectors, the backward recursion is only redone from the last stage whose data changed.

**Normalization.** The rows of $$M$$ are normalized, as for dense QPs, which DAQP's tolerances and soft-constraint weights assume. The norms follow without forming $$M$$: with $$v \sim \mathcal{N}(0, I)$$, $$\Delta x_t$$ has a covariance $$\Sigma_t$$ that follows from a forward recursion, and $$\|M_a\|^2$$ is the norm of the row's block at its own stage plus $$g^T \Sigma_t g$$ for its coefficients $$g$$ on $$\Delta x_t$$. The norms are kept across refactorizations as long as they stay within a factor 2 (as between SQP iterations).

**Equality elimination.** At stage $$t$$, the eliminated equalities of the stage and the state equalities passed on from stage $$t+1$$ give $$G_u u_t + G_x x_t = g$$. A QR factorization with column pivoting of $$G_u^T$$ splits $$u_t$$ into the part the equalities determine and free inputs; the remaining combinations are state equalities passed on to stage $$t-1$$. The Riccati recursion and $$v$$ are formed in the free inputs only. If the elimination is ill-conditioned (inputs that only weakly affect the equalities), the equalities are kept in the working set instead, and the elimination is tried again at the next factorization.

**Free initial state.** If $$x_0$$ is a variable, an extra stage is added in front, with no state, input $$x_0$$, $$B = I$$ and zero cost.

**Semidefinite problems.** For proximal-point iterations, $$\epsilon I$$ is added to every stage Hessian in the original variables (the classical recursion is used and equalities are kept in the working set), and the inner problems have the linear term $$[r_t; q_t] - \epsilon z_{c,t}$$ for the current center $$z_c$$.

**Mixed-integer problems.** Branch and bound works on the LDP in $$v$$ and only fixes constraints in the working set, so it runs unchanged on the implicit $$M$$.

The approach (rows of $$M$$ and products $$Mv$$ from Riccati recursions in a dual active-set method, and the elimination of equalities) follows DAOCP ([A. Zaupa, S. Erickson, M. Johansson, "DAOCP: a dual active set solver for optimal control problems", arXiv:2609.33451, 2026](https://arxiv.org/abs/2609.33451)). The coupling constraints follow the multistage formulation of [PIQP](https://github.com/PREDICT-EPFL/piqp).
