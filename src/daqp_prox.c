#include "daqp_prox.h"
#include "auxiliary.h"
#include "utils.h"

static int prox_step(DAQPWorkspace* work, c_float* s_prev);

/* --------------------------------------------------------------------------
 * daqp_prox  --  outer proximal-point / semi-proximal loop
 *
 * QP problems (Rinv or RinvD is set):
 *   Diagonal Hessians use a semi-proximal method, perturbing only singular
 *   coordinate directions. Dense Hessians use a full proximal shift when
 *   Cholesky detects numerical singularity; selectively shifting failed
 *   pivots is not a reliable nullspace regularization. If H is already
 *   positive definite (n_prox == 0), the inner QP equals the original and
 *   we exit after one solve.
 *
 * LP problems (Rinv == NULL && RinvD == NULL):
 *   Classical regularisation-based smoothing with adaptive eps.
 * --------------------------------------------------------------------------*/
int daqp_prox(DAQPWorkspace *work){
    int i, total_iter = 0;
    c_float s_prev = -1; // Exact step length of the latest prox_step (-1: none)
    int center_relaxed = 0;
    int nx, is_lp;
    const c_float relaxation = 1.5;
    int exitflag = DAQP_EXIT_ITERLIMIT; // If no iteration can be taken
    c_float *swp_ptr;
    c_float max_diff, tol_stat;
    c_float eta = work->settings->eta_prox;
    c_float eps;

    // nh act as a counter for outer iterations
    work->nh = 0;

    nx = work->n;
    is_lp = (work->Rinv == NULL && work->RinvD == NULL);
    eps = is_lp ? 1.0 : daqp_get_proximal_regularization(work);

    // For a QP whose Hessian is already positive definite (n_prox == 0),
    // no direction needs a proximal shift.  The inner QP equals the
    // original problem, so one solve gives the exact solution.
    const int all_pd = (!is_lp) && (work->n_prox == 0);

    // A negative eta selects an automatic tolerance. Preserve the established
    // default, but tighten it when the user requests a non-default dual
    // tolerance. Skip this entirely for positive-definite QPs, where no
    // proximal convergence test is needed.
    if(!all_pd && eta < 0.0){
        eta = DAQP_AUTO_ETA_CAP;
        if(work->settings->dual_tol != DAQP_DEFAULT_DUAL_TOL &&
           0.1 * work->settings->dual_tol < eta)
            eta = 0.1 * work->settings->dual_tol;
    }

    while(total_iter < work->settings->iter_limit){
        /* ----------------------------------------------------------------
         * Perturb the problem: form v = R'\(f - eps_mask * x_old)
         * ----------------------------------------------------------------*/
        if(is_lp){
            // No Hessian factor.  Adapt eps heuristically: grow when the
            // inner LP stalls (iterations==1), shrink otherwise to improve
            // accuracy. Keep the deterministic initial value for the first
            // solve; work->iterations has no current-loop value yet.
            if(total_iter > 0)
                eps *= (work->iterations == 1) ? 10.0 : 0.9;
            if(eps > 1e3) eps = 1e3;
            for(i = 0; i < nx; i++)
                work->v[i] = work->qp->f[i]*eps - work->x[i];
        }
        else{
            if(work->prox_mask == NULL || work->n_prox == nx){
                // Dense singular Hessians use a full shift. Avoid a mask
                // load and branch for every component on every outer step.
                if(work->qp->f != NULL)
                    for(i = 0; i < nx; i++)
                        work->v[i] = work->qp->f[i] - eps * work->x[i];
                else
                    for(i = 0; i < nx; i++)
                        work->v[i] = -eps * work->x[i];
            }
            else{
                // Diagonal Hessians can regularize only singular directions.
                if(work->qp->f != NULL)
                    for(i = 0; i < nx; i++)
                        work->v[i] = work->qp->f[i]
                                     - (work->prox_mask[i] ? eps : 0.0) * work->x[i];
                else
                    for(i = 0; i < nx; i++)
                        work->v[i] = -(work->prox_mask[i] ? eps : 0.0)
                                     * work->x[i];
            }
            daqp_update_v(work->v, work);
        }

        daqp_update_d(work, work->qp->bupper, work->qp->blower);

        // xold <-- x  (pointer swap avoids copying)
        swp_ptr = work->xold; work->xold = work->x; work->x = swp_ptr;

        /* ----------------------------------------------------------------
         * Solve the (regularised) least-distance problem
         * ----------------------------------------------------------------*/
        work->u = work->x;
        work->nh++;
        exitflag = daqp_ldp(work);

        total_iter += work->iterations;
        if(exitflag < 0)
            break;              // Inner solver failed -- propagate error
        ldp2qp_solution(work); // Recover QP primal from LDP dual

        if(eps == 0) break;     // No regularisation -> single outer step

        /* ----------------------------------------------------------------
         * If H is fully positive definite, the inner QP is the original
         * problem.  The first solve gives the exact solution.
         * ----------------------------------------------------------------*/
        if(all_pd){
            exitflag = DAQP_EXIT_OPTIMAL;
            break;
        }

        /* ----------------------------------------------------------------
         * Convergence check: fixed point  ||x - x_old||_inf < tol_stat.
         *
         * A fixed point is a valid stationarity certificate regardless of
         * how many active-set changes the inner solve needed.  Checking it
         * after every successful solve avoids extra outer iterations and,
         * unlike objective stagnation, cannot label a non-stationary point
         * as optimal.
         * ----------------------------------------------------------------*/
        tol_stat = is_lp ? eta*eps : eta/eps;
        for(i = 0; i < nx; i++){
            max_diff = work->x[i] - work->xold[i];
            if(max_diff > tol_stat || max_diff < -tol_stat) break;
        }
        if(i == nx){
            if(center_relaxed &&
                    total_iter < work->settings->iter_limit){
                center_relaxed = 0;
                continue; // Confirm convergence from the feasible iterate.
            }
            exitflag = DAQP_EXIT_OPTIMAL;
            break;
        }

        // With an unchanged working set the proximal map is locally affine,
        // and its steps slow down along directions with little curvature (or
        // none, for an LP). Move the center along the step (prox_step), or
        // relax the step for an AVI. The center is then confirmed by a plain
        // proximal step before convergence is declared.
        center_relaxed = 0;
        if(work->iterations != 1) s_prev = -1; // The working set has changed
        if(work->iterations == 1 && work->n_active < nx &&
                total_iter < work->settings->iter_limit){
            if(work->avi != NULL){
                for(i = 0; i < nx; i++)
                    work->x[i] = work->xold[i]
                        + relaxation*(work->x[i] - work->xold[i]);
                center_relaxed = 1;
            }
            else{
                const int step_flag = prox_step(work,&s_prev);
                if(step_flag == DAQP_EXIT_UNBOUNDED){
                    exitflag = DAQP_EXIT_UNBOUNDED;
                    break;
                }
                center_relaxed = step_flag;
            }
        }
    }

    // Finalize
    if(total_iter >= work->settings->iter_limit) exitflag = DAQP_EXIT_ITERLIMIT;
    // The final iterate is the solution of the latest inner problem (not
    // refined after a short solve, see DAQP_REFINE_MIN_ITER)
    if(exitflag > 0 && total_iter > DAQP_REFINE_MIN_ITER) daqp_refine_primal(work);
    if(is_lp){
        for(i = 0; i < work->n_active; i++)
            work->lam_star[i] /= eps; // Rescale dual variables
    }
    else{
        /*
         * daqp_extract_result forms 0.5*(fval - ||v||^2). Correct the
         * regularized objective here while the reconstructed eps is local,
         * avoiding a second reconstruction during result extraction.
         */
        c_float prox_norm = 0.0;
        for(i = 0; i < nx; i++){
            if(work->prox_mask == NULL || work->prox_mask[i])
                prox_norm += work->x[i]*work->x[i];
        }
        work->fval += eps*prox_norm;
    }
    work->iterations = total_iter;
    return exitflag;
}

/*
 * Step length along d = x - x_old that minimizes the objective 0.5x'Hx + f'x:
 * -g'd/d'Hd with g = Hx+f (d in work->xldl on return). Returns DAQP_INF for a
 * direction without curvature (an LP, or curvature at the level of rounding
 * errors), and -1 if d is not a descent direction.
 */
static c_float curvature_step(DAQPWorkspace* work){
    int i, j;
    const int n = work->n;
    const DAQPProblem* qp = work->qp;
    // Scratch: d and Hd (formed anew by the next inner solve, since
    // daqp_update_d resets reuse_ind before it)
    c_float *d = work->xldl, *hd = work->zldl;
    c_float gd = 0, dhd = 0, dd = 0, hmax = 0;

    for(i = 0; i < n; i++) d[i] = work->x[i] - work->xold[i];
    if(qp->f != NULL) for(i = 0; i < n; i++) gd += qp->f[i]*d[i];
    if(qp->H != NULL){
        for(i = 0; i < n; i++){
            const c_float* Hi = qp->H+(size_t)i*n;
            const c_float hii = Hi[i] < 0 ? -Hi[i] : Hi[i];
            c_float sum = 0;
            for(j = 0; j < n; j++) sum += Hi[j]*d[j];
            hd[i] = sum;
            if(hii > hmax) hmax = hii;
        }
        for(i = 0; i < n; i++){
            gd += work->x[i]*hd[i];
            dhd += d[i]*hd[i];
            dd += d[i]*d[i];
        }
    }
    if(gd >= 0) return -1;
    return (dhd > work->settings->zero_tol*hmax*dd) ? -gd/dhd : DAQP_INF;
}

/*
 * The first constraint that blocks the step x + s*(x - x_old) for s < *s,
 * among the ones that are not active, immutable, or set aside. Returns its
 * index (DAQP_EMPTY_IND if none), with *s shortened to the step that reaches
 * it and *lower marking whether it is its lower bound.
 */
static int blocking_constraint(DAQPWorkspace* work, c_float* s, int* lower){
    int i, j, ind = DAQP_EMPTY_IND;
    const int n = work->n, m = work->m, ms = work->ms;
    const DAQPProblem* qp = work->qp;
    c_float ad, ax, sb;
    for(i = 0; i < m; i++){
        if(work->sense[i] & (DAQP_ACTIVE + DAQP_IMMUTABLE + DAQP_SET_ASIDE)) continue;
        if(i < ms){ ax = work->x[i]; ad = ax - work->xold[i]; }
        else{
            const c_float* a = qp->A+(size_t)(i-ms)*n;
            for(j = 0, ad = 0, ax = 0; j < n; j++){
                ax += a[j]*work->x[j];
                ad += a[j]*(work->x[j] - work->xold[j]);
            }
        }
        if(ad > 0 && qp->bupper[i] < DAQP_INF) sb = (qp->bupper[i]-ax)/ad;
        else if(ad < 0 && qp->blower[i] > -DAQP_INF) sb = (qp->blower[i]-ax)/ad;
        else continue;
        if(sb < *s){
            *s = sb;
            *lower = ad < 0;
            ind = i;
        }
    }
    return ind;
}

/* --------------------------------------------------------------------------
 * prox_step  --  step along the latest proximal step d = x - x_old
 *
 * With an unchanged working set, d keeps the active constraints active, and
 * the proximal iterates contract along d only by eps/(eps+mu), where mu =
 * d'Hd/d'd is the curvature along d (for an LP, they do not converge along d
 * at all). x is therefore moved along d (curvature_step), or to the first
 * inactive constraint that blocks the step, which is then added to the working
 * set. A blocking constraint that is linearly dependent on the active ones
 * does not change the working set (and daqp_ldp would remove it again), so it
 * is set aside and the step continues to the next one.
 *
 * Returns 1 if the center was moved and 0 if not (d is not a descent
 * direction, or no constraint blocks a step without curvature in a QP, which
 * may be flat only up to rounding errors). An LP whose descent direction is
 * not blocked is unbounded: DAQP_EXIT_UNBOUNDED.
 * --------------------------------------------------------------------------*/
static int prox_step(DAQPWorkspace* work, c_float* s_prev){
    int i, k, ind, lower = 0, moved = 0, skipped = 0, first = 1;
    c_float s;
    while((s = curvature_step(work)) >= 0){
        const c_float* d = work->xldl;
        if(first){ // Lagged step length (see above)
            const c_float s_exact = s;
            if(s_exact < DAQP_INF && *s_prev >= 0)
                s = (*s_prev < 2*s_exact) ? *s_prev : 2*s_exact;
            *s_prev = (s_exact < DAQP_INF) ? s_exact : -1;
            first = 0;
        }
        ind = blocking_constraint(work,&s,&lower);
        if(ind == DAQP_EMPTY_IND){
            if(s < DAQP_INF){ // The minimizer along d
                for(k = 0; k < work->n; k++) work->x[k] += s*d[k];
                moved = 1;
            }
            else if(work->qp->H == NULL && !moved && !skipped)
                return DAQP_EXIT_UNBOUNDED;
            break;
        }
        *s_prev = -1; // Blocked: the working set changes
        // Advance to the blocking constraint and activate it. A constraint
        // that x already violates (within the tolerance) blocks at once: the
        // step is not reversed, which would undo proximal progress.
        if(s >= 0) for(k = 0; k < work->n; k++) work->x[k] += s*d[k];
        moved = 1;
        if(lower) DAQP_SET_LOWER(ind);
        else DAQP_SET_UPPER(ind);
        daqp_add_constraint(work, ind, lower ? -1.0 : 1.0);
        if(work->sing_ind == DAQP_EMPTY_IND) break;
        // Linearly dependent on the active constraints: set it aside
        work->n_active--;
        DAQP_SET_INACTIVE(ind);
        work->sense[ind] |= DAQP_SET_ASIDE;
        work->sing_ind = DAQP_EMPTY_IND;
        if(work->reuse_ind > work->n_active) work->reuse_ind = work->n_active;
        skipped = 1;
    }
    if(skipped)
        for(i = 0; i < work->m; i++) work->sense[i] &= ~DAQP_SET_ASIDE;
    return moved;
}
