#include "daqp_prox.h"
#include "utils.h"
#include <math.h>

static int gradient_step(DAQPWorkspace* work);

/*
 * Return the proximal shift used consistently by factorization, the outer
 * iteration, and objective-value recovery.  A relative floor prevents the
 * regularized LDP from inheriting a condition number that is too large for
 * the requested primal accuracy.
 */
c_float daqp_proximal_regularization_scaled(
        const DAQPWorkspace *work, c_float hessian_scale){
    c_float eps;
    if(work == NULL || work->settings == NULL) return 0.0;
    eps = work->settings->eps_prox;
    if(eps <= 0 || work->qp == NULL || work->qp->H == NULL ||
            work->qp->problem_type == 2)
        return eps;

    if(hessian_scale > 0){
        c_float floor = sqrt(work->settings->primal_tol) * hessian_scale;
        if(eps < floor) eps = floor;
    }
    return eps;
}

c_float daqp_proximal_regularization(const DAQPWorkspace *work){
    int i;
    c_float hessian_scale = 0.0;

    if(work == NULL || work->settings == NULL) return 0.0;
    if(work->settings->eps_prox <= 0 || work->qp == NULL ||
            work->qp->H == NULL || work->qp->problem_type == 2)
        return work->settings->eps_prox;

    for(i = 0; i < work->n; i++){
        c_float abs_diag = work->qp->H[i*work->n+i];
        if(abs_diag < 0) abs_diag = -abs_diag;
        if(abs_diag > hessian_scale) hessian_scale = abs_diag;
    }
    return daqp_proximal_regularization_scaled(work, hessian_scale);
}

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
    const int nx = work->n;
    int exitflag;
    c_float *swp_ptr;
    c_float max_diff, tol_stat;
    c_float eta = work->settings->eta_prox;
    int rhs_is_current = 1;

    const int is_lp = (work->Rinv == NULL && work->RinvD == NULL);
    c_float eps = is_lp ? 1.0 : daqp_proximal_regularization(work);

    // For a QP whose Hessian is already positive definite (n_prox == 0),
    // no direction needs a proximal shift.  The inner QP equals the
    // original problem, so one solve gives the exact solution.
    const int all_pd = (!is_lp) && (work->n_prox == 0);

    /*
     * Setup has already formed v and d for x=0. Preserve that work for the
     * common cold-start case instead of rebuilding the same transformed
     * gradient and every constraint RHS before the first inner solve.
     */
    for(i = 0; i < nx; i++){
        if(work->x[i] != 0.0){
            rhs_is_current = 0;
            break;
        }
    }

    while(total_iter < work->settings->iter_limit){

        /* ----------------------------------------------------------------
         * Perturb the problem: form v = R'\(f - eps_mask * x_old)
         * ----------------------------------------------------------------*/
        if(!rhs_is_current && is_lp){
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
        else if(!rhs_is_current){
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
            daqp_update_v(work->v, work, 0);
        }

        if(!rhs_is_current)
            daqp_update_d(work, work->qp->bupper, work->qp->blower);
        rhs_is_current = 0;

        // xold <-- x  (pointer swap avoids copying)
        swp_ptr = work->xold; work->xold = work->x; work->x = swp_ptr;

        /* ----------------------------------------------------------------
         * Solve the (regularised) least-distance problem
         * ----------------------------------------------------------------*/
        work->u = work->x;
        exitflag = daqp_ldp(work);

        total_iter += work->iterations;
        if(exitflag < 0)
            break;              // Inner solver failed -- propagate error
        else
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
            exitflag = DAQP_EXIT_OPTIMAL;
            break;
        }

        if(work->iterations == 1){
            // LP: when not at a vertex take a gradient step toward the
            // nearest constraint to escape from the interior.
            if(is_lp && work->n_active != nx){
                if(gradient_step(work) == DAQP_EMPTY_IND){
                    exitflag = DAQP_EXIT_UNBOUNDED;
                    break;
                }
            }
        }
    }

    // Finalize
    if(total_iter >= work->settings->iter_limit) exitflag = DAQP_EXIT_ITERLIMIT;
    if(is_lp){
        for(i = 0; i < work->n_active; i++)
            work->lam_star[i] /= eps; // Rescale dual variables
    }
    work->iterations = total_iter;
    return exitflag;
}

/* --------------------------------------------------------------------------
 * gradient_step  --  line-search toward the first blocking constraint
 *
 * Used for LP problems when the current iterate is not at a vertex (the
 * active set does not span all n variables).  Advances x along the
 * direction delta_x = x - x_old until the nearest constraint boundary is
 * hit, then returns the index of that constraint.
 * Returns DAQP_EMPTY_IND if the problem is unbounded in that direction.
 * --------------------------------------------------------------------------*/
static int gradient_step(DAQPWorkspace* work){
    int j, k, disp, add_ind = DAQP_EMPTY_IND;
    const int nx = work->n;
    const int m  = work->m;
    const int ms = work->ms;
    c_float Ax, delta_s, min_alpha = DAQP_INF;

    // Simple bounds: find first blocking constraint along delta_x
    for(j = 0; j < ms; j++){
        if(work->sense[j] & (DAQP_ACTIVE + DAQP_IMMUTABLE)) continue;
        delta_s = work->x[j] - work->xold[j];
        if(delta_s > 0 &&
                work->qp->bupper[j] < DAQP_INF &&
                work->qp->bupper[j] - work->x[j] < min_alpha*delta_s){
            add_ind   = j;
            min_alpha = (work->qp->bupper[j] - work->x[j]) / delta_s;
        }
        else if(delta_s < 0 &&
                work->qp->blower[j] > -DAQP_INF &&
                work->qp->blower[j] - work->x[j] > min_alpha*delta_s){
            add_ind   = j;
            min_alpha = (work->qp->blower[j] - work->x[j]) / delta_s;
        }
    }

    // General bounds
    for(j = ms, disp = 0; j < m; j++){
        if(work->sense[j] & (DAQP_ACTIVE + DAQP_IMMUTABLE)){
            disp += nx;
            continue;
        }
        if(work->Mu != NULL){
            /*
             * At a completed inner solve Mu=M*u and x=u-v, while
             * d=b_scaled+M*v. Hence the normalized upper slack is
             * dupper-Mu (similarly for the lower slack). Only M*(x-xold)
             * still needs a dot product; normalization cancels in the
             * step-length ratio.
             */
            for(k = 0, delta_s = 0; k < nx; k++)
                delta_s += work->M[disp++] * (work->x[k]-work->xold[k]);
            Ax = work->Mu[j-ms];
            if(delta_s > 0 &&
                    work->qp->bupper[j] < DAQP_INF &&
                    work->dupper[j] - Ax < delta_s*min_alpha){
                add_ind   = j;
                min_alpha = (work->dupper[j] - Ax) / delta_s;
            }
            else if(delta_s < 0 &&
                    work->qp->blower[j] > -DAQP_INF &&
                    work->dlower[j] - Ax > delta_s*min_alpha){
                add_ind   = j;
                min_alpha = (work->dlower[j] - Ax) / delta_s;
            }
            continue;
        }

        // Fallback for static workspaces without the batched Mu cache.
        for(k = 0, delta_s = 0, Ax = 0; k < nx; k++){
            Ax      += work->M[disp]   * work->x[k];
            delta_s -= work->M[disp++] * work->xold[k];
        }
        delta_s += Ax;
        if(work->scaling != NULL){
            Ax /= work->scaling[j];
            delta_s /= work->scaling[j];
        }
        if(delta_s > 0 &&
                work->qp->bupper[j] < DAQP_INF &&
                work->qp->bupper[j] - Ax < delta_s*min_alpha){
            add_ind   = j;
            min_alpha = (work->qp->bupper[j] - Ax) / delta_s;
        }
        else if(delta_s < 0 &&
                work->qp->blower[j] > -DAQP_INF &&
                work->qp->blower[j] - Ax > delta_s*min_alpha){
            add_ind   = j;
            min_alpha = (work->qp->blower[j] - Ax) / delta_s;
        }
    }

    // Advance: x <-- x + min_alpha * (x - x_old)
    if(add_ind != DAQP_EMPTY_IND)
        for(k = 0; k < nx; k++)
            work->x[k] += min_alpha * (work->x[k] - work->xold[k]);

    return add_ind;
}
