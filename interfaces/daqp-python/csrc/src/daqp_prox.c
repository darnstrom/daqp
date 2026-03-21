#include "daqp_prox.h"
#include "utils.h"

static int gradient_step(DAQPWorkspace* work);

/* --------------------------------------------------------------------------
 * daqp_prox  --  outer proximal-point / semi-proximal loop
 *
 * QP problems (Rinv or RinvD is set):
 *   Implements a Semi-Proximal Method.  Only the directions i where the
 *   Cholesky diagonal would have been non-positive without regularisation
 *   (prox_mask[i] == 1) receive the eps·x_old[i] perturbation in the
 *   right-hand side.  If H is already positive definite everywhere
 *   (n_prox == 0) the inner QP equals the original and we exit after one
 *   solve.
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
    c_float eps     = work->settings->eps_prox;
    c_float eta     = work->settings->eta_prox;
    int  cycle_counter = 0;
    c_float best_fval  = DAQP_INF;

    const int is_lp = (work->Rinv == NULL && work->RinvD == NULL);

    // For a QP whose Hessian is already positive definite (n_prox == 0),
    // no direction needs a proximal shift.  The inner QP equals the
    // original problem, so one solve gives the exact solution.
    const int all_pd = (!is_lp) && (work->n_prox == 0);

    while(total_iter < work->settings->iter_limit){

        /* ----------------------------------------------------------------
         * Perturb the problem: form v = R'\(f - eps_mask * x_old)
         * ----------------------------------------------------------------*/
        if(is_lp){
            // No Hessian factor.  Adapt eps heuristically: grow when the
            // inner LP stalls (iterations==1), shrink otherwise to improve
            // accuracy.
            eps *= (work->iterations == 1) ? 10.0 : 0.9;
            if(eps > 1e3) eps = 1e3;
            for(i = 0; i < nx; i++)
                work->v[i] = work->qp->f[i]*eps - work->x[i];
        }
        else{
            // Semi-proximal: shift f only for directions that needed
            // regularisation to make the Cholesky factor non-singular.
            if(work->prox_mask != NULL){
                for(i = 0; i < nx; i++)
                    work->v[i] = work->qp->f[i]
                                 - (work->prox_mask[i] ? eps : 0.0) * work->x[i];
            }
            else{
                // prox_mask unavailable -- fall back to full proximal
                for(i = 0; i < nx; i++)
                    work->v[i] = work->qp->f[i] - eps * work->x[i];
            }
            daqp_update_v(work->v, work, 0);
        }

        // Perturb RHS of constraints to match the shifted objective
        daqp_update_d(work, work->qp->bupper, work->qp->blower);

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
         * Convergence check: fixed point  ||x - x_old||_inf < tol_stat
         * Only test when the active set did not change (iterations == 1),
         * because that is the cheapest indicator that a stationary point
         * has been reached for the inner QP.
         * ----------------------------------------------------------------*/
        if(work->iterations == 1){
            tol_stat = is_lp ? eta*eps : eta/eps;
            for(i = 0; i < nx; i++){
                max_diff = work->x[i] - work->xold[i];
                if(max_diff > tol_stat || max_diff < -tol_stat) break;
            }
            if(i == nx){
                exitflag = DAQP_EXIT_OPTIMAL; // Fixed point reached
                break;
            }
            // LP: when not at a vertex take a gradient step toward the
            // nearest constraint to escape from the interior.
            if(is_lp && work->n_active != nx){
                if(gradient_step(work) == DAQP_EMPTY_IND){
                    exitflag = DAQP_EXIT_UNBOUNDED;
                    break;
                }
            }
        }

        /* ----------------------------------------------------------------
         * Stagnation / cycle detection
         * ----------------------------------------------------------------*/
        if(is_lp){
            // Track LP objective f'x; no improvement over several
            // consecutive iterations signals a fixed point.
            c_float lp_obj = 0.0;
            for(i = 0; i < nx; i++) lp_obj += work->qp->f[i] * work->x[i];
            max_diff = best_fval - lp_obj;
            if(max_diff < work->settings->progress_tol){
                if(++cycle_counter > work->settings->cycle_tol){
                    exitflag = DAQP_EXIT_OPTIMAL; // Stagnated at fixed point
                    break;
                }
            }
            else{
                best_fval    = lp_obj;
                cycle_counter = 0;
            }
        }
        else{
            // QP: track the inner LDP dual objective ||u||^2 (work->fval).
            // When the outer iterate converges, v stabilises and so does
            // work->fval; stagnation therefore indicates a fixed point.
            max_diff = best_fval - work->fval;
            if(max_diff < work->settings->progress_tol){
                if(++cycle_counter > work->settings->cycle_tol){
                    exitflag = DAQP_EXIT_OPTIMAL; // Stagnated at fixed point
                    break;
                }
            }
            else{
                best_fval    = work->fval;
                cycle_counter = 0;
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
        for(k = 0, delta_s = 0, Ax = 0; k < nx; k++){
            Ax      += work->M[disp]   * work->x[k];
            delta_s -= work->M[disp++] * work->xold[k];
        }
        delta_s += Ax;

        if(work->scaling != NULL){
            Ax      /= work->scaling[j];
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
