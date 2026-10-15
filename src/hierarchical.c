#include "hierarchical.h"
#include "auxiliary.h"
#include <string.h>
#include <math.h>

/*
 * Hierarchical QP: the violations of the constraints in each level are
 * minimized (in the least-squares sense) without increasing those of the
 * levels above it. The first level is hard.
 *
 * Each level is solved as an LDP with its constraints softened, and re-solved
 * centered at its solution until the slacks settle (proximal iterations), such
 * that the slacks and multipliers are those of the level without the soft
 * penalty. The bounds of the active constraints are then shifted by their
 * slacks. A constraint that stays active in the lower levels (an equality, a
 * violated constraint, or one with a nonzero multiplier) is eliminated:
 * u = Q*[y; z], with Q from Householder reflections and y fixed, and the lower
 * levels are solved as an LDP in z (since ||u||^2 = ||y||^2 + ||z||^2). The rows
 * of the constraints are copied at the first elimination (to work->hierarchy->M)
 * and reduced there, so the LDP of the workspace is unchanged.
 */

// An equality that is this close (relatively) to the span of the ones
// eliminated before it is implied by them
#define DAQP_HQ_DEPENDENT ((c_float)1e-10)
// Maximum number of proximal iterations in a level
#define DAQP_HQ_MAX_PROX 20
// The proximal iterations have converged when a step is this small
// (relatively, squared)
#define DAQP_HQ_PROX_STEP ((c_float)1e-12)
// Slacks and multipliers (times rho_soft) above this factor times primal_tol
// are taken as nonzero when freezing
#define DAQP_HQ_FREEZE_TOL 10

// x <-- (I - tau*v*v')*x
static void daqp_hq_reflect(const c_float *v, const c_float tau, c_float *x, const int len){
    int k;
    c_float dot = 0;
    for(k = 0; k < len; k++) dot += v[k]*x[k];
    dot *= tau;
    for(k = 0; k < len; k++) x[k] -= dot*v[k];
}

// Copy the rows of all constraints to h->M (row id = constraint id), which
// is done before the first elimination. A simple bound is its row of the
// (normalized) Rinv, or a unit vector if Rinv is diagonal
static void daqp_hq_copy_rows(DAQPWorkspace *work, DAQPHierarchy *h, const int m){
    int j;
    const int n = work->n, ms = work->ms;
    for(j = 0; j < m; j++){
        c_float *row = h->M+(size_t)j*n;
        if(j >= ms) memcpy(row,work->M+(size_t)(j-ms)*n,n*sizeof(c_float));
        else{
            memset(row,0,n*sizeof(c_float));
            if(work->Rinv == NULL) row[j] = 1;
            else memcpy(row+j,work->Rinv+DAQP_R_OFFSET(j,n)+j,(n-j)*sizeof(c_float));
        }
    }
    work->M = h->M; work->Mu = h->Mu; work->ms = 0;
}

// The value of constraint j at v
static c_float daqp_hq_row_value(DAQPWorkspace *work, const int j, const c_float *v){
    int k;
    const int n = work->n;
    c_float val = 0;
    if(j >= work->ms){
        const c_float *row = work->M+(size_t)(j-work->ms)*n;
        for(k = 0; k < n; k++) val += row[k]*v[k];
    }
    else if(work->Rinv == NULL) val = v[j];
    else{
        const c_float *row = work->Rinv+DAQP_R_OFFSET(j,n);
        for(k = j; k < n; k++) val += row[k]*v[k];
    }
    return val;
}

// d <-- d - sign*M*v for the first m constraints (the LDP centered at v)
static void daqp_hq_shift(DAQPWorkspace *work, const int m, const c_float *v, const c_float sign){
    int j;
    for(j = 0; j < m; j++){
        const c_float val = sign*daqp_hq_row_value(work,j,v);
        work->dupper[j] -= val;
        work->dlower[j] -= val;
    }
    work->reuse_ind = 0; // (the right-hand side of the working set changed)
}

// Reflections of at least this length are applied to four rows at a time
#define DAQP_HQ_BLOCK_LEN 24

// daqp_hq_reflect for four vectors
static void daqp_hq_reflect4(const c_float *v, const c_float tau, c_float *x0,
        c_float *x1, c_float *x2, c_float *x3, const int len){
    int k;
    c_float d0 = 0, d1 = 0, d2 = 0, d3 = 0;
    for(k = 0; k < len; k++){
        d0 += v[k]*x0[k]; d1 += v[k]*x1[k]; d2 += v[k]*x2[k]; d3 += v[k]*x3[k];
    }
    d0 *= tau; d1 *= tau; d2 *= tau; d3 *= tau;
    for(k = 0; k < len; k++){
        x0[k] -= d0*v[k]; x1[k] -= d1*v[k]; x2[k] -= d2*v[k]; x3[k] -= d3*v[k];
    }
}

// The reflection of coordinate k (acting on u[k:n], with n the original size)
#define DAQP_HQ_V(h,k,n) ((h)->V+(k)*(n)-(k)*((k)-1)/2)

// Eliminate the immutable constraints in the working set (see above). The
// coordinates that they fix are kept in front of work->u, which is moved past
// them. nfixed is the number of coordinates fixed so far (updated)
static void daqp_hq_eliminate(DAQPWorkspace *work, DAQPHierarchy *h, const int m, int *nfixed){
    int i, k, r, id, p = 0;
    const int nc = work->n, n0 = nc+*nfixed;
    c_float *M = h->M;
    for(i = 0; i < work->n_active; i++){
        id = work->WS[i];
        if(!DAQP_IS_IMMUTABLE(id)) continue;
        if(work->M != M) daqp_hq_copy_rows(work,h,m); // (first elimination)
        // Ignored from now on (IMMUTABLE and inactive), until it is restored
        DAQP_SET_INACTIVE(id);
        work->sense[id] |= DAQP_HQ_ELIMINATED;
        c_float *a = M+(size_t)id*nc;
        c_float s = 0, an = 0;
        for(k = 0; k < nc; k++) an += a[k]*a[k];
        for(k = p; k < nc; k++) s += a[k]*a[k];
        if(s <= DAQP_HQ_DEPENDENT*DAQP_HQ_DEPENDENT*an) continue; // Implied
        // Reflection that zeros a[p+1:nc], for coordinate g = nfixed+p
        const int g = *nfixed+p, len = nc-p;
        s = sqrt(s);
        c_float *v = DAQP_HQ_V(h,g,n0);
        v[0] = a[p] > 0 ? a[p]+s : a[p]-s;
        for(k = p+1; k < nc; k++) v[k-p] = a[k];
        const c_float tau = h->tau[g] = 1/(s*(s+fabs(a[p]))); // 2/||v||^2
        // Apply it to u and the rows of the other constraints
        daqp_hq_reflect(v,tau,work->u+p,len);
        int rows[4], nr = 0;
        for(r = 0; r < m; r++){
            if(work->sense[r] & DAQP_HQ_ELIMINATED) continue;
            if(len < DAQP_HQ_BLOCK_LEN){
                daqp_hq_reflect(v,tau,M+(size_t)r*nc+p,len);
                continue;
            }
            rows[nr++] = r;
            if(nr == 4){
                daqp_hq_reflect4(v,tau,M+(size_t)rows[0]*nc+p,M+(size_t)rows[1]*nc+p,
                        M+(size_t)rows[2]*nc+p,M+(size_t)rows[3]*nc+p,len);
                nr = 0;
            }
        }
        for(r = 0; r < nr; r++) daqp_hq_reflect(v,tau,M+(size_t)rows[r]*nc+p,len);
        p++;
    }
    if(p == 0) return;
    // Fix u[0:p] and shift the bounds by the fixed part of each row
    const int nz = nc-p;
    for(r = 0; r < m; r++){
        c_float *row = M+(size_t)r*nc, shift = 0;
        for(k = 0; k < p; k++) shift += row[k]*work->u[k];
        work->dupper[r] -= shift;
        work->dlower[r] -= shift;
        memmove(M+(size_t)r*nz,row+p,nz*sizeof(c_float));
    }
    work->u += p;
    work->n = nz;
    *nfixed += p;
}

// Re-solve the LDP of a level centered at its solution, while there are
// violated soft constraints and the steps are not small, accumulating the center
// in c. Returns 1 if this has converged (the slacks and multipliers are then
// those of the level without the soft penalty)
static int daqp_hq_proximal(DAQPWorkspace *work, const int m, c_float *c, int *iterations){
    int it, j, converged = 0;
    c_float cn = 0;
    for(j = 0; j < work->n; j++) c[j] = 0;
    for(it = 0; it < DAQP_HQ_MAX_PROX; it++){
        for(j = 0; j < work->n_active; j++)
            if(DAQP_IS_SOFT(work->WS[j]) &&
                    fabs(daqp_soft_slack(work,j)) > work->settings->primal_tol) break;
        if(j == work->n_active) break; // (no violations)
        // Move the center to the current point
        c_float un = 0;
        daqp_hq_shift(work,m,work->u,1);
        for(j = 0, cn = 0; j < work->n; j++){
            un += work->u[j]*work->u[j];
            c[j] += work->u[j];
            cn += c[j]*c[j];
            work->u[j] = 0;
        }
        if(it > 0 && un <= DAQP_HQ_PROX_STEP*cn){ converged = 1; break; }
        const int exitflag = daqp_ldp(work);
        *iterations += work->iterations;
        if(exitflag < 0){ // Keep the current center
            for(j = 0; j < work->n; j++) work->u[j] = 0;
            break;
        }
    }
    if(cn > 0){ // Back to the original coordinates
        daqp_hq_shift(work,m,c,-1);
        for(j = 0; j < work->n; j++) work->u[j] += c[j];
    }
    return converged;
}

int daqp_hiqp(DAQPWorkspace *work, c_float *lambda){
    int i, j, k, id, start, end, iterations = 0, exitflag = 0;
    const c_float freeze_tol = DAQP_HQ_FREEZE_TOL*work->settings->primal_tol;

    // If only one level -> just solve normal LDP
    if(!DAQP_IS_HIERARCHICAL(work)) return daqp_ldp(work);

    // A previous solve shifted d by the slacks of the soft levels, so it has
    // to be reformed by an update first (without a qp, the caller resets d)
    if((work->state & DAQP_UPDATE_d) && work->qp != NULL) return DAQP_EXIT_UNSUPPORTED;

    DAQPHierarchy *h = work->hierarchy;
    const int nh = h->nh, *break_points = h->break_points;
    const int n0 = work->n, ms0 = work->ms, m = break_points[nh-1];

    // Reset lambda for output
    if(lambda != NULL) for(i=0;i<m;i++) lambda[i]=0;

    // Start from the (hard) first level
    start=break_points[0];
    work->m = start;
    reset_daqp_workspace(work);
    exitflag = daqp_activate_constraints(work);
    if(exitflag < 0) return exitflag;

    c_float *M0 = work->M, *Mu0 = work->Mu, *u0 = work->u;
    int nfixed = 0;

    for(i =1; i < nh; i++){
        // initialize current level
        end=break_points[i];
        work->m = end;
        // Soften constraints and activate
        for(j =start;j<end;j++){
            DAQP_SET_SOFT(j);
            if(DAQP_IS_ACTIVE(j)){
                daqp_add_constraint(work,j,DAQP_IS_LOWER(j) ? -1.0 : 1.0);
                // A dependent constraint (e.g., from a warm start) is left out
                // (due to pivoting, it is not necessarily j)
                if(work->sing_ind != DAQP_EMPTY_IND) daqp_drop_singular_last(work);
            }
        }

        // Save the solution of the previous level in case daqp_ldp fails
        for(j = 0; j<work->n;j++) work->xold[j] = work->u[j];
        exitflag = daqp_ldp(work);
        iterations+=work->iterations;
        if(exitflag < 0) break;

        // (xold is not needed anymore, and holds the center)
        const int converged = daqp_hq_proximal(work,end,work->xold,&iterations);
        if(iterations >= work->settings->iter_limit){
            exitflag = DAQP_EXIT_CONSTRAINED_POINT;
            break;
        }

        // Shift the bounds of the active soft constraints by their slacks, and
        // freeze the equalities and, after converged proximal iterations, the
        // violated constraints and the ones with nonzero multipliers
        for(j=0; j<work->n_active;j++){
            id=work->WS[j];
            const int soft = DAQP_IS_SOFT(id);
            const c_float w = soft ? daqp_soft_slack(work,j) : 0;
            if(work->dupper[id]-work->dlower[id] < work->settings->zero_tol || (converged &&
                        (soft ? fabs(w) : work->settings->rho_soft*fabs(work->lam_star[j])) > freeze_tol)){
                const c_float b = (DAQP_IS_LOWER(id) ? work->dlower[id] : work->dupper[id]) + w;
                work->dupper[id] = b;
                work->dlower[id] = b;
                if(!DAQP_IS_IMMUTABLE(id)) work->sense[id] |= DAQP_IMMUTABLE|DAQP_HQ_FROZEN;
            }
            else if(w < -work->settings->primal_tol)
                work->dlower[id]+=w;
            else if(w > work->settings->primal_tol)
                work->dupper[id]+=w;
            if(soft && lambda != NULL) // (offset for weakly active)
                lambda[id] = w + (DAQP_IS_LOWER(id) ? -1e-14 : 1e-14);
        }

        // Make constraints in current level hard
        for(j=start; j<end;j++) DAQP_SET_HARD(j);

        if(i==nh-1) break;

        // Eliminate the frozen constraints, and refactorize the rest of the working set
        daqp_hq_eliminate(work,h,m,&nfixed);
        if(work->n == 0){ // No degrees of freedom left => the lower levels are fixed
            if(lambda != NULL){
                for(j = end; j < m; j++){
                    if(DAQP_IS_IMMUTABLE(j) && !DAQP_IS_ACTIVE(j)) continue;
                    // (the rows have no entries left: A_j x is in the shifted bounds)
                    if(work->dupper[j] < 0) lambda[j] = -work->dupper[j];
                    else if(work->dlower[j] > 0) lambda[j] = -work->dlower[j];
                }
            }
            break;
        }
        for(j = 0, k = 0; j < work->n_active; j++){
            if(DAQP_IS_IMMUTABLE(work->WS[j])) continue;
            work->WS[k] = work->WS[j];
            work->lam_star[k++] = work->lam_star[j];
        }
        reset_daqp_workspace(work);
        for(j = 0; j < k; j++){
            daqp_add_constraint(work,work->WS[j],work->lam_star[j]);
            if(work->sing_ind != DAQP_EMPTY_IND) daqp_drop_singular_last(work);
        }

        // Move down the hierarchy
        start = end;
    }
    if(exitflag < 0){ // Return the last point that satisfies the levels above the failed one
        for(j = 0; j<work->n;j++) work->u[j] = work->xold[j];
        exitflag = DAQP_EXIT_CONSTRAINED_POINT;
    }

    // Map u back: u = Q*[y; z], with the reflections in reverse
    for(k = nfixed-1; k >= 0; k--)
        daqp_hq_reflect(DAQP_HQ_V(h,k,n0),h->tau[k],u0+k,n0-k);

    // Restore the LDP (the bounds are reformed by the next update), and keep
    // the active set of the constraints for warm starting the next solve
    work->M = M0; work->Mu = Mu0; work->ms = ms0; work->m = m;
    work->u = u0; work->n = n0;
    for(j = 0; j < m; j++){
        if(work->sense[j] & DAQP_HQ_ELIMINATED) DAQP_SET_ACTIVE(j);
        if(work->sense[j] & DAQP_HQ_FROZEN) work->sense[j] &= ~DAQP_IMMUTABLE;
        work->sense[j] &= ~(DAQP_HQ_ELIMINATED|DAQP_HQ_FROZEN);
    }
    reset_daqp_workspace(work); // (the working set is formed again by the next solve)
    work->iterations = iterations; // Append total number of iterations
    work->state |= DAQP_UPDATE_d; // The levels have shifted d
    return exitflag;
}
