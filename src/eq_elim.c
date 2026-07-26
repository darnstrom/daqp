#include "api.h"
#include "utils.h"
#include "eq_elim.h"
#include "constants.h"
#include "factorization.h"
#include <stdlib.h>
#include <math.h>

/*
 * The equality constraints are eliminated from
 *   min 0.5||u||^2  s.t.  dlower <= M u <= dupper,   M = A*Rinv, u = R*x+v
 * through the QR factorization M_E' = [Q1 Q2][R; 0] and u = Q1*y1 + Q2*w,
 * where R'y1 = d_E. Since ||u||^2 = ||y1||^2+||w||^2, the reduced problem in w
 * is again a least-distance problem, so nothing has to be refactorized.
 *
 * The reduced constraints are M*Q2 = A*(Rinv*Q2), which is why W = Rinv*Q2 is
 * formed: the rows of the reduced problem are then A_i*W, the rows belonging
 * to simple bounds are rows of W, and the solution is x = xp + W*w. Only the
 * equality rows of M are needed, so the full M never has to be formed.
 */

// A constraint is eliminated if it is a general equality constraint
// (equalities are marked as active and immutable in sense)
#define IS_EQ_CAND(work,i) \
    (((work)->sense[i]&(DAQP_ACTIVE+DAQP_IMMUTABLE+DAQP_SOFT+DAQP_BINARY)) \
     == (DAQP_ACTIVE+DAQP_IMMUTABLE))

/* ---------------------------------------------------------------------------
 * Products with the Hessian factor
 *
 * The first ms rows of a dense Rinv are normalized in place, which is undone
 * here so that everything below refers to the same Rinv as the QP does.
 * -------------------------------------------------------------------------*/

// dst <-- (row i of A)*Rinv, for a general constraint
static void arinv_row(const DAQPWorkspace* work, const int i, c_float* dst){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n, ms = eq->ms;
    const c_float* a = work->qp->A+(size_t)(i-ms)*n;
    int j, k;
    for(j = 0; j < n; j++) dst[j] = 0;
    if(eq->Rinv_full != NULL){
        for(k = 0; k < n; k++){
            const c_float* Rk = eq->Rinv_full+DAQP_R_OFFSET(k,n);
            const c_float ak = (k < ms) ? a[k]/eq->scaling_full[k] : a[k];
            if(ak == 0) continue;
            for(j = k; j < n; j++) dst[j] += ak*Rk[j];
        }
    }
    else if(eq->RinvD_full != NULL)
        for(k = 0; k < n; k++) dst[k] = a[k]*eq->RinvD_full[k];
    else
        for(k = 0; k < n; k++) dst[k] = a[k];
}

// dst <-- Rinv*y
static void rinv_times(const DAQPWorkspace* work, const c_float* y, c_float* dst){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n, ms = eq->ms;
    int i, j;
    if(eq->Rinv_full != NULL){
        for(i = 0; i < n; i++){
            const c_float* Ri = eq->Rinv_full+DAQP_R_OFFSET(i,n);
            c_float sum = 0;
            for(j = i; j < n; j++) sum += Ri[j]*y[j];
            dst[i] = (i < ms) ? sum/eq->scaling_full[i] : sum;
        }
    }
    else if(eq->RinvD_full != NULL)
        for(i = 0; i < n; i++) dst[i] = eq->RinvD_full[i]*y[i];
    else
        for(i = 0; i < n; i++) dst[i] = y[i];
}

// dst <-- Rinv'*y
static void rinvT_times(const DAQPWorkspace* work, const c_float* y, c_float* dst){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n, ms = eq->ms;
    int i, j;
    if(eq->Rinv_full != NULL){
        for(j = 0; j < n; j++) dst[j] = 0;
        for(i = 0; i < n; i++){
            const c_float* Ri = eq->Rinv_full+DAQP_R_OFFSET(i,n);
            const c_float yi = (i < ms) ? y[i]/eq->scaling_full[i] : y[i];
            if(yi == 0) continue;
            for(j = i; j < n; j++) dst[j] += Ri[j]*yi;
        }
    }
    else if(eq->RinvD_full != NULL)
        for(i = 0; i < n; i++) dst[i] = eq->RinvD_full[i]*y[i];
    else
        for(i = 0; i < n; i++) dst[i] = y[i];
}

/* ---------------------------------------------------------------------------
 * Householder reflectors of the QR factorization
 * -------------------------------------------------------------------------*/

// y <-- (I - tau_k v_k v_k') y, with v_k[k] = 1 implicitly and v_k stored in Q
static void apply_reflector(const DAQPEqElim* eq, const int k, c_float* y){
    const int n = eq->n;
    const c_float* v = eq->Q+(size_t)k*n;
    int i;
    c_float w = y[k];
    for(i = k+1; i < n; i++) w += v[i]*y[i];
    w *= eq->tau[k];
    y[k] -= w;
    for(i = k+1; i < n; i++) y[i] -= w*v[i];
}

// y <-- Q'y, whose leading part is Q1'y
static void apply_QT(const DAQPEqElim* eq, c_float* y){
    int k;
    for(k = 0; k < eq->neq; k++) apply_reflector(eq,k,y);
}

// y <-- Qy
static void apply_Q(const DAQPEqElim* eq, c_float* y){
    int k;
    for(k = eq->neq-1; k >= 0; k--) apply_reflector(eq,k,y);
}

/* ---------------------------------------------------------------------------
 * Deciding whether to eliminate
 * -------------------------------------------------------------------------*/

static int is_eq_elim_eligible(const DAQPWorkspace* work){
    if(work->settings->eq_elim == 0) return 0;
    if(work->avi != NULL || work->bnb != NULL || work->nh > 1) return 0;
    if(work->qp == NULL || work->qp->A == NULL) return 0;
    /*
     * A singular Hessian is handled by the proximal method, which solves a
     * sequence of problems that only differ in their bounds. The reduction is
     * reused over those, but the multipliers of an LP are rescaled by the
     * (adapted) regularization, which the elimination does not account for.
     */
    if(work->n_prox > 0 && work->Rinv == NULL && work->RinvD == NULL) return 0;
    return 1;
}

// Number of equality constraints that are candidates for elimination
static int count_eq(const DAQPWorkspace* work){
    int i, n_eq = 0;
    for(i = work->ms; i < work->m; i++)
        if(IS_EQ_CAND(work,i)) n_eq++;
    return n_eq;
}

/*
 * Whether eliminating n_eq equalities is expected to pay off.
 *
 * Every iteration scans all constraints for the most violated one, which is
 * the dominant cost of a solve. The elimination reduces the length of the
 * rows from n to n-n_eq, but it also turns the simple bounds into general
 * constraints, which are cheaper to scan when they are kept as bounds.
 *
 * The elimination also removes the equalities from the working set, which
 * pays off when many active-set iterations are needed. How many that will be
 * is not known beforehand, so eq_elim > 1 skips the scan comparison.
 */
static int is_eq_elim_worthwhile(const DAQPWorkspace* work, const int n_eq){
    const int n = work->n, m = work->m, ms = work->ms;
    c_float scan_full, scan_reduced;
    if(n_eq <= DAQP_EQ_MIN_COUNT || DAQP_EQ_MIN_RATIO*n_eq <= n) return 0;
    if(work->settings->eq_elim > 1) return 1; // Always eliminate
    // A simple bound is a row of Rinv, or a unit vector for a diagonal Hessian
    scan_full = (work->Rinv != NULL) ? (c_float)ms*n/2 : (c_float)ms;
    scan_full += (c_float)(m-ms-n_eq)*n;
    scan_reduced = (c_float)(m-n_eq)*(n-n_eq);
    return scan_reduced < scan_full;
}

int daqp_eq_will_reduce(const DAQPWorkspace* work){
    if(!is_eq_elim_eligible(work)) return 0;
    return is_eq_elim_worthwhile(work,count_eq(work));
}

// Whether the equality constraints are the ones that the QR was formed for
static int is_qr_valid(const DAQPWorkspace* work){
    const DAQPEqElim* eq = work->eq;
    const int end = eq->neq+eq->nign;
    int i, a = 0, b = eq->neq, n_eq = 0;
    if(eq->neq == 0) return 0;
    for(i = eq->ms; i < eq->m; i++){
        if(!IS_EQ_CAND(work,i)) continue;
        n_eq++;
        // The eliminated and the ignored ids are both stored in ascending order
        if(a < eq->neq && eq->eq_ids[a] == i) a++;
        else if(b < end && eq->eq_ids[b] == i) b++;
        else return 0;
    }
    return n_eq == end;
}

/* ---------------------------------------------------------------------------
 * Forming the reduction
 * -------------------------------------------------------------------------*/

static void allocate_daqp_eq(DAQPWorkspace* work, const int n_eq){
    DAQPEqElim* eq = work->eq;
    const int n = work->n, m = work->m;
    if(eq == NULL){
        eq = calloc(1,sizeof(DAQPEqElim));
        work->eq = eq;
        eq->n = n;
        eq->m = m;
        eq->ms = work->ms;
        // The workspace holds the full problem whenever nothing is installed
        eq->M_full = work->M;
        eq->dupper_full = work->dupper;
        eq->dlower_full = work->dlower;
        eq->scaling_full = work->scaling;
        eq->sense_full = work->sense;
        eq->Mu_full = work->Mu;
        eq->Rinv_full = work->Rinv;
        eq->RinvD_full = work->RinvD;
        eq->v_full = work->v;

        eq->Q = malloc((size_t)n*n*sizeof(c_float));
        eq->xp = malloc(n*sizeof(c_float));
        eq->tmp = malloc(2*(size_t)n*sizeof(c_float));
        eq->map = malloc(m*sizeof(int));
        eq->drop_ids = malloc(m*sizeof(int));
        eq->sense = malloc(m*sizeof(int));
        eq->dupper = malloc(m*sizeof(c_float));
        eq->dlower = malloc(m*sizeof(c_float));
        eq->scaling = malloc(m*sizeof(c_float));
        eq->Mu = (work->Mu == NULL) ? NULL : malloc(m*sizeof(c_float));
    }
    if(eq->ncand != n_eq){
        free(eq->eq_ids);
        free(eq->R);
        free(eq->tau);
        free(eq->s_eq);
        free(eq->y1);
        free(eq->lam_eq);
        eq->eq_ids = malloc(n_eq*sizeof(int));
        eq->R = malloc(((size_t)n_eq*(n_eq+1)/2)*sizeof(c_float));
        eq->tau = malloc(n_eq*sizeof(c_float));
        eq->s_eq = malloc(n_eq*sizeof(c_float));
        eq->y1 = malloc(n_eq*sizeof(c_float));
        eq->lam_eq = calloc(n_eq,sizeof(c_float));
        eq->ncand = n_eq;
    }
}

/*
 * Householder QR of M_E' (one candidate column at a time). Equalities that are
 * linearly dependent on the already factorized ones are not eliminated; their
 * ids are stored after the eliminated ones in eq_ids.
 */
static void build_qr(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    const int n = eq->n;
    const c_float tol = sqrt(work->settings->zero_tol);
    int c, i, j, k = 0;
    for(c = 0; c < eq->ncand && k < n; c++){
        c_float alpha = 0, beta, d, nrm = 0;
        c_float* col = eq->Q+(size_t)k*n;
        const int id = eq->eq_ids[c]; // read before eq_ids[k] is written (k <= c)
        arinv_row(work,id,col);
        for(i = 0; i < n; i++) nrm += col[i]*col[i];
        if(nrm <= work->settings->zero_tol) continue; // Empty constraint
        nrm = 1/sqrt(nrm);
        for(i = 0; i < n; i++) col[i] *= nrm; // Normalized row of M
        for(j = 0; j < k; j++) apply_reflector(eq,j,col);
        for(i = k; i < n; i++) alpha += col[i]*col[i];
        alpha = sqrt(alpha);
        if(alpha <= tol) continue; // Linearly dependent
        beta = (col[k] > 0) ? -alpha : alpha;
        d = col[k]-beta;
        eq->tau[k] = -d/beta;
        for(i = k+1; i < n; i++) col[i] /= d;
        col[k] = beta;
        for(i = 0; i <= k; i++) eq->R[DAQP_ARSUM(k)+i] = col[i];
        eq->s_eq[k] = nrm;
        eq->eq_ids[k++] = id;
    }
    eq->neq = k;
    // Append the ids of the equalities that were not eliminated. Both lists
    // are ascending, since the candidates are factorized in order.
    for(i = eq->ms, j = 0; i < eq->m && k < eq->ncand; i++){
        if(!IS_EQ_CAND(work,i)) continue;
        if(j < eq->neq && eq->eq_ids[j] == i) j++;
        else eq->eq_ids[k++] = i;
    }
    eq->nign = k-eq->neq;

    /*
     * Accumulate Q2 = Q*[0;I] in the trailing columns. Q1 is never formed:
     * the leading columns keep the Householder vectors, which is all that the
     * products with Q1 need.
     */
    for(j = eq->neq; j < n; j++){
        c_float* col = eq->Q+(size_t)j*n;
        for(i = 0; i < n; i++) col[i] = 0;
        col[j] = 1;
    }
    for(k = eq->neq-1; k >= 0; k--){
        const c_float* v = eq->Q+(size_t)k*n;
        const c_float tau = eq->tau[k];
        // Reflect four columns at a time to reuse the Householder vector
        for(j = eq->neq; j+3 < n; j += 4){
            c_float *c0 = eq->Q+(size_t)j*n, *c1 = c0+n, *c2 = c1+n, *c3 = c2+n;
            c_float w0 = c0[k], w1 = c1[k], w2 = c2[k], w3 = c3[k];
            for(i = k+1; i < n; i++){
                const c_float vi = v[i];
                w0 += vi*c0[i]; w1 += vi*c1[i]; w2 += vi*c2[i]; w3 += vi*c3[i];
            }
            w0 *= tau; w1 *= tau; w2 *= tau; w3 *= tau;
            c0[k] -= w0; c1[k] -= w1; c2[k] -= w2; c3[k] -= w3;
            for(i = k+1; i < n; i++){
                const c_float vi = v[i];
                c0[i] -= w0*vi; c1[i] -= w1*vi; c2[i] -= w2*vi; c3[i] -= w3*vi;
            }
        }
        for(; j < n; j++){
            c_float* col = eq->Q+(size_t)j*n;
            c_float w = col[k];
            for(i = k+1; i < n; i++) w += v[i]*col[i];
            w *= tau;
            col[k] -= w;
            for(i = k+1; i < n; i++) col[i] -= w*v[i];
        }
    }
}

/*
 * W = Rinv*Q2, stored row by row. The rows of the reduced problem are A_i*W,
 * the ones of the simple bounds are rows of W, and the solution of the QP is
 * x = xp + W*w.
 */
static void form_W(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    const int n = eq->n, ms = eq->ms, nz = eq->nz;
    const c_float* Q2 = eq->Q+(size_t)eq->neq*n;
    int i, j, k;
    for(i = 0; i < n; i++){
        c_float* Wi = eq->W+(size_t)i*nz;
        if(eq->Rinv_full != NULL){
            const c_float* Ri = eq->Rinv_full+DAQP_R_OFFSET(i,n);
            const c_float s = (i < ms) ? 1/eq->scaling_full[i] : 1;
            for(j = 0; j+3 < nz; j += 4){
                const c_float *q0 = Q2+(size_t)j*n, *q1 = q0+n, *q2 = q1+n, *q3 = q2+n;
                c_float s0 = 0, s1 = 0, s2 = 0, s3 = 0;
                for(k = i; k < n; k++){
                    const c_float rk = Ri[k];
                    s0 += rk*q0[k]; s1 += rk*q1[k]; s2 += rk*q2[k]; s3 += rk*q3[k];
                }
                Wi[j] = s*s0; Wi[j+1] = s*s1; Wi[j+2] = s*s2; Wi[j+3] = s*s3;
            }
            for(; j < nz; j++)
                Wi[j] = s*daqp_dot_inline(Ri+i,Q2+(size_t)j*n+i,n-i);
        }
        else{
            const c_float s = (eq->RinvD_full != NULL) ? eq->RinvD_full[i] : 1;
            for(j = 0; j < nz; j++) Wi[j] = s*Q2[(size_t)j*n+i];
        }
    }
}

// Reduced row of constraint i: A_i*W, or row i of W for a simple bound
static void reduced_row(const DAQPWorkspace* work, const int i, c_float* row){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n, ms = eq->ms, nz = eq->nz;
    int j, k;
    if(i < ms){
        const c_float* Wi = eq->W+(size_t)i*nz;
        for(j = 0; j < nz; j++) row[j] = Wi[j];
        return;
    }
    {
        const c_float* a = work->qp->A+(size_t)(i-ms)*n;
        for(j = 0; j < nz; j++) row[j] = 0;
        for(k = 0; k < n; k++){
            const c_float ak = a[k];
            const c_float* Wk = eq->W+(size_t)k*nz;
            if(ak == 0) continue;
            for(j = 0; j < nz; j++) row[j] += ak*Wk[j];
        }
    }
}

// A_i*xp, or entry i of xp for a simple bound
static c_float row_dot_xp(const DAQPWorkspace* work, const int i){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n, ms = eq->ms;
    if(i < ms) return eq->xp[i];
    return daqp_dot_inline(work->qp->A+(size_t)(i-ms)*n,eq->xp,n);
}

/*
 * Solve R'y1 = d_E and form xp = Rinv*(up-v), which is the part of the
 * solution that the eliminated equality constraints determine.
 */
static void compute_particular(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    const int n = eq->n;
    c_float* up = eq->tmp;
    c_float* Rv = eq->tmp+n;
    int i, k;
    c_float norm2 = 0;
    if(eq->v_full != NULL) rinv_times(work,eq->v_full,Rv);
    else for(i = 0; i < n; i++) Rv[i] = 0;
    for(k = 0; k < eq->neq; k++){
        const int id = eq->eq_ids[k];
        const c_float* Rk = eq->R+DAQP_ARSUM(k); // Column k of R (row k of R')
        const c_float b = (work->sense[id]&DAQP_LOWER) ?
            work->qp->blower[id] : work->qp->bupper[id];
        const c_float* a = work->qp->A+(size_t)(id-eq->ms)*n;
        c_float sum = eq->s_eq[k]*(b+daqp_dot_inline(a,Rv,n));
        for(i = 0; i < k; i++) sum -= Rk[i]*eq->y1[i];
        eq->y1[k] = sum/Rk[k];
        norm2 += eq->y1[k]*eq->y1[k];
    }
    eq->up_norm2 = norm2; // ||up||^2 = ||y1||^2 since Q is orthogonal
    for(i = 0; i < eq->neq; i++) up[i] = eq->y1[i];
    for(; i < n; i++) up[i] = 0;
    apply_Q(eq,up);
    rinv_times(work,up,eq->xp);
    for(i = 0; i < n; i++) eq->xp[i] -= Rv[i]; // xp = Rinv*(up-v)
}

int daqp_eq_reduce(DAQPWorkspace* work, const int mask){
    DAQPEqElim* eq;
    const int n = work->n, m = work->m, ms = work->ms;
    const c_float primal_tol = work->settings->primal_tol;
    const c_float zero_tol = work->settings->zero_tol;
    int i, j, k, rebuild;

    if(!is_eq_elim_eligible(work)) return 0;

    // Start over if the workspace no longer holds the problem that was reduced
    if(work->eq != NULL && (work->eq->n != n || work->eq->m != m ||
                work->eq->ms != ms || work->eq->M_full != work->M))
        free_daqp_eq(work);

    rebuild = (work->eq == NULL) || (mask&(DAQP_UPDATE_Rinv+DAQP_UPDATE_M))
        || !is_qr_valid(work);

    if(rebuild){
        const int n_eq = count_eq(work);
        if(!is_eq_elim_worthwhile(work,n_eq)) return 0;
        allocate_daqp_eq(work,n_eq);
        eq = work->eq;
        for(i = ms, k = 0; i < m; i++)
            if(IS_EQ_CAND(work,i)) eq->eq_ids[k++] = i;
        build_qr(work);
        if(eq->neq == 0 || eq->neq == n){ // Nothing (or everything) eliminated
            eq->neq = 0; // Marks that the QR has to be redone
            return 0;
        }
        if(eq->nz != n-eq->neq){
            free(eq->W);
            free(eq->M);
            eq->nz = n-eq->neq;
            eq->W = malloc((size_t)n*eq->nz*sizeof(c_float));
            eq->M = malloc((size_t)eq->nz*(m-eq->neq)*sizeof(c_float));
        }
        form_W(work);
    }
    eq = work->eq;

    compute_particular(work);

    if(rebuild){
        /*
         * Sweep through the constraints and keep the ones that the reduced
         * variables can affect. The eliminated equalities, the ones that are
         * linearly dependent on them, and the constraints that they imply are
         * all left out, so that no rows without influence are scanned.
         */
        int i_elim = 0, i_ign = eq->neq;
        const int end = eq->neq+eq->nign;
        c_float* row = eq->M;
        j = 0;
        eq->ndrop = 0;
        for(i = 0; i < m; i++){
            c_float nrm2 = 0, c, shift;
            int is_ignored = 0;
            if(i_elim < eq->neq && eq->eq_ids[i_elim] == i){ i_elim++; continue; }
            if(i_ign < end && eq->eq_ids[i_ign] == i){ i_ign++; is_ignored = 1; }
            shift = -row_dot_xp(work,i);
            if(!is_ignored){
                reduced_row(work,i,row);
                for(k = 0; k < eq->nz; k++) nrm2 += row[k]*row[k];
            }
            if(nrm2 <= zero_tol && !DAQP_IS_SOFT(i)){
                /*
                 * The constraint is implied by the equalities: it has to be
                 * consistent, and can then never be violated.
                 */
                if(work->qp->bupper[i]+shift < -primal_tol ||
                        work->qp->blower[i]+shift > primal_tol){
                    eq->neq = 0; // The partial reduction cannot be reused
                    // An inconsistent dependent equality is reported the same
                    // way as when it is detected while forming the working set
                    return is_ignored ? DAQP_EXIT_OVERDETERMINED_INITIAL
                        : DAQP_EXIT_INFEASIBLE;
                }
                eq->drop_ids[eq->ndrop++] = i;
                continue;
            }
            if(DAQP_IS_SOFT(i)){
                /*
                 * The penalty of a soft constraint refers to the normalization
                 * of the full problem, so keep it rather than normalizing the
                 * reduced row (which would change the penalized slack).
                 */
                if(i < ms) c = eq->scaling_full[i];
                else{
                    c_float nrm = 0;
                    arinv_row(work,i,eq->tmp);
                    for(k = 0; k < n; k++) nrm += eq->tmp[k]*eq->tmp[k];
                    c = (nrm <= zero_tol) ? 1 : 1/sqrt(nrm);
                }
            }
            else c = 1/sqrt(nrm2);
            for(k = 0; k < eq->nz; k++) row[k] *= c;
            eq->map[j] = i;
            eq->sense[j] = work->sense[i];
            eq->scaling[j] = c;
            eq->dupper[j] = (work->qp->bupper[i] >= DAQP_INF) ?
                DAQP_INF : c*(work->qp->bupper[i]+shift);
            eq->dlower[j] = (work->qp->blower[i] <= -DAQP_INF) ?
                -DAQP_INF : c*(work->qp->blower[i]+shift);
            row += eq->nz;
            j++;
        }
        eq->m_r = j;
    }
    else{
        // Only the bounds depend on data that can have changed
        if(mask&DAQP_UPDATE_sense)
            for(j = 0; j < eq->m_r; j++) eq->sense[j] = work->sense[eq->map[j]];
        for(j = 0; j < eq->m_r; j++){
            const int id = eq->map[j];
            const c_float c = eq->scaling[j];
            const c_float shift = -row_dot_xp(work,id);
            eq->dupper[j] = (work->qp->bupper[id] >= DAQP_INF) ?
                DAQP_INF : c*(work->qp->bupper[id]+shift);
            eq->dlower[j] = (work->qp->blower[id] <= -DAQP_INF) ?
                -DAQP_INF : c*(work->qp->blower[id]+shift);
        }
        // The constraints that were left out still have to be consistent
        for(k = 0; k < eq->ndrop; k++){
            const int id = eq->drop_ids[k];
            const c_float shift = -row_dot_xp(work,id);
            if(work->qp->bupper[id]+shift < -primal_tol ||
                    work->qp->blower[id]+shift > primal_tol){
                eq->neq = 0;
                return DAQP_EXIT_INFEASIBLE;
            }
        }
    }
    work->reuse_ind = 0; // The right hand side changed

    // Install the reduced problem
    work->M = eq->M;
    work->dupper = eq->dupper;
    work->dlower = eq->dlower;
    work->scaling = eq->scaling;
    work->sense = eq->sense;
    work->Mu = eq->Mu;
    // The reduced problem is a plain least-distance problem
    work->Rinv = NULL;
    work->RinvD = NULL;
    work->v = NULL;
    work->n = eq->nz;
    work->m = eq->m_r;
    work->ms = 0;
    eq->installed = 1;
    eq->expanded = 0;

    return rebuild ? 2 : 1;
}

void daqp_eq_restore(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    if(eq == NULL || eq->installed == 0) return;
    work->M = eq->M_full;
    work->dupper = eq->dupper_full;
    work->dlower = eq->dlower_full;
    work->scaling = eq->scaling_full;
    work->sense = eq->sense_full;
    work->Mu = eq->Mu_full;
    work->Rinv = eq->Rinv_full;
    work->RinvD = eq->RinvD_full;
    work->v = eq->v_full;
    work->n = eq->n;
    work->m = eq->m;
    work->ms = eq->ms;
    eq->installed = 0;
}

/* ---------------------------------------------------------------------------
 * Retrieving the solution of the original problem
 * -------------------------------------------------------------------------*/

/*
 * Multipliers of the eliminated equalities. The reduced problem does not see
 * the part of the stationarity condition that lies in the range of Q1:
 *   R*mu = -(y1 + Q1'*Rinv'*A'*lam),   lam_eq = s_eq*mu,
 * where lam are the multipliers of the constraints that were kept.
 */
static void compute_lam_eq(DAQPWorkspace* work, const c_float* lam){
    DAQPEqElim* eq = work->eq;
    const int n = eq->n, ms = eq->ms, neq = eq->neq;
    c_float* z = eq->tmp;
    c_float* g = eq->tmp+n;
    int i, k;
    for(i = 0; i < n; i++) z[i] = 0;
    for(k = 0; k < eq->m_r; k++){
        const int id = eq->map[k];
        const c_float l = lam[id];
        if(l == 0) continue;
        if(id < ms) z[id] += l;
        else{
            const c_float* a = work->qp->A+(size_t)(id-ms)*n;
            for(i = 0; i < n; i++) z[i] += l*a[i];
        }
    }
    rinvT_times(work,z,g);
    apply_QT(eq,g);
    for(k = 0; k < neq; k++) eq->lam_eq[k] = -(eq->y1[k]+g[k]);
    for(k = neq-1; k >= 0; k--){ // Back substitution
        const c_float* Rk = eq->R+DAQP_ARSUM(k);
        const c_float mu = eq->lam_eq[k]/Rk[k];
        eq->lam_eq[k] = mu;
        for(i = 0; i < k; i++) eq->lam_eq[i] -= Rk[i]*mu;
    }
    for(k = 0; k < neq; k++) eq->lam_eq[k] *= eq->s_eq[k];
}

// x <-- xp + W*w, where w solves the reduced problem
static void expand_solution(const DAQPWorkspace* work, const c_float* w, c_float* x){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n, nz = eq->nz;
    int i, j;
    for(i = 0; i < n; i++){
        const c_float* Wi = eq->W+(size_t)i*nz;
        c_float sum = eq->xp[i];
        for(j = 0; j < nz; j++) sum += Wi[j]*w[j];
        x[i] = sum;
    }
}

/*
 * Turn the solution of the reduced problem into the solution of the full one,
 * and hand the full problem back to the workspace.
 */
void daqp_eq_expand(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    int i;
    for(i = 0; i < work->n_active; i++) // Scale to the original constraints
        work->lam_star[i] *= eq->scaling[work->WS[i]];
    work->fval += eq->up_norm2; // ||u||^2 = ||up||^2+||w||^2
    expand_solution(work,work->u,eq->tmp);
    daqp_eq_restore(work);
    for(i = 0; i < eq->n; i++) work->x[i] = eq->tmp[i];
    eq->expanded = 1;
}

/*
 * Eliminate the equality constraints of the problem that the workspace holds.
 * Returns the number of eliminated constraints (0 if it was left intact) or a
 * negative exit flag if the equality constraints cannot be satisfied.
 *
 * A singular Hessian is handled by the proximal method, which re-forms the
 * problem in the full space between its iterations. The reduction is then
 * prepared but not installed; daqp_prox applies it to each inner problem.
 */
/*
 * Form the constraints of the full problem, which the setup leaves to the
 * elimination whenever one is expected, and put it in the workspace.
 */
static int daqp_eq_form_full(DAQPWorkspace* work){
    int error_flag;
    daqp_eq_restore(work);
    if(work->eq != NULL) work->eq->neq = 0; // Nothing to retrieve
    if(work->qp == NULL || work->qp->A == NULL) return 0;
    error_flag = daqp_update_M(work,work->qp->A,0);
    if(error_flag < 0) return error_flag;
    daqp_update_d(work,work->qp->bupper,work->qp->blower);
    reset_daqp_workspace(work);
    return daqp_activate_constraints(work);
}

int daqp_eq_eliminate(DAQPWorkspace* work){
    int flag, error_flag;
    /*
     * Only the bounds are known to have changed here; daqp_update_ldp marks
     * the factorization as invalid when the data it is formed from changes.
     */
    flag = daqp_eq_reduce(work,DAQP_UPDATE_d);
    if(flag <= 0){
        // The setup leaves forming the constraints to the elimination when one
        // is expected, so form them here if no reduction was installed
        if(daqp_eq_will_reduce(work)){
            error_flag = daqp_eq_form_full(work);
            if(error_flag < 0) return error_flag;
        }
        else if(work->eq != NULL) work->eq->neq = 0;
        return (flag < 0) ? flag : 0;
    }
    // Form the working set of the reduced problem
    reset_daqp_workspace(work);
    error_flag = daqp_activate_constraints(work);
    if(error_flag < 0) return error_flag;
    if(work->n_prox > 0) daqp_eq_restore(work);
    return work->eq->neq;
}

/*
 * Retrieve the solution of the original problem from the solution of the
 * reduced one. Called after an ordinary solve, so that nothing in between has
 * to be aware of the elimination.
 */
void daqp_eq_retrieve(DAQPResult* res, DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    int i, j;
    if(eq == NULL || eq->neq == 0) return;
    if(eq->installed) expand_solution(work,work->u,eq->tmp);
    if((eq->installed || eq->expanded) && res->lam != NULL && work->nh < 2){
        // Scatter the multipliers onto the constraints of the original problem
        const int m_r = eq->m_r;
        for(i = eq->m; i > m_r; ) res->lam[--i] = 0;
        for(j = m_r-1; j >= 0; j--){
            const c_float l = res->lam[j];
            res->lam[j] = 0;
            res->lam[eq->map[j]] = l;
        }
    }
    if(eq->installed){
        // The solve ran on the reduced problem
        res->fval = work->fval+eq->up_norm2;
        daqp_eq_restore(work);
        for(i = 0; i < eq->n; i++){
            work->x[i] = eq->tmp[i];
            res->x[i] = eq->tmp[i];
        }
        if(work->v != NULL) // Shift back the objective function value
            for(i = 0; i < eq->n; i++) res->fval -= work->v[i]*work->v[i];
        res->fval *= 0.5;
    }
    eq->expanded = 0;
    if(res->lam != NULL && work->nh < 2){
        compute_lam_eq(work,res->lam);
        for(i = 0; i < eq->neq; i++) res->lam[eq->eq_ids[i]] = eq->lam_eq[i];
        for(i = 0; i < eq->ndrop; i++) res->lam[eq->drop_ids[i]] = 0;
    }
}

void free_daqp_eq(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    if(eq == NULL) return;
    daqp_eq_restore(work);
    free(eq->eq_ids);
    free(eq->drop_ids);
    free(eq->map);
    free(eq->Q);
    free(eq->R);
    free(eq->tau);
    free(eq->s_eq);
    free(eq->y1);
    free(eq->lam_eq);
    free(eq->W);
    free(eq->xp);
    free(eq->tmp);
    free(eq->M);
    free(eq->sense);
    free(eq->dupper);
    free(eq->dlower);
    free(eq->scaling);
    free(eq->Mu);
    free(eq);
    work->eq = NULL;
}
