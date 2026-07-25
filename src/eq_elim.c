#include "api.h"
#include "eq_elim.h"
#include "constants.h"
#include "factorization.h"
#include <stdlib.h>
#include <math.h>

// A constraint is eliminated if it is a general equality constraint
// (equalities are marked as active and immutable in sense)
#define IS_EQ_CAND(work,i) \
    (((work)->sense[i]&(DAQP_ACTIVE+DAQP_IMMUTABLE+DAQP_SOFT+DAQP_BINARY)) \
     == (DAQP_ACTIVE+DAQP_IMMUTABLE))

// Expand row i of the full (normalized) LDP constraint matrix into dst
static void ldp_row(const DAQPWorkspace* work, const int i, c_float* dst){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n;
    int j;
    if(i < eq->ms){
        for(j = 0; j < n; j++) dst[j] = 0;
        if(eq->Rinv_full != NULL){
            const c_float* R = eq->Rinv_full+DAQP_R_OFFSET(i,n);
            for(j = i; j < n; j++) dst[j] = R[j];
        }
        else dst[i] = 1; // Identity/diagonal Hessian => normalized row is a unit vector
    }
    else{
        const c_float* Mi = eq->M_full+n*(i-eq->ms);
        for(j = 0; j < n; j++) dst[j] = Mi[j];
    }
}

// Inner product between row i of the full LDP constraint matrix and v
static c_float ldp_row_dot(const DAQPWorkspace* work, const int i, const c_float* v){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n;
    if(i >= eq->ms) return daqp_dot_inline(eq->M_full+n*(i-eq->ms),v,n);
    if(eq->Rinv_full == NULL) return v[i];
    return daqp_dot_inline(eq->Rinv_full+DAQP_R_OFFSET(i,n)+i,v+i,n-i);
}

// y <-- y + alpha*(row i of the full LDP constraint matrix)
static void ldp_row_axpy(const DAQPWorkspace* work, const int i, const c_float alpha,
        c_float* y){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n;
    int j;
    if(i >= eq->ms){
        const c_float* Mi = eq->M_full+n*(i-eq->ms);
        for(j = 0; j < n; j++) y[j] += alpha*Mi[j];
    }
    else if(eq->Rinv_full == NULL) y[i] += alpha;
    else{
        const c_float* R = eq->Rinv_full+DAQP_R_OFFSET(i,n);
        for(j = i; j < n; j++) y[j] += alpha*R[j];
    }
}

// y <-- (I - tau_k v_k v_k') y, with v_k[k] = 1 implicitly and v_k stored in Q
static void apply_reflector(const DAQPEqElim* eq, const int k, c_float* y){
    const int n = eq->n;
    const c_float* v = eq->Q+k*n;
    int i;
    c_float w = y[k];
    for(i = k+1; i < n; i++) w += v[i]*y[i];
    w *= eq->tau[k];
    y[k] -= w;
    for(i = k+1; i < n; i++) y[i] -= w*v[i];
}

// y <-- Q'y, i.e., [Q1'y; Q2'y]
static void apply_QT(const DAQPEqElim* eq, c_float* y){
    int k;
    for(k = 0; k < eq->neq; k++) apply_reflector(eq,k,y);
}

// y <-- Qy
static void apply_Q(const DAQPEqElim* eq, c_float* y){
    int k;
    for(k = eq->neq-1; k >= 0; k--) apply_reflector(eq,k,y);
}

/*
 * Mrow <-- Q2'*(row i of the full LDP constraint matrix), returning ||Mrow||^2.
 * This projection dominates the setup, so the columns of Q2 are traversed in
 * blocks to keep the constraint row in registers and the accumulations
 * independent.
 */
static c_float project_row(const DAQPWorkspace* work, const int i, c_float* Mrow){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n, nz = n-eq->neq;
    const c_float* Q2 = eq->Q+eq->neq*n;
    const c_float* Mi;
    int j, k, kstart;
    c_float nrm2 = 0;
    if(i >= eq->ms){
        Mi = eq->M_full+n*(i-eq->ms);
        kstart = 0;
    }
    else if(eq->Rinv_full == NULL){
        // Identity/diagonal Hessian: the row is a unit vector
        for(j = 0; j < nz; j++){
            const c_float q = Q2[j*n+i];
            Mrow[j] = q;
            nrm2 += q*q;
        }
        return nrm2;
    }
    else{
        Mi = eq->Rinv_full+DAQP_R_OFFSET(i,n);
        kstart = i;
    }
    for(j = 0; j+3 < nz; j += 4){
        const c_float *q0 = Q2+j*n, *q1 = q0+n, *q2 = q1+n, *q3 = q2+n;
        c_float s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for(k = kstart; k < n; k++){
            const c_float mk = Mi[k];
            s0 += mk*q0[k];
            s1 += mk*q1[k];
            s2 += mk*q2[k];
            s3 += mk*q3[k];
        }
        Mrow[j] = s0; Mrow[j+1] = s1; Mrow[j+2] = s2; Mrow[j+3] = s3;
        nrm2 += s0*s0+s1*s1+s2*s2+s3*s3;
    }
    for(; j < nz; j++){
        const c_float s = daqp_dot_inline(Mi+kstart,Q2+j*n+kstart,n-kstart);
        Mrow[j] = s;
        nrm2 += s*s;
    }
    return nrm2;
}

/*
 * Project all general constraints onto the nullspace of the equalities.
 * Four constraints are projected against two columns of Q2 at a time, which
 * keeps Q2 from being streamed once per constraint. The eliminated rows are
 * projected as well (they are overwritten with zeros afterwards) so that the
 * blocking is not broken up.
 */
static void project_general(const DAQPWorkspace* work){
    const DAQPEqElim* eq = work->eq;
    const int n = eq->n, nz = n-eq->neq, ms = eq->ms, m = eq->m;
    const c_float* Q2 = eq->Q+eq->neq*n;
    int i, j, k;
    for(i = ms; i+3 < m; i += 4){
        const c_float *m0 = eq->M_full+n*(i-ms), *m1 = m0+n, *m2 = m1+n, *m3 = m2+n;
        c_float *r0 = eq->M+(size_t)nz*i, *r1 = r0+nz, *r2 = r1+nz, *r3 = r2+nz;
        for(j = 0; j+1 < nz; j += 2){
            const c_float *q0 = Q2+j*n, *q1 = q0+n;
            c_float s00 = 0, s01 = 0, s10 = 0, s11 = 0;
            c_float s20 = 0, s21 = 0, s30 = 0, s31 = 0;
            for(k = 0; k < n; k++){
                const c_float b0 = q0[k], b1 = q1[k];
                c_float a = m0[k];
                s00 += a*b0; s01 += a*b1;
                a = m1[k];
                s10 += a*b0; s11 += a*b1;
                a = m2[k];
                s20 += a*b0; s21 += a*b1;
                a = m3[k];
                s30 += a*b0; s31 += a*b1;
            }
            r0[j] = s00; r0[j+1] = s01;
            r1[j] = s10; r1[j+1] = s11;
            r2[j] = s20; r2[j+1] = s21;
            r3[j] = s30; r3[j+1] = s31;
        }
        for(; j < nz; j++){
            const c_float* q = Q2+j*n;
            r0[j] = daqp_dot_inline(m0,q,n);
            r1[j] = daqp_dot_inline(m1,q,n);
            r2[j] = daqp_dot_inline(m2,q,n);
            r3[j] = daqp_dot_inline(m3,q,n);
        }
    }
    for(; i < m; i++) project_row(work,i,eq->M+(size_t)nz*i);
}

static int is_eq_elim_eligible(const DAQPWorkspace* work){
    if(work->settings->eq_elim == 0) return 0;
    if(work->avi != NULL || work->bnb != NULL || work->nh > 1) return 0;
    /*
     * A singular Hessian is handled by the proximal method, which solves a
     * sequence of LDPs that only differ in their bounds. The reduction is
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
    scan_reduced = (c_float)m*(n-n_eq);
    return scan_reduced < scan_full;
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

static void allocate_daqp_eq(DAQPWorkspace* work, const int n_eq){
    DAQPEqElim* eq = work->eq;
    const int n = work->n, m = work->m;
    if(eq == NULL){
        eq = calloc(1,sizeof(DAQPEqElim));
        work->eq = eq;
        eq->n = n;
        eq->m = m;
        eq->ms = work->ms;
        // The workspace holds the full LDP whenever nothing is installed
        eq->M_full = work->M;
        eq->dupper_full = work->dupper;
        eq->dlower_full = work->dlower;
        eq->scaling_full = work->scaling;
        eq->Mu_full = work->Mu;
        eq->Rinv_full = work->Rinv;
        eq->RinvD_full = work->RinvD;
        eq->v_full = work->v;

        eq->Q = malloc(n*n*sizeof(c_float));
        eq->up = malloc(n*sizeof(c_float));
        eq->tmp = malloc(n*sizeof(c_float));
        eq->t = malloc(m*sizeof(c_float));
        eq->dupper = malloc(m*sizeof(c_float));
        eq->dlower = malloc(m*sizeof(c_float));
        eq->scaling = malloc(m*sizeof(c_float));
        eq->Mu = (work->Mu == NULL) ? NULL : malloc(m*sizeof(c_float));
    }
    if(eq->ncand != n_eq){
        free(eq->eq_ids);
        free(eq->R);
        free(eq->tau);
        free(eq->y1);
        free(eq->lam_eq);
        eq->eq_ids = malloc(n_eq*sizeof(int));
        eq->R = malloc((n_eq*(n_eq+1)/2)*sizeof(c_float));
        eq->tau = malloc(n_eq*sizeof(c_float));
        eq->y1 = malloc(n_eq*sizeof(c_float));
        eq->lam_eq = calloc(n_eq,sizeof(c_float));
        eq->ncand = n_eq;
    }
}

/*
 * Householder QR of M_E' (one candidate column at a time). Equalities that are
 * linearly dependent on the already factorized ones are not eliminated; their
 * ids are stored after the eliminated ones in eq_ids. The orthogonal factor is
 * accumulated explicitly, which makes the projection of the constraints and the
 * expansion of the solution ordinary matrix-vector products.
 */
static void build_qr(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    const int n = eq->n;
    const c_float tol = sqrt(work->settings->zero_tol);
    int c, i, j, k = 0;
    for(c = 0; c < eq->ncand && k < n; c++){
        c_float alpha = 0, beta, d;
        c_float* col = eq->Q+k*n;
        const int id = eq->eq_ids[c]; // read before eq_ids[k] is written (k <= c)
        ldp_row(work,id,col);
        for(j = 0; j < k; j++) apply_reflector(eq,j,col);
        for(i = k; i < n; i++) alpha += col[i]*col[i];
        alpha = sqrt(alpha);
        if(alpha <= tol) continue; // Linearly dependent
        beta = (col[k] > 0) ? -alpha : alpha;
        d = col[k]-beta;
        eq->tau[k] = -d/beta;
        for(i = k+1; i < n; i++) col[i] /= d;
        col[k] = beta;
        for(i = 0; i <= k; i++) eq->R[DAQP_ARSUM(k)+i] = col[i]; // R is overwritten below
        eq->eq_ids[k++] = id;
    }
    eq->neq = k;
    // Append the ids of the equalities that were not eliminated
    for(i = eq->ms, j = 0; i < eq->m; i++){
        if(!IS_EQ_CAND(work,i)) continue;
        if(j < eq->neq && eq->eq_ids[j] == i) j++;
        else eq->eq_ids[k++] = i;
    }
    eq->nign = k-eq->neq;

    /*
     * Accumulate Q2 = Q*[0;I] in the trailing columns. Q1 is never formed:
     * the leading columns keep the Householder vectors (and R above the
     * diagonal), which is all that the Q1-products in the pre/postprocessing
     * need, and only the projection requires an explicit factor.
     */
    for(j = eq->neq; j < n; j++){
        c_float* col = eq->Q+j*n;
        for(i = 0; i < n; i++) col[i] = 0;
        col[j] = 1;
    }
    for(k = eq->neq-1; k >= 0; k--){
        const c_float* v = eq->Q+k*n;
        const c_float tau = eq->tau[k];
        // Reflect four columns at a time to reuse the Householder vector
        for(j = eq->neq; j+3 < n; j += 4){
            c_float *c0 = eq->Q+j*n, *c1 = c0+n, *c2 = c1+n, *c3 = c2+n;
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
            c_float* col = eq->Q+j*n;
            c_float w = col[k];
            for(i = k+1; i < n; i++) w += v[i]*col[i];
            w *= tau;
            col[k] -= w;
            for(i = k+1; i < n; i++) col[i] -= w*v[i];
        }
    }
}

/*
 * Solve R'y1 = d_E and form the particular solution up = Q1*y1.
 * Only the equality rows of d are needed, so they are formed here
 * (d_i = scaling_i*b_i + m_i'v) rather than taken from the full LDP.
 */
static void compute_particular(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    const int n = eq->n;
    int i, k;
    c_float sum, norm2 = 0;
    for(k = 0; k < eq->neq; k++){
        const int id = eq->eq_ids[k];
        const c_float* Rk = eq->R+DAQP_ARSUM(k); // Column k of R (row k of R')
        const c_float b = (work->sense[id]&DAQP_LOWER) ?
            work->qp->blower[id] : work->qp->bupper[id];
        sum = b*eq->scaling_full[id];
        if(eq->v_full != NULL) sum += ldp_row_dot(work,id,eq->v_full);
        for(i = 0; i < k; i++) sum -= Rk[i]*eq->y1[i];
        eq->y1[k] = sum/Rk[k];
        norm2 += eq->y1[k]*eq->y1[k];
    }
    eq->up_norm2 = norm2; // ||up||^2 = ||y1||^2 since Q is orthogonal
    for(i = 0; i < eq->neq; i++) eq->up[i] = eq->y1[i];
    for(; i < n; i++) eq->up[i] = 0;
    apply_Q(eq,eq->up); // up = Q*[y1; 0]
}

/*
 * Whether the equality constraints will be eliminated. Used to skip work that
 * the elimination would make redundant, namely forming the bounds of the full
 * LDP and putting the equality constraints in the working set.
 */
int daqp_eq_will_reduce(const DAQPWorkspace* work){
    if(!is_eq_elim_eligible(work)) return 0;
    return is_eq_elim_worthwhile(work,count_eq(work));
}

int daqp_eq_reduce(DAQPWorkspace* work, const int mask){
    DAQPEqElim* eq;
    const int n = work->n, m = work->m;
    const c_float primal_tol = work->settings->primal_tol;
    const c_float zero_tol = work->settings->zero_tol;
    int i, j, k, nz, rebuild;
    int i_elim, i_ign;
    c_float* Mrow;

    if(!is_eq_elim_eligible(work)) return 0;

    // Start over if the workspace no longer holds the LDP that was reduced
    if(work->eq != NULL && (work->eq->n != n || work->eq->m != m ||
                work->eq->ms != work->ms || work->eq->M_full != work->M))
        free_daqp_eq(work);

    rebuild = (work->eq == NULL) || (mask&(DAQP_UPDATE_Rinv+DAQP_UPDATE_M))
        || !is_qr_valid(work);

    if(rebuild){
        const int n_eq = count_eq(work);
        if(!is_eq_elim_worthwhile(work,n_eq)) return 0;
        allocate_daqp_eq(work,n_eq);
        eq = work->eq;
        for(i = eq->ms, k = 0; i < m; i++)
            if(IS_EQ_CAND(work,i)) eq->eq_ids[k++] = i;
        build_qr(work);
        if(eq->neq == 0 || eq->neq == n){ // Nothing (or everything) eliminated
            eq->neq = 0; // Marks that the QR has to be redone
            return 0;
        }
        nz = n-eq->neq;
        if(eq->M == NULL || eq->nz != nz){
            free(eq->M);
            eq->nz = nz;
            eq->M = malloc(nz*m*sizeof(c_float));
        }
    }
    eq = work->eq;
    nz = n-eq->neq;

    // Nothing that the reduced LDP depends on changed: only reinstall it
    if(!rebuild && !(mask&(DAQP_UPDATE_Rinv+DAQP_UPDATE_M
                    +DAQP_UPDATE_v+DAQP_UPDATE_d)))
        goto install;

    compute_particular(work);

    /*
     * The reduced bounds only depend on v and up through v-up, which lets the
     * shift be folded into the (single) sweep through the constraints below.
     */
    if(eq->v_full != NULL)
        for(i = 0; i < n; i++) eq->tmp[i] = eq->v_full[i]-eq->up[i];
    else
        for(i = 0; i < n; i++) eq->tmp[i] = -eq->up[i];
    work->reuse_ind = 0; // The right hand side changed

    /*
     * Shift and project the remaining constraints. The reduced rows are
     * normalized, and the normalization is included in scaling to retrieve
     * multipliers in the original scaling.
     */
    if(rebuild) project_general(work);

    i_elim = 0;
    i_ign = eq->neq;
    Mrow = eq->M;
    for(i = 0; i < m; i++, Mrow += nz){
        c_float Mvu, nrm2;
        int is_eliminated = 0, is_ignored = 0;
        if(i_elim < eq->neq && eq->eq_ids[i_elim] == i){
            i_elim++;
            is_eliminated = 1;
        }
        else if(i_ign < eq->neq+eq->nign && eq->eq_ids[i_ign] == i){
            i_ign++;
            is_ignored = 1;
        }
        else if(rebuild == 0 && eq->t[i] != 0){
            // The projection is unchanged; only the bounds have to be redone
            Mvu = ldp_row_dot(work,i,eq->tmp);
            eq->dupper[i] = (work->qp->bupper[i] >= DAQP_INF) ? DAQP_INF :
                eq->t[i]*(work->qp->bupper[i]*eq->scaling_full[i]+Mvu);
            eq->dlower[i] = (work->qp->blower[i] <= -DAQP_INF) ? -DAQP_INF :
                eq->t[i]*(work->qp->blower[i]*eq->scaling_full[i]+Mvu);
            continue;
        }

        eq->scaling[i] = eq->scaling_full[i];
        if(is_eliminated){
            // The row is in the range of the eliminated equalities
            if(rebuild) for(j = 0; j < nz; j++) Mrow[j] = 0;
            eq->t[i] = 0;
            eq->dupper[i] = DAQP_INF;
            eq->dlower[i] = -DAQP_INF;
            continue;
        }

        Mvu = ldp_row_dot(work,i,eq->tmp);

        if(rebuild){
            // The general constraints were projected above
            if(i < eq->ms) project_row(work,i,Mrow);
            for(j = 0, nrm2 = 0; j < nz; j++) nrm2 += Mrow[j]*Mrow[j];
            // Soft constraints are not renormalized, to keep the slack (and
            // hence the penalty) identical to the one in the full problem
            if(DAQP_IS_SOFT(i)) eq->t[i] = 1;
            else if(nrm2 <= zero_tol) eq->t[i] = 0; // Implied by the equalities
            else{
                eq->t[i] = 1/sqrt(nrm2);
                for(j = 0; j < nz; j++) Mrow[j] *= eq->t[i];
            }
        }

        if(eq->t[i] == 0){
            // The constraint is implied by the equalities: it has to be
            // consistent, and can then never be violated
            eq->dupper[i] = DAQP_INF;
            eq->dlower[i] = -DAQP_INF;
            if(work->qp->bupper[i]*eq->scaling_full[i]+Mvu < -primal_tol ||
                    work->qp->blower[i]*eq->scaling_full[i]+Mvu > primal_tol){
                eq->neq = 0; // The partial reduction cannot be reused
                // An inconsistent dependent equality is reported the same way
                // as when it is detected while forming the working set
                return is_ignored ? DAQP_EXIT_OVERDETERMINED_INITIAL
                    : DAQP_EXIT_INFEASIBLE;
            }
            continue;
        }

        eq->scaling[i] = eq->t[i]*eq->scaling_full[i];
        eq->dupper[i] = (work->qp->bupper[i] >= DAQP_INF) ? DAQP_INF :
            eq->t[i]*(work->qp->bupper[i]*eq->scaling_full[i]+Mvu);
        eq->dlower[i] = (work->qp->blower[i] <= -DAQP_INF) ? -DAQP_INF :
            eq->t[i]*(work->qp->blower[i]*eq->scaling_full[i]+Mvu);
    }

install:
    // Install the reduced LDP (the constraint indexing is kept)
    work->M = eq->M;
    work->dupper = eq->dupper;
    work->dlower = eq->dlower;
    work->scaling = eq->scaling;
    work->Mu = eq->Mu;
    // The reduced problem is a plain least-distance problem
    work->Rinv = NULL;
    work->RinvD = NULL;
    work->v = NULL;
    work->n = nz;
    work->ms = 0;
    for(k = 0; k < eq->neq+eq->nign; k++)
        work->sense[eq->eq_ids[k]] = DAQP_IMMUTABLE; // Never activate
    eq->installed = 1;

    return rebuild ? 2 : 1;
}

void daqp_eq_restore(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    int k;
    if(eq == NULL || eq->installed == 0) return;
    work->M = eq->M_full;
    work->dupper = eq->dupper_full;
    work->dlower = eq->dlower_full;
    work->scaling = eq->scaling_full;
    work->Mu = eq->Mu_full;
    work->Rinv = eq->Rinv_full;
    work->RinvD = eq->RinvD_full;
    work->v = eq->v_full;
    work->n = eq->n;
    work->m = eq->m;
    work->ms = eq->ms;
    for(k = 0; k < eq->neq+eq->nign; k++)
        work->sense[eq->eq_ids[k]] |= DAQP_ACTIVE+DAQP_IMMUTABLE;
    eq->installed = 0;
}

/*
 * Multipliers of the eliminated equalities. They are determined by the part of
 * the stationarity condition that the reduced problem does not see:
 *   R*lam_eq = -Q1'*(u + sum_i lam_i*m_i).
 * The multipliers of the reduced constraints are given in the normalization of
 * the reduced problem, unless they have already been scaled back (is_scaled).
 */
// u <-- up + Q2*w, where w is the iterate of the reduced problem
static void expand_iterate(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    const int n = eq->n, neq = eq->neq, nz = n-neq;
    c_float* w = eq->tmp;
    int i, j;
    // ||u||^2 = ||up||^2 + ||w||^2 since up is orthogonal to the nullspace
    work->fval += eq->up_norm2;
    for(j = 0; j < nz; j++) w[j] = work->u[j];
    for(i = 0; i < n; i++) work->u[i] = eq->up[i];
    for(j = 0; j < nz; j++){
        const c_float* q = eq->Q+(neq+j)*n;
        const c_float wj = w[j];
        for(i = 0; i < n; i++) work->u[i] += wj*q[i];
    }
}

static void compute_lam_eq(DAQPWorkspace* work, const int is_scaled){
    DAQPEqElim* eq = work->eq;
    const int n = eq->n, neq = eq->neq;
    c_float* g = eq->tmp;
    int i, k;
    for(i = 0; i < n; i++) g[i] = work->u[i];
    for(k = 0; k < work->n_active; k++){
        const int id = work->WS[k];
        const c_float lam = is_scaled ? work->lam_star[k]/eq->scaling_full[id]
            : eq->t[id]*work->lam_star[k];
        ldp_row_axpy(work,id,lam,g);
    }
    apply_QT(eq,g); // The leading part of Q'g is Q1'g
    for(k = 0; k < neq; k++) eq->lam_eq[k] = -g[k];
    for(k = neq-1; k >= 0; k--){ // Back substitution
        const c_float* Rk = eq->R+DAQP_ARSUM(k);
        const c_float lam = eq->lam_eq[k]/Rk[k];
        eq->lam_eq[k] = lam;
        for(i = 0; i < k; i++) eq->lam_eq[i] -= Rk[i]*lam;
    }
    // Scale the multipliers to the original constraint scaling
    for(k = 0; k < neq; k++) eq->lam_eq[k] *= eq->scaling_full[eq->eq_ids[k]];
}

/*
 * Expand the reduced iterate w into the iterate u of the full LDP and put the
 * full problem back into the workspace. Everything downstream (retrieving the
 * QP solution in particular) is then oblivious to the elimination.
 */
void daqp_eq_expand(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    int k;
    expand_iterate(work);
    compute_lam_eq(work,0);
    /*
     * The multipliers are given in the normalization of the reduced
     * constraints; put them in the one of the full LDP, so that retrieving
     * the QP solution scales them like any other multiplier.
     */
    for(k = 0; k < work->n_active; k++)
        work->lam_star[k] *= eq->t[work->WS[k]];
    daqp_eq_restore(work);
}

/*
 * Eliminate the equality constraints of the LDP that the workspace holds.
 * Returns the number of eliminated constraints (0 if the LDP was left intact).
 *
 * A singular Hessian is handled by the proximal method, which re-forms the
 * problem in the full space between its iterations. The reduction is then
 * prepared but not installed; daqp_prox applies it to each inner problem.
 */
int daqp_eq_eliminate(DAQPWorkspace* work){
    const int flag = daqp_eq_reduce(work,DAQP_UPDATE_M+DAQP_UPDATE_d);
    int error_flag;
    if(flag < 0){
        // Leave the full problem behind, so that a workspace whose caller
        // ignores the error still describes the problem that was set up
        daqp_eq_restore(work);
        if(work->eq != NULL) work->eq->neq = 0;
        reset_daqp_workspace(work);
        daqp_activate_constraints(work);
        return flag; // The equality constraints are inconsistent
    }
    if(flag == 0){
        if(work->eq != NULL) work->eq->neq = 0; // Nothing to retrieve
        // The setup may have left the working set to the elimination
        if(work->n_active == 0){
            reset_daqp_workspace(work);
            error_flag = daqp_activate_constraints(work);
            if(error_flag < 0) return error_flag;
        }
        return 0;
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
 * reduced one. Called after an ordinary solve, so that nothing in between
 * has to be aware of the elimination.
 */
void daqp_eq_retrieve(DAQPResult* res, DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    int i;
    if(eq == NULL || eq->neq == 0) return;
    if(eq->installed){
        // The solve ran on the reduced LDP: expand it into the full problem
        // and transform it into the solution of the QP
        expand_iterate(work);
        compute_lam_eq(work,1);
        daqp_eq_restore(work);
        ldp2qp_solution(work);
        for(i = 0; i < eq->n; i++) res->x[i] = work->x[i];
        if(work->v != NULL && (work->Rinv != NULL || work->RinvD != NULL)){
            res->fval = work->fval;
            for(i = 0; i < eq->n; i++) res->fval -= work->v[i]*work->v[i];
            res->fval *= 0.5;
        }
    }
    if(res->lam != NULL && work->nh < 2)
        for(i = 0; i < eq->neq; i++) res->lam[eq->eq_ids[i]] = eq->lam_eq[i];
}

void free_daqp_eq(DAQPWorkspace* work){
    DAQPEqElim* eq = work->eq;
    if(eq == NULL) return;
    daqp_eq_restore(work);
    free(eq->eq_ids);
    free(eq->Q);
    free(eq->R);
    free(eq->tau);
    free(eq->y1);
    free(eq->lam_eq);
    free(eq->up);
    free(eq->tmp);
    free(eq->t);
    free(eq->M);
    free(eq->dupper);
    free(eq->dlower);
    free(eq->scaling);
    free(eq->Mu);
    free(eq);
    work->eq = NULL;
}
