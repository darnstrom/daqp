#include "presolve.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "auxiliary.h"
#include "utils.h"

/* --------------------------------------------------------------------------
 * presolve_problem  --  tighten a mixed-integer problem before branch and bound
 *
 * Bound propagation over the general constraints tightens the simple bounds
 * (and fixes binaries), and coefficient tightening strengthens the
 * coefficients of binaries in one-sided constraints (e.g., big-M constraints).
 * Neither changes the set of integer-feasible points, but both tighten the
 * relaxations, which reduces the number of nodes. qp is modified in place, and
 * the coefficients are tightened with the (unpropagated) bounds of bq.
 * Soft constraints and bounds are not relied on, and immutable ones (and the
 * Hessian) are not changed.
 *
 * Returns the number of tightenings, or DAQP_EXIT_INFEASIBLE if the problem is
 * infeasible.
 * --------------------------------------------------------------------------*/

#define DAQP_PRESOLVE_IS_INF(v) ((v) >= DAQP_INF || (v) <= -DAQP_INF)

// Bounds of x_j that hold (only the first ms variables have simple bounds, and
// soft bounds can be violated)
static c_float presolve_lower(const DAQPProblem* qp, const int j){
    return (j < qp->ms && !(qp->sense[j] & DAQP_SOFT)) ? qp->blower[j] : -DAQP_INF;
}
static c_float presolve_upper(const DAQPProblem* qp, const int j){
    return (j < qp->ms && !(qp->sense[j] & DAQP_SOFT)) ? qp->bupper[j] : DAQP_INF;
}

// Tighten the upper (or lower) bound of x_j to v. Returns 1 if it was
// tightened, 0 if not, and DAQP_EXIT_INFEASIBLE if the bounds cross
static int presolve_bound(DAQPProblem* qp, const int j, c_float v, const int upper){
    c_float *bu = qp->bupper+j, *bl = qp->blower+j, tol;
    if(j >= qp->ms || (qp->sense[j] & (DAQP_SOFT+DAQP_IMMUTABLE))) return 0;
    if(qp->sense[j] & DAQP_BINARY){ // x_j in {bl, bu}: fix it if one value is cut off
        tol = DAQP_PRESOLVE_FEAS_TOL*(1+fabs(*bu)+fabs(*bl));
        if(upper ? v < *bl - tol : v > *bu + tol) return DAQP_EXIT_INFEASIBLE;
        if(*bl == *bu || (upper ? v >= *bu - tol : v <= *bl + tol)) return 0;
        if(upper) *bu = *bl;
        else *bl = *bu;
        return 1;
    }
    v += (upper ? 1 : -1)*DAQP_PRESOLVE_MARGIN*(1+fabs(v)); // Safety margin
    tol = DAQP_PRESOLVE_MIN_IMPROVEMENT*((DAQP_PRESOLVE_IS_INF(*bu) || DAQP_PRESOLVE_IS_INF(*bl)) ?
            1+fabs(v) : *bu-*bl);
    if(upper){
        if(v >= *bu - tol) return 0;
        if(v < *bl - DAQP_PRESOLVE_FEAS_TOL*(1+fabs(*bl))) return DAQP_EXIT_INFEASIBLE;
        *bu = v < *bl ? *bl : v;
    }
    else{
        if(v <= *bl + tol) return 0;
        if(v > *bu + DAQP_PRESOLVE_FEAS_TOL*(1+fabs(*bu))) return DAQP_EXIT_INFEASIBLE;
        *bl = v > *bu ? *bu : v;
    }
    return 1;
}

// Bound propagation over row r of A: a'x <= ru bounds a_j*x_j by ru minus the
// smallest value of the other terms (and similarly for rl). Returns the number
// of tightened bounds or DAQP_EXIT_INFEASIBLE
static int presolve_propagate(DAQPProblem* qp, const int r){
    int j, flag, nmin = 0, nmax = 0, changed = 0;
    const int n = qp->n, id = qp->ms+r;
    const c_float* a = qp->A+(size_t)r*n;
    const c_float ru = qp->bupper[id], rl = qp->blower[id];
    c_float amin = 0, amax = 0, scale = 0, span = 0, lo, hi, rest;
    if(qp->sense[id] & DAQP_SOFT) return 0;
    for(j = 0; j < n; j++){ // Smallest and largest value of a'x (and number of infinite terms)
        if(a[j] == 0) continue;
        lo = a[j] > 0 ? presolve_lower(qp,j) : presolve_upper(qp,j);
        hi = a[j] > 0 ? presolve_upper(qp,j) : presolve_lower(qp,j);
        if(DAQP_PRESOLVE_IS_INF(lo)) nmin++; else amin += a[j]*lo;
        if(DAQP_PRESOLVE_IS_INF(hi)) nmax++; else amax += a[j]*hi;
        if(fabs(a[j]) > scale) scale = fabs(a[j]);
        if(a[j]*(hi-lo) > span) span = a[j]*(hi-lo); // Largest range of a term
    }
    if(nmin == 0 && ru < DAQP_INF && amin > ru + DAQP_PRESOLVE_FEAS_TOL*scale*(1+fabs(ru)))
        return DAQP_EXIT_INFEASIBLE;
    if(nmax == 0 && rl > -DAQP_INF && amax < rl - DAQP_PRESOLVE_FEAS_TOL*scale*(1+fabs(rl)))
        return DAQP_EXIT_INFEASIBLE;
    // A bound can only be tightened if the slack of the row is less than the range of
    // its term (the slack is infinite if more than one term is unbounded)
    if((ru >= DAQP_INF || nmin > 1 || (nmin == 0 && ru-amin >= span)) &&
            (rl <= -DAQP_INF || nmax > 1 || (nmax == 0 && amax-rl >= span))) return 0;
    for(j = 0; j < n; j++){
        // A bound from a negligible coefficient would mostly be amplified rounding
        if(fabs(a[j]) <= DAQP_PRESOLVE_MARGIN*scale) continue;
        lo = a[j] > 0 ? presolve_lower(qp,j) : presolve_upper(qp,j);
        hi = a[j] > 0 ? presolve_upper(qp,j) : presolve_lower(qp,j);
        if(ru < DAQP_INF && (nmin == 0 || (nmin == 1 && DAQP_PRESOLVE_IS_INF(lo)))){
            rest = nmin == 0 ? amin - a[j]*lo : amin;
            if((flag = presolve_bound(qp,j,(ru-rest)/a[j],a[j] > 0)) < 0) return flag;
            changed += flag;
        }
        if(rl > -DAQP_INF && (nmax == 0 || (nmax == 1 && DAQP_PRESOLVE_IS_INF(hi)))){
            rest = nmax == 0 ? amax - a[j]*hi : amax;
            if((flag = presolve_bound(qp,j,(rl-rest)/a[j],a[j] < 0)) < 0) return flag;
            changed += flag;
        }
    }
    return changed;
}

// Coefficient tightening of the binaries in a one-sided row r of A. With the
// row as rest + c*b <= beta, where x_k = l + (u-l)*b and b in {0,1}, the
// constraint is redundant for b = 1 (c < 0) or b = 0 (c > 0) if the largest
// value R of rest is small enough, and c (and beta) can then be reduced until
// the constraint holds with equality at R. R is formed from the bounds of bq (the
// problem before the bounds were propagated), since a row that holds with equality
// at the propagated bounds makes them degenerate (which can lead to cycling).
// Returns the number of changes
static int presolve_coefficients(DAQPProblem* qp, const DAQPProblem* bq, const int r){
    int j, k, changed = 0;
    const int n = qp->n, id = qp->ms+r;
    c_float* a = qp->A+(size_t)r*n;
    const int upper = qp->bupper[id] < DAQP_INF;
    const c_float s = upper ? 1 : -1; // The row is s*a'x <= s*bound
    c_float c, beta, R, v, l, w, tol, d;
    if((qp->sense[id] & (DAQP_SOFT+DAQP_IMMUTABLE+DAQP_BINARY)) ||
            upper == (qp->blower[id] > -DAQP_INF)) return 0; // Not one-sided
    for(k = 0; k < qp->ms && k < n; k++){
        if(a[k] == 0 || !(qp->sense[k] & DAQP_BINARY) ||
                (qp->sense[k] & (DAQP_SOFT+DAQP_IMMUTABLE))) continue;
        l = qp->blower[k];
        w = qp->bupper[k]-l;
        if(w <= 0) continue;
        c = s*a[k]*w;
        beta = (upper ? qp->bupper[id] : -qp->blower[id]) - s*a[k]*l;
        for(j = 0, R = 0; j < n; j++){ // Largest value of the other terms
            if(j == k || a[j] == 0) continue;
            v = s*a[j] > 0 ? presolve_upper(bq,j) : presolve_lower(bq,j);
            if(DAQP_PRESOLVE_IS_INF(v)) break;
            R += s*a[j]*v;
        }
        if(j < n) continue;
        tol = DAQP_PRESOLVE_MARGIN*(1+fabs(beta)+fabs(R));
        if(c < 0 && R - beta < -c - tol) c = beta - R; // Redundant for b = 1
        else if(c > 0 && R < beta - tol){              // Redundant for b = 0
            d = beta - R < c ? beta - R : c;
            c -= d;
            beta -= d;
        }
        else continue;
        if(fabs(c) <= tol) c = 0; // Cancelled up to rounding
        a[k] = s*c/w;
        beta += s*a[k]*l;
        if(upper) qp->bupper[id] = beta;
        else qp->blower[id] = -beta;
        changed++;
    }
    return changed;
}

static int presolve_problem(DAQPProblem* qp, const DAQPProblem* bq){
    int i, r, flag, total = 0, nb = 0;
    const int mA = qp->m-qp->ms;
    if(qp->sense == NULL || qp->A == NULL || mA <= 0 || qp->nh > 1 || qp->problem_type == 1)
        return 0;
    for(i = 0; i < qp->m; i++) if(qp->sense[i] & DAQP_BINARY) nb++;
    if(nb == 0) return 0; // Only mixed-integer problems are presolved

    // One pass of bound propagation, then one of coefficient tightening (further
    // passes mostly give small tightenings, which cost more than they save)
    for(r = 0; r < mA; r++){
        if((flag = presolve_propagate(qp,r)) < 0) return flag;
        total += flag;
    }
    for(r = 0; r < mA; r++) total += presolve_coefficients(qp,bq,r);
    return total;
}

/* --------------------------------------------------------------------------
 * daqp_presolve  --  presolve for the LDP (settings->presolve)
 *
 * qp is presolved into a copy (in scratch allocated once, since the dimensions
 * are fixed at setup), which the LDP is formed from instead of qp. qp itself is
 * left untouched (and the solution is checked against it). Returns the copy, or
 * qp if the presolve changed nothing or found the problem infeasible (BnB then
 * finds it infeasible).
 * --------------------------------------------------------------------------*/
DAQPProblem* daqp_presolve(DAQPWorkspace* work, DAQPProblem* qp){
    const int n = qp->n, m = qp->m, mA = m-qp->ms;
    DAQPProblem* cp = work->bnb->ps_qp;
    if(mA <= 0 || qp->A == NULL) return qp;
    if(cp == NULL)
        cp = work->bnb->ps_qp = malloc(sizeof(DAQPProblem)+((size_t)mA*n+2*m)*sizeof(c_float));
    *cp = *qp; // H, f and the sense are shared
    cp->A = (c_float*)(cp+1);
    cp->bupper = cp->A+(size_t)mA*n;
    cp->blower = cp->bupper+m;
    memcpy(cp->A,qp->A,(size_t)mA*n*sizeof(c_float));
    memcpy(cp->bupper,qp->bupper,m*sizeof(c_float));
    memcpy(cp->blower,qp->blower,m*sizeof(c_float));
    return presolve_problem(cp,qp) > 0 ? cp : qp;
}
