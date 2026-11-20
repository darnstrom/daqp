#ifdef DAQP_MULTISTAGE // (requires BLASFEO)
#include "multistage.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "auxiliary.h"
#include "daqp_prox.h"
#include "bnb.h"
#include "api.h"
#include <blasfeo.h>
#ifndef D_PS
#define D_PS 4
#endif

/*
 * Structure-exploiting multistage backend: DAQP's dual active-set method on the LDP in
 * v, with M accessed implicitly (DAQPMOperator) through Riccati recursions in
 * BLASFEO. The method (elimination of equalities, normalization, proximal-point
 * iterations, branch and bound, x_0) is described in docs/docs/multistage.md.
 *
 * Based on DAOCP: A. Zaupa, S. Erickson, M. Johansson, "DAOCP: a dual active set
 * solver for optimal control problems", arXiv:2609.33451, 2026
 * Coupling constraints follow PIQP (https://github.com/PREDICT-EPFL/piqp).
 */

// A stage (with the stage of x_0 first if x0free)
typedef struct{
    int nx, nu, nb, ng, nc; // nc: coupling constraints with stage t+1
    int nn; // Free inputs (nu minus those fixed by eliminated equalities)
    int cv; // Offset of v_t in v
    int c0; // Index of the first constraint of the stage
    int *idxb;
    struct blasfeo_dmat *BAbt, *RSQ, *DCt, *Et, *Ft; // The user's
    int ngb, *gblk; // Blocks of columns of DCt: [column, columns, first row, rows] (multistage_gblocks)
    struct blasfeo_dmat *BAr, *RSr; // Of the free inputs: BAbt, RSQ or T'BAbt, T'RSQ T
    c_float *lb, *ub, *lg, *ug, *lc, *uc;
    struct blasfeo_dvec rq, b; // [r; q] and b (copies)
    int has_b; // b is not zero
    struct blasfeo_dmat D; // [Luu; Lxu] and Lxx (sqrt) or P_t (classical)
    struct blasfeo_dvec z, dz; // [u; x] of z_lqr and dz
    struct blasfeo_dvec yp; // [y; p_t] of the backward recursion (kv_t = -y; p_t from nn)
    struct blasfeo_dvec Pb; // P_{t+1} b (kept until b or the factorization changes)
    c_float bPb; // 0.5 b'P_{t+1} b (with Pb)
    int Pb_valid;
    struct blasfeo_dvec zc; // [u; x] of the prox center
    c_float cacc, *h; // c_t of the backward recursion, h_t (elimination)
    int infb; // Inconsistent eliminated equalities at stages t..N
    // Elimination of equality constraints (see multistage_eq_stage)
    int ncand, *cand; // Candidates (equal bounds at setup)
    int q, r, rh, ke; // Stacked equalities, rank in u_t, rows of H_t, consistency conditions
    c_float *Yg, *Z, *Z2, *om; // Row major (Yg: R11^-T, r x r); om: multipliers of the stacked equalities
    struct blasfeo_dmat Ht; // H_t' (nx x rh)
    int *piv, *piv2; // Pivots of the QR of Gu' and of the dependent rows (Z = R11\R12 of each)
    struct blasfeo_dmat Qd, Kt, BAm; // Q, K1' (T = [Q2 -Q1 K1; 0 I]), reduced BAbt
    struct blasfeo_dvec u0; // Inputs fixed by the equalities
}Stage;

typedef struct{
    DAQPMultistage *user;
    DAQPWorkspace *work;
    int N, m, nv, nvmax, npi; // Of the problem that is solved (npi: number of costates)
    int nsmax; // Largest number of soft constraints (the core workspace is for nvmax + nsmax)
    int nxmax, numax, ngmax; // Largest stage dimensions
    Stage *st;
    int x0free; // x_0 is a variable (stage 0 is its stage)
    c_float *x0; // x_0 (if fixed)
    int *x0src; // Constraint that fixes each component of x_0 (if fixed)
    int nx0c; // Number of constraints on x_0 alone (if fixed)
    int *ct; // Stage of each constraint
    int *eqf; // DAQP_MULTISTAGE_EQ_* flags
    int neq_now; // Constraints with DAQP_MULTISTAGE_EQ_NOW (-1: not counted yet)
    int neqc; // Number of candidates for elimination
    int no_elim; // Equalities kept in the working set (elimination failed)
    int ill_cond; // Latest elimination ill-conditioned
    int elim_bad; // The dual pass found the elimination inaccurate: equalities in the working set until they change
    int elim_check; // The latest factorization eliminated inputs and no solution of it has been checked yet
    int infeasible; // The eliminated equalities are inconsistent
    c_float *eqval; // Values of the candidates at the latest update
    struct blasfeo_dmat B0, R0; // The stage of x_0: B = I and R = 0
    struct blasfeo_dmat sT, sH, sAt, sCw, sG, sC0, sC1, sC2, sV, sK, sX; // Scratch
    struct blasfeo_dmat sY, sV2, sB, sF; // Scratch of the elimination
    struct blasfeo_dmat sQv, sQw; // Scratch of multistage_householder_q
    struct blasfeo_dmat sRS; // Reduced RSQ of the stage being factored
    struct blasfeo_dvec ws, wg, wc0, wc1, wr, wt, wb;
    struct blasfeo_dmat sSig[2], *Sig; // Sig_t of the stage being normalized (rolling)
    int bw_from; // Redo the backward recursion from this stage (-1: not at all)
    c_float *sGe, *sW2, *sTau, *sg, *sy; // Scratch of the elimination
    c_float *sng; // Multipliers of the eliminated general constraints of a stage (multistage_duals)
    c_float *scale; // 1/|M_i|
    char *zrow; // Row i of M is zero (latest normalization)
    int norm_stale; // Recompute the normalization at the next factorization
    int ric_sqrt; // Square-root recursion (P_t = Lxx Lxx' in D)
    c_float eps; // Prox regularization of the Hessian (0: none)
    // Branch and bound: the bound each binary constraint is at in the latest solution (1: lower, 2: upper, 0: none),
    // the other bound while a candidate fixes it at that one, and the best candidate (multistage_bnb_candidate): its
    // solution v, objective (in the scale of work->fval) and soft slack
    unsigned char *bside;
    c_float *bsave, *vcand, cand_fval, cand_slack;
    c_float *pic, *lamv; // Costates and multipliers of the latest solve
    int pi_valid, lam_valid;
    int lam_done; // lamv is complete (else: the eliminated equalities and the bounds fixing x_0 are missing)
    c_float *vlast; // v that dz was last computed for
    int dz_valid;
    c_float Jlqr; // Objective at z_lqr
    // Cache of rows of M (length nvmax)
    c_float *rows;
    int *slot_id, *id_slot, *rows_len, nslots;
    const int *pin; int npin; // Constraints whose rows multistage_rows_block forms (kept in the cache meanwhile)
    int *ids, *smax; // Rows formed in blocks (multistage_rows_active)
    // Block activation of the warm start (multistage_activate): DAQP_MULTISTAGE_ACT_COLS columns of the rows of the
    // working set, and the reciprocals of the diagonal of the Cholesky factor
    struct blasfeo_dmat sAb;
    c_float *sic;
    c_float *dupper, *dlower, *Mu; int *sense; // The vectors of the LDP (of work)
    void *mem; // The block of all memory of the multistage problem
    int own_mem; // Allocated by setup_daqp_multistage (else given by the caller)
    DAQPMOperator op;
    // Warm start (multistage_ws_record): final working sets of the latest solves (newest first), as the state of each
    // constraint (wsc: its DAQP_MULTISTAGE_WS_BITS, 0 if inactive) and the list of the active ones (wsh, wsl of them)
    int ws_n, ws_alt, ws_two;
    // DAQP_MULTISTAGE_UPDATE_SHIFT: the constraint of the next stage whose state each constraint takes (-1: none; wsdst
    // is the inverse), and whether the shifted working set fit the latest solution better
    int *wssrc, *wsdst, ws_shift;
    unsigned char *wsc[3];
    int *wsh[3], wsl[3]; // (at most nvmax + nsmax + 1 constraints, as DAQP's WS)
}DAQPMultistageWork;

#define DAQP_MULTISTAGE_NB(prob,t) ((prob)->nb != NULL ? (prob)->nb[t] : 0)
#define DAQP_MULTISTAGE_NG(prob,t) ((prob)->ng != NULL ? (prob)->ng[t] : 0)
#define DAQP_MULTISTAGE_NC(prob,t) (((prob)->nc != NULL && (t) < (prob)->N) ? (prob)->nc[t] : 0)
#define DAQP_MULTISTAGE_NU(prob,t) ((t) < (prob)->N ? (prob)->nu[t] : 0)

// Kinds of constraints
#define DAQP_MULTISTAGE_BOUND 0
#define DAQP_MULTISTAGE_GENERAL 1
#define DAQP_MULTISTAGE_COUPLING 2

// Columns of the rows of the working set that the block activation packs at a time (multistage_activate)
#define DAQP_MULTISTAGE_ACT_COLS 64
// Smallest working set activated at once (fewer constraints are added one by one, which is then as fast)
#define DAQP_MULTISTAGE_ACT_MIN 16
// Rows stored one after the other start at multiples of this many entries (aligned loads)
#define DAQP_MULTISTAGE_ROW_ALIGN 4
#define DAQP_MULTISTAGE_STRIDE(len) (((len) + DAQP_MULTISTAGE_ROW_ALIGN-1) & ~(DAQP_MULTISTAGE_ROW_ALIGN-1))
// Rows of M formed together (multistage_rows_block)
#define DAQP_MULTISTAGE_BLOCK 32

// Stages with at most DAQP_MULTISTAGE_SMALL variables use loops instead of BLASFEO calls
#ifndef DAQP_MULTISTAGE_SMALL
#define DAQP_MULTISTAGE_SMALL 6
#endif

// Largest drift of the norm of a normalized row before it is renormalized
#define DAQP_MULTISTAGE_NORM_DRIFT 2

// Relative tolerance for the numerical rank of the equality constraints
#define DAQP_MULTISTAGE_RANK_TOL 1e-10
// Smallest relative input pivot of the elimination (else ill-conditioned)
#define DAQP_MULTISTAGE_ELIM_TOL 1e-4
// Largest stationarity residual in the free inputs of an elimination, at a stage and relative to its terms (checked
// in the dual pass; else the equalities are kept in the working set from the next factorization)
#ifndef DAQP_MULTISTAGE_ELIM_CHECK
#define DAQP_MULTISTAGE_ELIM_CHECK 1e-5
#endif

// Flags of equality constraints
#define DAQP_MULTISTAGE_EQ_CAND 1 // Equal bounds at setup (can be eliminated)
#define DAQP_MULTISTAGE_EQ_NOW 2 // Equal bounds in the latest update (marked as an equality in sense)
#define DAQP_MULTISTAGE_EQ_ELIM 4 // Eliminated
#define DAQP_MULTISTAGE_EQ_X0 8 // Bound on a fixed x_0 alone (not in the LDP)

int daqp_multistage_nz(const DAQPMultistage *prob){
    int t, nz = 0;
    for(t = 0; t <= prob->N; t++) nz += prob->nx[t] + DAQP_MULTISTAGE_NU(prob,t);
    return nz;
}

int daqp_multistage_m(const DAQPMultistage *prob){
    int t, m = 0;
    for(t = 0; t <= prob->N; t++)
        m += DAQP_MULTISTAGE_NB(prob,t) + DAQP_MULTISTAGE_NG(prob,t) + DAQP_MULTISTAGE_NC(prob,t);
    return m;
}

/* BLASFEO */

#define DAQP_MULTISTAGE_MAX_GBLOCKS 4

// Consecutive columns of DCt (general constraints) on inputs only, on states only, or on both, so that
// products with DCt skip zero blocks. Too many blocks: one block of all
static void multistage_gblocks(Stage *S){
    const int nu = S->nu, nx = S->nx;
    int k, i, prev = -1;
    S->ngb = 0;
    for(k = 0; k < S->ng; k++){
        int hu = 0, hx = 0, cls;
        for(i = 0; i < nu && !hu; i++) hu = BLASFEO_DMATEL(S->DCt, i, k) != 0;
        for(i = nu; i < nu+nx && !hx; i++) hx = BLASFEO_DMATEL(S->DCt, i, k) != 0;
        cls = hu && hx ? 2 : (hx ? 1 : 0); // (a zero column is on inputs)
        if(cls == prev) S->gblk[4*(S->ngb-1)+1]++;
        else{
            int *b = S->gblk + 4*S->ngb++;
            b[0] = k; b[1] = 1; b[2] = cls == 1 ? nu : 0; b[3] = cls == 0 ? nu : (cls == 1 ? nx : nu+nx);
            prev = cls;
        }
    }
    if(S->ngb > DAQP_MULTISTAGE_MAX_GBLOCKS){
        S->ngb = 1;
        S->gblk[0] = 0; S->gblk[1] = S->ng; S->gblk[2] = 0; S->gblk[3] = nu+nx;
    }
}

// z[zi..zi+ng) = DCt'x (x: [u; x] of the stage)
static void multistage_dct_t(Stage *S, struct blasfeo_dvec *x, struct blasfeo_dvec *z, const int zi){
    int j;
    for(j = 0; j < S->ngb; j++){
        const int *b = S->gblk + 4*j;
        if(b[3] > 0) blasfeo_dgemv_t(b[3], b[1], 1.0, S->DCt, b[2], b[0], x, b[2], 0.0, z, zi+b[0], z, zi+b[0]);
        else blasfeo_dvecse(b[1], 0.0, z, zi+b[0]);
    }
}

// g += DCt y[yi..yi+ng) in rows r0..r0+nr of g (the rows of [u; x] of the stage; skipped blocks are zero there)
static void multistage_dct_n(Stage *S, struct blasfeo_dvec *y, const int yi, const int r0, const int nr,
        struct blasfeo_dvec *g){
    int j;
    for(j = 0; j < S->ngb; j++){
        const int *b = S->gblk + 4*j;
        const int lo = b[2] > r0 ? b[2] : r0, hi = b[2]+b[3] < r0+nr ? b[2]+b[3] : r0+nr;
        if(hi > lo) blasfeo_dgemv_n(hi-lo, b[1], 1.0, S->DCt, lo, b[0], y, yi+b[0], 1.0, g, lo, g, lo);
    }
}

// x[0:n] <-- L(0:n,0:n)'\(x[0:n] - L(n:m,0:n)' x[n:m])
// (the C kernel of trsv_ltn_mn, used by GENERIC and ARM targets, is wrong for m = n+1, n%4 = 3)
static void bf_trsv_ltn_mn(const int m, const int n, struct blasfeo_dmat *L, struct blasfeo_dvec *x){
    if(m == n+1 && n%4 == 3){
        blasfeo_dgemv_t(1, n, -1.0, L, n, 0, x, n, 1.0, x, 0, x, 0);
        blasfeo_dtrsv_ltn(n, L, 0, 0, x, 0, x, 0);
    }
    else blasfeo_dtrsv_ltn_mn(m, n, L, 0, 0, x, 0, x, 0);
}

// Memory: one aligned block (only the size is counted if base is NULL). Only the BLASFEO matrices need the
// alignment of the block (DAQP_MULTISTAGE_ALIGN); vectors (also created within plain arrays) and plain arrays less
#define DAQP_MULTISTAGE_ALIGN 64
#define DAQP_MULTISTAGE_ALIGN_VEC 32
#define DAQP_MULTISTAGE_ALIGN_PLAIN 8
typedef struct{ char *base; size_t n; }MultistageMem;
static void* multistage_take_aligned(MultistageMem *M, const size_t bytes, const size_t align){
    void *p;
    M->n = (M->n + align-1) & ~(align-1);
    p = M->base != NULL ? M->base + M->n : NULL;
    M->n += bytes;
    return p;
}
static void* multistage_take(MultistageMem *M, const size_t bytes){
    return multistage_take_aligned(M, bytes, DAQP_MULTISTAGE_ALIGN_PLAIN);
}
static void multistage_new_mat(MultistageMem *M, struct blasfeo_dmat *A, int m, int n){
    void *p;
    if(m < 1) m = 1;
    if(n < 1) n = 1;
    if((p = multistage_take_aligned(M, blasfeo_memsize_dmat(m, n), DAQP_MULTISTAGE_ALIGN)) != NULL)
        blasfeo_create_dmat(m, n, A, p);
}
static void multistage_new_vec(MultistageMem *M, struct blasfeo_dvec *v, int m){
    void *p;
    if(m < 1) m = 1;
    if((p = multistage_take_aligned(M, blasfeo_memsize_dvec(m), DAQP_MULTISTAGE_ALIGN_VEC)) != NULL)
        blasfeo_create_dvec(m, v, p);
}

/* Constraints */

// Kind of constraint id, its index within that kind, and its stage (*t)
static int multistage_decode(const DAQPMultistageWork *ow, const int id, int *t, int *k){
    const Stage *S = ow->st + (*t = ow->ct[id]);
    *k = id - S->c0;
    if(*k < S->nb) return DAQP_MULTISTAGE_BOUND;
    *k -= S->nb;
    if(*k < S->ng) return DAQP_MULTISTAGE_GENERAL;
    *k -= S->ng;
    return DAQP_MULTISTAGE_COUPLING;
}

// Last stage with coefficients of constraint id (that of the next stage for a coupling constraint)
static int multistage_last_stage(const DAQPMultistageWork *ow, const int id){
    int t, k;
    return multistage_decode(ow, id, &t, &k) == DAQP_MULTISTAGE_COUPLING ? t+1 : t;
}

// Coefficients [cu; cx] of constraint id on stage s; returns whether it has u terms
static int multistage_coeffs(DAQPMultistageWork *ow, const int id, const int s, struct blasfeo_dvec *c){
    const Stage *S = ow->st+s;
    const int n = S->nu + S->nx;
    struct blasfeo_dmat *A;
    int i, t, k;
    const int kind = multistage_decode(ow, id, &t, &k);
    blasfeo_dvecse(n, 0.0, c, 0);
    if(kind == DAQP_MULTISTAGE_BOUND){
        if(s != t) return 0;
        c->pa[S->idxb[k]] = 1; // idxb indexes [u; x]
        return S->idxb[k] < S->nu;
    }
    if(kind == DAQP_MULTISTAGE_GENERAL) A = s == t ? S->DCt : NULL;
    else A = s == t ? ow->st[t].Et : (s == t+1 ? ow->st[t].Ft : NULL);
    if(A == NULL || n == 0) return 0;
    blasfeo_dcolex(n, A, 0, k, c, 0); // Column of [D'; C'], [Eu'; Ex'] or [Fu'; Fx']
    for(i = 0; i < S->nu; i++) if(c->pa[i] != 0) return 1;
    return 0;
}

// Bounds into work->dlower/dupper (shifted to the LDP by multistage_lqr_forward)
static void multistage_fetch_bounds(const DAQPMultistageWork *ow){
    c_float *lb = ow->work->dlower, *ub = ow->work->dupper;
    int t;
#define DAQP_MULTISTAGE_FETCH(n, l, u) do{ if((n) > 0){ \
        memcpy(lb+c, l, (n)*sizeof(c_float)); memcpy(ub+c, u, (n)*sizeof(c_float)); c += (n); } }while(0)
    for(t = 0; t <= ow->N; t++){
        const Stage *S = ow->st+t;
        int c = S->c0;
        DAQP_MULTISTAGE_FETCH(S->nb, S->lb, S->ub);
        DAQP_MULTISTAGE_FETCH(S->ng, S->lg, S->ug);
        DAQP_MULTISTAGE_FETCH(S->nc, S->lc, S->uc);
    }
#undef DAQP_MULTISTAGE_FETCH
}

// Constraint values of z_lqr, or of dz (dz nonzero)
static void multistage_eval(DAQPMultistageWork *ow, const int dz, c_float *val){
    int t, i;
    for(t = 0; t <= ow->N; t++){
        Stage *S = ow->st+t;
        const int n = S->nu + S->nx;
        struct blasfeo_dvec *ux = dz ? &S->dz : &S->z, sv;
        blasfeo_create_dvec(S->nb+S->ng+S->nc, &sv, val + S->c0);
        for(i = 0; i < S->nb; i++) sv.pa[i] = ux->pa[S->idxb[i]];
        multistage_dct_t(S, ux, &sv, S->nb);
        if(S->nc > 0){
            Stage *Sn = S+1;
            const int i0 = S->nb + S->ng;
            blasfeo_dgemv_t(n, S->nc, 1.0, S->Et, 0, 0, ux, 0, 0.0, &sv, i0, &sv, i0);
            blasfeo_dgemv_t(Sn->nu+Sn->nx, S->nc, 1.0, S->Ft, 0, 0, dz ? &Sn->dz : &Sn->z, 0, 1.0, &sv, i0, &sv, i0);
        }
    }
}

/* Elimination of equality constraints */

// Householder QR with column pivoting of W (m x n, column major; swaps also apply to rows
// m..ms-1). Stops at a column norm below DAQP_MULTISTAGE_RANK_TOL*max(ref, largest norm). Returns the rank
static int multistage_qrp(c_float *W, const int ld, const int m, const int n, const int ms,
        int *piv, c_float *tau, const c_float ref){
    int i, j, k;
    const int kmax = m < n ? m : n;
    c_float r0 = 0;
    for(j = 0; j < n; j++) piv[j] = j;
    for(k = 0; k < kmax; k++){
        // Pivot: the remaining column with the largest norm
        int best = k;
        c_float bn = -1, nrm, x0, beta, scal, *wk = W + (size_t)k*ld;
        for(j = k; j < n; j++){
            const c_float *wj = W + (size_t)j*ld;
            c_float sq = 0;
            for(i = k; i < m; i++) sq += wj[i]*wj[i];
            if(sq > bn){ bn = sq; best = j; }
        }
        if(best != k){
            c_float *wb = W + (size_t)best*ld;
            for(i = 0; i < ms; i++){ const c_float tmp = wk[i]; wk[i] = wb[i]; wb[i] = tmp; }
            j = piv[k]; piv[k] = piv[best]; piv[best] = j;
        }
        nrm = sqrt(bn);
        if(k == 0) r0 = nrm > ref ? nrm : ref;
        if(nrm <= DAQP_MULTISTAGE_RANK_TOL*r0 || nrm == 0) return k;
        // Householder reflection that maps column k to beta*e_k
        x0 = wk[k];
        beta = x0 >= 0 ? -nrm : nrm;
        tau[k] = (beta-x0)/beta;
        scal = 1/(x0-beta);
        for(i = k+1; i < m; i++) wk[i] *= scal;
        wk[k] = beta;
        for(j = k+1; j < n; j++){
            c_float *wj = W + (size_t)j*ld, sq = wj[k];
            for(i = k+1; i < m; i++) sq += wk[i]*wj[i];
            sq *= tau[k];
            wj[k] -= sq;
            for(i = k+1; i < m; i++) wj[i] -= sq*wk[i];
        }
    }
    return kmax;
}

// Z = R11\R12 (r x (q-r), row major): column piv[r+j] = sum_k Z[k,j] column piv[k]
static void multistage_qrp_z(const c_float *W, const int ld, const int q, const int r, c_float *Z){
    int i, j, l;
    for(j = 0; j < q-r; j++)
        for(i = r-1; i >= 0; i--){
            c_float sq = W[(size_t)(r+j)*ld+i];
            for(l = i+1; l < r; l++) sq -= W[(size_t)l*ld+i]*Z[l*(q-r)+j];
            Z[i*(q-r)+j] = sq/W[(size_t)i*ld+i];
        }
}

// Q = H_0 ... H_{r-1} = I - V T V' in Qd, from the reflections of multistage_qrp in G (v_k = [e_k; G[k+1:nu, k]];
// T upper triangular, as in LAPACK's dlarft). Scratch: Vs (nu x r), Ts (r x r), Ms (nu x r)
static void multistage_householder_q(const c_float *G, const int ld, const int nu, const int r, const c_float *tau,
        struct blasfeo_dmat *Vs, struct blasfeo_dmat *Ts, struct blasfeo_dmat *Ms, struct blasfeo_dmat *Qd){
    int i, j, l;
    blasfeo_pack_dmat(nu, r, (c_float*)G, ld, Vs, 0, 0);
    for(j = 0; j < r; j++){
        for(i = 0; i < j; i++) BLASFEO_DMATEL(Vs, i, j) = 0;
        BLASFEO_DMATEL(Vs, j, j) = 1;
    }
    // V'V in Ts; its strictly lower part (z = V(:, 0:i)'v_i in row i) is read while T is formed above it
    blasfeo_dgemm_tn(r, r, nu, 1.0, Vs, 0, 0, Vs, 0, 0, 0.0, Ts, 0, 0, Ts, 0, 0);
    for(i = 0; i < r; i++){ // T(0:i, i) = -tau_i T(0:i, 0:i) z
        for(j = 0; j < i; j++){
            c_float s = 0;
            for(l = j; l < i; l++) s += BLASFEO_DMATEL(Ts, j, l)*BLASFEO_DMATEL(Ts, i, l);
            BLASFEO_DMATEL(Ts, j, i) = -tau[i]*s;
        }
        BLASFEO_DMATEL(Ts, i, i) = tau[i];
    }
    for(j = 0; j < r; j++) for(i = j+1; i < r; i++) BLASFEO_DMATEL(Ts, i, j) = 0;
    // Q = I - (V T) V'
    blasfeo_dgemm_nn(nu, r, r, 1.0, Vs, 0, 0, Ts, 0, 0, 0.0, Ms, 0, 0, Ms, 0, 0);
    blasfeo_dgese(nu, nu, 0.0, Qd, 0, 0);
    blasfeo_ddiare(nu, 1.0, Qd, 0, 0);
    blasfeo_dgemm_nt(nu, nu, r, -1.0, Ms, 0, 0, Vs, 0, 0, 1.0, Qd, 0, 0, Qd, 0, 0);
}

// Reduced stage t, T = [Q2 -Q1 K1; 0 I]: with X = RSQ_{:,u} Q, Hq = Q'X_u, V = Xx1 - 0.5 K1'Hq11,
// RSr = [Hq22 .; Xx2 - K1'Hq12  Hxx - V K1 - K1'V'] and BAr = [Bq2; Ax - K1'Bq1], Bq = Q'B
static void multistage_reduce_stage(DAQPMultistageWork *ow, const int t){
    Stage *S = ow->st+t;
    const int nx = S->nx, nu = S->nu, n = nu+nx, r = S->r, nn = S->nn, nx1 = t < ow->N ? S[1].nx : 0;
    blasfeo_dgecp(n, nu, S->RSQ, 0, 0, &ow->sH, 0, 0);
    blasfeo_dtrtr_l(nu, &ow->sH, 0, 0, &ow->sH, 0, 0);
    blasfeo_dgemm_nn(n, nu, nu, 1.0, &ow->sH, 0, 0, &S->Qd, 0, 0, 0.0, &ow->sX, 0, 0, &ow->sX, 0, 0);
    blasfeo_dgemm_tn(nu, nu, nu, 1.0, &S->Qd, 0, 0, &ow->sX, 0, 0, 0.0, &ow->sY, 0, 0, &ow->sY, 0, 0);
    blasfeo_dgecp(nn, nn, &ow->sY, r, r, &ow->sRS, 0, 0);
    if(nx > 0){
        if(nn > 0){
            blasfeo_dgecp(nx, nn, &ow->sX, nu, r, &ow->sRS, nn, 0);
            blasfeo_dgemm_nn(nx, nn, r, -1.0, &S->Kt, 0, 0, &ow->sY, 0, r, 1.0, &ow->sRS, nn, 0, &ow->sRS, nn, 0);
        }
        blasfeo_dgemm_nn(nx, r, r, -0.5, &S->Kt, 0, 0, &ow->sY, 0, 0, 1.0, &ow->sX, nu, 0, &ow->sV2, 0, 0);
        blasfeo_dsyr2k_ln(nx, r, -1.0, &ow->sV2, 0, 0, &S->Kt, 0, 0, 1.0, S->RSQ, nu, nu, &ow->sRS, nn, nn);
    }
    if(nx1 > 0){
        blasfeo_dgemm_tn(nu, nx1, nu, 1.0, &S->Qd, 0, 0, S->BAbt, 0, 0, 0.0, &ow->sB, 0, 0, &ow->sB, 0, 0);
        blasfeo_dgecp(nn, nx1, &ow->sB, r, 0, &S->BAm, 0, 0);
        blasfeo_dgemm_nn(nx, nx1, r, -1.0, &S->Kt, 0, 0, &ow->sB, 0, 0, 1.0, S->BAbt, nu, 0, &S->BAm, nn, 0);
    }
    S->RSr = &ow->sRS; S->BAr = &S->BAm;
}

// Dependent rows of stage t (after multistage_qrp, Z = R11\R12 in S->Z): with Craw g = g_dep - Z'g_ind (the
// rows of Craw: piv), (Craw Gx) x_t = Craw g; the independent ones are passed on (H_t, h_t = rows piv2 of Craw g),
// the others require Craw g to be consistent (Z2 = R11\R12 of the second QR)
static void multistage_eq_pass_on(DAQPMultistageWork *ow, const int t, const c_float gscale){
    Stage *S = ow->st+t;
    const int nx = S->nx, nu = S->nu, n = nu+nx, q = S->q, r = S->r, nr = q-r;
    c_float *G = ow->sGe, *W2 = ow->sW2, cscale = 0;
    int i, j, k;
    for(j = 0; j < nr; j++){
        c_float *hj = G + (size_t)(r+j)*n + nu, s1 = 1;
        for(k = 0; k < r; k++){
            const c_float rkj = G[(size_t)(r+j)*n+k], *xk = G + (size_t)k*n + nu;
            for(i = 0; i < nx; i++) hj[i] -= rkj*xk[i];
            s1 += fabs(S->Z[k*nr+j]);
        }
        memcpy(W2 + (size_t)j*nx, hj, nx*sizeof(c_float));
        if(s1 > cscale) cscale = s1;
    }
    S->rh = nx > 0 ? multistage_qrp(W2, nx, nx, nr, nx, S->piv2, ow->sTau, gscale*cscale) : 0;
    S->ke = nr-S->rh;
    if(S->rh > 0 && S->ke > 0) multistage_qrp_z(W2, nx, nr, S->rh, S->Z2);
    for(j = 0; nx == 0 && j < nr; j++) S->piv2[j] = j;
    for(k = 0; k < S->rh; k++)
        blasfeo_pack_dmat(nx, 1, G + (size_t)(r+S->piv2[k])*n + nu, nx, &S->Ht, 0, k);
}

// Eliminate the equalities of stage t (its own and H_{t+1}): QR of Gu' gives Q, Yg, H_t (passed
// on), Nc (consistency), and the reduced stage BAr = T'BAbt, RSr = T'RSQ T
static void multistage_eq_stage(DAQPMultistageWork *ow, const int t){
    Stage *S = ow->st+t;
    const int nx = S->nx, nu = S->nu, n = nu+nx, nx1 = t < ow->N ? S[1].nx : 0;
    const int mh = t < ow->N ? S[1].rh : 0;
    int i, j, k, q = 0, r;
    c_float gscale = 0, *G = ow->sGe;
    struct blasfeo_dvec *c = &ow->wc0;

    S->BAr = S->BAbt; S->RSr = S->RSQ; S->nn = nu;
    S->q = S->r = S->rh = S->ke = 0;
    for(i = 0; i < S->ncand; i++){
        if(!(ow->eqf[S->cand[i]] & DAQP_MULTISTAGE_EQ_ELIM)) continue;
        multistage_coeffs(ow, S->cand[i], t, c);
        memcpy(G + (size_t)q*n, c->pa, n*sizeof(c_float));
        q++;
    }
    if(mh > 0){ // BAbt H_{t+1}'
        blasfeo_dgemm_nn(n, mh, nx1, 1.0, S->BAbt, 0, 0, &S[1].Ht, 0, 0, 0.0, &ow->sF, 0, 0, &ow->sF, 0, 0);
        blasfeo_unpack_dmat(n, mh, &ow->sF, 0, 0, G + (size_t)q*n, n);
        q += mh;
    }
    if((S->q = q) == 0) return;
    // Scale: the largest row norm of [Gu Gx]
    for(j = 0; j < q; j++){
        c_float sq = 0;
        for(i = 0; i < n; i++) sq += G[(size_t)j*n+i]*G[(size_t)j*n+i];
        if(sq > gscale) gscale = sq;
    }
    gscale = sqrt(gscale);

    // The independent rows of Gu determine u_t in the directions Q1
    r = S->r = nu > 0 ? multistage_qrp(G, n, nu, q, n, S->piv, ow->sTau, gscale) : 0;
    for(j = 0; nu == 0 && j < q; j++) S->piv[j] = j;
    if(r > 0 && fabs(G[(size_t)(r-1)*n+r-1]) < DAQP_MULTISTAGE_ELIM_TOL*gscale) ow->ill_cond = 1; // (|R11| <= gscale)
    if(r > 0 && r < q) multistage_qrp_z(G, n, q, r, S->Z);
    if(r > 0){
        multistage_householder_q(G, n, nu, r, ow->sTau, &ow->sQv, &ow->sRS, &ow->sQw, &S->Qd);
        // Yg = R11^{-T}: mu = Yg (g - Gx x) on the rows piv[0..r)
        for(k = 0; k < r; k++)
            for(i = 0; i < r; i++){ // Column k of R11^{-T}
                c_float sq = (i == k) ? 1 : 0;
                for(j = 0; j < i; j++) sq -= G[(size_t)i*n+j]*S->Yg[j*r+k];
                S->Yg[i*r+k] = sq/G[(size_t)i*n+i];
            }
        // K1' = Gx_P1' R11^{-1} (in place of Gx_P1')
        for(k = 0; k < r; k++){
            c_float *xk = G + (size_t)k*n + nu;
            for(j = 0; j < k; j++){
                const c_float rjk = G[(size_t)k*n+j], *xj = G + (size_t)j*n + nu;
                for(i = 0; i < nx; i++) xk[i] -= rjk*xj[i];
            }
            for(i = 0; i < nx; i++) xk[i] /= G[(size_t)k*n+k];
        }
        blasfeo_pack_dmat(nx, r, G + nu, n, &S->Kt, 0, 0);
        S->nn = nu-r;
        multistage_reduce_stage(ow, t);
    }
    if(r < q) multistage_eq_pass_on(ow, t, gscale);
}

// Right-hand sides g of the equalities of stage t: u0_t = Q1 Yg g, h_t (rows piv2 of Craw g). Returns 1 if they are
// inconsistent
static int multistage_eq_vectors(DAQPMultistageWork *ow, const int t, const c_float *hn, c_float *h){
    Stage *S = ow->st+t;
    const int q = S->q, r = S->r, nu = S->nu, mh = t < ow->N ? S[1].rh : 0;
    c_float *g = ow->sg, *y = ow->sy, gmax = 0;
    int i, j = 0, k, inconsistent = 0;
    blasfeo_dvecse(nu, 0.0, &S->u0, 0);
    if(q == 0) return 0;
    for(i = 0; i < S->ncand; i++)
        if(ow->eqf[S->cand[i]] & DAQP_MULTISTAGE_EQ_ELIM) g[j++] = ow->work->dlower[S->cand[i]]; // Not shifted yet
    for(k = 0; k < mh; k++, j++){
        g[j] = hn[k];
        if(S->has_b) for(i = 0; i < S[1].nx; i++) g[j] -= BLASFEO_DMATEL(&S[1].Ht, i, k)*S->b.pa[i];
    }
    for(i = 0; i < q; i++) if(fabs(g[i]) > gmax) gmax = fabs(g[i]);
    for(j = 0; j < q-r; j++){ // y = Craw g
        c_float s = g[S->piv[r+j]];
        for(k = 0; k < r; k++) s -= S->Z[k*(q-r)+j]*g[S->piv[k]];
        y[j] = s;
    }
    for(j = 0; j < S->ke; j++){ // Consistency of the dependent rows
        c_float s = y[S->piv2[S->rh+j]];
        for(k = 0; k < S->rh; k++) s -= S->Z2[k*S->ke+j]*y[S->piv2[k]];
        if(fabs(s) > ow->work->settings->primal_tol*(1+gmax)) inconsistent = 1;
    }
    if(r > 0){
        for(i = 0; i < r; i++){ // tau = Yg g (rows piv[0..r) of g)
            c_float sq = 0;
            for(k = 0; k < r; k++) sq += S->Yg[i*r+k]*g[S->piv[k]];
            ow->sTau[i] = sq;
        }
        struct blasfeo_dvec tv;
        blasfeo_create_dvec(r, &tv, ow->sTau);
        blasfeo_dgemv_n(nu, r, 1.0, &S->Qd, 0, 0, &tv, 0, 0.0, &S->u0, 0, &S->u0, 0); // u0 = Q1 tau
    }
    for(k = 0; k < S->rh; k++) h[k] = y[S->piv2[k]];
    return inconsistent;
}

/* Riccati recursion (on the reduced stages) */

// Input pivots of the Cholesky factor L positive and, without prox, above zero_tol relative to
// the largest input diagonal of H (Hd, or from L) and of R (not P_t, which can be much larger)
static int multistage_pivots_ok(DAQPMultistageWork *ow, Stage *S, struct blasfeo_dmat *L,
        struct blasfeo_dmat *Hd, const int k){
    c_float hmax = 0, h;
    int i, j;
    for(i = 0; i < k; i++) if(!(BLASFEO_DMATEL(L, i, i) > 0)) return 0;
    if(ow->eps > 0) return 1; // H + eps I is positive definite
    for(i = 0; i < S->nu; i++) if((h = fabs(BLASFEO_DMATEL(S->RSQ, i, i))) > hmax) hmax = h;
    for(i = 0; i < k; i++){
        if(Hd != NULL) h = BLASFEO_DMATEL(Hd, i, i);
        else for(j = 0, h = 0; j <= i; j++) h += BLASFEO_DMATEL(L, i, j)*BLASFEO_DMATEL(L, i, j);
        if(h > hmax) hmax = h;
    }
    for(i = 0; i < k; i++){
        const c_float d = BLASFEO_DMATEL(L, i, i);
        if(!(d*d > ow->work->settings->zero_tol*hmax)) return 0;
    }
    return 1;
}

// Riccati recursion: square-root (ric_sqrt, as in HPIPM: D_t = chol(RSr + BAr P BAr'), P = Lxx Lxx') or classical
// (D_t = [Luu; Lxu] and P_t, with eps I added). Returns -1 if P_t is not positive definite (square-root) or an input
// pivot is negligible
static int multistage_factor(DAQPMultistageWork *ow){
    Stage *SN = ow->st + ow->N;
    int t, i;
    if(ow->neqc > 0) multistage_eq_stage(ow, ow->N);
    // P_N = Q_N (+ eps I)
    if(SN->nx > 0 && ow->ric_sqrt){
        blasfeo_dpotrf_l_mn(SN->nx, SN->nx, SN->RSQ, 0, 0, &SN->D, 0, 0);
        for(i = 0; i < SN->nx; i++) if(!(BLASFEO_DMATEL(&SN->D, i, i) > 0)) return -1;
    }
    else if(SN->nx > 0){
        blasfeo_dgecp(SN->nx, SN->nx, SN->RSQ, 0, 0, &SN->D, 0, 0);
        if(ow->eps > 0) blasfeo_ddiare(SN->nx, ow->eps, &SN->D, 0, 0);
        blasfeo_dtrtr_l(SN->nx, &SN->D, 0, 0, &SN->D, 0, 0);
    }
    for(t = ow->N-1; t >= 0; t--){
        Stage *S = ow->st+t;
        const int nx = S->nx, nx1 = S[1].nx;
        int nu, n;
        if(ow->neqc > 0) multistage_eq_stage(ow, t);
        nu = S->nn; n = nu+nx;
        if(ow->ric_sqrt){
            if(nx1 > 0){
                blasfeo_dtrmm_rlnn(n, nx1, 1.0, &S[1].D, S[1].nn, S[1].nn, S->BAr, 0, 0, &ow->sT, 0, 0);
                blasfeo_dsyrk_dpotrf_ln_mn(n, n, nx1, &ow->sT, 0, 0, &ow->sT, 0, 0, S->RSr, 0, 0, &S->D, 0, 0);
            }
            else blasfeo_dpotrf_l_mn(n, n, S->RSr, 0, 0, &S->D, 0, 0);
            if(!multistage_pivots_ok(ow, S, &S->D, NULL, nu)) return -1;
            for(i = nu; i < n; i++) if(!(BLASFEO_DMATEL(&S->D, i, i) > 0)) return -1;
            continue;
        }
        // H = RSr + BAr P BAr' (lower) (+ eps I)
        if(nx1 > 0){
            blasfeo_dgemm_nt(n, nx1, nx1, 1.0, S->BAr, 0, 0, &S[1].D, S[1].nn, S[1].nn, 0.0, &ow->sT, 0, 0, &ow->sT, 0, 0);
            blasfeo_dsyrk_ln(n, nx1, 1.0, &ow->sT, 0, 0, S->BAr, 0, 0, 1.0, S->RSr, 0, 0, &ow->sH, 0, 0);
        }
        else blasfeo_dgecp(n, n, S->RSr, 0, 0, &ow->sH, 0, 0);
        if(ow->eps > 0) blasfeo_ddiare(n, ow->eps, &ow->sH, 0, 0);
        // [Luu; Lxu]: Cholesky of the first nu columns
        if(nu > 0){
            blasfeo_dpotrf_l_mn(n, nu, &ow->sH, 0, 0, &S->D, 0, 0);
            if(!multistage_pivots_ok(ow, S, &S->D, &ow->sH, nu)) return -1;
        }
        // P = Hxx - Lxu Lxu'
        if(nx > 0){
            if(nu > 0){ // Lxu copied: syrk_ln with row offsets is wrong in some BLASFEO targets
                blasfeo_dgecp(nx, nu, &S->D, nu, 0, &ow->sAt, 0, 0);
                blasfeo_dsyrk_ln(nx, nu, -1.0, &ow->sAt, 0, 0, &ow->sAt, 0, 0, 1.0, &ow->sH, nu, nu, &S->D, nu, nu);
            }
            else blasfeo_dgecp(nx, nx, &ow->sH, 0, 0, &S->D, 0, 0);
            blasfeo_dtrtr_l(nx, &S->D, nu, nu, &S->D, nu, nu);
        }
    }
    return 0;
}

// y <-- P_t x (tmp: scratch of length nx)
static void multistage_Pmul(DAQPMultistageWork *ow, Stage *S, struct blasfeo_dvec *x, const int xi,
        struct blasfeo_dvec *tmp, struct blasfeo_dvec *y){
    if(S->nx == 0) return;
    if(ow->ric_sqrt){
        blasfeo_dtrmv_ltn(S->nx, &S->D, S->nn, S->nn, x, xi, tmp, 0);
        blasfeo_dtrmv_lnn(S->nx, &S->D, S->nn, S->nn, tmp, 0, y, 0);
    }
    else blasfeo_dgemv_n(S->nx, S->nx, 1.0, &S->D, S->nn, S->nn, x, xi, 0.0, y, 0, y, 0); // P_t, symmetric
}

// c <-- T'c = [Q2'c_u; c_x - K1'Q1'c_u] (length nu+nx to nn+nx; w: scratch)
static void multistage_tt(Stage *S, struct blasfeo_dvec *c, struct blasfeo_dvec *w){
    const int nu = S->nu, nx = S->nx, nn = S->nn;
    blasfeo_dgemv_t(nu, nu, 1.0, &S->Qd, 0, 0, c, 0, 0.0, w, 0, w, 0);
    if(nx > 0){
        blasfeo_dgemv_n(nx, S->r, -1.0, &S->Kt, 0, 0, w, 0, 1.0, c, nu, c, nu);
        memmove(c->pa + nn, c->pa + nu, nx*sizeof(c_float));
    }
    blasfeo_dveccp(nn, w, S->r, c, 0);
}

// Backward recursion for P_t, p_t, c_t and kv_t from stage `from` (later stages are kept)
static void multistage_lqr_backward(DAQPMultistageWork *ow, int from){
    struct blasfeo_dvec *s = &ow->ws, *g = &ow->wg, *e = &ow->wt;
    Stage *SN = ow->st + ow->N;
    int t;
    if(from >= ow->N){
        SN->infb = ow->neqc > 0 && multistage_eq_vectors(ow, ow->N, NULL, SN->h);
        if(ow->eps > 0) blasfeo_daxpy(SN->nx, -ow->eps, &SN->zc, 0, &SN->rq, 0, &SN->yp, SN->nn);
        else blasfeo_dveccp(SN->nx, &SN->rq, 0, &SN->yp, SN->nn);
        SN->cacc = 0;
        from = ow->N-1;
    }
    for(t = from; t >= 0; t--){
        Stage *S = ow->st+t;
        struct blasfeo_dvec *p = &S[1].yp, *gy = &S->rq, *bb = &S->b, *sp = p; // (p, and s at sp, from pi = S[1].nn)
        const int nx = S->nx, nu = S->nu, n = nu+nx, nn = S->nn, nx1 = S[1].nx, pi = S[1].nn;
        int has_bb = S->has_b, si = pi;
        c_float c = S[1].cacc;
        S->infb = S[1].infb | (ow->neqc > 0 && multistage_eq_vectors(ow, t, S[1].h, S->h));
        // g = [r; q] (- eps zc) (+ RSQ [u0; 0]), bb = b (+ B u0), c += 0.5 u0'R u0 + r'u0 (gy: g so far, rq if unchanged)
        if(ow->eps > 0){ blasfeo_daxpy(n, -ow->eps, &S->zc, 0, &S->rq, 0, g, 0); gy = g; }
        if(S->r > 0){
            if(gy != g){ blasfeo_dveccp(n, &S->rq, 0, g, 0); gy = g; }
            blasfeo_dvecse(n, 0.0, e, 0);
            blasfeo_dveccp(nu, &S->u0, 0, e, 0);
            c += 0.5*blasfeo_ddot(nu, &S->u0, 0, g, 0);
            blasfeo_dsymv_l(n, 1.0, S->RSQ, 0, 0, e, 0, 1.0, g, 0, g, 0);
            c += 0.5*blasfeo_ddot(nu, &S->u0, 0, g, 0);
            blasfeo_dgemv_t(nu, nx1, 1.0, S->BAbt, 0, 0, &S->u0, 0, S->has_b ? 1.0 : 0.0, &S->b, 0, &ow->wb, 0);
            bb = &ow->wb;
            has_bb = 1;
        }
        // s = P_{t+1} bb + p_{t+1}, c += 0.5 bb'P bb + p'bb
        if(has_bb && nx1 > 0){
            struct blasfeo_dvec *Pbb = S->r > 0 ? &ow->wc0 : &S->Pb; // (b + B u0 changes with the equalities)
            if(S->r > 0 || !S->Pb_valid){
                multistage_Pmul(ow, S+1, bb, 0, &ow->wr, Pbb);
                S->bPb = 0.5*blasfeo_ddot(nx1, bb, 0, Pbb, 0);
            }
            S->Pb_valid = S->r == 0;
            blasfeo_daxpy(nx1, 1.0, Pbb, 0, p, pi, s, 0);
            c += S->bPb + blasfeo_ddot(nx1, bb, 0, p, pi);
            sp = s; si = 0;
        }
        // g = gy + BAbt s, and of the free inputs: T'g
        blasfeo_dgemv_n(n, nx1, 1.0, S->BAbt, 0, 0, sp, si, 1.0, gy, 0, g, 0);
        if(S->r > 0) multistage_tt(S, g, &ow->wc1);
        // y = Luu\g_u, p_t = g_x - Lxu y (in yp), c -= 0.5|y|^2
        if(nn > 0){
            blasfeo_dtrsv_lnn_mn(nn+nx, nn, &S->D, 0, 0, g, 0, &S->yp, 0);
            c -= 0.5*blasfeo_ddot(nn, &S->yp, 0, &S->yp, 0);
        }
        else blasfeo_dveccp(nx, g, 0, &S->yp, 0);
        S->cacc = c;
    }
}

// Trajectory with nu_t = Luu'\(w_t - Lxu'x_t): v NULL gives z_lqr (w = kv, with x0, u0, b),
// otherwise the perturbation dz (w = v, without affine terms)
static void multistage_simulate(DAQPMultistageWork *ow, c_float *v){
    const int aff = v == NULL;
    Stage *S0 = ow->st;
    struct blasfeo_dvec *xi = &ow->wc1;
    int t;
    if(S0->nx > 0){
        if(aff) blasfeo_pack_dvec(S0->nx, ow->x0, 1, &S0->z, S0->nu);
        else blasfeo_dvecse(S0->nx, 0.0, &S0->dz, S0->nu);
    }
    for(t = 0; t < ow->N; t++){
        Stage *S = ow->st+t;
        const int nx = S->nx, nu = S->nu, nn = S->nn;
        struct blasfeo_dvec *ux = aff ? &S->z : &S->dz, *uxn = aff ? &S[1].z : &S[1].dz;
        struct blasfeo_dvec sv, *w = aff ? &S->yp : &sv; // (kv = -y)
        blasfeo_create_dvec(nn, &sv, aff ? NULL : v + S->cv);
        if(S->r == 0 && nu+nx <= DAQP_MULTISTAGE_SMALL){
            // Small stages: the same in loops
            c_float *z = ux->pa, *zn = uxn->pa + S[1].nu, sum;
            int i, j;
            for(j = 0; j < nu; j++){
                for(i = 0, sum = aff ? -w->pa[j] : w->pa[j]; i < nx; i++) sum -= BLASFEO_DMATEL(&S->D, nu+i, j)*z[nu+i];
                z[j] = sum;
            }
            for(j = nu-1; j >= 0; j--){
                for(i = j+1, sum = z[j]; i < nu; i++) sum -= BLASFEO_DMATEL(&S->D, i, j)*z[i];
                z[j] = sum/BLASFEO_DMATEL(&S->D, j, j);
            }
            for(j = 0; j < S[1].nx; j++){
                for(i = 0, sum = aff && S->has_b ? S->b.pa[j] : 0; i < nu+nx; i++) sum += BLASFEO_DMATEL(S->BAbt, i, j)*z[i];
                zn[j] = sum;
            }
            continue;
        }
        // [nu; x] in ux, or in xi with eliminated inputs: then u = T [nu; x] (+ u0) = Q [-K1 x; nu] (+ u0)
        struct blasfeo_dvec *y = S->r > 0 ? xi : ux;
        if(S->r > 0) blasfeo_dveccp(nx, ux, nu, xi, nn);
        if(nn > 0){
            blasfeo_dveccpsc(nn, aff ? -1.0 : 1.0, w, 0, y, 0);
            bf_trsv_ltn_mn(nn+nx, nn, &S->D, y);
        }
        if(S->r > 0){
            blasfeo_dveccp(nn, xi, 0, &ow->wt, S->r);
            blasfeo_dgemv_t(nx, S->r, -1.0, &S->Kt, 0, 0, xi, nn, 0.0, &ow->wt, 0, &ow->wt, 0);
            blasfeo_dgemv_n(nu, nu, 1.0, &S->Qd, 0, 0, &ow->wt, 0, aff ? 1.0 : 0.0, &S->u0, 0, ux, 0);
        }
        blasfeo_dgemv_t(nu+nx, S[1].nx, 1.0, S->BAbt, 0, 0, ux, 0, aff && S->has_b ? 1.0 : 0.0, &S->b, 0, uxn, S[1].nu);
    }
}

// Forward recursion: z_lqr (depends on x0) and the LDP bounds (from multistage_fetch_bounds)
static void multistage_lqr_forward(DAQPMultistageWork *ow){
    DAQPWorkspace *work = ow->work;
    Stage *S0 = ow->st;
    int i, k;
    // The state equalities that remain for x0
    ow->infeasible = S0->infb;
    if(ow->neqc > 0 && S0->rh > 0 && !ow->x0free){
        c_float xmax = 0;
        for(i = 0; i < S0->nx; i++) if(fabs(ow->x0[i]) > xmax) xmax = fabs(ow->x0[i]);
        for(k = 0; k < S0->rh; k++){
            c_float s = -S0->h[k];
            for(i = 0; i < S0->nx; i++) s += BLASFEO_DMATEL(&S0->Ht, i, k)*ow->x0[i];
            if(fabs(s) > work->settings->primal_tol*(1+xmax+fabs(S0->h[k]))) ow->infeasible = 1;
        }
    }
    multistage_simulate(ow, NULL);
    // J(z_lqr) = 0.5 x0'P_0 x0 + p_0'x0 + c_0
    ow->Jlqr = S0->cacc;
    if(S0->nx > 0){
        multistage_Pmul(ow, S0, &S0->z, S0->nu, &ow->wc0, &ow->wr);
        blasfeo_daxpy(S0->nx, 2.0, &S0->yp, S0->nn, &ow->wr, 0, &ow->wr, 0);
        ow->Jlqr += 0.5*blasfeo_ddot(S0->nx, &ow->wr, 0, &S0->z, S0->nu);
    }
    work->fval_offset = ow->Jlqr;
    // Bounds of the LDP: d = b - a'z_lqr
    multistage_eval(ow, 0, work->Mu);
    for(i = 0; i < ow->m; i++){
        const c_float u = work->dupper[i], l = work->dlower[i], mu = work->Mu[i], s = ow->scale[i];
        work->dupper[i] = u < DAQP_INF ? (u-mu)*s : u;
        work->dlower[i] = l > -DAQP_INF ? (l-mu)*s : l;
    }
}

/* The operator */

// Mu <-- M*u for all constraints
static void multistage_product(void *data, const c_float *v, c_float *Mu){
    DAQPMultistageWork *ow = (DAQPMultistageWork*)data;
    int i, t;
    for(i = 0; i < ow->nv && v[i] == 0; i++);
    memcpy(ow->vlast, v, ow->nv*sizeof(c_float));
    ow->dz_valid = 1;
    if(i == ow->nv){ // v = 0 (e.g. an empty working set, z_lqr is the candidate): dz = 0, Mu = 0
        for(t = 0; t <= ow->N; t++)
            blasfeo_dvecse(ow->st[t].nu + ow->st[t].nx, 0.0, &ow->st[t].dz, 0);
        memset(Mu, 0, ow->m*sizeof(c_float));
        return;
    }
    multistage_simulate(ow, (c_float*)v);
    multistage_eval(ow, 1, Mu);
    for(i = 0; i < ow->m; i++) Mu[i] *= ow->scale[i];
}

// Prox: the LDP of the inner problem around the center
static void multistage_prox_center(void *data){
    DAQPMultistageWork *ow = (DAQPMultistageWork*)data;
    multistage_fetch_bounds(ow);
    multistage_lqr_backward(ow, ow->N);
    multistage_lqr_forward(ow);
    ow->work->reuse_ind = 0;
    ow->dz_valid = 0;
}

// The latest inner solution becomes the center; returns the largest change
static c_float multistage_prox_primal(void *data){
    DAQPMultistageWork *ow = (DAQPMultistageWork*)data;
    c_float diff = 0;
    int t, i;
    multistage_simulate(ow, ow->work->u);
    memcpy(ow->vlast, ow->work->u, ow->nv*sizeof(c_float));
    ow->dz_valid = 1;
    for(t = 0; t <= ow->N; t++){
        Stage *S = ow->st+t;
        for(i = 0; i < S->nu + S->nx; i++){
            const c_float z = S->z.pa[i] + S->dz.pa[i];
            if(fabs(z - S->zc.pa[i]) > diff) diff = fabs(z - S->zc.pa[i]);
            S->zc.pa[i] = z;
        }
    }
    return diff;
}

// Backward recursion for the row of id down to stage smin, by vectors (gemv/trsv: blocks of one row in
// multistage_rows_block would use matrix-matrix kernels for a single row). slot >= 0: the blocks go to that
// slot of the cache (smin = 0); else returns their squared norm, and *g gets [cu; cx] of stage smin (cx from nn)
static c_float multistage_row_vec(DAQPMultistageWork *ow, const int id, const int smin, const int slot,
        struct blasfeo_dvec **g){
    const int smax = multistage_last_stage(ow, id);
    int s, has_u;
    struct blasfeo_dvec *c = &ow->wc0, *cn = &ow->wc1, *swp;
    c_float n2 = 0, *row = slot >= 0 ? ow->rows + (size_t)slot*ow->nvmax : NULL;

    has_u = multistage_coeffs(ow, id, smax, c);
    if(row != NULL) ow->rows_len[slot] = ow->st[smax].cv + (has_u ? ow->st[smax].nn : 0); // (all of it is written)
    for(s = smax; s >= smin; s--){
        Stage *S = ow->st+s;
        const int nx = S->nx, nu = S->nu, nn = S->nn;
        if(s < smax){
            // c = [cu; cx] of stage s + BAbt_s p (p: cx of stage s+1)
            multistage_coeffs(ow, id, s, c);
            blasfeo_dgemv_n(nu+nx, S[1].nx, 1.0, S->BAbt, 0, 0, cn, S[1].nn, 1.0, c, 0, c, 0);
        }
        if(S->r > 0 && (s < smax || has_u)) multistage_tt(S, c, &ow->wt); // T'c
        else if(S->r > 0 && nx > 0) memmove(c->pa + nn, c->pa + nu, nx*sizeof(c_float)); // No input terms
        if(nn > 0 && (s < smax || has_u)){ // No input terms: zero block
            // v_s block: Luu\cu, and cx -= Lxu Luu\cu (in c)
            blasfeo_dtrsv_lnn_mn(nn+nx, nn, &S->D, 0, 0, c, 0, c, 0);
            if(row != NULL) memcpy(row + S->cv, c->pa, nn*sizeof(c_float));
            else n2 += blasfeo_ddot(nn, c, 0, c, 0);
        }
        if(s > smin){ swp = c; c = cn; cn = swp; }
    }
    if(g != NULL) *g = c;
    return n2;
}

// |c|^2 of the coefficients of constraint id
static c_float multistage_coef_norm2(DAQPMultistageWork *ow, const int id){
    int s, k;
    c_float c2 = 0;
    if(multistage_decode(ow, id, &s, &k) == DAQP_MULTISTAGE_BOUND) return 1;
    for(s = ow->ct[id]; s <= multistage_last_stage(ow, id); s++){
        multistage_coeffs(ow, id, s, &ow->wc0);
        c2 += blasfeo_ddot(ow->st[s].nu + ow->st[s].nx, &ow->wc0, 0, &ow->wc0, 0);
    }
    return c2;
}

// Whether row id of M (|M_id|^2 = n2) is zero relative to the coefficients of the constraint
static int multistage_zero_row(DAQPMultistageWork *ow, const int id, const c_float n2){
    return !(n2 > ow->work->settings->zero_tol*multistage_coef_norm2(ow, id));
}

// The normalization is recomputed at the next factorization if row id (|M_id|^2 = n2, scaled) has drifted, or if it
// was zero and no longer is
static void multistage_norm_check(DAQPMultistageWork *ow, const int id, const c_float n2){
    const c_float drift2 = DAQP_MULTISTAGE_NORM_DRIFT*DAQP_MULTISTAGE_NORM_DRIFT;
    if(ow->zrow[id] ? !multistage_zero_row(ow, id, n2) : n2*drift2 < 1 || n2 > drift2)
        ow->norm_stale = 1;
}

// Scale the row in slot s and check its normalization
static void multistage_row_finish(DAQPMultistageWork *ow, const int id, const int s){
    c_float *row = ow->rows + (size_t)s*ow->nvmax, n2 = 0;
    int i;
    for(i = 0; i < ow->rows_len[s]; i++){ row[i] *= ow->scale[id]; n2 += row[i]*row[i]; }
    multistage_norm_check(ow, id, n2);
}

// Empty the cache of rows
static void multistage_cache_clear(DAQPMultistageWork *ow){
    int i;
    for(i = 0; i < ow->nslots; i++) ow->slot_id[i] = -1;
    for(i = 0; i < ow->m; i++) ow->id_slot[i] = -1;
}

// A free slot (its constraint is no longer active) for constraint id
static int multistage_new_slot(DAQPMultistageWork *ow, const int id){
    int s;
    for(s = 0; s < ow->nslots; s++){
        const int old = ow->slot_id[s];
        if(old < 0) break;
        if(!(ow->work->sense[old] & DAQP_ACTIVE) || ow->id_slot[old] != s){
            ow->id_slot[old] = -1;
            break;
        }
    }
    // All slots taken by constraints marked active (possible during an activation): one not in the working set nor pinned
    for(s = s < ow->nslots ? s : 0; s < ow->nslots-1 && ow->slot_id[s] >= 0 && ow->id_slot[ow->slot_id[s]] == s; s++){
        const int old = ow->slot_id[s];
        int j, used = 0;
        for(j = 0; j < ow->work->n_active && !used; j++) used = ow->work->WS[j] == old;
        for(j = 0; j < ow->npin && !used; j++) used = ow->pin[j] == old;
        if(!used){
            ow->id_slot[old] = -1;
            break;
        }
    }
    ow->slot_id[s] = id;
    return ow->id_slot[id] = s;
}

// One stage of the backward recursion for rows C = [Cu Cx]: C T, V = Cu Luu^{-T}, Cx -= V Lxu'
static void multistage_block_stage(DAQPMultistageWork *ow, Stage *S, struct blasfeo_dmat *C, const int kk){
    const int nu = S->nu, nx = S->nx, nn = S->nn;
    if(kk == 0) return;
    if(S->r > 0){ // Cu Q in sC2
        blasfeo_dgemm_nn(kk, nu, nu, 1.0, C, 0, 0, &S->Qd, 0, 0, 0.0, &ow->sC2, 0, 0, &ow->sC2, 0, 0);
        if(nx > 0){
            blasfeo_dgemm_nt(kk, nx, S->r, -1.0, &ow->sC2, 0, 0, &S->Kt, 0, 0, 1.0, C, 0, nu, &ow->sC2, 0, nu);
            blasfeo_dgecp(kk, nx, &ow->sC2, 0, nu, C, 0, nn);
        }
        blasfeo_dgecp(kk, nn, &ow->sC2, 0, S->r, C, 0, 0);
    }
    if(nn > 0){
        blasfeo_dtrsm_rltn(kk, nn, 1.0, &S->D, 0, 0, C, 0, 0, &ow->sV, 0, 0);
        blasfeo_dgemm_nt(kk, nx, nn, -1.0, &ow->sV, 0, 0, &S->D, nn, 0, 1.0, C, 0, nn, C, 0, nn);
    }
}

// Rows of (at most DAQP_MULTISTAGE_BLOCK) constraints sorted by last stage (decreasing), by one backward
// recursion for all (Luu\cu per stage, with cx propagated by BAbt). If their slots are consecutive (as in an empty
// cache), the blocks of a stage are copied to all rows at once
static void multistage_rows_block(DAQPMultistageWork *ow, const int *ids, const int *smax, const int k){
    struct blasfeo_dmat *C = &ow->sC0, *Cn = &ow->sC1, *swp;
    struct blasfeo_dvec *c = &ow->wc0;
    int j, s, kk = 0, slot[DAQP_MULTISTAGE_BLOCK] = {0}, has_u[DAQP_MULTISTAGE_BLOCK] = {0}, consec = 1;
    ow->pin = ids; ow->npin = k;
    for(s = smax[0]; s >= 0; s--){
        Stage *S = ow->st+s;
        const int nu = S->nu, nx = S->nx, n = nu+nx, nn = S->nn;
        // The rows that have joined: C = Cx_{s+1} BAbt'
        if(kk > 0 && n > 0){
            if(s < ow->N && S[1].nx > 0)
                blasfeo_dgemm_nt(kk, n, S[1].nx, 1.0, Cn, 0, S[1].nn, S->BAbt, 0, 0, 0.0, C, 0, 0, C, 0, 0);
            else blasfeo_dgese(kk, n, 0.0, C, 0, 0);
        }
        // + the coefficients of stage s (of joining constraints, and of coupling ones of stage s)
        for(j = 0; j < k && smax[j] >= s; j++){
            if(smax[j] == s){
                has_u[j] = multistage_coeffs(ow, ids[j], s, c);
                blasfeo_drowin(n, 1.0, c, 0, C, j, 0);
                slot[j] = multistage_new_slot(ow, ids[j]);
                consec &= slot[j] == slot[0] + j;
                ow->rows_len[slot[j]] = S->cv + (has_u[j] ? nn : 0);
                memset(ow->rows + (size_t)slot[j]*ow->nvmax, 0, ow->rows_len[slot[j]]*sizeof(c_float));
                kk = j+1;
            }
            else if(ow->ct[ids[j]] == s){ // (a coupling constraint, which ends at stage s+1)
                multistage_coeffs(ow, ids[j], s, c);
                blasfeo_drowad(n, 1.0, c, 0, C, j, 0);
            }
        }
        multistage_block_stage(ow, S, C, kk);
        // (the block of a row without input terms at its last stage is zero there: also when copied, beyond its length)
        if(consec && kk > 0 && nn > 0)
            blasfeo_unpack_tran_dmat(kk, nn, &ow->sV, 0, 0, ow->rows + (size_t)slot[0]*ow->nvmax + S->cv, ow->nvmax);
        else for(j = 0; j < kk && nn > 0; j++){
            struct blasfeo_dvec r;
            blasfeo_create_dvec(nn, &r, ow->rows + (size_t)slot[j]*ow->nvmax + S->cv);
            if(s < smax[j] || has_u[j]) blasfeo_drowex(nn, 1.0, &ow->sV, j, 0, &r, 0);
        }
        swp = C; C = Cn; Cn = swp;
    }
    for(j = 0; j < k; j++) multistage_row_finish(ow, ids[j], slot[j]);
    ow->npin = 0;
}

// Row of constraint id, computed into the cache if needed
static const c_float* multistage_row(void *data, const int id, int *len){
    DAQPMultistageWork *ow = (DAQPMultistageWork*)data;
    int s = ow->id_slot[id];
    if(s < 0){
        multistage_row_vec(ow, id, 0, s = multistage_new_slot(ow, id), NULL);
        multistage_row_finish(ow, id, s);
    }
    *len = ow->rows_len[s];
    return ow->rows + (size_t)s*ow->nvmax;
}

// Rows of the constraints marked active (and not cached), in blocks
static void multistage_rows_active(DAQPMultistageWork *ow){
    int *ids = ow->ids, *smax = ow->smax;
    int i, j, k = 0;
    for(i = 0; i < ow->m && k < ow->nslots-1; i++){
        if(!(ow->work->sense[i] & DAQP_ACTIVE) || ow->id_slot[i] >= 0 ||
                (ow->eqf[i] & (DAQP_MULTISTAGE_EQ_X0 | DAQP_MULTISTAGE_EQ_ELIM))) continue;
        // Insertion sort by last stage, decreasing
        const int sm = multistage_last_stage(ow, i);
        for(j = k++; j > 0 && smax[j-1] < sm; j--){ ids[j] = ids[j-1]; smax[j] = smax[j-1]; }
        ids[j] = i; smax[j] = sm;
    }
    for(i = 0; i < k; i += DAQP_MULTISTAGE_BLOCK)
        multistage_rows_block(ow, ids+i, smax+i, k-i < DAQP_MULTISTAGE_BLOCK ? k-i : DAQP_MULTISTAGE_BLOCK);
}

// After multistage_activate formed L at the start of the cache, with the rows of the working set in L (len[i] entries
// at L + off[i], one after the other): the two swap places (through the scratch sAb), and the rows go to slots 0..k-1
// (backwards: row i moves to i*nvmax >= off[i], beyond its successor's source)
static void multistage_rows_back(DAQPMultistageWork *ow, const int k, const int *len, const int *off){
    DAQPWorkspace *work = ow->work;
    const size_t nl = DAQP_ARSUM((size_t)k), nr = (size_t)off[k-1] + len[k-1], n = nl > nr ? nl : nr;
    const size_t chunk = ow->sAb.memsize/sizeof(c_float);
    c_float *buf = ow->sAb.pA;
    size_t o, c;
    int i;
    for(o = 0; o < n; o += c){
        c = n - o < chunk ? n - o : chunk;
        memcpy(buf, work->L + o, c*sizeof(c_float));
        memcpy(work->L + o, ow->rows + o, c*sizeof(c_float));
        memcpy(ow->rows + o, buf, c*sizeof(c_float));
    }
    multistage_cache_clear(ow);
    for(i = k-1; i >= 0; i--){
        memmove(ow->rows + (size_t)i*ow->nvmax, ow->rows + off[i], len[i]*sizeof(c_float));
        ow->rows_len[i] = len[i];
        ow->slot_id[i] = work->WS[i];
        ow->id_slot[work->WS[i]] = i;
    }
}

// Warm start after an update (as daqp_activate_constraints): the LDL factorization of the working set at once, as the
// Cholesky factorization of M_A M_A' + diag(rho) (rho: the weight of a soft constraint with a free slack, else 0;
// equalities first). Falls back to one by one additions in branch and bound, for a small pivot, and for few
// constraints. The factor is formed in DAQP's L, else in the row cache (the rows then moved to L), else for the longest
// leading part of the working set that fits, with the rest added one by one
static int multistage_activate(DAQPMultistageWork *ow){
    DAQPWorkspace *work = ow->work;
    const int n = work->n;
    int *len = ow->smax, *off = ow->ids; // Rows of the working set: len[i] entries at rows + off[i]
    c_float *rho = ow->sic; // (until the reciprocals of the diagonal replace it)
    int i, j, k = 0, ks = 0, kb, kl, kc, pass, j0, b, kr, rs, c, fused = 0, moved = 0;
    const size_t lsize = DAQP_ARSUM((size_t)ow->nvmax+ow->nsmax+1)*sizeof(c_float);
    const c_float *rows = ow->rows;
    char *la;
    struct blasfeo_dmat G;
    if(work->bnb != NULL || work->n_active > 0) return daqp_activate_constraints(work);
    for(i = 0; i < work->m; i++){
        if(!DAQP_IS_ACTIVE(i)) continue;
        k++;
        if(DAQP_IS_SOFT(i) && (DAQP_IS_SLACK_FREE(i) || daqp_soft_w(work, i) == 0)) ks++;
    }
    if(k < DAQP_MULTISTAGE_ACT_MIN || k > n + ks) return daqp_activate_constraints(work);
    for(pass = 0, k = 0; pass < 2; pass++)
        for(i = 0; i < work->m; i++)
            if(DAQP_IS_ACTIVE(i) && (DAQP_IS_IMMUTABLE(i) != 0) == (pass == 0)) work->WS[k++] = i;
    for(i = 0; i < k; i++){
        const int id = work->WS[i];
        rho[i] = DAQP_IS_SOFT(id) && (DAQP_IS_SLACK_FREE(id) || daqp_soft_w(work, id) == 0) ? daqp_soft_rho(work, id) : 0;
    }
    // G in L, which holds DAQP_ARSUM(nvmax+nsmax+1) entries (the core workspace for nvmax variables and nsmax soft
    // constraints)
    for(i = 0; i < k; i++) off[i] = (int)(multistage_row(ow, work->WS[i], len+i) - ow->rows); // (all in the cache)
    // The longest leading part of the working set whose factor fits in L (kl), or in the row cache with its rows in L
    // (kc)
    la = (char*)(((uintptr_t)work->L + DAQP_MULTISTAGE_ALIGN-1) & ~(uintptr_t)(DAQP_MULTISTAGE_ALIGN-1));
    for(kl = k; kl > 0 && blasfeo_memsize_dmat(kl, kl) + (size_t)(la - (char*)work->L) > lsize; kl--);
    kb = k;
    if(kl < k){
        char *lc = (char*)(((uintptr_t)ow->rows + DAQP_MULTISTAGE_ALIGN-1) & ~(uintptr_t)(DAQP_MULTISTAGE_ALIGN-1));
        size_t tot = 0;
        for(kc = 0; kc < k && (tot + len[kc])*sizeof(c_float) <= lsize &&
                blasfeo_memsize_dmat(kc+1, kc+1) + (size_t)(lc - (char*)ow->rows) <=
                (size_t)ow->nslots*ow->nvmax*sizeof(c_float); tot += len[kc++]);
        if(kc > kl){
            for(i = 0, tot = 0; i < kc; tot += len[i++]){
                memcpy(work->L + tot, ow->rows + off[i], len[i]*sizeof(c_float));
                off[i] = (int)tot;
            }
            rows = work->L;
            la = lc;
            moved = 1;
        }
        kb = kc > kl ? kc : kl;
        if(kb < DAQP_MULTISTAGE_ACT_MIN) return daqp_activate_constraints(work);
    }
    blasfeo_create_dmat(kb, kb, &G, la);
    // M_A M_A' = C C' (lower Cholesky factor in G), accumulated over blocks of DAQP_MULTISTAGE_ACT_COLS columns of M_A
    // (packed into sAb; zero beyond the end of the rows). The rows that end before a block (first in the working set,
    // as they are ordered by stage) do not change G there: a block updates G from the first row that reaches into it
    // (rounded down to a panel; staircase)
    blasfeo_dgese(kb, kb, 0.0, &G, 0, 0);
    for(i = 0; i < kb; i++) if(rho[i] != 0) BLASFEO_DMATEL(&G, i, i) = rho[i];
    for(j0 = 0; j0 < n; j0 += DAQP_MULTISTAGE_ACT_COLS){
        b = n - j0 < DAQP_MULTISTAGE_ACT_COLS ? n - j0 : DAQP_MULTISTAGE_ACT_COLS;
        for(kr = 0; kr < kb && len[kr] <= j0; kr++);
        if(kr == kb) continue; // (no row reaches into the block)
        rs = kr/D_PS*D_PS;
        for(i = rs; i < kb; i++){
            c = len[i] - j0 < b ? len[i] - j0 : b;
            if(c < 0) c = 0;
            blasfeo_pack_tran_dmat(c, 1, (c_float*)rows + off[i] + j0, len[i], &ow->sAb, i, 0);
            if(c < b) blasfeo_dgese(1, b-c, 0.0, &ow->sAb, i, c);
        }
        if(j0 + b < n || rs > 0) blasfeo_dsyrk_ln(kb-rs, b, 1.0, &ow->sAb, rs, 0, &ow->sAb, rs, 0, 1.0, &G, rs, rs, &G, rs, rs);
        else{ blasfeo_dsyrk_dpotrf_ln(kb, b, &ow->sAb, 0, 0, &ow->sAb, 0, 0, &G, 0, 0, &G, 0, 0); fused = 1; }
    }
    if(!fused) blasfeo_dpotrf_l(kb, &G, 0, 0, &G, 0, 0);
    for(i = 0; i < kb; i++){
        const c_float ci = BLASFEO_DMATEL(&G, i, i), d = ci*ci;
        if(!(ci > 0) || d < work->settings->sing_tol || d < work->settings->pivot_tol) break;
        work->D[i] = d;
        ow->sic[i] = 1/ci;
    }
    if(i < kb){ // (the additions one by one decide)
        if(moved) multistage_cache_clear(ow); // (C overwrote the rows)
        return daqp_activate_constraints(work);
    }
    // L (unit lower, packed by rows) and D of DAQP: l_ij = c_ij/c_jj, d_i = c_ii^2. By panels of rows of C, copied
    // first (C is in the memory it is converted into: the rows of L up to a panel end before it)
    for(b = 0; b < kb; b += D_PS){
        const int r = kb - b < D_PS ? kb - b : D_PS;
        blasfeo_unpack_dmat(r, b+r, &G, b, 0, ow->sAb.pA, r);
        for(c = 0; c < r; c++){
            c_float *Li = (moved ? ow->rows : work->L) + DAQP_ARSUM(b+c);
            for(j = 0; j < b+c; j++) Li[j] = ow->sAb.pA[j*r+c]*ow->sic[j];
            Li[b+c] = 1;
        }
    }
    if(moved) multistage_rows_back(ow, kb, len, off);
    for(i = 0; i < kb; i++){ // (multipliers and slack states as daqp_activate_constraint)
        const int id = work->WS[i];
        c_float lam = 1.0;
        work->sense[id] &= ~DAQP_SLACK_SWITCHED;
        if(DAQP_IS_SOFT(id)){
            const c_float w = daqp_soft_w(work, id);
            DAQP_SET_MUTABLE(id);
            if(w > 0) lam = DAQP_IS_SLACK_FREE(id) ? w : 0;
            else DAQP_SET_SLACK_FREE(id);
        }
        work->lam[i] = DAQP_IS_LOWER(id) ? -lam : lam;
    }
    work->n_active = kb;
    work->sing_ind = DAQP_EMPTY_IND;
    work->reuse_ind = 0;
    if(kb == k) return 1;
    // The rest one by one: daqp_activate_constraints adds those marked active, so the ones added are unmarked meanwhile
    for(i = 0; i < kb; i++){ off[i] = work->WS[i]; DAQP_SET_INACTIVE(off[i]); }
    j = daqp_activate_constraints(work);
    for(i = 0; i < kb; i++) DAQP_SET_ACTIVE(off[i]);
    return j;
}

/* Normalization of the rows of M */

// Sig <-- Sig_{t+1} = W'W + Acl Sig_t Acl' (W' = Br Luu^{-T}, Acl = Ar - W'Lxu')
static void multistage_covariance_next(DAQPMultistageWork *ow, const int t){
    Stage *S = ow->st+t;
    const int nx = S->nx, nu = S->nn, nx1 = S[1].nx;
    struct blasfeo_dmat *Sn = ow->Sig == ow->sSig ? ow->sSig+1 : ow->sSig;
    if(nx1 > 0){
        blasfeo_dgetr(nu+nx, nx1, S->BAr, 0, 0, &ow->sAt, 0, 0);
        if(nu > 0){
            blasfeo_dtrsm_rltn(nx1, nu, 1.0, &S->D, 0, 0, &ow->sAt, 0, 0, &ow->sCw, 0, 0);
            blasfeo_dsyrk_ln(nx1, nu, 1.0, &ow->sCw, 0, 0, &ow->sCw, 0, 0, 0.0, Sn, 0, 0, Sn, 0, 0);
            blasfeo_dgemm_nt(nx1, nx, nu, -1.0, &ow->sCw, 0, 0, &S->D, nu, 0, 1.0, &ow->sAt, 0, nu, &ow->sAt, 0, nu);
        }
        else blasfeo_dgese(nx1, nx1, 0.0, Sn, 0, 0);
        if(nx > 0){
            blasfeo_dgemm_nn(nx1, nx, nx, 1.0, &ow->sAt, 0, nu, ow->Sig, 0, 0, 0.0, &ow->sG, 0, 0, &ow->sG, 0, 0);
            blasfeo_dsyrk_ln(nx1, nx, 1.0, &ow->sG, 0, 0, &ow->sAt, 0, nu, 1.0, Sn, 0, 0, Sn, 0, 0);
        }
        blasfeo_dtrtr_l(nx1, Sn, 0, 0, Sn, 0, 0);
    }
    ow->Sig = Sn;
}

// Scale from |M_id|^2 = n2; a row that is zero relative to the coefficients is not scaled
static void multistage_set_scale(DAQPMultistageWork *ow, const int id, const c_float n2){
    ow->zrow[id] = multistage_zero_row(ow, id, n2);
    if(!ow->zrow[id]) ow->scale[id] = 1/sqrt(n2);
}

// Scales of k bounds and general constraints of stage t at once: |M_i|^2 = diag(V V' + Cx Sig Cx')
static void multistage_scale_block(DAQPMultistageWork *ow, const int t, const int *ids, const int k){
    Stage *S = ow->st+t;
    const int n = S->nu + S->nx, nx = S->nx, nn = S->nn;
    struct blasfeo_dmat *C = &ow->sC0, *K = &ow->sK;
    int j;
    for(j = 0; j < k; j++){
        multistage_coeffs(ow, ids[j], t, &ow->wc0);
        blasfeo_drowin(n, 1.0, &ow->wc0, 0, C, j, 0);
    }
    multistage_block_stage(ow, S, C, k);
    if(nn > 0) blasfeo_dsyrk_ln(k, nn, 1.0, &ow->sV, 0, 0, &ow->sV, 0, 0, 0.0, K, 0, 0, K, 0, 0);
    else blasfeo_dgese(k, k, 0.0, K, 0, 0);
    if(nx > 0){
        blasfeo_dgemm_nn(k, nx, nx, 1.0, C, 0, nn, ow->Sig, 0, 0, 0.0, &ow->sC1, 0, 0, &ow->sC1, 0, 0);
        blasfeo_dsyrk_ln(k, nx, 1.0, &ow->sC1, 0, 0, C, 0, nn, 1.0, K, 0, 0, K, 0, 0);
    }
    for(j = 0; j < k; j++) multistage_set_scale(ow, ids[j], BLASFEO_DMATEL(K, j, j));
}

// Normalize the rows of M (after factorization): |M_i|^2 = |blocks|^2 + g'Sig_t g
static void multistage_normalize(DAQPMultistageWork *ow){
    int i, t, ids[DAQP_MULTISTAGE_BLOCK], k;
    ow->norm_stale = 0;
    ow->Sig = ow->sSig; // Sig_0 = 0
    blasfeo_dgese(ow->st[0].nx, ow->st[0].nx, 0.0, ow->Sig, 0, 0);
    for(i = 0; i < ow->m; i++){ ow->scale[i] = 1; ow->zrow[i] = 0; }
    for(t = 0; t <= ow->N; t++){
        Stage *S = ow->st+t;
        for(i = S->c0, k = 0; i < S->c0 + S->nb + S->ng + S->nc; i++){
            int tt, kc, j;
            const int kind = multistage_decode(ow, i, &tt, &kc);
            if(ow->eqf[i] & (DAQP_MULTISTAGE_EQ_X0 | DAQP_MULTISTAGE_EQ_ELIM)) continue;
            if(kind == DAQP_MULTISTAGE_BOUND && (j = S->idxb[kc] - S->nu) >= 0)
                multistage_set_scale(ow, i, blasfeo_dgeex1(ow->Sig, j, j)); // Bound on dx_t: its variance
            else if(kind == DAQP_MULTISTAGE_COUPLING){
                struct blasfeo_dvec *g;
                c_float m2 = multistage_row_vec(ow, i, t, -1, &g);
                blasfeo_dgemv_n(S->nx, S->nx, 1.0, ow->Sig, 0, 0, g, S->nn, 0.0, &ow->wr, 0, &ow->wr, 0);
                multistage_set_scale(ow, i, m2 + blasfeo_ddot(S->nx, g, S->nn, &ow->wr, 0));
            }
            else{
                ids[k++] = i;
                if(k == DAQP_MULTISTAGE_BLOCK){ multistage_scale_block(ow, t, ids, k); k = 0; }
            }
        }
        if(k > 0) multistage_scale_block(ow, t, ids, k);
        if(t < ow->N) multistage_covariance_next(ow, t);
    }
}

/* Warm start */

// The state of a constraint in the working set that a warm start reproduces
#define DAQP_MULTISTAGE_WS_BITS (DAQP_ACTIVE | DAQP_LOWER | DAQP_SLACK_FIXED)

// Whether constraint i takes part in warm starts (mutable, and in the LDP)
static int multistage_ws_mutable(const DAQPMultistageWork *ow, const int i){
    return !(ow->work->sense[i] & DAQP_IMMUTABLE) && !(ow->eqf[i] & (DAQP_MULTISTAGE_EQ_X0 | DAQP_MULTISTAGE_EQ_ELIM));
}
// DAQP_MULTISTAGE_UPDATE_SHIFT, for a problem shifted by one stage: constraint d of stage t takes the state of the
// constraint wssrc[d] of stage t+1 on the same variables (a bound on the same entry of [u; x], one to one; the general
// and coupling constraint with the same index, if the stages have as many); -1 for the others, and those of the last
// stage. wsdst is the inverse
static void multistage_ws_shift_map(DAQPMultistageWork *ow){
    int t, j, k;
    for(j = 0; j < ow->m; j++) ow->wssrc[j] = ow->wsdst[j] = -1;
    for(t = ow->x0free; t < ow->N; t++){ // (not the stage of a free x_0: its variables are x_0)
        const Stage *S = ow->st+t, *Sn = S+1;
        const int *ib = S->nb > 0 ? ow->user->idxb[t - ow->x0free] : NULL, *ibn = Sn->nb > 0 ? ow->user->idxb[t+1 - ow->x0free] : NULL;
        if(S->nu != Sn->nu || S->nx != Sn->nx) continue;
        for(j = 0; j < S->nb; j++){ // (idxb of the user: the stages are bound at updates)
            for(k = 0; k < Sn->nb && (ibn[k] != ib[j] || ow->wsdst[Sn->c0 + k] >= 0); k++);
            if(k < Sn->nb){ ow->wssrc[S->c0 + j] = Sn->c0 + k; ow->wsdst[Sn->c0 + k] = S->c0 + j; }
        }
        for(j = 0; S->ng == Sn->ng && j < S->ng; j++){
            ow->wssrc[S->c0 + S->nb + j] = Sn->c0 + Sn->nb + j;
            ow->wsdst[Sn->c0 + Sn->nb + j] = S->c0 + S->nb + j;
        }
        for(j = 0; t+1 < ow->N && S->nc == Sn->nc && j < S->nc; j++){
            ow->wssrc[S->c0 + S->nb + S->ng + j] = Sn->c0 + Sn->nb + Sn->ng + j;
            ow->wsdst[Sn->c0 + Sn->nb + Sn->ng + j] = S->c0 + S->nb + S->ng + j;
        }
    }
}
// The constraint whose state constraint d takes in the shift, or -1 (not to a zero row, as of a constraint on a fixed
// x_0 alone: it cannot be in the working set)
static int multistage_ws_src(const DAQPMultistageWork *ow, const int d){
    const int s = ow->wssrc[d];
    return s >= 0 && multistage_ws_mutable(ow, d) && multistage_ws_mutable(ow, s) && !ow->zrow[d] ? s : -1;
}
// Drop the history of working sets
static void multistage_ws_forget(DAQPMultistageWork *ow){
    ow->ws_n = ow->ws_alt = ow->ws_two = ow->ws_shift = 0;
}
// Number of constraints whose state differs between the working sets a and b of the history
static int multistage_ws_diff(const DAQPMultistageWork *ow, const int a, const int b){
    const unsigned char *ca = ow->wsc[a], *cb = ow->wsc[b];
    int i, d = 0;
    for(i = 0; i < ow->wsl[a]; i++) d += ca[ow->wsh[a][i]] != cb[ow->wsh[a][i]];
    for(i = 0; i < ow->wsl[b]; i++) d += ca[ow->wsh[b][i]] == 0;
    return d;
}

// After a solve: its final working set (ok: solved, else the history is dropped). The next warm start is the working
// set of two solves back if the working sets alternate (as when a closed loop chatters at saturated inputs): if it fit
// clearly better than the previous one (at most half as many differences, and the previous one at least 2) for this
// solve and the one before
static void multistage_ws_record(DAQPMultistageWork *ow, const int ok){
    const DAQPWorkspace *work = ow->work;
    int i, n = 0, d1, d2, ds, alt, *l;
    unsigned char *c;
    if(!ok){ multistage_ws_forget(ow); return; }
    // The oldest working set (cleared) becomes the newest
    c = ow->wsc[2]; l = ow->wsh[2];
    for(i = 0; i < ow->wsl[2]; i++) c[l[i]] = 0;
    ow->wsc[2] = ow->wsc[1]; ow->wsh[2] = ow->wsh[1]; ow->wsl[2] = ow->wsl[1];
    ow->wsc[1] = ow->wsc[0]; ow->wsh[1] = ow->wsh[0]; ow->wsl[1] = ow->wsl[0];
    ow->wsc[0] = c; ow->wsh[0] = l;
    for(i = 0; i < work->n_active; i++){
        const int id = work->WS[i];
        if(id < 0 || id >= ow->m || !multistage_ws_mutable(ow, id)) continue;
        c[id] = (unsigned char)(work->sense[id] & DAQP_MULTISTAGE_WS_BITS);
        l[n++] = id;
    }
    ow->wsl[0] = n;
    if(ow->ws_n < 3) ow->ws_n++;
    d1 = ow->ws_n >= 2 ? multistage_ws_diff(ow, 0, 1) : 0;
    d2 = ow->ws_n == 3 ? multistage_ws_diff(ow, 0, 2) : 0;
    // The working set of the previous solve shifted by one stage (multistage_ws_shift): used at a shift if it fit the
    // latest solution better than the previous one did (fewer differing constraints)
    // (over the constraints in either working set: the others are inactive in both)
    for(i = 0, ds = 0; ow->ws_n >= 2 && i < n; i++){
        const int src = multistage_ws_src(ow, l[i]);
        ds += c[l[i]] != ow->wsc[1][src >= 0 ? src : l[i]];
    }
    for(i = 0; ow->ws_n >= 2 && i < ow->wsl[1]; i++){ // (inactive now, and shifted from, or kept at, one active before)
        const int j = ow->wsh[1][i], d = ow->wsdst[j];
        ds += c[j] == 0 && multistage_ws_src(ow, j) < 0;
        ds += d >= 0 && c[d] == 0 && multistage_ws_src(ow, d) == j;
    }
    ow->ws_shift = ow->ws_n >= 2 && ds < d1;
    alt = ow->ws_n == 3 && 2*d2 <= d1;
    ow->ws_two = alt && ow->ws_alt && d1 >= 2;
    ow->ws_alt = alt;
}

// At an update: the working set of two solves back is marked if chosen (returns 1), unless the working set is no longer
// that of the latest solve (changed by the caller, e.g. a cold start)
static int multistage_ws_apply(DAQPMultistageWork *ow){
    DAQPWorkspace *work = ow->work;
    int i, n = 0;
    if(!ow->ws_two || work->bnb != NULL) return 0;
    ow->ws_two = 0;
    for(i = 0; i < ow->m; i++) n += multistage_ws_mutable(ow, i) && (work->sense[i] & DAQP_ACTIVE);
    if(n != ow->wsl[0]) return 0;
    for(i = 0; i < ow->wsl[0]; i++){
        const int id = ow->wsh[0][i];
        if((work->sense[id] & DAQP_MULTISTAGE_WS_BITS) != ow->wsc[0][id]) return 0;
    }
    for(i = 0; i < ow->wsl[0]; i++) work->sense[ow->wsh[0][i]] &= ~DAQP_MULTISTAGE_WS_BITS;
    for(i = 0; i < ow->wsl[1]; i++){
        const int id = ow->wsh[1][i];
        if(multistage_ws_mutable(ow, id)) work->sense[id] = (work->sense[id] & ~DAQP_MULTISTAGE_WS_BITS) | ow->wsc[1][id];
    }
    return 1;
}

// The working set and the multipliers (lamv) of the latest solve one stage earlier (in place: stage t+1 is read before
// it is set)
static void multistage_ws_shift(DAQPMultistageWork *ow){
    int *sense = ow->work->sense;
    int d;
    for(d = 0; d < ow->m; d++){
        const int src = multistage_ws_src(ow, d);
        if(src < 0) continue;
        sense[d] = (sense[d] & ~DAQP_MULTISTAGE_WS_BITS) | (sense[src] & DAQP_MULTISTAGE_WS_BITS);
        if(!(sense[d] & DAQP_SOFT)) sense[d] &= ~DAQP_SLACK_FIXED;
        ow->lamv[d] = ow->lamv[src];
    }
}

// After a warm start is activated: the multipliers of the latest solve (lamv) instead of the ones that the
// activation chose (1, 0 or w), where they have the sign of the side of the constraint, in the range of its slack
// state. The dual iteration then starts where the latest solve ended, instead of at the bounds of these ranges
// (where it would leave a free slack, or a constraint, at a step of length zero)
static void multistage_ws_multipliers(DAQPMultistageWork *ow){
    DAQPWorkspace *work = ow->work;
    int i;
    for(i = 0; i < work->n_active; i++){
        const int id = work->WS[i], lower = DAQP_IS_LOWER(id);
        c_float l = ow->lamv[id];
        if(!multistage_ws_mutable(ow, id) || (lower ? l >= 0 : l <= 0)) continue;
        l = fabs(l)/ow->scale[id]; // (internal multipliers are of normalized rows)
        if(DAQP_IS_SOFT(id)){
            const c_float w = daqp_soft_w(work, id);
            if(DAQP_IS_SLACK_FIXED(id) ? l > w : l < w) l = w;
        }
        work->lam[i] = lower ? -l : l;
    }
}

/* Setup and solve */

// Whether bounds l and u are equal and finite
static int multistage_eq_bounds(const c_float l, const c_float u){
    return l == u && l < DAQP_INF && l > -DAQP_INF;
}

// Whether constraint id is a hard equality constraint (after multistage_fetch_bounds)
static int multistage_is_equality(const DAQPMultistageWork *ow, const int id){
    return !(ow->work->sense[id] & DAQP_SOFT) && multistage_eq_bounds(ow->work->dlower[id], ow->work->dupper[id]);
}

// Mark equalities (eliminated: immutable, inactive; others: immutable, active). Returns 1 if changed
static int multistage_mark_equalities(DAQPMultistageWork *ow){
    const c_float *lo = ow->work->dlower, *up = ow->work->dupper;
    const int i0 = ow->x0free ? 0 : ow->st[0].nb; // (the bounds on a fixed x_0 alone are among the first ones)
    int *sense = ow->work->sense;
    int i, changed = 0, any = 0;
    // Without equalities (as is typical), nothing changes unless bounds are now equal
    if(ow->neq_now == 0){
        for(i = 0; i < i0; i++) any |= !(ow->eqf[i] & DAQP_MULTISTAGE_EQ_X0) && lo[i] == up[i];
        for(i = i0; i < ow->m; i++) any |= lo[i] == up[i];
        if(!any) return 0;
    }
    ow->neq_now = 0;
    for(i = 0; i < ow->m; i++){
        int f = ow->eqf[i] & DAQP_MULTISTAGE_EQ_CAND;
        if(ow->eqf[i] & DAQP_MULTISTAGE_EQ_X0) continue;
        if(multistage_is_equality(ow, i)){
            f |= DAQP_MULTISTAGE_EQ_NOW;
            if((f & DAQP_MULTISTAGE_EQ_CAND) && !ow->no_elim){
                f |= DAQP_MULTISTAGE_EQ_ELIM;
                sense[i] = (sense[i] & ~DAQP_ACTIVE) | DAQP_IMMUTABLE;
            }
            // (one the solver has made mutable, as dependent on the others, is left so: only a new
            // equality is marked active, since only a changed set of equalities rebuilds the working set)
            else if(!(ow->eqf[i] & DAQP_MULTISTAGE_EQ_NOW) || (sense[i] & DAQP_IMMUTABLE))
                sense[i] |= DAQP_ACTIVE | DAQP_IMMUTABLE;
        }
        else if(ow->eqf[i] & DAQP_MULTISTAGE_EQ_NOW) // No longer an equality
            sense[i] &= ~(DAQP_ACTIVE | DAQP_IMMUTABLE);
        if((f ^ ow->eqf[i]) & DAQP_MULTISTAGE_EQ_NOW) changed = 1;
        ow->neq_now += (f & DAQP_MULTISTAGE_EQ_NOW) != 0;
        ow->eqf[i] = f;
    }
    return changed;
}

// The data of stage t of the backward recursion has changed
static void multistage_dirty(DAQPMultistageWork *ow, const int t){
    if(t > ow->bw_from) ow->bw_from = t;
}

// Mark the stages where the value of an eliminated equality has changed (eqval: the values of the latest update)
static void multistage_eq_values_changed(DAQPMultistageWork *ow){
    int t, k, i = 0;
    for(t = 0; t <= ow->N; t++)
        for(k = 0; k < ow->st[t].ncand; k++, i++){
            const c_float val = ow->work->dlower[ow->st[t].cand[k]];
            if((ow->eqf[ow->st[t].cand[k]] & DAQP_MULTISTAGE_EQ_ELIM) && val != ow->eqval[i]) multistage_dirty(ow, t);
            ow->eqval[i] = val;
        }
}

// Keep the equalities in the working set (no_elim), or eliminate them
static void multistage_set_elimination(DAQPMultistageWork *ow, const int no_elim){
    ow->no_elim = no_elim;
    multistage_mark_equalities(ow);
}

// Classical recursion with eps I added (prox): eps = max(|eps_prox|, sqrt(zero_tol)*hmax)
static int multistage_factor_prox(DAQPMultistageWork *ow){
    const DAQPSettings *set = ow->work->settings;
    c_float hmax = 0, eps = fabs(set->eps_prox);
    int t, i;
    for(t = 0; t <= ow->N; t++)
        for(i = 0; i < ow->st[t].nu + ow->st[t].nx; i++){
            const c_float h = fabs(BLASFEO_DMATEL(ow->st[t].RSQ, i, i));
            if(h > hmax) hmax = h;
        }
    if(eps < sqrt(set->zero_tol)*hmax) eps = sqrt(set->zero_tol)*hmax;
    ow->eps = eps;
    ow->ric_sqrt = 0;
    if(ow->neqc > 0 && !ow->no_elim) multistage_set_elimination(ow, 1);
    return multistage_factor(ow);
}

// Factorization: the elimination is retried each time (else equalities stay in the working set);
// prox if the recursion breaks down (eps_prox < 0) or always (eps_prox > 0). -1 on breakdown
static int multistage_factorize(DAQPMultistageWork *ow){
    const c_float eps_prox = ow->work->settings->eps_prox, eps_old = ow->eps;
    const int no_elim_old = ow->no_elim;
    int flag = -1, t;
    ow->eps = 0;
    ow->ric_sqrt = 1; // (the classical recursion only where the square-root one breaks down, at each factorization)
    for(t = 0; t <= ow->N; t++){
        if(ow->st[t].ng > 0) multistage_gblocks(ow->st+t); // (DCt may have changed)
        ow->st[t].Pb_valid = 0;
    }
    if(eps_prox > 0) flag = multistage_factor_prox(ow);
    else{
        if(ow->neqc > 0 && ow->no_elim != ow->elim_bad) multistage_set_elimination(ow, ow->elim_bad); // (retry, or not)
        for(;;){
            ow->ill_cond = 0;
            if((flag = multistage_factor(ow)) < 0 && ow->ric_sqrt){ ow->ric_sqrt = 0; flag = multistage_factor(ow); }
            if(ow->neqc == 0 || ow->no_elim || (flag >= 0 && !ow->ill_cond)) break;
            multistage_set_elimination(ow, 1);
        }
        if(flag < 0 && eps_prox < 0) flag = multistage_factor_prox(ow);
    }
    // Other rows: renormalize
    if((ow->eps > 0) != (eps_old > 0) || ow->no_elim != no_elim_old) ow->norm_stale = 1;
    // Offsets of the free inputs in v
    for(t = 0, ow->nv = 0; t <= ow->N; t++){ ow->st[t].cv = ow->nv; ow->nv += ow->st[t].nn; }
    ow->work->n = ow->nv;
    ow->work->n_prox = ow->eps > 0 ? ow->nv : 0;
    ow->op.prox_eps = ow->eps;
    if(flag >= 0 && ow->norm_stale) multistage_normalize(ow);
    // (checked by the dual pass of the next solve: daqp_multistage_solve)
    for(t = 0, ow->elim_check = 0; flag >= 0 && !ow->no_elim && ow->eps == 0 && t <= ow->N; t++) ow->elim_check |= ow->st[t].r > 0;
    return flag;
}

// v[vi..vi+n) <-- src (NULL: zero). Returns 1 if it changed
static int multistage_sync_vec(const int n, const c_float *src, struct blasfeo_dvec *v, const int vi){
    c_float *d = v->pa + vi;
    int i, changed = 0;
    if(src != NULL){
        for(i = 0; i < n; i++) changed |= d[i] != src[i]; // (no early exit: vectorized)
        if(changed) memcpy(d, src, n*sizeof(c_float));
    }
    else{
        for(i = 0; i < n; i++) changed |= d[i] != 0;
        if(changed) memset(d, 0, n*sizeof(c_float));
    }
    return changed;
}

// Stage data from the user's problem: pointers, and copies of b, q, r
static void multistage_bind(DAQPMultistageWork *ow){
    const DAQPMultistage *u = ow->user;
    int t;
    for(t = ow->x0free; t <= ow->N; t++){
        Stage *S = ow->st+t;
        const int tu = t - ow->x0free; // Stage of the user
        S->BAbt = t < ow->N ? u->BAbt + tu : NULL;
        S->RSQ = u->RSQrq + tu;
        S->DCt = S->ng > 0 ? u->DCt + tu : NULL;
        S->Et = S->nc > 0 ? u->Et + tu : NULL;
        S->Ft = S->nc > 0 ? u->Ft + tu : NULL;
        S->idxb = S->nb > 0 ? u->idxb[tu] : NULL;
        if(S->nb > 0){ S->lb = u->lb[tu]; S->ub = u->ub[tu]; }
        if(S->ng > 0){ S->lg = u->lg[tu]; S->ug = u->ug[tu]; }
        if(S->nc > 0){ S->lc = u->lc[tu]; S->uc = u->uc[tu]; }
        // NULL is zero; b = 0 is skipped
        if(multistage_sync_vec(S->nu, u->r != NULL ? u->r[tu] : NULL, &S->rq, 0) |
                multistage_sync_vec(S->nx, u->q != NULL ? u->q[tu] : NULL, &S->rq, S->nu)) multistage_dirty(ow, t);
        if(t < ow->N){
            const c_float *b = u->b != NULL ? u->b[tu] : NULL;
            int i, has_b = 0;
            for(i = 0; b != NULL && i < S[1].nx; i++) if(b[i] != 0) has_b = 1;
            int changed = has_b != S->has_b;
            if(has_b) changed |= multistage_sync_vec(S[1].nx, b, &S->b, 0);
            S->has_b = has_b;
            if(changed){ multistage_dirty(ow, t); S->Pb_valid = 0; }
        }
    }
    if(ow->x0free){ // The stage of x_0
        ow->st[0].BAbt = &ow->B0;
        ow->st[0].RSQ = &ow->R0;
    }
    // The reduced stages are formed by multistage_eq_stage
    for(t = 0; t <= ow->N && ow->neqc == 0; t++){ ow->st[t].BAr = ow->st[t].BAbt; ow->st[t].RSr = ow->st[t].RSQ; }
}

// Assign the memory (or only count it if M->base is NULL)
static void multistage_memory(DAQPMultistageWork *ow, MultistageMem *M){
    const int m = ow->m, nv = ow->nvmax;
    const int nxmax = ow->nxmax, numax = ow->numax, nmax = nxmax + numax;
    int t, qall = 0;
    for(t = 0; t <= ow->N && ow->neqc > 0; t++){
        const int qm = ow->st[t].ncand + (t < ow->N ? ow->st[t+1].nx : 0);
        if(qm > qall) qall = qm;
    }
    for(t = 0; t <= ow->N; t++){
        Stage *S = ow->st+t;
        const int n = S->nu + S->nx, nx1 = t < ow->N ? S[1].nx : 0, qm = S->ncand + nx1;
        multistage_new_mat(M, &S->D, n, n);
        multistage_new_vec(M, &S->z, n);
        multistage_new_vec(M, &S->dz, n);
        multistage_new_vec(M, &S->yp, n);
        multistage_new_vec(M, &S->zc, n);
        multistage_new_vec(M, &S->rq, n);
        multistage_new_vec(M, &S->b, nx1);
        multistage_new_vec(M, &S->Pb, nx1);
        S->cand = multistage_take(M, S->ncand*sizeof(int));
        S->gblk = multistage_take(M, 4*S->ng*sizeof(int));
        if(ow->neqc == 0) continue;
        S->h = multistage_take(M, S->nx*sizeof(c_float));
        multistage_new_mat(M, &S->Qd, S->nu, S->nu);
        multistage_new_mat(M, &S->Kt, S->nx, S->nu);
        multistage_new_mat(M, &S->BAm, n, nx1);
        multistage_new_vec(M, &S->u0, S->nu);
        S->Yg = multistage_take(M, S->nu*S->nu*sizeof(c_float));
        S->Z = multistage_take(M, S->nu*qm*sizeof(c_float));
        S->Z2 = multistage_take(M, S->nx*qm*sizeof(c_float));
        multistage_new_mat(M, &S->Ht, S->nx, S->nx);
        S->piv = multistage_take(M, qm*sizeof(int));
        S->piv2 = multistage_take(M, qm*sizeof(int));
        S->om = multistage_take(M, qm*sizeof(c_float));
    }
    if(ow->x0free){
        multistage_new_mat(M, &ow->B0, ow->st[0].nu, ow->st[1].nx);
        multistage_new_mat(M, &ow->R0, ow->st[0].nu, ow->st[0].nu);
    }
    multistage_new_mat(M, &ow->sT, nmax, nxmax);
    multistage_new_mat(M, &ow->sH, nmax, nmax);
    multistage_new_mat(M, &ow->sAt, nxmax, nmax);
    // Rows of M (multistage_rows_block) and their normalization (multistage_normalize)
    multistage_new_mat(M, &ow->sCw, nxmax, numax);
    multistage_new_mat(M, &ow->sG, nxmax, nxmax);
    multistage_new_mat(M, &ow->sSig[0], nxmax, nxmax);
    multistage_new_mat(M, &ow->sSig[1], nxmax, nxmax);
    multistage_new_mat(M, &ow->sC0, DAQP_MULTISTAGE_BLOCK, nmax);
    multistage_new_mat(M, &ow->sC1, DAQP_MULTISTAGE_BLOCK, nmax);
    multistage_new_mat(M, &ow->sV, DAQP_MULTISTAGE_BLOCK, numax);
    multistage_new_mat(M, &ow->sK, DAQP_MULTISTAGE_BLOCK, DAQP_MULTISTAGE_BLOCK);
    if(ow->neqc > 0){
        multistage_new_mat(M, &ow->sC2, DAQP_MULTISTAGE_BLOCK, nmax);
        multistage_new_mat(M, &ow->sX, nmax, nmax);
        multistage_new_mat(M, &ow->sRS, nmax, nmax);
        multistage_new_mat(M, &ow->sY, numax, numax);
        multistage_new_mat(M, &ow->sQv, numax, numax);
        multistage_new_mat(M, &ow->sQw, numax, numax);
        multistage_new_mat(M, &ow->sV2, nxmax, numax);
        multistage_new_mat(M, &ow->sB, numax, nxmax);
        multistage_new_mat(M, &ow->sF, nmax, qall);
        ow->sGe = multistage_take(M, nmax*qall*sizeof(c_float));
        ow->sW2 = multistage_take(M, nxmax*qall*sizeof(c_float));
        ow->sTau = multistage_take(M, (qall+nmax)*sizeof(c_float));
        ow->sg = multistage_take(M, qall*sizeof(c_float));
        ow->sy = multistage_take(M, qall*sizeof(c_float));
        ow->eqval = multistage_take(M, ow->neqc*sizeof(c_float));
        ow->sng = multistage_take(M, ow->ngmax*sizeof(c_float));
    }
    multistage_new_vec(M, &ow->ws, nxmax); multistage_new_vec(M, &ow->wb, nxmax);
    multistage_new_vec(M, &ow->wg, nmax); multistage_new_vec(M, &ow->wc0, nmax);
    multistage_new_vec(M, &ow->wc1, nmax); multistage_new_vec(M, &ow->wr, nmax); multistage_new_vec(M, &ow->wt, nmax);
    ow->ct = multistage_take(M, m*sizeof(int));
    ow->wssrc = multistage_take(M, m*sizeof(int));
    ow->wsdst = multistage_take(M, m*sizeof(int));
    ow->eqf = multistage_take(M, m*sizeof(int));
    for(t = 0; t < 3; t++){
        ow->wsc[t] = multistage_take(M, m*sizeof(unsigned char));
        ow->wsh[t] = multistage_take(M, (m < nv + ow->nsmax + 1 ? m : nv + ow->nsmax + 1)*sizeof(int));
    }
    ow->id_slot = multistage_take(M, m*sizeof(int));
    ow->slot_id = multistage_take(M, ow->nslots*sizeof(int));
    ow->rows_len = multistage_take(M, ow->nslots*sizeof(int));
    ow->ids = multistage_take(M, ow->nslots*sizeof(int));
    ow->smax = multistage_take(M, ow->nslots*sizeof(int));
    ow->scale = multistage_take(M, m*sizeof(c_float));
    ow->zrow = multistage_take(M, m*sizeof(char));
    ow->lamv = multistage_take(M, m*sizeof(c_float));
    ow->pic = multistage_take(M, ow->npi*sizeof(c_float));
    ow->vlast = multistage_take(M, nv*sizeof(c_float));
    ow->bside = multistage_take(M, m*sizeof(unsigned char));
    ow->bsave = multistage_take(M, m*sizeof(c_float));
    ow->vcand = multistage_take(M, nv*sizeof(c_float));
    ow->dupper = multistage_take(M, m*sizeof(c_float));
    ow->dlower = multistage_take(M, m*sizeof(c_float));
    ow->Mu = multistage_take(M, m*sizeof(c_float));
    ow->sense = multistage_take(M, m*sizeof(int));
    // (multistage_activate packs blocks of at most min(nv, DAQP_MULTISTAGE_ACT_COLS) columns, and panels of rows of L)
    multistage_new_mat(M, &ow->sAb, nv + ow->nsmax,
            nv < DAQP_MULTISTAGE_ACT_COLS ? (nv > D_PS ? nv : D_PS) : DAQP_MULTISTAGE_ACT_COLS);
    ow->sic = multistage_take(M, (nv + ow->nsmax + 1)*sizeof(c_float));
    ow->rows = multistage_take(M, (size_t)ow->nslots*nv*sizeof(c_float));
}

// x_0 from the bounds that fix it; -1 if a component is no longer fixed
static int multistage_fetch_x0(DAQPMultistageWork *ow){
    int i;
    if(ow->x0free) return 0;
    for(i = 0; i < ow->st[0].nx; i++){
        const int id = ow->x0src[i];
        if(!multistage_is_equality(ow, id)) return -1;
        ow->x0[i] = ow->work->dlower[id];
    }
    return 0;
}

// Whether a bound on the fixed x_0 is violated (infeasible)
static int multistage_x0_infeasible(DAQPMultistageWork *ow){
    const Stage *S0 = ow->st;
    int i;
    for(i = 0; i < S0->nb; i++){
        c_float x, tol;
        if(!(ow->eqf[i] & DAQP_MULTISTAGE_EQ_X0)) continue;
        x = ow->x0[S0->idxb[i]-S0->nu];
        tol = ow->work->settings->primal_tol*(1+fabs(x));
        if(x > S0->ub[i]+tol || x < S0->lb[i]-tol) return 1;
    }
    return 0;
}

// Whether user constraint j of stage tu (id overall) can be eliminated: hard, equal finite
// bounds at setup, and not a bound on a fixed x_0
static int multistage_is_candidate(const DAQPMultistage *user, const int x0free, const int tu, const int j, const int id){
    const int nb = DAQP_MULTISTAGE_NB(user,tu);
    if(!multistage_eq_bounds(j < nb ? user->lb[tu][j] : user->lg[tu][j-nb], j < nb ? user->ub[tu][j] : user->ug[tu][j-nb]))
        return 0;
    if(user->sense != NULL && (user->sense[id] & DAQP_SOFT)) return 0;
    return !(tu == 0 && !x0free && j < nb && user->idxb[0][j] >= DAQP_MULTISTAGE_NU(user,0));
}

// x_0 is fixed if the hard bounds at stage 0 fix all of its components (x0src: the bound that fixes each)
static int multistage_x0free(const DAQPMultistage *user, int *x0src){
    const int nu0 = DAQP_MULTISTAGE_NU(user,0);
    int i, j, x0free = 0;
    for(j = 0; j < user->nx[0]; j++){
        int src = -1;
        for(i = 0; i < DAQP_MULTISTAGE_NB(user,0) && src < 0; i++)
            if(user->idxb[0][i]-nu0 == j && !(user->sense != NULL && (user->sense[i] & DAQP_SOFT)) &&
                    multistage_eq_bounds(user->lb[0][i], user->ub[0][i])) src = i;
        if(x0src != NULL) x0src[j] = src;
        x0free |= src < 0;
    }
    return x0free;
}

// Stages (with that of x_0 first if x0free) and candidates for elimination; worst: from the dimensions
// alone (every bound and general constraint a candidate). ns: number of soft constraints
static void multistage_dims(DAQPMultistageWork *ow, const DAQPMultistage *user, const int worst, const int ns){
    int t, i;
    ow->N = user->N + ow->x0free;
    if(ow->x0free) ow->st[0].nu = ow->st[0].nn = user->nx[0];
    for(t = 0; t <= ow->N; t++){
        Stage *S = ow->st+t;
        const int tu = t - ow->x0free;
        if(tu >= 0){
            S->nx = user->nx[tu]; S->nu = S->nn = DAQP_MULTISTAGE_NU(user,tu);
            S->nb = DAQP_MULTISTAGE_NB(user,tu); S->ng = DAQP_MULTISTAGE_NG(user,tu);
            S->nc = DAQP_MULTISTAGE_NC(user,tu);
            for(i = 0; i < S->nb+S->ng; i++) S->ncand += worst || multistage_is_candidate(user, ow->x0free, tu, i, ow->m+i);
        }
        S->cv = ow->nvmax; S->c0 = ow->m;
        if(S->nx > ow->nxmax) ow->nxmax = S->nx;
        if(S->nu > ow->numax) ow->numax = S->nu;
        if(S->ng > ow->ngmax) ow->ngmax = S->ng;
        ow->nvmax += S->nu;
        ow->m += S->nb + S->ng + S->nc;
        if(t > ow->x0free) ow->npi += S->nx;
        ow->neqc += S->ncand;
    }
    ow->nv = ow->nvmax;
    ow->nsmax = ns;
    ow->nslots = ow->nvmax + ns + 2; // One more than the working set can hold
}

// The head of the block: the workspace, x_0 and N+2 stages (zeroed; only counted if M->base is NULL)
static DAQPMultistageWork* multistage_head(MultistageMem *M, const DAQPMultistage *user){
    const int nx0 = user->nx[0];
    DAQPMultistageWork *ow = multistage_take(M, sizeof(DAQPMultistageWork));
    int *x0src = multistage_take(M, (nx0+1)*sizeof(int));
    c_float *x0 = multistage_take(M, (nx0+1)*sizeof(c_float));
    Stage *st = multistage_take(M, (user->N+2)*sizeof(Stage));
    if(ow == NULL) return NULL;
    memset(ow, 0, sizeof(*ow));
    memset(st, 0, (user->N+2)*sizeof(Stage));
    ow->user = (DAQPMultistage*)user; ow->x0src = x0src; ow->x0 = x0; ow->st = st;
    return ow;
}

// Bytes of the block: the head and the data
static size_t multistage_block_size(DAQPMultistageWork *ow){
    MultistageMem M = {NULL, 0};
    multistage_head(&M, ow->user);
    multistage_memory(ow, &M);
    return M.n + DAQP_MULTISTAGE_ALIGN;
}

// Bytes of the block of prob, from its dimensions (0: out of memory)
static size_t multistage_size(const DAQPMultistage *prob, const int x0free, const int worst, const int ns){
    DAQPMultistageWork ow;
    size_t size;
    memset(&ow, 0, sizeof(ow));
    ow.user = (DAQPMultistage*)prob;
    ow.x0free = x0free;
    if((ow.st = calloc(prob->N+2, sizeof(Stage))) == NULL) return 0;
    multistage_dims(&ow, prob, worst, ns);
    size = multistage_block_size(&ow);
    free(ow.st);
    return size;
}

size_t daqp_multistage_memsize(const DAQPMultistage *prob, const int ns){
    return multistage_size(prob, prob->nx[0] > 0, 1, ns); // (x_0 may be free)
}

int daqp_multistage_nv(const DAQPMultistage *prob){
    return daqp_multistage_nvmax(prob) - (multistage_x0free(prob, NULL) ? 0 : prob->nx[0]);
}

int daqp_multistage_nvmax(const DAQPMultistage *prob){
    int t, n = prob->nx[0];
    for(t = 0; t < prob->N; t++) n += prob->nu[t];
    return n;
}

// Set up in mem (size bytes), or in a block that is allocated (mem NULL)
static int multistage_setup(DAQPMultistage *user, DAQPWorkspace *work, void *mem, size_t size){
    const int nx0 = user->nx[0], nu0 = DAQP_MULTISTAGE_NU(user,0), own = mem == NULL;
    const int own_settings = work->settings == NULL, m = daqp_multistage_m(user);
    int t, i, k, ns = 0;
    DAQPMultistageWork *ow;
    MultistageMem M;

    // The BLASFEO matrices are required (daqp_multistage_pack forms them)
    if(user->RSQrq == NULL || (user->N > 0 && user->BAbt == NULL)) return DAQP_EXIT_UNSUPPORTED;
    for(t = 0; t <= user->N; t++)
        if((DAQP_MULTISTAGE_NG(user,t) > 0 && user->DCt == NULL) ||
                (DAQP_MULTISTAGE_NC(user,t) > 0 && (user->Et == NULL || user->Ft == NULL))) return DAQP_EXIT_UNSUPPORTED;
    // (given memory: also the settings and the arrays of the core workspace)
    if(!own && work->settings == NULL) return DAQP_EXIT_UNSUPPORTED;
    if(user->sense != NULL)
        for(i = 0; i < m; i++) if(user->sense[i] & DAQP_SOFT) ns++;

    // A block for the dimensions of the problem
    if(own && ((size = multistage_size(user, multistage_x0free(user, NULL), 0, ns)) == 0 || (mem = malloc(size)) == NULL))
        return DAQP_EXIT_UNSUPPORTED;
    // The head first (its size from the dimensions)
    M.base = NULL; M.n = 0;
    multistage_head(&M, user);
    if(M.n + DAQP_MULTISTAGE_ALIGN > size) return DAQP_EXIT_UNSUPPORTED; // (given memory too small)
    M.base = (char*)(((size_t)mem + DAQP_MULTISTAGE_ALIGN-1) & ~(size_t)(DAQP_MULTISTAGE_ALIGN-1));
    M.n = 0;
    ow = multistage_head(&M, user);
    ow->work = work;
    ow->mem = mem;
    ow->own_mem = own;
    ow->x0free = multistage_x0free(user, ow->x0src);
    multistage_dims(ow, user, 0, ns);
    if(!own && multistage_block_size(ow) > size) return DAQP_EXIT_UNSUPPORTED; // (given memory too small)
    // The data. Zeroed (BLASFEO reads y for beta = 0; R0 = 0; the history of working sets; bside), except the row cache
    {
        char *data = M.base + M.n;
        multistage_memory(ow, &M);
        memset(data, 0, (char*)ow->rows - data);
    }
    for(t = 0; t <= ow->N; t++){
        Stage *S = ow->st+t;
        const int tu = t - ow->x0free;
        for(i = S->c0, k = 0; i < S->c0 + S->nb + S->ng + S->nc; i++){
            ow->ct[i] = t;
            if(tu < 0 || i - S->c0 >= S->nb + S->ng || !multistage_is_candidate(user, ow->x0free, tu, i - S->c0, i)) continue;
            S->cand[k++] = i;
            ow->eqf[i] = DAQP_MULTISTAGE_EQ_CAND;
        }
    }
    for(i = 0; ow->x0free && i < nx0; i++) blasfeo_dgein1(1.0, &ow->B0, i, i);
    multistage_ws_shift_map(ow);

    // DAQP workspace for the LDP in v (M is implicit)
    if(own){
        allocate_daqp_settings(work);
        allocate_daqp_workspace(work, ow->nvmax, ns);
        free(work->prox_mask); // (unused: the proximal iterations of Mop shift all directions; xold holds the incumbent of BnB)
        work->prox_mask = NULL;
    }
    else{ // (the arrays lam, lam_star, WS, L, D, xldl, zldl, u, the soft weights, and xold with binaries are given)
        init_daqp_workspace(work, ow->nvmax);
        work->prox_mask = NULL;
    }
    work->m = ow->m;
    work->ms = 0;
    work->n = ow->nv;
    work->qp = NULL;
    work->M = NULL;
    work->dupper = ow->dupper;
    work->dlower = ow->dlower;
    work->Mu = ow->Mu;
    work->sense = ow->sense;
    for(i = 0; i < ow->m; i++){
        work->sense[i] = user->sense != NULL ? user->sense[i] : 0;
        if(user->binary_horizon > 0 && ow->ct[i] - ow->x0free >= user->binary_horizon) work->sense[i] &= ~DAQP_BINARY;
    }
    work->scaling = ow->scale; // Row normalization
    ow->op.data = ow;
    ow->op.row = multistage_row;
    ow->op.product = multistage_product;
    ow->op.prox_center = multistage_prox_center;
    ow->op.prox_primal = multistage_prox_primal;
    work->Mop = &ow->op;
    // Binary constraints: branch and bound (sized for nvmax; allocated, also with given memory)
    for(i = 0, k = 0; i < ow->m; i++) k += (work->sense[i] & DAQP_BINARY) != 0;
    if(!own && k == 0) work->xold = NULL;
    k = k == 0 ? 1 : !own && work->xold == NULL ? DAQP_EXIT_UNSUPPORTED : setup_daqp_bnb(work, work->sense, k, ns, work->m);

    // The bounds on a fixed x_0 alone are left out of the LDP
    for(i = 0; !ow->x0free && i < DAQP_MULTISTAGE_NB(user,0); i++)
        if(user->idxb[0][i] >= nu0){
            ow->nx0c++;
            ow->eqf[i] = DAQP_MULTISTAGE_EQ_X0;
            work->sense[i] = (work->sense[i] & ~DAQP_ACTIVE) | DAQP_IMMUTABLE;
        }

    // The data: a full update
    ow->neq_now = -1;
    ow->norm_stale = 1; // First factorization
    i = k < 0 ? k : daqp_multistage_update(work, DAQP_MULTISTAGE_UPDATE_MATRICES);
    // As setup_daqp: freed on failure (but not settings given by the caller)
    if(i < 0){
        if(!own_settings) work->settings = NULL;
        free_daqp_multistage(work);
    }
    return i;
}
int setup_daqp_multistage(DAQPMultistage *user, DAQPWorkspace *work){
    return multistage_setup(user, work, NULL, 0);
}
int setup_daqp_multistage_mem(DAQPMultistage *user, DAQPWorkspace *work, void *mem, size_t size){
    return mem == NULL ? DAQP_EXIT_UNSUPPORTED : multistage_setup(user, work, mem, size);
}

int daqp_multistage_update(DAQPWorkspace *work, const int mask){
    DAQPMultistageWork *ow = (DAQPMultistageWork*)work->Mop->data;
    int refactor = mask & DAQP_MULTISTAGE_UPDATE_MATRICES, reset, ws_switch, flag;
    ow->pi_valid = ow->lam_valid = 0; // Of the previous problem
    ow->bw_from = -1;
    multistage_bind(ow);
    multistage_fetch_bounds(ow);
    if(multistage_fetch_x0(ow) < 0) return DAQP_EXIT_UNSUPPORTED; // Requires a new setup
    // A changed set of equalities resets the working set (and the parametrization)
    reset = multistage_mark_equalities(ow);
    if(reset) ow->elim_bad = 0; // (other equalities: eliminate them again)
    if(reset && ow->neqc > 0) refactor = ow->norm_stale = 1;
    if(ow->elim_bad && !ow->no_elim) refactor = 1; // Found inaccurate by the dual pass: keep the equalities in the working set
    // Warm start: the working set of the latest solve, or of two solves back (multistage_ws_record); other equalities
    // make the history of other problems
    if(reset){
        multistage_ws_forget(ow);
        memset(ow->bside, 0, ow->m); // (and the binary assignment of branch and bound)
    }
    ws_switch = multistage_ws_apply(ow);
    // A shifted problem: the working set and the multipliers one stage earlier (if that fit the latest solutions
    // better, and unless they alternate)
    if(!ws_switch && (mask & DAQP_MULTISTAGE_UPDATE_SHIFT) && !reset && work->bnb == NULL && ow->ws_shift){
        multistage_ws_shift(ow);
        ws_switch = 1;
    }
    if(refactor){
        if(multistage_factorize(ow) < 0) return DAQP_EXIT_NONCONVEX;
        multistage_cache_clear(ow); // All rows of M change
    }
    if(refactor || reset || ws_switch) reset_daqp_workspace(work);
    // Backward recursion from the last stage where b, q, r or eliminated values changed
    if(ow->neqc > 0) multistage_eq_values_changed(ow);
    if(refactor) ow->bw_from = ow->N;
    if(ow->bw_from >= 0) multistage_lqr_backward(ow, ow->bw_from);
    multistage_lqr_forward(ow);
    work->reuse_ind = 0; // The right-hand side of the working set changed
    if(!(refactor || reset || ws_switch)) return 1;
    // Warm start: the constraints marked active (soft ones with their slack state), at once from their rows (formed
    // in blocks), with the multipliers of the latest solve
    multistage_rows_active(ow);
    flag = multistage_activate(ow);
    if(work->bnb == NULL) multistage_ws_multipliers(ow);
    return flag;
}

// Whether any of v[0..n) is nonzero
static int multistage_any(const c_float *v, const int n){
    int i;
    for(i = 0; i < n; i++) if(v[i] != 0) return 1;
    return 0;
}

// g = RSQ [u; x] + [r; q] + sum lam_i a_i at stage s for the solution (in dz)
static void multistage_grad(DAQPMultistageWork *ow, const int s, c_float *lam, struct blasfeo_dvec *g){
    Stage *S = ow->st+s;
    const int n = S->nu + S->nx;
    struct blasfeo_dvec sl;
    int i;
    if(n == 0) return;
    blasfeo_create_dvec(ow->m, &sl, lam);
    blasfeo_dsymv_l(n, 1.0, S->RSQ, 0, 0, &S->dz, 0, 1.0, &S->rq, 0, g, 0);
    for(i = 0; i < S->nb; i++) g->pa[S->idxb[i]] += lam[S->c0+i];
    // Skip blocks without nonzero multipliers
    if(multistage_any(lam + S->c0+S->nb, S->ng)) multistage_dct_n(S, &sl, S->c0+S->nb, 0, n, g);
    if(multistage_any(lam + S->c0+S->nb+S->ng, S->nc))
        blasfeo_dgemv_n(n, S->nc, 1.0, S->Et, 0, 0, &sl, S->c0+S->nb+S->ng, 1.0, g, 0, g, 0);
    if(s > 0 && multistage_any(lam + S[-1].c0+S[-1].nb+S[-1].ng, S[-1].nc))
        blasfeo_dgemv_n(n, S[-1].nc, 1.0, S[-1].Ft, 0, 0, &sl, S[-1].c0+S[-1].nb+S[-1].ng, 1.0, g, 0, g, 0);
}

// Gives the costates in pic (pi_valid) and, if full, completes lamv (lam_done) (the eliminated equalities and the bounds that fix x_0)
// and checks the elimination, in one backward pass over the gradient g of the Lagrangian of the working set:
//   om_t = Yg'Q1'(-g_u - B'l_{t+1}),  l_t = g_x + A'(l_{t+1} + H_{t+1}'om_H) + Gx'om_own,
// then forwards: the multipliers w_t of H_t to the rows it combines, and pi_{t-1} = l_t + H_t'w_t.
// x_0 fixed: pi_{-1} = l_0 and the multipliers of the bounds that fix it are -pi_{-1}
static void multistage_duals(DAQPMultistageWork *ow, const int full){
    const int N = ow->N, s0 = ow->x0free, x0mult = full && ow->nx0c > 0; // Stage s is stage s-s0 of the user
    c_float *lam = ow->lamv, *p = ow->pic, *pn = NULL;
    struct blasfeo_dvec *g = &ow->wg, *dl = &ow->ws, ln;
    int t, i, k, t0, tlow;
    for(t0 = 0; t0 <= N && ow->st[t0].q == 0; t0++);
    tlow = x0mult ? 0 : s0+1;
    if(t0 < tlow) tlow = t0;
    for(t = s0+1; t <= N; t++) p += ow->st[t].nx;
    for(t = N; t >= tlow; t--){
        Stage *S = ow->st+t;
        const int nx = S->nx, nu = S->nu, q = S->q, r = S->r;
        const int mh = t < N ? S[1].rh : 0, ne = q-mh;
        multistage_grad(ow, t, lam, g);
        blasfeo_create_dvec(t < N ? S[1].nx : 0, &ln, pn);
        const c_float gu2 = r > 0 ? blasfeo_ddot(nu, g, 0, g, 0) : 0; // (|g_u|^2 before B'l, for the scale of the check)
        if(t < N && r > 0) blasfeo_dgemv_n(nu, S[1].nx, 1.0, S->BAbt, 0, 0, &ln, 0, 1.0, g, 0, g, 0); // B'l_{t+1} (for om)
        if(q > 0) memset(S->om, 0, q*sizeof(c_float));
        if(r > 0){
            // om = Yg'(Q1'(-g_u))
            struct blasfeo_dvec tv;
            blasfeo_create_dvec(r, &tv, ow->sTau);
            blasfeo_dgemv_t(nu, r, -1.0, &S->Qd, 0, 0, g, 0, 0.0, &tv, 0, &tv, 0);
            { // Stationarity in the free inputs: |Q2's|^2 = |s|^2 - |Q1's|^2 for s = g_u + B'l, relative to |g_u| + |B'l|
                const c_float ss = blasfeo_ddot(nu, g, 0, g, 0), res2 = ss - blasfeo_ddot(r, &tv, 0, &tv, 0);
                const c_float scale = 1 + 2*sqrt(gu2) + sqrt(ss); // (|B'l| <= |s| + |g_u|)
                if(full && res2 > DAQP_MULTISTAGE_ELIM_CHECK*DAQP_MULTISTAGE_ELIM_CHECK*scale*scale) ow->elim_bad = 1;
            }
            for(k = 0; k < r; k++){ // (to the rows piv[0..r))
                c_float sq = 0;
                for(i = 0; i < r; i++) sq += S->Yg[i*r+k]*ow->sTau[i];
                S->om[S->piv[k]] = sq;
            }
        }
        // l_t = g_x + A'(l_{t+1} + H_{t+1}'om_H) + Gx'om_own (one product with A')
        if(t < N){
            const int nx1 = S[1].nx;
            struct blasfeo_dvec *w = &ln;
            if(mh > 0){ // dl = l_{t+1} + H_{t+1}'om_H
                struct blasfeo_dvec oh;
                blasfeo_create_dvec(mh, &oh, S->om + ne);
                blasfeo_dgemv_n(nx1, mh, 1.0, &S[1].Ht, 0, 0, &oh, 0, 1.0, &ln, 0, dl, 0);
                w = dl;
            }
            blasfeo_dgemv_n(nx, nx1, 1.0, S->BAbt, nu, 0, w, 0, 1.0, g, nu, g, nu);
        }
        if(q > 0 && nx > 0){ // (bounds directly, general constraints by one product with DCt)
            int ngen = 0, tc, kc;
            struct blasfeo_dvec sg;
            blasfeo_create_dvec(S->ng, &sg, ow->sng);
            for(i = 0, k = 0; i < S->ncand; i++){
                const int id = S->cand[i];
                if(!(ow->eqf[id] & DAQP_MULTISTAGE_EQ_ELIM)) continue;
                if(multistage_decode(ow, id, &tc, &kc) == DAQP_MULTISTAGE_BOUND){
                    if(S->idxb[kc] >= nu) g->pa[S->idxb[kc]] += S->om[k];
                }
                else{
                    if(!ngen) memset(ow->sng, 0, S->ng*sizeof(c_float));
                    ow->sng[kc] = S->om[k];
                    ngen = 1;
                }
                k++;
            }
            if(ngen) multistage_dct_n(S, &sg, 0, nu, nx, g);
        }
        p = t > s0 ? p-nx : ow->wr.pa; // l_0 (pi_{-1}) in scratch
        if(nx > 0) memcpy(p, g->pa + nu, nx*sizeof(c_float));
        pn = p;
    }
    // Forwards: multipliers of H_{t+1} to the rows of stage t+1 and the eliminated equalities
    for(t = t0; t <= N; t++){
        Stage *S = ow->st+t;
        const int q = S->q, mh = t < N ? S[1].rh : 0;
        for(k = 0; t > t0 && k < S->rh; k++){ // om += Craw'(w in rows piv2), w: multipliers of H_t
            const c_float w = S[-1].om[S[-1].q - S->rh + k];
            const int j = S->piv2[k], r = S->r;
            S->om[S->piv[r+j]] += w;
            for(i = 0; i < r; i++) S->om[S->piv[i]] -= S->Z[i*(q-r)+j]*w;
        }
        for(i = 0, k = 0; full && i < S->ncand && k < q-mh; i++)
            if(ow->eqf[S->cand[i]] & DAQP_MULTISTAGE_EQ_ELIM) lam[S->cand[i]] = S->om[k++];
    }
    // pi_{t-1} = l_t + H_t'w_t (w_t final after the forward pass)
    for(t = s0+1, p = ow->pic; t <= N; p += ow->st[t].nx, t++){
        Stage *S = ow->st+t;
        struct blasfeo_dvec w, pv;
        if(t <= t0 || S->rh == 0) continue;
        blasfeo_create_dvec(S->rh, &w, S[-1].om + S[-1].q - S->rh);
        blasfeo_create_dvec(S->nx, &pv, p);
        blasfeo_dgemv_n(S->nx, S->rh, 1.0, &S->Ht, 0, 0, &w, 0, 1.0, &pv, 0, &pv, 0);
    }
    if(x0mult) for(i = 0; i < ow->st[0].nx; i++) lam[ow->x0src[i]] = -ow->wr.pa[i];
    ow->pi_valid = 1;
    if(full) ow->lam_done = 1;
}

// Branch and bound: the bound binary constraint i is at in the latest solution (shift: that of the same constraint one
// stage later), or 0 if none
static int multistage_bnb_side(const DAQPMultistageWork *ow, const int i, const int shift){
    const int s = shift ? multistage_ws_src(ow, i) : -1;
    return s >= 0 && ow->bside[s] ? ow->bside[s] : ow->bside[i];
}

// Branch and bound: solve the problem with the binary constraints fixed at their bounds in the latest solution (shift:
// one stage later), bounded by fb (as fval_bound). Returns whether it is feasible; it is then the best candidate
static int multistage_bnb_candidate(DAQPMultistageWork *ow, const int shift, const c_float fb){
    DAQPWorkspace *work = ow->work;
    DAQPBnB *b = work->bnb;
    int i, j, side, exitflag = -1;
    // (the other bound of each binary is saved and restored below; daqp_bnb leaves none of them in the working set)
    for(j = 0; j < b->nb; j++){
        i = b->bin_ids[j];
        if(!(side = multistage_bnb_side(ow, i, shift))) break;
        if(side == 1){ ow->bsave[i] = work->dupper[i]; work->dupper[i] = work->dlower[i]; }
        else{ ow->bsave[i] = work->dlower[i]; work->dlower[i] = work->dupper[i]; }
    }
    if(j == b->nb){ // (nothing to branch over)
        work->settings->fval_bound = fb;
        exitflag = daqp_bnb(work);
    }
    for(j--; j >= 0; j--){
        i = b->bin_ids[j];
        if(multistage_bnb_side(ow, i, shift) == 1) work->dupper[i] = ow->bsave[i]; else work->dlower[i] = ow->bsave[i];
    }
    if(exitflag <= 0) return 0;
    memcpy(ow->vcand, work->u, ow->nv*sizeof(c_float));
    ow->cand_fval = work->fval;
    ow->cand_slack = work->soft_slack;
    return 1;
}

// Branch and bound: record the bound each binary constraint is at in the solution (work->u); none if it is strictly
// between them (as a soft binary constraint can be)
static void multistage_bnb_record(DAQPMultistageWork *ow){
    DAQPWorkspace *work = ow->work;
    daqp_compute_Mu(work);
    for(int j = 0; j < work->bnb->nb; j++){
        const int i = work->bnb->bin_ids[j];
        const c_float lo = work->Mu[i] - work->dlower[i], up = work->dupper[i] - work->Mu[i];
        const c_float tol = work->settings->primal_tol*work->scaling[i];
        ow->bside[i] = lo <= tol ? 1 : up <= tol ? 2 : 0;
    }
}

// Whether the solution violates a hard constraint by more than primal_tol (in the units of the constraint, as
// daqp_violates_hard for dense problems; constraints that the LDP ignores are skipped)
static int multistage_violates_hard(DAQPWorkspace *work){
    int i;
    daqp_compute_Mu(work);
    for(i = 0; i < work->m; i++){
        if(DAQP_IS_SOFT(i) || (DAQP_IS_IMMUTABLE(i) && !DAQP_IS_ACTIVE(i))) continue;
        const c_float tol = work->settings->primal_tol*work->scaling[i];
        if(work->Mu[i] > work->dupper[i] + tol || work->Mu[i] < work->dlower[i] - tol) return 1;
    }
    return 0;
}

// Prox for a regularized Hessian, BnB for binaries (not combined)
static int multistage_solve_ldp(DAQPMultistageWork *ow){
    DAQPWorkspace *work = ow->work;
    int exitflag;
    if(work->bnb != NULL && ow->eps > 0) return DAQP_EXIT_NONCONVEX;
    if(work->bnb != NULL){
        // Candidates: the binary assignment of the latest solution, and that one stage later (if it differs). The better
        // feasible one is the incumbent of the tree
        const c_float fb = work->settings->fval_bound;
        int iter = 0, nodes = 0, has_cand = 0, shift, j, n;
        for(shift = 0; shift < 2; shift++){
            for(j = 0; shift && j < work->bnb->nb; j++){
                const int i = work->bnb->bin_ids[j];
                if(multistage_bnb_side(ow, i, 1) != multistage_bnb_side(ow, i, 0)) break;
            }
            if(shift && j == work->bnb->nb) break; // (the same assignment)
            work->iterations = work->bnb->nodecount = 0;
            has_cand |= multistage_bnb_candidate(ow, shift, has_cand ? 0.5*ow->cand_fval - daqp_bnb_offset(work) : fb);
            iter += work->iterations;
            nodes += work->bnb->nodecount;
        }
        // The root of the tree is warm started by the working set of the candidate (for the current problem), without
        // the binary constraints that it fixed
        for(j = 0, n = 0; has_cand && j < work->bnb->n_root_WS; j++)
            if(!(work->sense[DAQP_REMOVE_LOWER_FLAG(work->bnb->root_WS[j])] & DAQP_BINARY))
                work->bnb->root_WS[n++] = work->bnb->root_WS[j];
        if(has_cand) work->bnb->n_root_WS = n;
        work->settings->fval_bound = fb;
        exitflag = has_cand ? daqp_bnb_from(work, ow->vcand, ow->cand_fval, ow->cand_slack) : daqp_bnb(work);
        work->iterations += iter;
        work->bnb->nodecount += nodes;
        if(exitflag > 0) multistage_bnb_record(ow);
        else memset(ow->bside, 0, ow->m); // (no assignment older than the latest solve)
        return exitflag;
    }
    exitflag = ow->eps > 0 ? daqp_prox(work) : daqp_ldp(work);
    // Refine the active constraints if the working set is nearly singular (rounding in the
    // recursions, amplified by large multipliers), as for dense problems
    if(exitflag > 0 && work->n_active > 0 && work->sing_ind == DAQP_EMPTY_IND){
        c_float min_D = work->D[0];
        int i;
        for(i = 1; i < work->n_active; i++) if(work->D[i] < min_D) min_D = work->D[i];
        if(min_D < DAQP_REFINE_PIVOT){
            daqp_refine_active(work);
            ow->dz_valid = 0;
        }
    }
    // Constraints were only added above the rounding level (cycling persisted): check the solution, as daqp_solve
    if(exitflag > 0 && (work->state & DAQP_STATE_NOISE_FLOOR) && multistage_violates_hard(work))
        exitflag = DAQP_EXIT_OPTIMAL_INEXACT;
    return exitflag;
}

// Solution z_lqr + dz (in dz; to x if not NULL) and its objective
static c_float multistage_solution(DAQPMultistageWork *ow, c_float *x){
    DAQPWorkspace *work = ow->work;
    c_float fval = 0;
    int i, t;
    if(!ow->dz_valid || memcmp(ow->vlast, work->u, ow->nv*sizeof(c_float)) != 0)
        multistage_simulate(ow, work->u);
    ow->dz_valid = 0;
    for(t = 0; t <= ow->N; t++){
        Stage *S = ow->st+t;
        const int n = S->nu + S->nx;
        if(n == 0) continue;
        blasfeo_daxpy(n, 1.0, &S->z, 0, &S->dz, 0, &S->dz, 0);
        if(ow->eps > 0){ // (Jlqr is of the regularized problem)
            blasfeo_dsymv_l(n, 0.5, S->RSQ, 0, 0, &S->dz, 0, 1.0, &S->rq, 0, &ow->wg, 0);
            fval += blasfeo_ddot(n, &ow->wg, 0, &S->dz, 0);
        }
    }
    if(ow->eps == 0){
        fval = ow->Jlqr;
        for(i = 0; i < ow->nv; i++) fval += 0.5*work->u[i]*work->u[i];
    }
    for(t = ow->x0free; x != NULL && t <= ow->N; t++){ // [x_0; u_0; x_1; ...] of the user
        Stage *S = ow->st+t;
        memcpy(x, S->dz.pa + S->nu, S->nx*sizeof(c_float));
        memcpy(x + S->nx, S->dz.pa, S->nu*sizeof(c_float));
        x += S->nx + S->nu;
    }
    return fval;
}

// Multipliers of the working set (in lamv). The eliminated equalities, the bounds that fix x_0 and the
// costates follow on demand (daqp_multistage_duals). After BnB, the working set is the root's: none
static void multistage_multipliers(DAQPMultistageWork *ow, const int exitflag, const int infeasible){
    DAQPWorkspace *work = ow->work;
    const int mult = !infeasible && work->bnb == NULL;
    int i;
    for(i = 0; i < work->m; i++) ow->lamv[i] = 0;
    for(i = 0; mult && i < work->n_active; i++)
        ow->lamv[work->WS[i]] = work->lam_star[i]*ow->scale[work->WS[i]];
    ow->lam_valid = work->bnb == NULL;
    ow->lam_done = !(mult && exitflag > 0 && (ow->neqc > 0 || ow->nx0c > 0));
    ow->pi_valid = 0;
}

void daqp_multistage_solve(DAQPResult *res, DAQPWorkspace *work){
    DAQPMultistageWork *ow = (DAQPMultistageWork*)work->Mop->data;
    const int infeasible = ow->infeasible || multistage_x0_infeasible(ow);
    int i;
#ifdef PROFILING
    DAQPtimer timer;
    tic(&timer);
    if(work->settings->time_limit > 0) work->timer = &timer;
#endif
    if(infeasible){ // Inconsistent equalities, or a violated bound on the fixed x_0
        res->exitflag = DAQP_EXIT_INFEASIBLE;
        work->iterations = 0;
        work->soft_slack = 0;
        for(i = 0; i < ow->nv; i++) work->u[i] = 0;
    }
    else res->exitflag = multistage_solve_ldp(ow);
#ifdef PROFILING
    work->timer = NULL;
#endif
    res->fval = multistage_solution(ow, res->x);
    multistage_multipliers(ow, res->exitflag, infeasible);
    if(res->lam != NULL) memcpy(res->lam, ow->lamv, work->m*sizeof(c_float));
    res->soft_slack = work->soft_slack;
    res->iter = work->iterations;
    res->nodes = 0;
    // The first solution after a factorization that eliminated equalities is checked by the dual pass (kept for
    // daqp_multistage_duals). If the elimination was inaccurate, the equalities go to the working set (the update
    // refactorizes) and the problem is solved again
    int resolved = 0;
    if(ow->elim_check && res->exitflag > 0 && !ow->lam_done){
        ow->elim_check = 0;
        multistage_duals(ow, 1);
        if(ow->elim_bad && !ow->no_elim && daqp_multistage_update(work, 0) >= 0){
            const int iter = res->iter;
            daqp_multistage_solve(res, work); // (records its working set)
            res->iter += iter;
            resolved = 1;
        }
    }
    if(!resolved) multistage_ws_record(ow, res->exitflag > 0 && !infeasible && work->bnb == NULL);
#ifdef PROFILING
    toc(&timer);
    res->solve_time = get_time(&timer);
#else
    res->solve_time = 0;
#endif
}

int daqp_multistage_duals(DAQPWorkspace *work, c_float *lam, c_float *pi){
    DAQPMultistageWork *ow = (DAQPMultistageWork*)work->Mop->data;
    if(!ow->lam_valid) return 0;
    if(!ow->lam_done) multistage_duals(ow, 1); // One pass: the missing multipliers and the costates
    if(pi != NULL && !ow->pi_valid) multistage_duals(ow, 0); // (lamv complete)
    if(lam != NULL) memcpy(lam, ow->lamv, work->m*sizeof(c_float));
    if(pi != NULL) memcpy(pi, ow->pic, ow->npi*sizeof(c_float));
    return ow->elim_bad ? 2 : 1; // (2: the elimination was found inaccurate, see DAQP_MULTISTAGE_ELIM_CHECK)
}

int daqp_multistage_fixed_row(const DAQPWorkspace *work, const int id){
    const DAQPMultistageWork *ow = (const DAQPMultistageWork*)work->Mop->data;
    return ow->zrow[id] || (ow->eqf[id] & DAQP_MULTISTAGE_EQ_X0);
}

void free_daqp_multistage(DAQPWorkspace *work){
    DAQPMultistageWork *ow;
    void *mem;
    int own;
    if(work->Mop == NULL) return;
    ow = (DAQPMultistageWork*)work->Mop->data;
    mem = ow->mem; own = ow->own_mem; // (ow is in mem)
    work->Mop = NULL;
    work->scaling = work->dupper = work->dlower = work->Mu = work->M = NULL; // In mem
    work->sense = NULL;
    if(own){
        free(mem);
        free(work->rho_ls); // (the soft weights, if daqp_set_soft_weights allocated them: one block)
        work->rho_ls = work->rho_us = work->w_ls = work->w_us = NULL;
        free_daqp_workspace(work);
    }
    else free_daqp_bnb(work); // (allocated, also with given memory)
}

/* BLASFEO matrices of DAQPMultistage from its row-major matrices, in one block starting with BAbt.
 * A stage whose row-major matrices are those of the previous stage (the same pointers, and dimensions) shares its
 * BLASFEO matrix (e.g. all stages of a time-invariant problem): packed once, and read from cache by the recursions */

#define DAQP_MULTISTAGE_MAT(M,t) ((M) != NULL ? (M)[t] : NULL)
#define DAQP_MULTISTAGE_NKINDS 5 // BAbt, RSQrq, DCt, Et, Ft

// Dimensions of matrix kind k of stage t (0 if none) and the row-major matrices it is formed from
static void multistage_pack_src(const DAQPMultistage *prob, const int k, const int t, int dims[2], const c_float *src[3]){
    const int N = prob->N, nx = prob->nx[t], nu = DAQP_MULTISTAGE_NU(prob,t), nc = DAQP_MULTISTAGE_NC(prob,t);
    const int nx1 = t < N ? prob->nx[t+1] : 0, nu1 = t < N ? DAQP_MULTISTAGE_NU(prob,t+1) : 0;
    src[0] = src[1] = src[2] = NULL;
    dims[0] = nu+nx;
    switch(k){
        case 0: dims[1] = nx1; if(t < N){ src[0] = DAQP_MULTISTAGE_MAT(prob->B,t); src[1] = DAQP_MULTISTAGE_MAT(prob->A,t); } break;
        case 1: dims[1] = nu+nx;
            src[0] = DAQP_MULTISTAGE_MAT(prob->R,t); src[1] = DAQP_MULTISTAGE_MAT(prob->S,t); src[2] = DAQP_MULTISTAGE_MAT(prob->Q,t); break;
        case 2: dims[1] = DAQP_MULTISTAGE_NG(prob,t); src[0] = DAQP_MULTISTAGE_MAT(prob->D,t); src[1] = DAQP_MULTISTAGE_MAT(prob->C,t); break;
        case 3: dims[1] = nc; if(nc > 0){ src[0] = DAQP_MULTISTAGE_MAT(prob->Eu,t); src[1] = DAQP_MULTISTAGE_MAT(prob->Ex,t); } break;
        default: dims[0] = nu1+nx1; dims[1] = nc; if(nc > 0){ src[0] = DAQP_MULTISTAGE_MAT(prob->Fu,t); src[1] = DAQP_MULTISTAGE_MAT(prob->Fx,t); }
    }
}

// Whether matrix kind k of stage t is formed from the same data as that of stage t-1
static int multistage_pack_same(const DAQPMultistage *prob, const int k, const int t){
    int d0[2], d1[2];
    const c_float *s0[3], *s1[3];
    if(t == 0) return 0;
    multistage_pack_src(prob, k, t-1, d0, s0);
    multistage_pack_src(prob, k, t, d1, s1);
    // (the same split of the rows into inputs and states: of stage t+1 for Ft)
    const int s = k == 4 ? t+1 : t;
    if(s > prob->N || prob->nx[s] != prob->nx[s-1] || DAQP_MULTISTAGE_NU(prob,s) != DAQP_MULTISTAGE_NU(prob,s-1)) return 0;
    return d0[0] == d1[0] && d0[1] == d1[1] && d1[1] > 0 && s0[0] == s1[0] && s0[1] == s1[1] && s0[2] == s1[2];
}

// The BLASFEO matrices of kind k
static struct blasfeo_dmat* multistage_pack_mats(const DAQPMultistage *prob, const int k){
    struct blasfeo_dmat *mats[DAQP_MULTISTAGE_NKINDS] = {prob->BAbt, prob->RSQrq, prob->DCt, prob->Et, prob->Ft};
    return mats[k];
}

// Bytes of the matrices of stage t (in the block), or create them (at p); shared ones are copies of the previous
static size_t multistage_pack_stage(DAQPMultistage *prob, const int t, char *p){
    size_t size = 0;
    int k, dims[2];
    const c_float *src[3];
    for(k = 0; k < DAQP_MULTISTAGE_NKINDS; k++){
        struct blasfeo_dmat *A = p != NULL ? multistage_pack_mats(prob, k) + t : NULL;
        if(multistage_pack_same(prob, k, t)){
            if(A != NULL) *A = A[-1];
            continue;
        }
        multistage_pack_src(prob, k, t, dims, src);
        const int m = dims[0] > 0 ? dims[0] : 1, c = dims[1] > 0 ? dims[1] : 1;
        if(A != NULL) blasfeo_create_dmat(m, c, A, p + size);
        size += blasfeo_memsize_dmat(m, c);
    }
    return size;
}

// Whether the block of prob was formed with the sharing of matrices that its data has now
static int multistage_pack_fits(const DAQPMultistage *prob){
    int k, t;
    for(k = 0; k < DAQP_MULTISTAGE_NKINDS; k++){
        const struct blasfeo_dmat *A = multistage_pack_mats(prob, k);
        for(t = 1; t <= prob->N; t++)
            if((A[t].pA == A[t-1].pA) != multistage_pack_same(prob, k, t)) return 0;
    }
    return 1;
}

// Row major is the transpose in column major; NULL is zero (lower: only the lower triangle of a square block)
static void multistage_pack_block(const int m, const int n, const c_float *A, struct blasfeo_dmat *sB, const int bi,
        const int bj, const int lower){
    if(m <= 0 || n <= 0) return;
    if(A == NULL) blasfeo_dgese(m, n, 0.0, sB, bi, bj);
    else if(lower) blasfeo_pack_l_dmat(m, n, (c_float*)A, m, sB, bi, bj);
    else blasfeo_pack_dmat(m, n, (c_float*)A, m, sB, bi, bj);
}

// Pack matrix kind k of stage t (every block of it is written)
static void multistage_pack_one(DAQPMultistage *prob, const int k, const int t){
    const int nx = prob->nx[t], nu = DAQP_MULTISTAGE_NU(prob,t);
    struct blasfeo_dmat *A = multistage_pack_mats(prob, k) + t;
    int dims[2], i;
    const c_float *src[3];
    multistage_pack_src(prob, k, t, dims, src);
    if(dims[1] == 0) return;
    switch(k){
        case 0: // BAbt = [B'; A'] (NULL B is the identity)
            multistage_pack_block(nu, dims[1], src[0], A, 0, 0, 0);
            if(src[0] == NULL) for(i = 0; i < nu && i < dims[1]; i++) blasfeo_dgein1(1.0, A, i, i);
            multistage_pack_block(nx, dims[1], src[1], A, nu, 0, 0);
            break;
        case 1: // RSQ = [R S'; S Q] (lower triangle)
            multistage_pack_block(nu, nu, src[0], A, 0, 0, 1);
            multistage_pack_block(nx, nu, src[1], A, nu, 0, 0);
            multistage_pack_block(nx, nx, src[2], A, nu, nu, 1);
            break;
        case 4: // Ft = [Fu'; Fx'] (of stage t+1)
            multistage_pack_block(DAQP_MULTISTAGE_NU(prob,t+1), dims[1], src[0], A, 0, 0, 0);
            multistage_pack_block(prob->nx[t+1], dims[1], src[1], A, DAQP_MULTISTAGE_NU(prob,t+1), 0, 0);
            break;
        default: // DCt = [D'; C'], Et = [Eu'; Ex']
            multistage_pack_block(nu, dims[1], src[0], A, 0, 0, 0);
            multistage_pack_block(nx, dims[1], src[1], A, nu, 0, 0);
    }
}

// The block for the current sharing of matrices (freed and formed again if it has changed). -1 on failure
static int multistage_pack_alloc(DAQPMultistage *prob){
    const int N = prob->N;
    const size_t head = (DAQP_MULTISTAGE_NKINDS*(N+1)*sizeof(struct blasfeo_dmat) + DAQP_MULTISTAGE_ALIGN-1)
        & ~(size_t)(DAQP_MULTISTAGE_ALIGN-1);
    size_t size = 0;
    char *block, *p;
    int t;
    if(prob->BAbt != NULL){
        if(multistage_pack_fits(prob)) return 0;
        daqp_multistage_pack_free(prob);
    }
    for(t = 0; t <= N; t++) size += multistage_pack_stage(prob, t, NULL);
    if((block = calloc(head + size + DAQP_MULTISTAGE_ALIGN, 1)) == NULL) return -1;
    prob->BAbt = (struct blasfeo_dmat*)block;
    prob->RSQrq = prob->BAbt+N+1; prob->DCt = prob->RSQrq+N+1;
    prob->Et = prob->DCt+N+1; prob->Ft = prob->Et+N+1;
    p = (char*)(((size_t)block + head + DAQP_MULTISTAGE_ALIGN-1) & ~(size_t)(DAQP_MULTISTAGE_ALIGN-1));
    for(t = 0; t <= N; t++) p += multistage_pack_stage(prob, t, p);
    return 0;
}

int daqp_multistage_pack(DAQPMultistage *prob){
    int k, t;
    if(multistage_pack_alloc(prob) < 0) return -1;
    for(k = 0; k < DAQP_MULTISTAGE_NKINDS; k++){
        const struct blasfeo_dmat *A = multistage_pack_mats(prob, k);
        for(t = 0; t <= prob->N; t++)
            if(t == 0 || A[t].pA != A[t-1].pA) multistage_pack_one(prob, k, t); // (shared: once)
    }
    return 0;
}

void daqp_multistage_pack_free(DAQPMultistage *prob){
    free(prob->BAbt); // The block
    prob->BAbt = prob->RSQrq = prob->DCt = prob->Et = prob->Ft = NULL;
}

#endif // DAQP_MULTISTAGE
