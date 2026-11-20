// Compares the structure-exploiting multistage solver with dense DAQP on the
// condensed QP, for random problems (including warm starts).
#include "multistage.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double randu(void){ return 2.0*rand()/RAND_MAX - 1.0; }

typedef struct{
    DAQPMultistage prob;
    int nz, m, nv;
    int *ubase, *xbase; // Position of the first input and state bound of each stage
    double *rl, *ru, *wl, *wu; // Weights of the soft constraints (NULL: none)
    int *soft; // The soft constraints (if rl)
    int nbin, bin[8]; // The binary constraints
}Gen;

static double* rmat(int r, int c, double s){
    double* M = malloc((r*c > 0 ? r*c : 1)*sizeof(double));
    for(int i = 0; i < r*c; i++) M[i] = s*randu();
    return M;
}

// Random multistage problem with stage dimensions that vary if tv
// Features of the random problems
#define F_COUPLING 1 // Constraints between consecutive stages (input rates and mixed)
#define F_EQ 2 // Equality constraints of all kinds (eliminated)
#define F_EQ_BAD 4 // An inconsistent duplicate equality (infeasible)
#define F_EQ_DUP 8 // A consistent duplicate equality
#define F_EQ_COUPLING 16 // An equality between stages (kept in the working set)
#define F_X0FREE 32 // x_0 is a variable (with some bounds)
#define F_X0PART 64 // Some components of x_0 are fixed
#define F_X0BAD 128 // x_0 is fixed, and violates another bound on x_0
#define F_XFIRST 256 // The state bounds come before the input bounds
#define F_SOFT 512 // Soft inequality constraints (with tightened bounds)
#define F_PSD 1024 // Positive semidefinite stage Hessians (proximal-point iterations)
#define F_BIN 2048 // Binary input bounds (branch and bound)
#define F_GSPLIT 4096 // General constraints on inputs only, on states only, mixed, and zero (blocks of DCt)
#define F_SHARED 16384 // Consecutive stages with the same matrices (pointers: shared BLASFEO matrices), and a change of
                       // that sharing (the matrices are allocated again)
#define F_EQ_GAIN 8192 // A mixed equality with a large state coupling at every stage (loose bounds): an inaccurate
                       // elimination that daqp_multistage_duals detects, recovered by keeping the equalities in the working set

// Pointers to the bounds of constraint id (in the order of the constraints);
// *x0 is set if it is a bound on x_0 alone

// The solver reads the BLASFEO matrices (formed from the row-major data, again
// after it has changed)
static DAQPMultistage* packed;
// Given memory (ext_mem): sized from the dimensions, with the settings and core arrays of the caller
static int ext_mem = 0;
static void *ext_block = NULL;
static double *ext_d = NULL;
static int *ext_i = NULL;
static DAQPSettings ext_settings;
static int test_setup(DAQPMultistage* o, DAQPWorkspace* w){
    packed = o;
    daqp_multistage_pack(o);
    if(!ext_mem) return setup_daqp_multistage(o, w);
    // (core arrays for ext_mem 1: daqp_multistage_nvmax variables; 2: daqp_multistage_nv)
    const int n = ext_mem == 2 ? daqp_multistage_nv(o) : daqp_multistage_nvmax(o), m = daqp_multistage_m(o);
    int ns = 0;
    for(int i = 0; o->sense != NULL && i < m; i++) ns += (o->sense[i] & DAQP_SOFT) != 0;
    const int nn = n + ns + 1;
    const size_t size = daqp_multistage_memsize(o, ns);
    ext_block = malloc(size);
    double *d = ext_d = calloc(5*nn + nn*(nn+1)/2 + 2*n + 4*m, sizeof(double));
    w->WS = ext_i = calloc(nn, sizeof(int));
    w->lam = d; d += nn; w->lam_star = d; d += nn; w->D = d; d += nn; w->xldl = d; d += nn; w->zldl = d; d += nn;
    w->L = d; d += nn*(nn+1)/2; w->u = d; d += n; w->xold = d; d += n;
    w->rho_ls = d; d += m; w->rho_us = d; d += m; w->w_ls = d; d += m; w->w_us = d;
    for(int i = 0; i < m; i++) w->rho_ls[i] = w->rho_us[i] = DAQP_DEFAULT_RHO_SOFT;
    daqp_default_settings(&ext_settings);
    w->settings = &ext_settings;
    return setup_daqp_multistage_mem(o, w, ext_block, size);
}
static void ext_release(void){
    free(ext_block); free(ext_d); free(ext_i);
    ext_block = NULL; ext_d = NULL; ext_i = NULL;
}
static int test_update(DAQPWorkspace* w, const int mask){
    if(mask & DAQP_MULTISTAGE_UPDATE_MATRICES) daqp_multistage_pack(packed);
    return daqp_multistage_update(w, mask);
}
#define setup_daqp_multistage test_setup
#define daqp_multistage_update test_update

// F_SHARED: the matrices of stage t are those of stage t-1 where the dimensions agree (leaked, as all of the generator)
static void share_matrices(DAQPMultistage* o, const int feat){
    for(int t = 1; t <= o->N; t++){
        const int same = o->nx[t] == o->nx[t-1] && (t < o->N ? o->nu[t] == o->nu[t-1] : 1);
        if(!same) continue;
        o->Q[t] = o->Q[t-1];
        if(t == o->N) continue;
        o->A[t] = o->A[t-1]; o->B[t] = o->B[t-1]; o->S[t] = o->S[t-1]; o->R[t] = o->R[t-1];
        if(o->ng[t] == o->ng[t-1] && !(feat & F_GSPLIT)){ o->C[t] = o->C[t-1]; o->D[t] = o->D[t-1]; }
        if(o->nc[t] > 0 && o->nc[t] == o->nc[t-1] && t+1 < o->N){
            o->Ex[t] = o->Ex[t-1]; o->Eu[t] = o->Eu[t-1]; o->Fx[t] = o->Fx[t-1]; o->Fu[t] = o->Fu[t-1];
        }
    }
}

static void bound_ptr(DAQPMultistage* o, int id, double** lo, double** hi, int* x0){
    for(int t = 0; t <= o->N; t++){
        int nu = t < o->N ? o->nu[t] : 0, nc = (o->nc != NULL && t < o->N) ? o->nc[t] : 0;
        *x0 = 0;
        if(id < o->nb[t]){ *lo = o->lb[t]+id; *hi = o->ub[t]+id; *x0 = t == 0 && o->idxb[t][id] >= nu; return; }
        id -= o->nb[t];
        if(id < o->ng[t]){ *lo = o->lg[t]+id; *hi = o->ug[t]+id; return; }
        id -= o->ng[t];
        if(id < nc){ *lo = o->lc[t]+id; *hi = o->uc[t]+id; return; }
        id -= nc;
    }
}

static void generate(Gen* g, int N, int nxb, int nub, int tv, int feat){
    DAQPMultistage* o = &g->prob;
    memset(o, 0, sizeof(*o));
    o->N = N;
    o->nx = malloc((N+1)*sizeof(int)); o->nu = malloc((N+1)*sizeof(int));
    o->nb = malloc((N+1)*sizeof(int)); o->ng = malloc((N+1)*sizeof(int));
    g->ubase = malloc((N+1)*sizeof(int)); g->xbase = malloc((N+1)*sizeof(int));
#define ALLOCP(f) o->f = calloc(N+1,sizeof(double*))
    ALLOCP(A); ALLOCP(B); ALLOCP(b); ALLOCP(Q); ALLOCP(S); ALLOCP(R); ALLOCP(q); ALLOCP(r);
    ALLOCP(lb); ALLOCP(ub); ALLOCP(C); ALLOCP(D); ALLOCP(lg); ALLOCP(ug);
    o->idxb = calloc(N+1,sizeof(int*));
    for(int t = 0; t <= N; t++){
        o->nx[t] = nxb + (tv ? rand()%2 : 0);
        o->nu[t] = t < N ? nub + (tv ? rand()%2 : 0) : 0;
    }
    for(int t = 0; t <= N; t++){
        int nx = o->nx[t], nu = o->nu[t];
        // Stage Hessian W W' + 0.1 I in [x;u]. F_PSD: no curvature in the
        // states and in the last input (the rows of W, and of I, are zero), as
        // in time-optimal control: the Riccati recursion breaks down
        int nn = nx+nu;
        double* W = rmat(nn,nn,1.0), *H = calloc(nn*nn,sizeof(double));
        if(feat & F_PSD)
            for(int i = 0; i < nn; i++) if(i < nx || i == nn-1) for(int k = 0; k < nn; k++) W[i*nn+k] = 0;
        for(int i = 0; i < nn; i++) for(int j = 0; j < nn; j++){
            double s = (i == j && !((feat & F_PSD) && (i < nx || i == nn-1))) ? 0.1 : 0;
            for(int k = 0; k < nn; k++) s += W[i*nn+k]*W[j*nn+k];
            H[i*nn+j] = s/nn;
        }
        o->Q[t] = malloc(nx*nx*sizeof(double));
        for(int i = 0; i < nx; i++) for(int j = 0; j < nx; j++) o->Q[t][i*nx+j] = H[i*nn+j];
        if(t < N){
            o->R[t] = malloc(nu*nu*sizeof(double));
            o->S[t] = malloc(nu*nx*sizeof(double));
            for(int i = 0; i < nu; i++){
                for(int j = 0; j < nu; j++) o->R[t][i*nu+j] = H[(nx+i)*nn+nx+j];
                for(int j = 0; j < nx; j++) o->S[t][i*nx+j] = H[(nx+i)*nn+j];
            }
            int nx1 = o->nx[t+1];
            o->A[t] = rmat(nx1,nx,0.5);
            for(int i = 0; i < nx1 && i < nx; i++) o->A[t][i*nx+i] += 1.0;
            o->B[t] = rmat(nx1,nu,1.0);
            o->b[t] = rmat(nx1,1,0.05);
            o->r[t] = rmat(nu,1,0.5);
        }
        o->q[t] = rmat(nx,1,0.5);
        // Bounds: on all inputs, and on every other state (on x_0 depending on feat)
        int nbu = nu, nbx = 0, xi[16]; double xl[16], xu[16];
        if(t > 0) for(int i = 0; 2*i < nx; i++, nbx++){
            xi[nbx] = 2*i;
            xu[nbx] = 1.0 + fabs(randu());
            xl[nbx] = -1.0 - fabs(randu());
            if(i == 0 && t % 3 == 0) xl[nbx] = -DAQP_INF; // one-sided
        }
        else for(int i = 0; i < nx; i++){
            const int fix = (feat & F_X0FREE) ? 0 : (feat & F_X0PART) ? i % 2 == 0 : 1;
            xi[nbx] = i;
            if(fix){ xl[nbx] = xu[nbx] = 2.0*randu(); nbx++; }
            else if(i % 2 == 1){ xu[nbx] = 1.5 + fabs(randu()); xl[nbx] = -1.5 - fabs(randu()); nbx++; }
        }
        if(t == 0 && nx > 0 && !(feat & (F_X0FREE|F_X0PART))){
            // Another bound on the fixed x_0 (only checked)
            xi[nbx] = 0;
            xl[nbx] = (feat & F_X0BAD) ? -DAQP_INF : xl[0] - 1.0;
            xu[nbx] = (feat & F_X0BAD) ? xl[0] - 0.5 : xl[0] + 1.0;
            nbx++;
        }
        o->nb[t] = nbu + nbx;
        g->ubase[t] = (feat & F_XFIRST) ? nbx : 0;
        g->xbase[t] = (feat & F_XFIRST) ? 0 : nbu;
        o->idxb[t] = malloc((o->nb[t]+1)*sizeof(int));
        o->lb[t] = malloc((o->nb[t]+1)*sizeof(double)); o->ub[t] = malloc((o->nb[t]+1)*sizeof(double));
        for(int i = 0; i < nbu; i++){
            const int k = g->ubase[t]+i;
            o->idxb[t][k] = i;
            o->ub[t][k] = 0.3 + 0.5*fabs(randu());
            o->lb[t][k] = -0.3 - 0.5*fabs(randu());
        }
        for(int i = 0; i < nbx; i++){
            const int k = g->xbase[t]+i;
            o->idxb[t][k] = nu + xi[i]; o->lb[t][k] = xl[i]; o->ub[t][k] = xu[i];
        }
        // General constraints (mixed, x-only at N)
        o->ng[t] = (t == 0 && nu == 0) ? 0 : ((feat & F_GSPLIT) ? 2 + (t % 4) : 1 + (t % 2));
        o->C[t] = rmat(o->ng[t],nx,1.0);
        o->D[t] = t < N ? rmat(o->ng[t],nu,1.0) : NULL;
        if(feat & F_GSPLIT){ // Rows by t%4: [u x], [u x u], [x 0 mixed x] (u at t=0), [u x u x u] (more blocks than kept)
            for(int i = 0; i < o->ng[t]; i++){
                const int pat = t % 4;
                int cls = pat == 2 ? (i == 1 ? 3 : (i == 2 ? 2 : 1)) : (pat == 0 ? (i < 1 ? 0 : 1) : i % 2);
                if(cls == 1 && t == 0) cls = 2; // (not on the fixed x_0 alone)
                if((cls == 0 || cls == 3) && nx > 0) memset(o->C[t]+i*nx, 0, nx*sizeof(double));
                if((cls == 1 || cls == 3) && t < N && nu > 0) memset(o->D[t]+i*nu, 0, nu*sizeof(double));
            }
        }
        o->lg[t] = malloc(o->ng[t]*sizeof(double)); o->ug[t] = malloc(o->ng[t]*sizeof(double));
        for(int i = 0; i < o->ng[t]; i++){
            o->ug[t][i] = 1.0 + fabs(randu());
            o->lg[t][i] = -1.0 - fabs(randu());
        }
        if((feat & F_EQ) && t == N/2 && o->ng[t] > 0){
            o->lg[t][0] = o->ug[t][0] = 0.1; // mixed
            if((feat & (F_EQ_DUP|F_EQ_BAD)) && o->ng[t] > 1){
                memcpy(o->C[t]+nx, o->C[t], nx*sizeof(double));
                if(t < N) memcpy(o->D[t]+nu, o->D[t], nu*sizeof(double));
                o->lg[t][1] = o->ug[t][1] = (feat & F_EQ_BAD) ? 0.3 : 0.1;
            }
        }
        if((feat & F_EQ) && t == N && o->ng[t] > 0) o->lg[t][0] = o->ug[t][0] = 0; // terminal
        if((feat & F_EQ) && t == N/3 && t < N && nu > 0) // input
            o->lb[t][g->ubase[t]] = o->ub[t][g->ubase[t]] = 0.1*randu();
        if((feat & F_EQ) && t == (2*N)/3 && t > 1 && nbx > 0) // state
            o->lb[t][g->xbase[t]] = o->ub[t][g->xbase[t]] = 0.2;
        if((feat & F_EQ_GAIN) && t < N && nu > 0 && o->ng[t] > 0){
            for(int i = 0; i < nx; i++) o->C[t][i] *= 30;
            o->lg[t][0] = o->ug[t][0] = 0.1;
            for(int k = 0; k < o->nb[t]; k++) if(o->lb[t][k] != o->ub[t][k]){ o->lb[t][k] = -1e3; o->ub[t][k] = 1e3; }
        }
        free(W); free(H);
    }
    // Constraints between stage t and t+1
    o->nc = calloc(N+1,sizeof(int));
    o->Ex = calloc(N+1,sizeof(double*)); o->Eu = calloc(N+1,sizeof(double*));
    o->Fx = calloc(N+1,sizeof(double*)); o->Fu = calloc(N+1,sizeof(double*));
    o->lc = calloc(N+1,sizeof(double*)); o->uc = calloc(N+1,sizeof(double*));
    if(feat & (F_COUPLING|F_EQ_COUPLING)) for(int t = 0; t < N; t++){
        int nx = o->nx[t], nu = o->nu[t], nx1 = o->nx[t+1], nu1 = t+1 < N ? o->nu[t+1] : 0;
        int nc = 2;
        o->nc[t] = nc;
        o->Ex[t] = calloc(nc*nx,sizeof(double)); o->Eu[t] = calloc(nc*nu,sizeof(double));
        o->Fx[t] = calloc(nc*nx1,sizeof(double)); o->Fu[t] = nu1 > 0 ? calloc(nc*nu1,sizeof(double)) : NULL;
        o->lc[t] = malloc(nc*sizeof(double)); o->uc[t] = malloc(nc*sizeof(double));
        // Row 0: rate of the first input (or a random state difference at the end)
        if(nu1 > 0){ o->Eu[t][0] = -1; o->Fu[t][0] = 1; o->lc[t][0] = -0.3; o->uc[t][0] = 0.3; }
        else{ for(int j = 0; j < nx1; j++) o->Fx[t][j] = randu(); o->lc[t][0] = -2; o->uc[t][0] = 2; }
        // Row 1: random mixed
        for(int j = 0; j < nx; j++) o->Ex[t][nx+j] = randu();
        for(int j = 0; j < nu; j++) o->Eu[t][nu+j] = randu();
        for(int j = 0; j < nx1; j++) o->Fx[t][nx1+j] = randu();
        o->lc[t][1] = -1.5 - fabs(randu()); o->uc[t][1] = 1.5 + fabs(randu());
        if((feat & F_EQ_COUPLING) && t == 1) o->lc[t][1] = o->uc[t][1] = 0.05;
    }
    g->nz = daqp_multistage_nz(o);
    g->m = daqp_multistage_m(o);
    g->nv = o->nx[0];
    g->rl = g->ru = g->wl = g->wu = NULL;
    if(feat & F_SOFT){
        o->sense = calloc(g->m+1,sizeof(int));
        g->rl = malloc((g->m+1)*sizeof(double)); g->ru = malloc((g->m+1)*sizeof(double));
        g->wl = malloc((g->m+1)*sizeof(double)); g->wu = malloc((g->m+1)*sizeof(double));
        g->soft = calloc(g->m+1,sizeof(int));
        for(int k = 0; k < g->m; k++){
            double *lo = NULL, *hi = NULL; int x0;
            bound_ptr(o, k, &lo, &hi, &x0);
            g->rl[k] = 0.2 + fabs(randu()); g->ru[k] = 0.2 + fabs(randu());
            g->wl[k] = rand()%2 ? 0.5*fabs(randu()) : 0; g->wu[k] = rand()%2 ? 0.5*fabs(randu()) : 0;
            if(x0 || *lo == *hi || rand()%10 >= 3) continue;
            o->sense[k] = DAQP_SOFT;
            g->soft[k] = 1;
            // Tightened (so that some are violated)
            if(*lo > -DAQP_INF) *lo *= 0.2;
            if(*hi < DAQP_INF) *hi *= 0.2;
        }
    }
    // Binary constraints: the first input bound of some stages (not those that
    // other features make equalities)
    g->nbin = 0;
    if(feat & F_BIN){
        if(o->sense == NULL) o->sense = calloc(g->m+1,sizeof(int));
        for(int t = 1; t < N-1 && g->nbin < 5; t++){
            if(t == N/3 || o->nu[t] == 0) continue;
            for(int k = 0; k < g->m; k++){
                double *lo = NULL, *hi = NULL; int x0;
                bound_ptr(o, k, &lo, &hi, &x0);
                if(lo != o->lb[t] + g->ubase[t]) continue;
                o->sense[k] = DAQP_BINARY; // (not soft)
                if(g->rl != NULL) g->soft[k] = 0;
                g->bin[g->nbin++] = k;
            }
        }
    }
    for(int t = 0; t < N; t++) g->nv += o->nu[t];
}

static void constraints(const DAQPMultistage* o, const double* z, double* c, double* lb, double* ub);

static int solve_dense(const Gen* g, double* z, double* fval, double* time);

// Reference for binary constraints without branch and bound: the best of the
// QPs with each binary constraint fixed at one of its bounds (all 2^nbin).
// Returns whether one is feasible (with its objective in fbest)
static double objective_soft(const Gen* g, const double* z);
static int binaries_brute_force(Gen* g, double* fbest){
    DAQPMultistage* o = &g->prob;
    double *lo[8], *hi[8], l0[8], h0[8], *z = malloc(g->nz*sizeof(double)), f, t;
    int found = 0, x0;
    for(int i = 0; i < g->nbin; i++){
        bound_ptr(o, g->bin[i], lo+i, hi+i, &x0);
        l0[i] = *lo[i]; h0[i] = *hi[i];
        o->sense[g->bin[i]] &= ~DAQP_BINARY;
    }
    for(int c = 0; c < (1 << g->nbin); c++){
        for(int i = 0; i < g->nbin; i++) *lo[i] = *hi[i] = (c >> i) & 1 ? h0[i] : l0[i];
        if(solve_dense(g, z, &f, &t) > 0 && (f = objective_soft(g, z), !found || f < *fbest)){ *fbest = f; found = 1; }
    }
    for(int i = 0; i < g->nbin; i++){ *lo[i] = l0[i]; *hi[i] = h0[i]; o->sense[g->bin[i]] |= DAQP_BINARY; }
    free(z);
    return found;
}

// The binary constraints are at one of their bounds in z
static int binaries_ok(const Gen* g, const double* z){
    double *c = malloc(g->m*sizeof(double)), *lb = malloc(g->m*sizeof(double)), *ub = malloc(g->m*sizeof(double));
    int ok = 1;
    constraints(&g->prob, z, c, lb, ub);
    for(int i = 0; i < g->nbin; i++){
        const int k = g->bin[i];
        if(fabs(c[k]-lb[k]) > 1e-6 && fabs(c[k]-ub[k]) > 1e-6) ok = 0;
    }
    free(c); free(lb); free(ub);
    return ok;
}

/* Reference: condensed dense QP in x_0 and the stacked inputs, formed
 * independently of multistage.c */
static void simulate(const DAQPMultistage* o, const double* v, int use_b, double* z){
    int off = 0, vo = o->nx[0];
    memcpy(z, v, o->nx[0]*sizeof(double));
    for(int t = 0; t < o->N; t++){
        int nx = o->nx[t], nu = o->nu[t], nx1 = o->nx[t+1];
        double *x = z+off, *u = x+nx, *xn = u+nu;
        memcpy(u, v+vo, nu*sizeof(double));
        for(int i = 0; i < nx1; i++){
            double s = use_b ? o->b[t][i] : 0;
            for(int j = 0; j < nx; j++) s += o->A[t][i*nx+j]*x[j];
            for(int j = 0; j < nu; j++) s += o->B[t][i*nu+j]*u[j];
            xn[i] = s;
        }
        off += nx+nu; vo += nu;
    }
}
static void constraints(const DAQPMultistage* o, const double* z, double* c, double* lb, double* ub){
    int off = 0, k = 0;
    for(int t = 0; t <= o->N; t++){
        int nx = o->nx[t], nu = t < o->N ? o->nu[t] : 0;
        const double *x = z+off, *u = x+nx;
        for(int i = 0; i < o->nb[t]; i++, k++){
            const int j = o->idxb[t][i];
            c[k] = j < nu ? u[j] : x[j-nu];
            if(lb){ lb[k] = o->lb[t][i]; ub[k] = o->ub[t][i]; }
        }
        for(int i = 0; i < o->ng[t]; i++, k++){
            double s = 0;
            for(int j = 0; j < nx; j++) s += o->C[t][i*nx+j]*x[j];
            if(t < o->N) for(int j = 0; j < nu; j++) s += o->D[t][i*nu+j]*u[j];
            c[k] = s;
            if(lb){ lb[k] = o->lg[t][i]; ub[k] = o->ug[t][i]; }
        }
        if(o->nc != NULL && t < o->N){
            int nx1 = o->nx[t+1], nu1 = t+1 < o->N ? o->nu[t+1] : 0;
            const double *xn = z+off+nx+nu, *un = xn+nx1;
            for(int i = 0; i < o->nc[t]; i++, k++){
                double s = 0;
                for(int j = 0; j < nx; j++) s += o->Ex[t][i*nx+j]*x[j];
                for(int j = 0; j < nu; j++) s += o->Eu[t][i*nu+j]*u[j];
                for(int j = 0; j < nx1; j++) s += o->Fx[t][i*nx1+j]*xn[j];
                if(o->Fu[t] != NULL) for(int j = 0; j < nu1; j++) s += o->Fu[t][i*nu1+j]*un[j];
                c[k] = s;
                if(lb){ lb[k] = o->lc[t][i]; ub[k] = o->uc[t][i]; }
            }
        }
        off += nx+nu;
    }
}
static double objective(const DAQPMultistage* o, const double* z){
    double f = 0; int off = 0;
    for(int t = 0; t <= o->N; t++){
        int nx = o->nx[t], nu = t < o->N ? o->nu[t] : 0;
        const double *x = z+off, *u = x+nx;
        for(int i = 0; i < nx; i++){
            for(int j = 0; j < nx; j++) f += 0.5*x[i]*o->Q[t][i*nx+j]*x[j];
            f += o->q[t][i]*x[i];
        }
        for(int i = 0; i < nu; i++){
            for(int j = 0; j < nu; j++) f += 0.5*u[i]*o->R[t][i*nu+j]*u[j];
            for(int j = 0; j < nx; j++) f += u[i]*o->S[t][i*nx+j]*x[j];
            f += o->r[t][i]*u[i];
        }
        off += nx+nu;
    }
    return f;
}

// The latest condensed QP (kept for the KKT check of the multipliers)
static double *KH = NULL, *Kf = NULL, *KA = NULL;
static int Kn = 0, Km = 0;

// Solve the condensed QP with dense DAQP; z is the resulting trajectory
static int solve_dense(const Gen* g, double* z, double* fval, double* time){
    const DAQPMultistage* o = &g->prob;
    int nv = g->nv, nz = g->nz, m = g->m;
    double *G = calloc((size_t)nz*nv,sizeof(double)), *z0 = malloc(nz*sizeof(double));
    double *e = calloc(nv,sizeof(double)), *col = malloc(nz*sizeof(double));
    simulate(o, e, 1, z0);
    for(int j = 0; j < nv; j++){
        e[j] = 1; simulate(o, e, 0, col); e[j] = 0;
        for(int i = 0; i < nz; i++) G[(size_t)i*nv+j] = col[i];
    }
    // Objective: J(z0 + G v) = J(z0) + grad(z0)'Gv + 0.5 v'G'HzGv. Use finite
    // structure: obj is quadratic, so get Hz*w by polarization on columns.
    double *H = calloc((size_t)nv*nv,sizeof(double)), *f = calloc(nv,sizeof(double));
    double J0 = objective(o, z0);
    double *zc = malloc(nz*sizeof(double)), *zd = malloc(nz*sizeof(double));
    double *Jcol = malloc(nv*sizeof(double)); // J(z0 + G e_j)
    for(int j = 0; j < nv; j++){
        for(int i = 0; i < nz; i++) zc[i] = z0[i] + G[(size_t)i*nv+j];
        Jcol[j] = objective(o, zc);
        for(int i = 0; i < nz; i++) zd[i] = z0[i] - G[(size_t)i*nv+j];
        double Jm = objective(o, zd);
        H[(size_t)j*nv+j] = Jcol[j] + Jm - 2*J0;
        f[j] = 0.5*(Jcol[j] - Jm);
    }
    for(int j = 0; j < nv; j++) for(int k = 0; k < j; k++){
        for(int i = 0; i < nz; i++) zc[i] = z0[i] + G[(size_t)i*nv+j] + G[(size_t)i*nv+k];
        double Jjk = objective(o, zc);
        // J(z0+Gj+Gk) = J0 + fj + fk + 0.5(Hjj + Hkk) + Hjk
        H[(size_t)j*nv+k] = H[(size_t)k*nv+j] =
            Jjk - J0 - f[j] - f[k] - 0.5*(H[(size_t)j*nv+j] + H[(size_t)k*nv+k]);
    }
    // Constraints
    double *A = malloc((size_t)m*nv*sizeof(double)), *c0 = malloc(m*sizeof(double));
    double *bu = malloc(m*sizeof(double)), *bl = malloc(m*sizeof(double)), *cc = malloc(m*sizeof(double));
    constraints(o, z0, c0, bl, bu);
    for(int j = 0; j < nv; j++){
        for(int i = 0; i < nz; i++) col[i] = G[(size_t)i*nv+j];
        constraints(o, col, cc, NULL, NULL);
        for(int i = 0; i < m; i++) A[(size_t)i*nv+j] = cc[i];
    }
    for(int i = 0; i < m; i++){
        if(bu[i] < DAQP_INF) bu[i] -= c0[i];
        if(bl[i] > -DAQP_INF) bl[i] -= c0[i];
    }
    int* sense = calloc(m,sizeof(int));
    if(o->sense != NULL) memcpy(sense, o->sense, m*sizeof(int));
    DAQPProblem qp = {nv, m, 0, H, f, A, bu, bl, sense, NULL, 0, 0};
    double* v = malloc(nv*sizeof(double));
    DAQPResult res; res.x = v; res.lam = NULL;
    clock_t t0 = clock();
    if(g->rl == NULL) daqp_quadprog(&res, &qp, NULL);
    else{
        DAQPWorkspace w; memset(&w, 0, sizeof(w));
        double st;
        res.exitflag = setup_daqp(&qp, &w, &st);
        if(res.exitflag >= 0){
            daqp_set_soft_weights(&w, g->rl, g->ru, g->wl, g->wu);
            daqp_solve(&res, &w);
            free_daqp_workspace(&w);
            free_daqp_ldp(&w);
        }
    }
    *time = (double)(clock()-t0)/CLOCKS_PER_SEC;
    for(int i = 0; i < nz; i++){
        double s = z0[i];
        for(int j = 0; j < nv; j++) s += G[(size_t)i*nv+j]*v[j];
        z[i] = s;
    }
    *fval = objective(o, z);
    free(KH); free(Kf); free(KA);
    KH = H; Kf = f; KA = A; Kn = nv; Km = m;
    free(G); free(z0); free(e); free(col); free(zc); free(zd);
    free(Jcol); free(c0); free(bu); free(bl); free(cc); free(sense); free(v);
    return res.exitflag;
}

// Whether the hard constraints are feasible: the condensed QP without the soft
// constraints (which cannot make an infeasible problem feasible). Its
// solution replaces the latest condensed QP, so this is called before solve_dense
static int hard_feasible(Gen* g){
    DAQPMultistage* o = &g->prob;
    double *lo, *hi, *rl = g->rl, *z = malloc(g->nz*sizeof(double)), f, t;
    double *saved = malloc(2*(g->m+1)*sizeof(double));
    int x0, e;
    for(int k = 0; k < g->m; k++) if(g->soft[k]){
        bound_ptr(o, k, &lo, &hi, &x0);
        saved[2*k] = *lo; saved[2*k+1] = *hi;
        *lo = -DAQP_INF; *hi = DAQP_INF; o->sense[k] = 0;
    }
    g->rl = NULL;
    e = solve_dense(g, z, &f, &t);
    g->rl = rl;
    for(int k = 0; k < g->m; k++) if(g->soft[k]){
        bound_ptr(o, k, &lo, &hi, &x0);
        *lo = saved[2*k]; *hi = saved[2*k+1]; o->sense[k] = DAQP_SOFT;
    }
    free(z); free(saved);
    return e > 0;
}

// Relative KKT residual H v + f + A'lam of the latest condensed QP at the
// inputs of z (lam with the multipliers of all constraints)
static double kkt_residual(const Gen* g, const double* z, const double* lam){
    const DAQPMultistage* o = &g->prob;
    double* v = malloc(Kn*sizeof(double)), res = 0, scale = 1;
    int off = 0, vo = o->nx[0];
    memcpy(v, z, o->nx[0]*sizeof(double));
    for(int t = 0; t < o->N; t++){
        memcpy(v+vo, z+off+o->nx[t], o->nu[t]*sizeof(double));
        off += o->nx[t]+o->nu[t]; vo += o->nu[t];
    }
    for(int i = 0; i < Kn; i++){
        double r = Kf[i];
        for(int j = 0; j < Kn; j++) r += KH[(size_t)i*Kn+j]*v[j];
        for(int k = 0; k < Km; k++) r += KA[(size_t)k*Kn+i]*lam[k];
        if(fabs(Kf[i]) > scale) scale = fabs(Kf[i]);
        if(fabs(r) > res) res = fabs(r);
    }
    for(int k = 0; k < Km; k++) if(fabs(lam[k]) > scale) scale = fabs(lam[k]);
    free(v);
    return res/scale;
}

// Relative residual of the stationarity of the Lagrangian in all of z, with
// the multipliers lam of the constraints and pi of the dynamics
// (x_{t+1} = A_t x_t + B_t u_t + b_t, with grad_{x_{t+1}} = pi_t)
static double kkt_full(const Gen* g, const double* z, const double* lam, const double* pi){
    const DAQPMultistage* o = &g->prob;
    double *r = calloc(g->nz,sizeof(double)), *e = calloc(g->nz,sizeof(double));
    double *c = malloc((g->m+1)*sizeof(double)), res = 0, scale = 1;
    int off = 0, po = 0;
    for(int t = 0; t <= o->N; t++){
        int nx = o->nx[t], nu = t < o->N ? o->nu[t] : 0;
        const double *x = z+off, *u = x+nx;
        double *rx = r+off, *ru = rx+nx;
        for(int i = 0; i < nx; i++){
            rx[i] += o->q[t][i];
            for(int j = 0; j < nx; j++) rx[i] += o->Q[t][i*nx+j]*x[j];
            for(int j = 0; j < nu; j++) rx[i] += o->S[t][j*nx+i]*u[j];
        }
        for(int i = 0; i < nu; i++){
            ru[i] += o->r[t][i];
            for(int j = 0; j < nu; j++) ru[i] += o->R[t][i*nu+j]*u[j];
            for(int j = 0; j < nx; j++) ru[i] += o->S[t][i*nx+j]*x[j];
        }
        if(t < o->N){
            int nx1 = o->nx[t+1];
            const double *p = pi+po;
            for(int i = 0; i < nx1; i++){
                for(int j = 0; j < nx; j++) rx[j] += o->A[t][i*nx+j]*p[i];
                for(int j = 0; j < nu; j++) ru[j] += o->B[t][i*nu+j]*p[i];
                r[off+nx+nu+i] -= p[i];
                if(fabs(p[i]) > scale) scale = fabs(p[i]);
            }
            po += nx1;
        }
        off += nx+nu;
    }
    // Constraint terms: column j of the (linear) constraints is c(e_j)
    for(int j = 0; j < g->nz; j++){
        e[j] = 1; constraints(o, e, c, NULL, NULL); e[j] = 0;
        for(int k = 0; k < g->m; k++) r[j] += c[k]*lam[k];
    }
    for(int k = 0; k < g->m; k++) if(fabs(lam[k]) > scale) scale = fabs(lam[k]);
    for(int j = 0; j < g->nz; j++) if(fabs(r[j]) > res) res = fabs(r[j]);
    free(r); free(e); free(c);
    return res/scale;
}

// The objective with the penalties w*v + v^2/(2*rho) of the violations v of the soft constraints
static double objective_soft(const Gen* g, const double* z){
    if(g->rl == NULL) return objective(&g->prob, z);
    double *c = malloc(g->m*sizeof(double)), *lb = malloc(g->m*sizeof(double)), *ub = malloc(g->m*sizeof(double));
    double f = objective(&g->prob, z);
    constraints(&g->prob, z, c, lb, ub);
    for(int i = 0; i < g->m; i++){
        if(!(g->prob.sense[i] & DAQP_SOFT)) continue;
        const double v = c[i] > ub[i] ? c[i]-ub[i] : c[i] < lb[i] ? lb[i]-c[i] : 0;
        const double w = c[i] > ub[i] ? g->wu[i] : g->wl[i], rho = c[i] > ub[i] ? g->ru[i] : g->rl[i];
        f += w*v + v*v/(2*rho);
    }
    free(c); free(lb); free(ub);
    return f;
}

// The solutions agree: in the objective (with soft constraints, including the penalties of the slacks), or zo is
// better (the dense reference of an ill-conditioned condensed QP can be suboptimal; zo is checked to be feasible)
static int agree(const Gen* g, double fd, double fo, const double* zd, const double* zo){
    if(g->rl != NULL){ fd = objective_soft(g, zd); fo = objective_soft(g, zo); }
    return fabs(fd-fo) <= 1e-6*(1+fabs(fd)) || (g->rl != NULL ? fo : objective(&g->prob, zo)) < fd;
}

static double max_violation(const Gen* g, const double* z){
    double *c = malloc(g->m*sizeof(double)), *lb = malloc(g->m*sizeof(double)), *ub = malloc(g->m*sizeof(double));
    constraints(&g->prob, z, c, lb, ub);
    double v = 0;
    for(int i = 0; i < g->m; i++){
        if(g->prob.sense != NULL && (g->prob.sense[i] & DAQP_SOFT)) continue;
        if(isnan(c[i])) v = INFINITY;
        if(c[i]-ub[i] > v) v = c[i]-ub[i];
        if(lb[i]-c[i] > v) v = lb[i]-c[i];
    }
    free(c); free(lb); free(ub);
    return v;
}

// The KKT conditions of the latest condensed QP at z with the multipliers lam:
// feasibility, stationarity, and signs and complementarity of lam (which certify
// the solution of a convex QP without a reference; dense DAQP is unreliable for
// F_PSD, where the condensed QP is often an LP)
static int certified(const Gen* g, const double* z, const double* lam){
    double *c = malloc(g->m*sizeof(double)), *lb = malloc(g->m*sizeof(double)), *ub = malloc(g->m*sizeof(double));
    double comp = 0, scale = 1;
    int sign_ok = 1;
    constraints(&g->prob, z, c, lb, ub);
    for(int k = 0; k < g->m; k++){
        if(fabs(lam[k]) > scale) scale = fabs(lam[k]);
        if(lam[k] > 0){
            if(ub[k] >= DAQP_INF) sign_ok = 0;
            else if(lam[k]*fabs(ub[k]-c[k]) > comp) comp = lam[k]*fabs(ub[k]-c[k]);
        }
        if(lam[k] < 0){
            if(lb[k] <= -DAQP_INF) sign_ok = 0;
            else if(-lam[k]*fabs(c[k]-lb[k]) > comp) comp = -lam[k]*fabs(c[k]-lb[k]);
        }
    }
    free(c); free(lb); free(ub);
    // (the proximal-point iterations stop at a stationarity of 1e-6 in each
    // component of z, which the condensed residual sums over the horizon)
    return sign_ok && comp < 1e-6*scale && max_violation(g, z) < 1e-6 && kkt_residual(g, z, lam) < 1e-4;
}

// F_PSD: a solution of the multistage solver is judged by the KKT conditions alone
// (dense DAQP can return an infeasible point as optimal); a failure only if
// dense DAQP has found a feasible solution. Returns 0 on a mismatch
static int psd_ok(const Gen* g, int ed, int eo, const double* zd, const double* zo, const double* lam){
    if(eo > 0) return certified(g, zo, lam);
    return !(ed > 0 && max_violation(g, zd) < 1e-6);
}

// Solve, then complete the multipliers (daqp_multistage_duals). The multipliers of the solve
// (working set) have to be unchanged by the completion
static int dual_mismatch = 0;
static int elim_detected = 0;
// Solve and complete the multipliers (counts the solves whose elimination was found inaccurate: the solver then
// keeps the equalities in the working set and solves again itself)
static void ms_solve(DAQPResult* res, DAQPWorkspace* w){
    daqp_multistage_solve(res, w);
    if(res->lam == NULL) return;
    const int m = w->m;
    double* ws = malloc((m+1)*sizeof(double));
    memcpy(ws, res->lam, m*sizeof(double));
    const int d = daqp_multistage_duals(w, res->lam, NULL);
    if(d) for(int i = 0; i < m; i++) if(ws[i] != 0 && ws[i] != res->lam[i]) dual_mismatch++;
    free(ws);
    elim_detected += d == 2;
}

// fval_bound and rel_subopt refer to the objective of the QP: the input v, at its bound 0 with a large linear term,
// makes the objective of the LDP much larger than that of the QP
static int test_bnb_objective(void){
    enum { N = 20 };
    double A[1] = {1}, B[2] = {1, 0}, b[1] = {-0.4}, Q[1] = {1}, R[4] = {1e-3, 0, 0, 1}, r[2] = {0, -100};
    int nxs[N+1], nus[N+1], nbs[N+1], *idxb[N+1];
    double *As[N+1], *Bs[N+1], *bs[N+1], *Qs[N+1], *Rs[N+1], *rs[N+1], *lb[N+1], *ub[N+1];
    DAQPMultistage o; memset(&o, 0, sizeof(o));
    for(int t = 0; t <= N; t++){ // (bounds: the binary input and v, at stage 0 also x_0)
        nxs[t] = 1; nus[t] = t < N ? 2 : 0; nbs[t] = (t < N ? 2 : 0) + (t == 0);
        As[t] = A; Bs[t] = B; bs[t] = b; Qs[t] = Q; Rs[t] = R; rs[t] = r;
        idxb[t] = malloc(3*sizeof(int)); lb[t] = malloc(3*sizeof(double)); ub[t] = malloc(3*sizeof(double));
        for(int i = 0; i < nbs[t]; i++){ idxb[t][i] = i; lb[t][i] = i == 1 ? -1 : 0; ub[t][i] = i == 0 ? 1 : 0; }
    }
    lb[0][2] = ub[0][2] = 1.3;
    o.N = N; o.nx = nxs; o.nu = nus; o.nb = nbs; o.A = As; o.B = Bs; o.b = bs; o.Q = Qs; o.R = Rs; o.r = rs;
    o.idxb = idxb; o.lb = lb; o.ub = ub;
    const int m = daqp_multistage_m(&o);
    int *sense = calloc(m, sizeof(int));
    for(int t = 0, c = 0; t < N; c += nbs[t], t++) sense[c] = DAQP_BINARY;
    o.sense = sense;
    double *z = malloc(daqp_multistage_nz(&o)*sizeof(double)), f0 = 0, f = 0;
    int ok = 1;
    ext_mem = 0;
    for(int k = 0; k < 3 && ok; k++){
        DAQPWorkspace w; w.settings = NULL;
        DAQPResult res; memset(&res, 0, sizeof(res)); res.x = z;
        if(setup_daqp_multistage(&o, &w) < 0){ ok = 0; break; }
        w.settings->rel_subopt = k == 1 ? 0.05 : 0;
        if(k == 2) w.settings->fval_bound = f0 + 0.1;
        daqp_multistage_solve(&res, &w);
        free_daqp_multistage(&w);
        packed = &o;
        ok = res.exitflag > 0;
        if(k == 0) f0 = res.fval;
        else if(k == 1) f = res.fval;
        else ok = fabs(res.fval - f0) <= 1e-9*(1 + fabs(f0));
    }
    ok &= f >= f0 - 1e-9 && f - f0 <= 0.05*fabs(f) + 1e-9;
    if(!ok) printf("branch and bound with fval_bound or rel_subopt: failed\n");
    daqp_multistage_pack_free(&o);
    for(int t = 0; t <= N; t++){ free(idxb[t]); free(lb[t]); free(ub[t]); }
    free(sense); free(z);
    return ok;
}

// Chattering (working sets that alternate): a double integrator whose input saturates at its upper bound for one
// initial state and at its lower bound for the next. The solutions alternate, and once the alternation is detected,
// the warm start (the working set of two solves back) needs at most one iteration per solve
// Stage of constraint id
static int ms_stage(const DAQPMultistage* o, int id){
    int t = 0;
    for(; t < o->N; t++){
        const int nt = o->nb[t] + (o->ng != NULL ? o->ng[t] : 0) + (o->nc != NULL ? o->nc[t] : 0);
        if(id < nt) break;
        id -= nt;
    }
    return t;
}

// binary_horizon: the same solution as with DAQP_BINARY removed from the constraints of the later stages
static int test_binary_horizon(Gen* g){
    DAQPMultistage* o = &g->prob;
    const int H = ms_stage(o, g->bin[g->nbin/2]); // (the binaries from the middle one on relaxed)
    int *s0 = o->sense, *s1 = malloc(g->m*sizeof(int)), ok = 1, nh = 0;
    for(int i = 0; i < g->m; i++){
        s1[i] = ms_stage(o, i) >= H ? s0[i] & ~DAQP_BINARY : s0[i];
        nh += (s1[i] & DAQP_BINARY) != 0;
    }
    double *z = malloc(g->nz*sizeof(double));
    DAQPWorkspace wh, wr; wh.settings = wr.settings = NULL;
    DAQPResult rh, rr; memset(&rh, 0, sizeof(rh)); memset(&rr, 0, sizeof(rr)); rh.x = rr.x = z;
    ext_mem = 0;
    o->binary_horizon = H;
    const int eh = setup_daqp_multistage(o, &wh);
    o->binary_horizon = 0;
    o->sense = s1;
    const int er = setup_daqp_multistage(o, &wr);
    o->sense = s0;
    if(eh < 0 || er < 0) ok = eh == er;
    else{
        ok = (wh.bnb != NULL ? wh.bnb->nb : 0) == nh && (wr.bnb != NULL ? wr.bnb->nb : 0) == nh;
        daqp_multistage_solve(&rh, &wh);
        daqp_multistage_solve(&rr, &wr);
        ok &= (rh.exitflag > 0) == (rr.exitflag > 0) && (rh.exitflag <= 0 || fabs(rh.fval - rr.fval) <= 1e-6*(1 + fabs(rr.fval)));
    }
    if(eh >= 0) free_daqp_multistage(&wh);
    if(er >= 0) free_daqp_multistage(&wr);
    free(s1); free(z);
    return ok;
}

static int test_chattering(const int refactor, const int extra){
    enum { N = 15, nx = 2, nu = 1, steps = 16 };
    double Ad[4] = {1, 0.1, 0, 1}, Bd[2] = {0.005, 0.1}, Qd[4] = {1, 0, 0, 1}, Rd[1] = {0.01};
    int nxs[N+1], nus[N+1], nbs[N+1], *idxb[N+1];
    double *A[N+1], *B[N+1], *Q[N+1], *R[N+1], *lb[N+1], *ub[N+1];
    DAQPMultistage o; memset(&o, 0, sizeof(o));
    for(int t = 0; t <= N; t++){
        nxs[t] = nx; nus[t] = t < N ? nu : 0; nbs[t] = nus[t] + (t == 0 ? nx : 0);
        A[t] = Ad; B[t] = Bd; Q[t] = Qd; R[t] = Rd;
        idxb[t] = malloc((nbs[t]+1)*sizeof(int)); lb[t] = malloc((nbs[t]+1)*sizeof(double)); ub[t] = malloc((nbs[t]+1)*sizeof(double));
        for(int i = 0; i < nbs[t]; i++){ idxb[t][i] = i; lb[t][i] = -1; ub[t][i] = 1; }
    }
    o.N = N; o.nx = nxs; o.nu = nus; o.nb = nbs; o.A = A; o.B = B; o.Q = Q; o.R = R; o.idxb = idxb; o.lb = lb; o.ub = ub;
    const int nz = daqp_multistage_nz(&o);
    double *z = malloc(nz*sizeof(double)), *zs = malloc(2*nz*sizeof(double));
    long it = 0;
    int ok = 1;
    ext_mem = 0;
    DAQPWorkspace work; work.settings = NULL;
    lb[0][nu] = ub[0][nu] = 5; lb[0][nu+1] = ub[0][nu+1] = 0;
    const int es = setup_daqp_multistage(&o, &work);
    if(es < 0) ok = 0;
    for(int st = 0; st < steps && ok; st++){
        lb[0][nu] = ub[0][nu] = st%2 ? -5 : 5;
        daqp_multistage_update(&work, (refactor ? DAQP_MULTISTAGE_UPDATE_MATRICES : 0) | extra);
        DAQPResult res; memset(&res, 0, sizeof(res)); res.x = z; res.lam = NULL;
        daqp_multistage_solve(&res, &work);
        if(res.exitflag <= 0) ok = 0;
        if(st < 2) memcpy(zs + st*nz, z, nz*sizeof(double));
        else for(int i = 0; i < nz; i++) if(fabs(z[i] - zs[(st%2)*nz+i]) > 1e-6) ok = 0;
        if(st >= steps/2) it += res.iter;
    }
    if(es >= 0) free_daqp_multistage(&work);
    daqp_multistage_pack_free(&o);
    printf("chattering (%s%s): iterations of the last %d solves %ld\n",
            refactor ? "refactorizations" : "x0 updates", extra ? ", flagged as shifts" : "", steps/2, it);
    for(int t = 0; t <= N; t++){ free(idxb[t]); free(lb[t]); free(ub[t]); }
    free(z); free(zs);
    return ok && it <= steps/2;
}

// Soft constraints in the warm start of a refactorization: a double integrator that starts outside its soft position
// bounds (linear weight w > 0), with slacks both fixed and free at the solution. Solved again after a refactorization,
// the warm start must keep their states (one iteration, the same solution); with w then set to 0 (all slacks free),
// the solution must be that of a new setup
static int test_soft_warm(void){
    enum { N = 30, nx = 2, nu = 1 };
    double Ad[4] = {1, 0.1, 0, 1}, Bd[2] = {0.005, 0.1}, Qd[4] = {1, 0, 0, 1}, Rd[1] = {0.01};
    int nxs[N+1], nus[N+1], nbs[N+1], *idxb[N+1];
    double *A[N+1], *B[N+1], *Q[N+1], *R[N+1], *lb[N+1], *ub[N+1];
    DAQPMultistage o; memset(&o, 0, sizeof(o));
    for(int t = 0; t <= N; t++){ // (bounds: u, then the position; at stage 0 also x_0)
        nxs[t] = nx; nus[t] = t < N ? nu : 0; nbs[t] = nus[t] + (t == 0 ? nx : 1);
        A[t] = Ad; B[t] = Bd; Q[t] = Qd; R[t] = Rd;
        idxb[t] = malloc(nbs[t]*sizeof(int)); lb[t] = malloc(nbs[t]*sizeof(double)); ub[t] = malloc(nbs[t]*sizeof(double));
        for(int i = 0; i < nbs[t]; i++){ idxb[t][i] = i; lb[t][i] = i < nus[t] ? -1 : -0.5; ub[t][i] = i < nus[t] ? 1 : 0.5; }
    }
    lb[0][nu] = ub[0][nu] = 2; lb[0][nu+1] = ub[0][nu+1] = 0; // x_0
    o.N = N; o.nx = nxs; o.nu = nus; o.nb = nbs; o.A = A; o.B = B; o.Q = Q; o.R = R; o.idxb = idxb; o.lb = lb; o.ub = ub;
    const int m = daqp_multistage_m(&o), nz = daqp_multistage_nz(&o);
    int *sense = calloc(m, sizeof(int));
    double *rho = malloc(m*sizeof(double)), *w = malloc(m*sizeof(double)), *z = malloc(nz*sizeof(double));
    double *z1 = malloc(nz*sizeof(double)), *z2 = malloc(nz*sizeof(double));
    for(int t = 1, c = nbs[0]; t <= N; c += nbs[t], t++) sense[c + nus[t]] = DAQP_SOFT;
    for(int i = 0; i < m; i++){ rho[i] = 0.1; w[i] = 2; }
    o.sense = sense;
    const int mm = DAQP_MULTISTAGE_UPDATE_MATRICES;
    int ok = 1, nfixed = 0, nfree = 0, k = 0, it_warm = 0;
    ext_mem = 0;
    DAQPWorkspace work; work.settings = NULL;
    DAQPResult res; memset(&res, 0, sizeof(res)); res.x = z; res.lam = NULL;
    if(setup_daqp_multistage(&o, &work) < 0) return 0;
    daqp_set_soft_weights(&work, rho, rho, w, w);
    daqp_multistage_update(&work, mm);
    daqp_multistage_solve(&res, &work);
    ok &= res.exitflag > 0;
    memcpy(z1, z, nz*sizeof(double));
    for(int i = 0; i < m; i++){
        if(!(work.sense[i] & DAQP_ACTIVE) || !(work.sense[i] & DAQP_SOFT)) continue;
        if(work.sense[i] & DAQP_SLACK_FIXED) nfixed++; else nfree++;
    }
    k = work.n_active;
    // The same problem, refactorized
    daqp_multistage_update(&work, mm);
    daqp_multistage_solve(&res, &work);
    it_warm = res.iter;
    ok &= res.exitflag > 0 && it_warm == 1;
    for(int i = 0; i < nz; i++) if(fabs(z[i] - z1[i]) > 1e-8) ok = 0;
    // Linear weights 0, refactorized, and in a new setup
    for(int i = 0; i < m; i++) w[i] = 0;
    daqp_set_soft_weights(&work, rho, rho, w, w);
    daqp_multistage_update(&work, mm);
    daqp_multistage_solve(&res, &work);
    ok &= res.exitflag > 0;
    memcpy(z2, z, nz*sizeof(double));
    free_daqp_multistage(&work);
    daqp_multistage_pack_free(&o);
    work.settings = NULL;
    if(setup_daqp_multistage(&o, &work) < 0) return 0;
    daqp_set_soft_weights(&work, rho, rho, w, w);
    daqp_multistage_update(&work, mm);
    daqp_multistage_solve(&res, &work);
    ok &= res.exitflag > 0;
    for(int i = 0; i < nz; i++) if(fabs(z[i] - z2[i]) > 1e-6) ok = 0;
    free_daqp_multistage(&work);
    daqp_multistage_pack_free(&o);
    printf("soft warm start: %d active (%d soft with fixed slacks, %d with free), %d iteration(s) after a refactorization\n",
            k, nfixed, nfree, it_warm);
    for(int t = 0; t <= N; t++){ free(idxb[t]); free(lb[t]); free(ub[t]); }
    free(sense); free(rho); free(w); free(z); free(z1); free(z2);
    return ok && nfixed > 0 && nfree > 0 && k >= 16; // (else it would not test the activation at once)
}

// DAQP_MULTISTAGE_UPDATE_SHIFT: a closed loop of a double integrator (input bounds, soft position bounds with a linear
// weight, a general constraint at stages 1..N, no input at the last stage) whose initial state is the next state of the
// latest solution (and disturbed at every 5th step), solved with and without the shift of the warm start. The
// solutions must agree, and the shifted warm starts need fewer iterations
static int test_shift(void){
    enum { N = 25, nx = 2, nu = 1, steps = 40 };
    double Ad[4] = {1, 0.1, 0, 1}, Bd[2] = {0.005, 0.1}, Qd[4] = {1, 0, 0, 1}, Rd[1] = {0.01}, Cd[2] = {1, 0.5}, Dd[1] = {0.2};
    int nxs[N+1], nus[N+1], nbs[N+1], ngs[N+1], *idxb[N+1];
    double *A[N+1], *B[N+1], *Q[N+1], *R[N+1], *C[N+1], *Dm[N+1], *lb[N+1], *ub[N+1], *lg[N+1], *ug[N+1];
    DAQPMultistage o; memset(&o, 0, sizeof(o));
    for(int t = 0; t <= N; t++){ // (bounds: u, then the position; at stage 0 also x_0)
        nxs[t] = nx; nus[t] = t < N ? nu : 0; nbs[t] = nus[t] + (t == 0 ? nx : 1); ngs[t] = t > 0;
        A[t] = Ad; B[t] = Bd; Q[t] = Qd; R[t] = Rd; C[t] = Cd; Dm[t] = Dd;
        idxb[t] = malloc(nbs[t]*sizeof(int)); lb[t] = malloc(nbs[t]*sizeof(double)); ub[t] = malloc(nbs[t]*sizeof(double));
        lg[t] = malloc(sizeof(double)); ug[t] = malloc(sizeof(double));
        for(int i = 0; i < nbs[t]; i++){ idxb[t][i] = i; lb[t][i] = i < nus[t] ? -1 : -0.5; ub[t][i] = i < nus[t] ? 1 : 0.5; }
        lg[t][0] = -2; ug[t][0] = 2;
    }
    double x[nx]; x[0] = 1.5; x[1] = 0;
    for(int i = 0; i < nx; i++) lb[0][nu+i] = ub[0][nu+i] = x[i];
    o.N = N; o.nx = nxs; o.nu = nus; o.nb = nbs; o.ng = ngs; o.A = A; o.B = B; o.Q = Q; o.R = R; o.C = C; o.D = Dm;
    o.idxb = idxb; o.lb = lb; o.ub = ub; o.lg = lg; o.ug = ug;
    const int m = daqp_multistage_m(&o), nz = daqp_multistage_nz(&o);
    int *sense = calloc(m, sizeof(int));
    double *rho = malloc(m*sizeof(double)), *w = malloc(m*sizeof(double)), *z = malloc(nz*sizeof(double)), *zs = malloc(nz*sizeof(double));
    for(int t = 1, c = nbs[0] + ngs[0]; t <= N; c += nbs[t] + ngs[t], t++) sense[c + nus[t]] = DAQP_SOFT;
    for(int i = 0; i < m; i++){ rho[i] = 0.1; w[i] = 1; }
    o.sense = sense;
    ext_mem = 0;
    DAQPWorkspace ws[2]; long it[2] = {0, 0};
    int ok = 1;
    for(int v = 0; v < 2; v++){
        ws[v].settings = NULL;
        if(setup_daqp_multistage(&o, ws+v) < 0) return 0;
        daqp_set_soft_weights(ws+v, rho, rho, w, w);
    }
    for(int k = 0; k < steps && ok; k++){
        for(int i = 0; i < nx; i++) lb[0][nu+i] = ub[0][nu+i] = x[i];
        for(int v = 0; v < 2; v++){
            DAQPResult res; memset(&res, 0, sizeof(res)); res.x = z; res.lam = NULL;
            daqp_multistage_update(ws+v, v && k > 0 ? DAQP_MULTISTAGE_UPDATE_SHIFT : 0);
            daqp_multistage_solve(&res, ws+v);
            if(res.exitflag <= 0) ok = 0;
            if(k > 0) it[v] += res.iter;
            if(v == 0) memcpy(zs, z, nz*sizeof(double));
            else for(int i = 0; i < nz; i++) if(fabs(z[i] - zs[i]) > 1e-6) ok = 0;
        }
        for(int i = 0; i < nx; i++) x[i] = zs[nx+nu+i] + (k % 5 == 4 ? 0.3 : 0); // (x_1, disturbed at every 5th step)
    }
    for(int v = 0; v < 2; v++) free_daqp_multistage(ws+v);
    daqp_multistage_pack_free(&o);
    printf("shifted warm starts (%d solves): iterations without the shift %ld, with DAQP_MULTISTAGE_UPDATE_SHIFT %ld\n", steps-1, it[0], it[1]);
    for(int t = 0; t <= N; t++){ free(idxb[t]); free(lb[t]); free(ub[t]); free(lg[t]); free(ug[t]); }
    free(sense); free(rho); free(w); free(z); free(zs);
    return ok && it[1] < it[0];
}

// Branch and bound with the candidates of the latest solution (its binary assignment, and that one stage later): a closed
// loop of a double integrator with a forward and a backward thruster with a minimum thrust when on (u in {0} U [0.6, 1]:
// binary d, 0.6 d <= u <= d, a cost on d), disturbed at every 5th step, updated with DAQP_MULTISTAGE_UPDATE_SHIFT at
// every other step. Against a new setup at every step (no candidate), the objectives must agree, and the candidates
// must need fewer nodes, counting those of the candidates themselves
static int test_bnb_candidate(void){
    enum { N = 30, nx = 2, nu = 4, steps = 30 };
    const double dt = 0.2;
    double A[4] = {1, dt, 0, 1}, B[8] = {0.5*dt*dt, -0.5*dt*dt, 0, 0, dt, -dt, 0, 0}, Q[4] = {1, 0, 0, 0.1};
    double R[16] = {0.01,0,0,0, 0,0.01,0,0, 0,0,1e-4,0, 0,0,0,1e-4}, r[4] = {0, 0, 0.005, 0.005};
    double C[8] = {0}, Dm[16] = {1,0,-1,0, 0,1,0,-1, -1,0,0.6,0, 0,-1,0,0.6}; // u - d <= 0, 0.6 d - u <= 0
    int nxs[N+1], nus[N+1], nbs[N+1], ngs[N+1], *idxb[N+1];
    double *As[N+1], *Bs[N+1], *Qs[N+1], *Rs[N+1], *rs[N+1], *Cs[N+1], *Ds[N+1], *lb[N+1], *ub[N+1], *lg[N+1], *ug[N+1];
    DAQPMultistage o; memset(&o, 0, sizeof(o));
    for(int t = 0; t <= N; t++){ // (bounds: u and d, then the position; at stage 0 also x_0)
        nxs[t] = nx; nus[t] = t < N ? nu : 0; nbs[t] = t == 0 ? nu + nx : nus[t] + 1; ngs[t] = t < N ? 4 : 0;
        As[t] = A; Bs[t] = B; Qs[t] = Q; Rs[t] = R; rs[t] = r; Cs[t] = C; Ds[t] = Dm;
        idxb[t] = malloc(nbs[t]*sizeof(int)); lb[t] = malloc(nbs[t]*sizeof(double)); ub[t] = malloc(nbs[t]*sizeof(double));
        lg[t] = malloc(4*sizeof(double)); ug[t] = malloc(4*sizeof(double));
        for(int i = 0; i < 4; i++){ lg[t][i] = -DAQP_INF; ug[t][i] = 0; }
        for(int i = 0; i < nbs[t]; i++){
            idxb[t][i] = t == 0 || i < nus[t] ? i : nus[t];
            lb[t][i] = i < nus[t] ? 0 : -5; ub[t][i] = i < nus[t] ? 1 : 5;
        }
    }
    double x[nx]; x[0] = 3; x[1] = 0;
    for(int i = 0; i < nx; i++) lb[0][nu+i] = ub[0][nu+i] = x[i];
    o.N = N; o.nx = nxs; o.nu = nus; o.nb = nbs; o.ng = ngs; o.A = As; o.B = Bs; o.Q = Qs; o.R = Rs; o.r = rs;
    o.C = Cs; o.D = Ds; o.lg = lg; o.ug = ug; o.idxb = idxb; o.lb = lb; o.ub = ub;
    const int m = daqp_multistage_m(&o), nz = daqp_multistage_nz(&o);
    int *sense = calloc(m, sizeof(int));
    for(int t = 0, c = 0; t < N; c += nbs[t] + ngs[t], t++) sense[c+2] = sense[c+3] = DAQP_BINARY;
    o.sense = sense;
    double *z = malloc(nz*sizeof(double)), *zs = malloc(nz*sizeof(double));
    long nodes[2] = {0, 0};
    int ok = 1;
    ext_mem = 0;
    DAQPWorkspace work; work.settings = NULL;
    if(setup_daqp_multistage(&o, &work) < 0) return 0;
    for(int k = 0; k < steps && ok; k++){
        DAQPResult res; memset(&res, 0, sizeof(res)); res.x = z;
        for(int i = 0; i < nx; i++) lb[0][nu+i] = ub[0][nu+i] = x[i];
        daqp_multistage_update(&work, k % 2 ? DAQP_MULTISTAGE_UPDATE_SHIFT : 0); // (the flag does not apply to BnB)
        daqp_multistage_solve(&res, &work);
        const double f = res.fval;
        ok &= res.exitflag > 0;
        nodes[0] += work.bnb->nodecount;
        memcpy(zs, z, nz*sizeof(double));
        DAQPWorkspace ref; ref.settings = NULL; // (a new setup: no candidate)
        if(setup_daqp_multistage(&o, &ref) < 0){ ok = 0; break; }
        daqp_multistage_update(&ref, 0);
        daqp_multistage_solve(&res, &ref);
        ok &= res.exitflag > 0 && fabs(res.fval - f) <= 1e-6*(1 + fabs(f));
        nodes[1] += ref.bnb->nodecount;
        free_daqp_multistage(&ref);
        packed = &o; // (test_setup repacked o for ref)
        for(int i = 0; i < nx; i++) x[i] = zs[nx+nu+i] + (k % 5 == 4 ? (i == 0 ? 1.5 : -0.5) : 0);
    }
    free_daqp_multistage(&work);
    daqp_multistage_pack_free(&o);
    printf("branch and bound (%d solves): nodes with the candidates of the latest solution %ld, without %ld\n", steps, nodes[0], nodes[1]);
    for(int t = 0; t <= N; t++){ free(idxb[t]); free(lb[t]); free(ub[t]); free(lg[t]); free(ug[t]); }
    free(sense); free(z); free(zs);
    return ok && nodes[0] < nodes[1];
}

int main(int argc, char** argv){
    int ntests = argc > 1 ? atoi(argv[1]) : 200, fails = 0, compared = 0;
    // Problems with soft constraints whose hard constraints are infeasible are
    // not compared: the infeasibility is then not reliably detected (by dense
    // DAQP either), so only how often it is not is counted
    int infeas = 0, wrong_ocp = 0, wrong_dense = 0, cert = 0, ngain = 0, ngain_det = 0;
#define SKIP_INFEASIBLE(ed, eo) (g.rl != NULL && !hf ? (infeas++, wrong_dense += (ed) > 0, wrong_ocp += (eo) > 0, 1) : \
        DENSE_BIN_WRONG(ed) || DENSE_FAILED(ed, eo))
    // Dense DAQP's branch and bound can return a binary constraint strictly between its bounds when equalities fix
    // it (a known issue of the dense reference, fixed separately): not compared
    int dense_bin_wrong = 0, nshared = 0, nshared_realloc = 0;
#define DENSE_BIN_WRONG(ed) ((feat & F_BIN) && (ed) > 0 && !binaries_ok(&g, zd) ? (dense_bin_wrong++, 1) : 0)
    // Dense DAQP can reach its iteration limit or cycle on an ill-conditioned condensed QP: not compared
    int dense_failed = 0;
#define DENSE_FAILED(ed, eo) (((ed) == DAQP_EXIT_ITERLIMIT || (ed) == DAQP_EXIT_CYCLE) && (eo) > 0 ? (dense_failed++, 1) : 0)
    srand(argc > 2 ? atoi(argv[2]) : 1);
    for(int k = 0; k < ntests; k++){
        Gen g;
        int N = 3 + rand()%20, nx = 1 + rand()%5, nu = 1 + rand()%3;
        int feat = 0;
        if(k%2 == 1) feat |= F_COUPLING;
        if(k%5 == 2 || k%5 == 4) feat |= F_EQ;
        if(k%10 == 7) feat |= F_EQ | F_EQ_DUP;
        if(k%20 == 9) feat |= F_EQ | F_EQ_BAD;
        if(k%7 == 3) feat |= F_EQ_COUPLING;
        if(k%4 == 1) feat |= F_X0FREE;
        if(k%4 == 3) feat |= F_X0PART;
        if(k%13 == 5 && !(feat & (F_X0FREE|F_X0PART))) feat |= F_X0BAD;
        if(k%3 == 2) feat |= F_XFIRST;
        if(k%3 == 0) feat |= F_SOFT;
        if(k%12 == 4 || k%12 == 8) feat |= F_PSD; // (without soft constraints and a free x_0)
        if(k%12 == 2 || k%12 == 10) feat |= F_BIN; // (also)
        if(k%24 == 2 || k%24 == 10) feat |= F_SOFT; // (with binaries)
        if(k%6 == 1 || k%6 == 4) feat |= F_GSPLIT;
        if(k%11 == 6) feat = F_EQ_GAIN; // (alone)
        if(k%8 == 6 && k%3 != 1 && !(feat & F_EQ_GAIN)) feat |= F_SHARED; // (equal dimensions)
        generate(&g, N, nx, nu, k%3 == 1, feat);
        if(feat & F_SHARED){ share_matrices(&g.prob, feat); nshared++; }
        double *zd = malloc(g.nz*sizeof(double)), *zo = malloc(g.nz*sizeof(double)), fd, td;
        int hf = g.rl == NULL || hard_feasible(&g);
        int ed = solve_dense(&g, zd, &fd, &td);
        // (binary constraints: also against all fixings of them, without branch and bound)
        double fbrute = 0;
        const int bin_feasible = (feat & F_BIN) && binaries_brute_force(&g, &fbrute);
        if(feat & F_BIN) solve_dense(&g, zd, &fd, &td); // (the condensed QP of the problem again)

        DAQPWorkspace work; work.settings = NULL;
        ext_mem = (feat & F_BIN) ? (k/12)%2 : k%2 ? 1 + (k/2)%2 : 0; // (every other test in given memory)
        int es = setup_daqp_multistage(&g.prob, &work);
        if(es >= 0 && g.rl != NULL) daqp_set_soft_weights(&work, g.rl, g.ru, g.wl, g.wu);
        double* lam = malloc((g.m+1)*sizeof(double));
        DAQPResult res; res.x = zo; res.lam = lam;
        const int det0 = elim_detected;
        if(es >= 0) ms_solve(&res, &work); else res.exitflag = es;
        if(elim_detected > det0 && !(feat & F_EQ_GAIN)) printf("test %d: inaccurate elimination detected (feat %d)\n", k, feat);
        ngain_det += (feat & F_EQ_GAIN) && elim_detected > det0;
        ngain += (feat & F_EQ_GAIN) != 0;
        int eo = res.exitflag;

        int ok = 1;
        if(SKIP_INFEASIBLE(ed, eo));
        else if(feat & F_PSD){
            ok = psd_ok(&g, ed, eo, zd, zo, lam);
            cert++;
            if(!ok) printf("test %d (N=%d nx=%d nu=%d): F_PSD exitflags dense %d prob %d, not certified\n", k, N, nx, nu, ed, eo);
        }
        else if(ed > 0 && eo > 0){
            double dz = 0;
            for(int i = 0; i < g.nz; i++) if(fabs(zd[i]-zo[i]) > dz) dz = fabs(zd[i]-zo[i]);
            double viol = max_violation(&g, zo), kkt = kkt_residual(&g, zo, lam);
            double* pi = malloc((g.nz+1)*sizeof(double));
            daqp_multistage_duals(&work, NULL, pi);
            double kktf = kkt_full(&g, zo, lam, pi);
            free(pi);
            if(feat & F_BIN){ // (no multipliers after branch and bound)
                const double fo = objective_soft(&g, zo); // (the enumeration is of this objective)
                ok = agree(&g, fd, res.fval, zd, zo) && viol < 1e-6 && binaries_ok(&g, zo)
                    && bin_feasible && fabs(fbrute - fo) <= 1e-6*(1+fabs(fbrute));
                if(!binaries_ok(&g, zo)) printf("test %d: binary constraints not at a bound\n", k);
                if(!bin_feasible || fabs(fbrute - fo) > 1e-6*(1+fabs(fbrute)))
                    printf("test %d: branch and bound %.10e, enumeration %.10e (%s)\n", k, fo, fbrute, bin_feasible ? "feasible" : "infeasible");
            }
            else{
                // (F_EQ_GAIN: the large multipliers of the equalities are accurate to a few digits less)
                const double ktol = (feat & F_EQ_GAIN) ? 1e-5 : 1e-6;
                ok = agree(&g, fd, res.fval, zd, zo) && viol < 1e-6 && kkt < ktol && kktf < ktol;
                if(kkt >= ktol) printf("test %d: KKT residual of the multipliers %.2e (feat %d)\n", k, kkt, feat);
                if(kktf >= ktol) printf("test %d: full KKT residual with the costates %.2e (feat %d)\n", k, kktf, feat);
            }
            compared++;
            if(!ok) printf("test %d (N=%d nx=%d nu=%d): |dz|=%.2e df=%.2e viol=%.2e (dense viol %.2e, feat %d)\n",
                    k, N, nx, nu, dz, fd-res.fval, viol, max_violation(&g, zd), feat);
        }
        else{
            ok = (ed > 0) == (eo > 0) && !((feat & F_BIN) && hf && bin_feasible);
            if(!ok) printf("test %d (N=%d nx=%d nu=%d): exitflags dense %d prob %d%s\n", k, N, nx, nu, ed, eo,
                    (feat & F_BIN) && bin_feasible ? " (enumeration of the binaries: feasible)" : "");
        }

        // Warm start: perturb the fixed components of x_0 and resolve from the
        // previous working set
        if(es >= 0){
            DAQPMultistage* o = &g.prob;
            for(int i = 0; i < o->nb[0]; i++)
                if(o->idxb[0][i] >= (o->N > 0 ? o->nu[0] : 0) && o->lb[0][i] == o->ub[0][i])
                    o->lb[0][i] = o->ub[0][i] = o->lb[0][i] + 0.3*randu();
            // A new value of the eliminated input equality (which the x0-only
            // update has to detect)
            if(feat & F_EQ){
                const int t3 = o->N/3, i3 = g.ubase[t3];
                o->lb[t3][i3] = o->ub[t3][i3] = o->lb[t3][i3] + 0.05;
            }
            daqp_multistage_update(&work, 0);
            ms_solve(&res, &work);
            hf = g.rl == NULL || hard_feasible(&g);
            int ed2 = solve_dense(&g, zd, &fd, &td);
            if(SKIP_INFEASIBLE(ed2, res.exitflag));
            else if(feat & F_PSD){
                if(!psd_ok(&g, ed2, res.exitflag, zd, zo, lam)){ ok = 0; printf("test %d warm start: F_PSD exitflags dense %d prob %d, not certified\n", k, ed2, res.exitflag); }
            }
            else if(ed2 > 0 && res.exitflag > 0){
                int ok2 = agree(&g, fd, res.fval, zd, zo) && max_violation(&g, zo) < 1e-6;
                if(!ok2){ ok = 0; printf("test %d warm start: df=%.2e\n", k, fd-res.fval); }
            }
            else if((ed2 > 0) != (res.exitflag > 0)){
                ok = 0; printf("test %d warm start: exitflags dense %d prob %d\n", k, ed2, res.exitflag);
            }
            // New q, r, and b at one stage (which the update has to detect also
            // with the mask of x_0 only)
            {
                DAQPMultistage* o = &g.prob;
                const int t2 = o->N/2;
                for(int i = 0; i < o->nx[t2]; i++) o->q[t2][i] += 0.2*randu();
                for(int i = 0; t2 < o->N && i < o->nu[t2]; i++) o->r[t2][i] += 0.2*randu();
                for(int i = 0; t2 < o->N && i < o->nx[t2+1]; i++) o->b[t2][i] += 0.05*randu();
                daqp_multistage_update(&work, 0);
                ms_solve(&res, &work);
                hf = g.rl == NULL || hard_feasible(&g);
                ed2 = solve_dense(&g, zd, &fd, &td);
                if(SKIP_INFEASIBLE(ed2, res.exitflag));
                else if(feat & F_PSD){
                    if(!psd_ok(&g, ed2, res.exitflag, zd, zo, lam)){ ok = 0; printf("test %d vector update: F_PSD exitflags dense %d prob %d, not certified\n", k, ed2, res.exitflag); }
                }
                else if(ed2 > 0 && res.exitflag > 0){
                    int ok2 = agree(&g, fd, res.fval, zd, zo) && max_violation(&g, zo) < 1e-6;
                    if(!ok2){ ok = 0; printf("test %d vector update at stage %d: df=%.2e\n", k, t2, fd-res.fval); }
                }
                else if((ed2 > 0) != (res.exitflag > 0)){
                    ok = 0; printf("test %d vector update: exitflags dense %d prob %d\n", k, ed2, res.exitflag);
                }
            }
            // New matrices: scale the input weights and resolve (warm started)
            for(int t = 0; t < g.prob.N; t++)
                if(t == 0 || g.prob.R[t] != g.prob.R[t-1]) // (shared: once)
                    for(int i = 0; i < g.prob.nu[t]*g.prob.nu[t]; i++) g.prob.R[t][i] *= 2.0;
            if((feat & F_SHARED) && g.prob.A[1] == g.prob.A[0]){ // Stage 1 no longer shares A (a new block)
                const int na = g.prob.nx[2]*g.prob.nx[1];
                double* A1 = malloc(na*sizeof(double));
                for(int i = 0; i < na; i++) A1[i] = 1.01*g.prob.A[1][i];
                g.prob.A[1] = A1;
                nshared_realloc++;
            }
            int eu = daqp_multistage_update(&work, DAQP_MULTISTAGE_UPDATE_MATRICES);
            if(eu >= 0) ms_solve(&res, &work); else res.exitflag = eu;
            hf = g.rl == NULL || hard_feasible(&g);
            ed2 = solve_dense(&g, zd, &fd, &td);
            if(SKIP_INFEASIBLE(ed2, res.exitflag));
            else if(feat & F_PSD){
                if(!psd_ok(&g, ed2, res.exitflag, zd, zo, lam)){ ok = 0; printf("test %d matrix update: F_PSD exitflags dense %d prob %d, not certified\n", k, ed2, res.exitflag); }
            }
            else if(ed2 > 0 && res.exitflag > 0){
                int ok2 = agree(&g, fd, res.fval, zd, zo) && max_violation(&g, zo) < 1e-6;
                if(!ok2){ ok = 0; printf("test %d matrix update: df=%.2e\n", k, fd-res.fval); }
            }
            else if((ed2 > 0) != (res.exitflag > 0)){
                ok = 0; printf("test %d matrix update: exitflags dense %d prob %d\n", k, ed2, res.exitflag);
            }
            // Relaxed inequalities (equalities and a fixed x_0 kept): the working set empties and
            // z_lqr becomes the solution, after a solve with a nonzero v (product with v = 0).
            // (Not for binaries, whose meaning changes, nor F_PSD)
            if(!(feat & (F_BIN|F_PSD))){
                DAQPMultistage* o = &g.prob;
                for(int t = 0; t <= o->N; t++){
                    for(int i = 0; i < o->nb[t]; i++)
                        if(o->lb[t][i] != o->ub[t][i]){ o->lb[t][i] -= 1e3; o->ub[t][i] += 1e3; }
                    for(int i = 0; i < o->ng[t]; i++)
                        if(o->lg[t][i] != o->ug[t][i]){ o->lg[t][i] -= 1e3; o->ug[t][i] += 1e3; }
                    for(int i = 0; o->nc != NULL && t < o->N && i < o->nc[t]; i++)
                        if(o->lc[t][i] != o->uc[t][i]){ o->lc[t][i] -= 1e3; o->uc[t][i] += 1e3; }
                }
                int er = daqp_multistage_update(&work, 0);
                if(er >= 0) ms_solve(&res, &work); else res.exitflag = er;
                hf = g.rl == NULL || hard_feasible(&g);
                ed2 = solve_dense(&g, zd, &fd, &td);
                if(SKIP_INFEASIBLE(ed2, res.exitflag));
                else if((ed2 == DAQP_EXIT_OVERDETERMINED_INITIAL || ((feat & F_EQ_GAIN) && ed2 == DAQP_EXIT_INFEASIBLE))
                        && res.exitflag > 0){ // (dense reference fails, also on the ill-conditioned condensed F_EQ_GAIN QP: certify by KKT)
                    const double viol = max_violation(&g, zo), kkt = kkt_residual(&g, zo, lam);
                    if(!(viol < 1e-6 && kkt < 1e-6)){ ok = 0; printf("test %d relaxed (KKT): viol %.2e kkt %.2e\n", k, viol, kkt); }
                }
                else if(ed2 > 0 && res.exitflag > 0){
                    int ok2 = agree(&g, fd, res.fval, zd, zo) && max_violation(&g, zo) < 1e-6;
                    // (an ill-conditioned condensed QP can make the dense reference inaccurate: certify by KKT)
                    if(!ok2 && res.soft_slack == 0) ok2 = certified(&g, zo, lam);
                    if(!ok2){ ok = 0; printf("test %d relaxed: df=%.2e\n", k, fd-res.fval); }
                }
                else if((ed2 > 0) != (res.exitflag > 0)){
                    ok = 0; printf("test %d relaxed: exitflags dense %d prob %d\n", k, ed2, res.exitflag);
                }
            }
            // Changed equality pattern: a new equality after setup (kept in the working
            // set), and the eliminated input equality becomes an inequality
            {
                int tN = g.prob.N-1, t3 = g.prob.N/3, iN = g.ubase[tN], i3 = g.ubase[t3];
                g.prob.lb[tN][iN] = g.prob.ub[tN][iN] = 0.5*(g.prob.lb[tN][iN]+g.prob.ub[tN][iN]);
                if((feat & F_EQ) && t3 < g.prob.N) g.prob.ub[t3][i3] = g.prob.lb[t3][i3] + 0.2;
                int ev = daqp_multistage_update(&work, 0);
                if(ev >= 0) ms_solve(&res, &work); else res.exitflag = ev;
                hf = g.rl == NULL || hard_feasible(&g);
                ed2 = solve_dense(&g, zd, &fd, &td);
                if(SKIP_INFEASIBLE(ed2, res.exitflag));
                else if(feat & F_PSD){
                    if(!psd_ok(&g, ed2, res.exitflag, zd, zo, lam)){ ok = 0; printf("test %d equality change: F_PSD exitflags dense %d prob %d, not certified\n", k, ed2, res.exitflag); }
                }
                else if(ed2 > 0 && res.exitflag > 0){
                    int ok2 = agree(&g, fd, res.fval, zd, zo) && max_violation(&g, zo) < 1e-6;
                    // (an ill-conditioned condensed QP can make the dense reference inaccurate: certify by KKT)
                    if(!ok2 && res.soft_slack == 0) ok2 = certified(&g, zo, lam);
                    if(!ok2){ ok = 0; printf("test %d equality change: df=%.2e\n", k, fd-res.fval); }
                }
                else if((ed2 > 0) != (res.exitflag > 0)){
                    ok = 0; printf("test %d equality change: exitflags dense %d prob %d (N=%d, prob viol %.2e, dense viol %.2e)\n",
                            k, ed2, res.exitflag, g.prob.N, max_violation(&g, zo), max_violation(&g, zd));
                }
            }
            // Alternating data (as a closed loop that chatters at saturated inputs): x_0 alternates between two
            // values, with refactorizations (odd tests) or x_0 updates. The solutions agree with dense DAQP
            if(ok && !(feat & (F_BIN|F_PSD|F_SOFT|F_EQ_BAD|F_X0FREE|F_X0BAD|F_EQ_GAIN))){
                DAQPMultistage* o = &g.prob;
                const int nb0 = o->nb[0], nu0 = o->N > 0 ? o->nu[0] : 0;
                double *xa = malloc((nb0+1)*sizeof(double)), *xb = malloc((nb0+1)*sizeof(double));
                double *za = malloc(g.nz*sizeof(double)), *zb = malloc(g.nz*sizeof(double)), fa = 0, fb = 0;
                int fixed = 0, ea, eb;
                for(int i = 0; i < nb0; i++){
                    xa[i] = xb[i] = o->lb[0][i];
                    if(o->idxb[0][i] >= nu0 && o->lb[0][i] == o->ub[0][i]){ xb[i] = xa[i] + 4.0*randu(); fixed = 1; }
                }
#define SET_X0(x) for(int i = 0; i < nb0; i++) if(o->idxb[0][i] >= nu0 && o->lb[0][i] == o->ub[0][i]) o->lb[0][i] = o->ub[0][i] = (x)[i]
                SET_X0(xb); eb = solve_dense(&g, zb, &fb, &td);
                SET_X0(xa); ea = solve_dense(&g, za, &fa, &td);
                for(int st = 0; fixed && ea > 0 && eb > 0 && st < 8; st++){
                    SET_X0(st%2 ? xb : xa);
                    daqp_multistage_update(&work, k%2 ? DAQP_MULTISTAGE_UPDATE_MATRICES : 0);
                    ms_solve(&res, &work);
                    if(res.exitflag <= 0 || !agree(&g, st%2 ? fb : fa, res.fval, st%2 ? zb : za, zo)){
                        ok = 0; printf("test %d alternating x0, solve %d: exitflag %d, df=%.2e\n",
                                k, st, res.exitflag, (st%2 ? fb : fa) - res.fval);
                        break;
                    }
                }
#undef SET_X0
                free(xa); free(xb); free(za); free(zb);
            }
            // A fixed x_0 that is no longer fixed requires a new setup
            if(!(feat & (F_X0FREE|F_X0PART)) && g.prob.nx[0] > 0){
                DAQPMultistage* o = &g.prob;
                for(int i = 0; i < o->nb[0]; i++)
                    if(o->idxb[0][i] >= (o->N > 0 ? o->nu[0] : 0)) o->ub[0][i] = o->lb[0][i] + 1.0;
                int ex = daqp_multistage_update(&work, 0);
                if(ex != DAQP_EXIT_UNSUPPORTED){ ok = 0; printf("test %d: unfixed x0 gives %d\n", k, ex); }
            }
            free_daqp_multistage(&work);
        }
        ext_release();
        if((feat & F_BIN) && g.nbin > 1 && !test_binary_horizon(&g)){ ok = 0; printf("test %d: binary_horizon\n", k); }
        daqp_multistage_pack_free(&g.prob);
        fails += !ok;
        free(zd); free(zo); free(lam);
    }
    // (with DAQP_MULTISTAGE_UPDATE_SHIFT: the alternation still decides, the shifted working sets do not fit)
    const int chat_ok = test_chattering(0, 0) & test_chattering(1, 0) & test_chattering(0, DAQP_MULTISTAGE_UPDATE_SHIFT);
    const int soft_ok = test_soft_warm();
    const int shift_ok = test_shift();
    const int bnb_ok = test_bnb_candidate();
    const int bnbobj_ok = test_bnb_objective();
    printf("%d/%d passed (%d with optimal solutions compared)\n", ntests-fails, ntests, compared);
    if(infeas > 0)
        printf("not compared: %d solves with soft constraints and infeasible hard constraints"
                " (reported as solved: prob %d, dense %d)\n", infeas, wrong_ocp, wrong_dense);
    if(cert > 0) printf("F_PSD solves certified by the KKT conditions (instead of compared): %d\n", cert);
    if(dual_mismatch > 0) printf("working-set multipliers changed by daqp_multistage_duals: %d\n", dual_mismatch);
    if(ngain > 0) printf("F_EQ_GAIN: inaccurate elimination detected (and recovered) in %d of %d problems\n", ngain_det, ngain);
    if(dense_bin_wrong > 0) printf("not compared: %d solves where dense branch and bound left a binary between its bounds\n", dense_bin_wrong);
    if(dense_failed > 0) printf("not compared: %d solves where dense DAQP reached its iteration limit or cycled\n", dense_failed);
    printf("F_SHARED: %d problems, %d with a change of the sharing\n", nshared, nshared_realloc);
    return fails > 0 || dual_mismatch > 0 || (ngain >= 20 && ngain_det == 0) || !chat_ok || !soft_ok || !shift_ok || !bnb_ok || !bnbobj_ok;
}
