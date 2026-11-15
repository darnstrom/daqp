/*
 * Tests of the solver logic of the S-function (without Simulink). Build from
 * this directory with the sources in ../../../src and the includes in
 * ../../../include and .. (see .github/workflows/ci_simulink.yml)
 */
#include "daqp_sfunc_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int fails = 0;
#define CHECK(c, ...) do{ if(!(c)){ fails++; printf("FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } }while(0)

static double rnd(){ return 2.0*rand()/RAND_MAX-1.0; }

// Reference: one-shot solve with daqp_quadprog
static int ref_solve(DAQPProblem* qp, c_float* x, c_float* lam, c_float* fval){
    DAQPResult r = {0};
    r.x = x; r.lam = lam;
    DAQPProblem q = *qp; // daqp_quadprog may modify sense
    int* sense = malloc(qp->m*sizeof(int));
    for(int i = 0; i < qp->m; i++) sense[i] = qp->sense[i];
    q.sense = sense;
    daqp_quadprog(&r,&q,NULL);
    free(sense);
    *fval = r.fval;
    return r.exitflag;
}

static double maxdiff(c_float* a, c_float* b, int n){
    double d = 0;
    for(int i = 0; i < n; i++) d = fmax(d, fabs(a[i]-b[i]));
    return d;
}

// Random QP data in place (fixed H/A if keep_HA)
static void random_qp(DAQPProblem* qp, int keep_HA, double scale_b){
    int n = qp->n, m = qp->m, ms = qp->ms, mg = m-ms;
    if(!keep_HA){
        // H = M'M + I
        double* M = malloc(n*n*sizeof(double));
        for(int i = 0; i < n*n; i++) M[i] = rnd();
        for(int i = 0; i < n; i++) for(int j = 0; j < n; j++){
            double s = (i==j);
            for(int k = 0; k < n; k++) s += M[k*n+i]*M[k*n+j];
            qp->H[i*n+j] = s;
        }
        free(M);
        for(int i = 0; i < mg*n; i++) qp->A[i] = rnd();
    }
    for(int i = 0; i < n; i++) qp->f[i] = 5*rnd();
    for(int i = 0; i < m; i++){
        qp->bupper[i] = scale_b*(1+rnd());
        qp->blower[i] = -scale_b*(1+rnd());
    }
}

static void test_random(int warm_start){
    int n = 20, m = 30, ms = 10, T = 30;
    srand(1);
    DAQPSfunc* s = daqp_sfunc_new(n,m,ms,NULL,0,warm_start,1000);
    c_float x[20], lam[30], fval;
    int iters = 0;
    for(int t = 0; t < T; t++){
        // H/A change every 10 steps (except WARM_FIXED, where they never change)
        random_qp(&s->qp, warm_start == DAQP_SFUNC_WARM_FIXED ? t > 0 : t % 10 != 0, 0.5);
        if(t == 15){ s->qp.sense[ms+3] = 5; s->qp.bupper[ms+3] = s->qp.blower[ms+3] = 0.1; }
        if(t == 20) s->qp.sense[ms+3] = 0;
        int flag = daqp_sfunc_solve(s);
        int rflag = ref_solve(&s->qp,x,lam,&fval);
        iters += s->res.iter;
        CHECK(flag == 1 && rflag == 1, "ws=%d t=%d flag %d ref %d", warm_start, t, flag, rflag);
        CHECK(maxdiff(s->res.x,x,n) < 1e-6, "ws=%d t=%d x diff %g", warm_start, t, maxdiff(s->res.x,x,n));
        CHECK(maxdiff(s->res.lam,lam,m) < 1e-6, "ws=%d t=%d lam diff %g", warm_start, t, maxdiff(s->res.lam,lam,m));
        CHECK(fabs(s->res.fval-fval) < 1e-6, "ws=%d t=%d fval %g vs %g", warm_start, t, s->res.fval, fval);
    }
    printf("random ws=%d: total iterations %d\n", warm_start, iters);
    daqp_sfunc_free(s);
}

// MPC-like: fixed H and A, slowly varying f and bounds => warm start should help
static int test_slow(int warm_start){
    int n = 30, m = 60, ms = 30, T = 40, iters = 0;
    srand(2);
    DAQPSfunc* s = daqp_sfunc_new(n,m,ms,NULL,0,warm_start,1000);
    random_qp(&s->qp,0,0.3);
    c_float f0[30], x[30], lam[60], fval;
    for(int i = 0; i < n; i++) f0[i] = s->qp.f[i];
    for(int t = 0; t < T; t++){
        for(int i = 0; i < n; i++) s->qp.f[i] = f0[i]*cos(0.05*t);
        int flag = daqp_sfunc_solve(s);
        int rflag = ref_solve(&s->qp,x,lam,&fval);
        iters += s->res.iter;
        CHECK(flag == 1 && rflag == 1, "slow ws=%d t=%d flag %d ref %d", warm_start, t, flag, rflag);
        CHECK(maxdiff(s->res.x,x,n) < 1e-6, "slow ws=%d t=%d x diff %g", warm_start, t, maxdiff(s->res.x,x,n));
    }
    printf("slow ws=%d: total iterations %d\n", warm_start, iters);
    daqp_sfunc_free(s);
    return iters;
}

// From issue #161
static void test_hierarchical(int warm_start){
    int n = 9, m = 12, ms = 9, bp[4] = {9,10,11,12};
    double A[27] = {1,-1,0,-5,-5,0,0,0,0, -1,1,0,0,0,0,0,0,0, -0.5,-0.5,0,0,0,0,0,0,0};
    double lb[12] = {0,0,0,-1.8,-1.8,-1.8,-1.8,-1.8,-1.8,0,0,0};
    double ub[12] = {4,4,4,1.8,1.8,1.8,1.8,1.8,1.8,0,0,0};
    double rhs[] = {0, 10, 20, 26, 30, 26, 10, 0, -15, 26};
    int T = sizeof(rhs)/sizeof(rhs[0]);
    DAQPSfunc* s = daqp_sfunc_new(n,m,ms,bp,4,warm_start,1000);
    for(int t = 0; t < T; t++){
        for(int i = 0; i < n*n; i++) s->qp.H[i] = (i%(n+1)==0);
        for(int i = 0; i < n; i++) s->qp.f[i] = 0;
        for(int i = 0; i < 27; i++) s->qp.A[i] = A[i];
        for(int i = 0; i < m; i++){ s->qp.blower[i] = lb[i]; s->qp.bupper[i] = ub[i]; s->qp.sense[i] = i >= 9 ? 5 : 0; }
        s->qp.blower[9] = s->qp.bupper[9] = rhs[t];
        daqp_sfunc_solve(s);
        // Reference: fresh solver
        DAQPSfunc* r = daqp_sfunc_new(n,m,ms,bp,4,DAQP_SFUNC_COLD,1000);
        for(int i = 0; i < n*n; i++) r->qp.H[i] = s->qp.H[i];
        for(int i = 0; i < 27; i++) r->qp.A[i] = A[i];
        for(int i = 0; i < m; i++){ r->qp.blower[i] = s->qp.blower[i]; r->qp.bupper[i] = s->qp.bupper[i]; r->qp.sense[i] = s->qp.sense[i]; }
        daqp_sfunc_solve(r);
        double Ax[3];
        for(int k = 0; k < 3; k++){ Ax[k] = 0; for(int j = 0; j < n; j++) Ax[k] += A[k*n+j]*s->res.x[j]; }
        // Level 1: x1-x2-5x4-5x5 can reach [-4-18, 4+18] = [-22, 22]
        double best1 = fmax(-22, fmin(22, rhs[t]));
        CHECK(s->res.exitflag > 0, "hqp ws=%d t=%d flag %d", warm_start, t, s->res.exitflag);
        CHECK(fabs(Ax[0]-best1) < 1e-6, "hqp ws=%d t=%d level 1: %g (best %g)", warm_start, t, Ax[0], best1);
        CHECK(maxdiff(s->res.x,r->res.x,n) < 1e-6, "hqp ws=%d t=%d x diff %g", warm_start, t, maxdiff(s->res.x,r->res.x,n));
        daqp_sfunc_free(r);
    }
    daqp_sfunc_free(s);
}

// Infeasible step followed by feasible steps, soft constraints toggled, dummy bound
static void test_misc(int warm_start){
    int n = 3, m = 4, ms = 3;
    DAQPSfunc* s = daqp_sfunc_new(n,m,ms,NULL,0,warm_start,1000);
    for(int t = 0; t < 6; t++){
        for(int i = 0; i < 9; i++) s->qp.H[i] = (i%4 == 0) ? 2 : 0;
        for(int i = 0; i < 3; i++) s->qp.A[i] = 1;
        for(int i = 0; i < n; i++) s->qp.f[i] = -4+i;
        for(int i = 0; i < m; i++){ s->qp.bupper[i] = 1; s->qp.blower[i] = -1; s->qp.sense[i] = 0; }
        s->qp.bupper[3] = 1; s->qp.blower[3] = -INFINITY; // x1+x2+x3 <= 1
        if(t == 2 || t == 3){ s->qp.blower[3] = 3.5; s->qp.bupper[3] = 5; } // infeasible with |x| <= 1
        if(t == 3) s->qp.sense[3] = 8; // ... unless soft
        int flag = daqp_sfunc_solve(s);
        c_float x[3], lam[4], fval;
        int rflag = ref_solve(&s->qp,x,lam,&fval);
        CHECK(flag == rflag, "misc ws=%d t=%d flag %d ref %d", warm_start, t, flag, rflag);
        CHECK(flag == (t == 2 ? -1 : t == 3 ? 2 : 1), "misc ws=%d t=%d flag %d", warm_start, t, flag);
        if(flag > 0) CHECK(maxdiff(s->res.x,x,n) < 1e-6, "misc ws=%d t=%d x diff %g", warm_start, t, maxdiff(s->res.x,x,n));
    }
    daqp_sfunc_free(s);

    // No constraints: a dummy bound that is infinite
    s = daqp_sfunc_new(2,1,1,NULL,0,warm_start,1000);
    for(int t = 0; t < 3; t++){
        s->qp.H[0] = 1; s->qp.H[3] = 2; s->qp.f[0] = 1+t; s->qp.f[1] = -2;
        s->qp.bupper[0] = INFINITY; s->qp.blower[0] = -INFINITY;
        int flag = daqp_sfunc_solve(s);
        CHECK(flag == 1 && fabs(s->res.x[0]+1+t) < 1e-9 && fabs(s->res.x[1]-1) < 1e-9,
                "unconstrained ws=%d t=%d flag %d x [%g %g]", warm_start, t, flag, s->res.x[0], s->res.x[1]);
    }
    daqp_sfunc_free(s);
}

int main(){
    for(int ws = 0; ws <= 2; ws++) test_random(ws);
    int it_cold = test_slow(0), it_warm = test_slow(1), it_fixed = test_slow(2);
    CHECK(it_warm < it_cold && it_fixed < it_cold, "warm start did not reduce iterations");
    for(int ws = 0; ws <= 2; ws++) test_hierarchical(ws);
    for(int ws = 0; ws <= 2; ws++) test_misc(ws);
    printf(fails ? "%d FAILURES\n" : "ALL PASSED\n", fails);
    return fails != 0;
}
