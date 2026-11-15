// Duplicate soft constraints (two identical rows a'x <= 0, each with the
// linear weight w and reciprocal quadratic weight rho) in random dense QPs,
// with many other active constraints and data scaled up to 1e5. The
// multipliers of the duplicates sum to lam in (w, 2w) at the optimum, so one
// of them is at the threshold w with a zero slack, where rounding made DAQP
// switch its slack back and forth until the iteration limit (fixed by
// DAQP_SLACK_SWITCHED in daqp_remove_blocking). Since 2w > lam, the duplicates
// are an exact penalty, and the solution is that of the QP with a'x <= 0 hard.
//
// Build: gcc -O2 -Iinclude -Ibuild/include test/soft_duplicates_test.c build/libdaqpstat.a -lm
#include "api.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static double randu(void){ return 2.0*rand()/RAND_MAX-1; }

int main(int argc, char **argv){
    int ntests = argc > 1 ? atoi(argv[1]) : 1000, fails = 0, iterlim = 0;
    srand(argc > 2 ? atoi(argv[2]) : 1);
    for(int k = 0; k < ntests; k++){
        const int n = 2 + rand()%39, nd = 2, mo = n/2 + rand()%(n+1), m = nd + mo;
        static double H[40*40], f[40], A[200*40], bu[200], bl[200], x[40], xr[40], lam[200];
        static double rl[200], ru[200], wl[200], wu[200];
        double scale = pow(10, 5*fabs(randu()));
        static int sense[200];
        memset(H, 0, sizeof(H));
        // H = I + sum of n random rank-one terms
        for(int i = 0; i < n; i++) H[i*n+i] = 1;
        for(int l = 0; l < n; l++){ double wv[40]; for(int i = 0; i < n; i++) wv[i] = randu();
            for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) H[i*n+j] += wv[i]*wv[j]; }
        for(int i = 0; i < n; i++) f[i] = scale*randu();
        // the duplicates (rows 0 and 1), and other inactive-ish soft/hard constraints
        double a[40]; for(int j = 0; j < n; j++) a[j] = randu();
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++) A[i*n+j] = i < nd ? a[j] : randu();
            bu[i] = i < nd ? 0 : scale*0.3*randu(); bl[i] = -DAQP_INF;
            sense[i] = i < nd ? DAQP_SOFT : 0;
            rl[i] = ru[i] = 1 + fabs(randu());
            wl[i] = wu[i] = 0;
        }
        // the multiplier of the hard version a'x <= 0, to put w in (lam/2, lam)
        DAQPProblem qr = {n, m-nd+1, 0, H, f, A+(nd-1)*n, bu+nd-1, bl+nd-1, NULL, NULL, 0, 0};
        static int sr[200]; memset(sr, 0, sizeof(sr)); qr.sense = sr;
        DAQPResult rr; rr.x = xr; rr.lam = lam;
        daqp_quadprog(&rr, &qr, NULL);
        if(rr.exitflag <= 0 || lam[0] < 1e-6){ k--; continue; } // (want an active a'x <= 0)
        const double lam_hard = lam[0], w = lam_hard*(0.55 + 0.4*fabs(randu()));
        for(int i = 0; i < nd; i++) wl[i] = wu[i] = w;
        // the soft version (cold, and warm from both duplicates active with fixed slacks)
        for(int warm = 0; warm < 2; warm++){
            for(int i = 0; i < m; i++) sense[i] = (i < nd ? DAQP_SOFT : 0) | (warm && i < nd ? DAQP_ACTIVE | DAQP_SLACK_FIXED : 0);
            DAQPProblem qp = {n, m, 0, H, f, A, bu, bl, sense, NULL, 0, 0};
            DAQPWorkspace work; memset(&work, 0, sizeof(work));
            double st;
            if(setup_daqp(&qp, &work, &st) < 0){ fails++; continue; }
            work.settings->iter_limit = 2000;
            daqp_set_soft_weights(&work, rl, ru, wl, wu);
            DAQPResult res; res.x = x; res.lam = lam;
            daqp_solve(&res, &work);
            double dx = 0;
            for(int i = 0; i < n; i++) if(fabs(x[i]-xr[i]) > dx) dx = fabs(x[i]-xr[i]);
            int ok = res.exitflag > 0 && dx < 1e-6*(1+scale);
            if(!ok){
                fails++; iterlim += res.exitflag == DAQP_EXIT_ITERLIMIT;
                if(fails <= 5) printf("test %d (%s, n=%d, scale %.1e, w/lam %.2f): exitflag %d, %d iterations, |dx| %.2e\n",
                        k, warm ? "warm" : "cold", n, scale, w/lam_hard, res.exitflag, res.iter, dx);
            }
            free_daqp_workspace(&work); free_daqp_ldp(&work);
        }
    }
    printf("%d of %d solves wrong (%d at the iteration limit)\n", fails, 2*ntests, iterlim);
    return fails > 0;
}
