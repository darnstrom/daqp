#include "hierarchical.h"
#include "types.h"
#include <stdlib.h>

int daqp_hiqp(DAQPWorkspace *work){
    int i,j,id;
    int start,end;
    int iterations = 0;
    int exitflag=0;

    // If only one hiearchy -> just solve normal LDP
    if( work->nh < 2) return daqp_ldp(work);

    // Start moving down the hierarchy
    c_float w;
    start=work->break_points[0];
    int nfree = work->n;
    for(i =1; i < work->nh; i++){
        // initialize current level
        end=work->break_points[i];
        work->m = end;
        // Soften constraints and activate 
        for(j =start;j<end;j++){
            SET_SOFT(j);
            if(IS_ACTIVE(j)){
                if(IS_LOWER(j))
                    add_constraint(work,j, -1.0);
                else
                    add_constraint(work,j, 1.0);
                nfree--;
                if(work->sing_ind != EMPTY_IND) return EXIT_OVERDETERMINED_INITIAL;
                }
        }

        // Solve best solution in case daqp_ldp fails
        for(j = 0; j<work->n;j++) work->xold[j] = work->x[j];
        // Solve LDP
        exitflag = daqp_ldp(work);
        iterations+=work->iterations;
        if(exitflag < 0) break;

        if(iterations >= work->settings->iter_limit){
            exitflag = EXIT_ITERLIMIT;
            break;
        }

        // Perturb rhs with slacks in level 
        for(j=0; j<work->n_active;j++){
            id=work->WS[j];
            if(IS_SOFT(id)){ 
                w = work->lam_star[j]*work->settings->rho_soft;
                if(IS_LOWER(id))
                    work->dlower[id]+=w;
                else
                    work->dupper[id]+=w;
                if(IS_IMMUTABLE(id)) continue; // Already an equality
                nfree--;
            }
        }
        if(nfree <= 0 ) break;  // No degrees of freedom left 

        // Make constraints in current level hard
        for(j=start; j<end;j++) SET_HARD(j);
        
        // find first active constraint in current level 
        for(j=0;j < work->n_active; j++) if(work->WS[j]>=start) break;

        // reactive constraint in level current (to addresss soft->hard)
        // TODO: can factorization be directly reused?
        int n_active_old = (work->n_active < work->n) ? work->n_active : work->n;
        for(int jj=n_active_old; jj < work->n_active ;jj++) SET_INACTIVE(work->WS[jj]);
        work->n_active =j;
        work->reuse_ind=j;
        work->sing_ind = EMPTY_IND;
        for(; j<n_active_old ;j++){
            add_constraint(work,work->WS[j],work->lam_star[j]);
            // Abort if WS becomes overdtermined
            if(work->sing_ind != EMPTY_IND) return EXIT_OVERDETERMINED_INITIAL;
        }

        // Move up hierarchy
        start = end;
    }
    // Finalize
    if(exitflag < 0){ // Restore a point that was good before it failed
        for(j = 0; j<work->n;j++) work->x[j] = work->xold[j];
        exitflag = 3; // signify no degrees of freedoom left
    }
    work->iterations = iterations; // Append total number of iterations
    return exitflag;
}
