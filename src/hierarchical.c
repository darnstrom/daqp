#include "hierarchical.h"
#include "types.h"


int daqp_hiqp(DAQPWorkspace *work){
    int i,j,id;
    int start,end;
    int iterations = 0;
    int exitflag=0;
    start=0;
    work->settings->rho_soft = 1e-10;
    for(i =0; i < work->hier->nh; i++){
        // initialize current level
        end=work->hier->break_points[i];
        work->m = end;

        // Solve LDP
        exitflag = daqp_ldp(work);
        if(exitflag < 0) return exitflag;

        iterations+=work->iterations;
        // Perturb rhs with slacks in level 
        for(j=0; j<work->n_active;j++){
            id=work->WS[j];
            if(IS_SOFT(id)){ 
                if(IS_LOWER(id))
                    work->dlower[id]+=work->lam_star[j]/work->scaling[id];
                else
                    work->dupper[id]+=work->lam_star[j]/work->scaling[id];
            }
        }

        // Make constraints in current level hard
        for(j=start; j<end;j++) SET_HARD(j);
        
        // find first active constraint in current level 
        for(j=0;j < work->n_active; j++) 
            if(work->WS[j]>=start) break;
        // reactive constraint in level current (to addresss soft->hard) 
        // TODO: can factorization be directly reused? 
        int n_active_old = work->n_active;
        work->n_active =j;
        work->reuse_ind=j;
        work->sing_ind = EMPTY_IND;
        for(; j<n_active_old ;j++){
            add_constraint(work,work->WS[j],work->lam_star[j]);
            // Abort reactivation if WS becomes overdtermined
            if(work->sing_ind != EMPTY_IND){
                for(; j<n_active_old ;j++) SET_INACTIVE(work->WS[j]);
                break;
            }
        }
        // Move up hierarchy
        start = end;
    }
    work->iterations = iterations; // Append total number of iterations
    return exitflag;
}
