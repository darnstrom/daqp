#include "hierarchical.h"
#include "types.h"

#define DAQP_HIQP_SOFT ((c_float)1e-5)

int daqp_hiqp(DAQPWorkspace *work){
    int i,j,id;
    int start,end;
    int iterations = 0;
    int exitflag=0;
    start=0;
    work->settings->rho_soft = DAQP_HIQP_SOFT*DAQP_HIQP_SOFT;
    for(i =0; i < work->hier->nh; i++){
        // initialize current level
        end=work->hier->break_points[i];
        work->m = end;
        //printf("H-level: %d, start:%d, end:%d\n",i,start,end);
        //printf("sense: ");
        //for(int ii = 0; ii < work->m; ii++)
        //    printf(" %d",work->sense[ii]);
        //printf("\n");

        // Solve LDP
        exitflag = daqp_ldp(work);
        if(exitflag < 0) return exitflag;
        //printf("WS: {");
        //for(int ii = 0; ii < work->n_active; ii++)
        //    printf(" %d",work->WS[ii]);
        //printf(" }\n");

        iterations+=work->iterations;
        // Perturb rhs with slacks in level 
        for(j=0; j<work->n_active;j++){
            id=work->WS[j];
            if(IS_SOFT(id)){ 
                if(IS_LOWER(id))
                    work->dlower[id]+=work->lam_star[j]*DAQP_HIQP_SOFT;
                else
                    work->dupper[id]+=work->lam_star[j]*DAQP_HIQP_SOFT;
            }
        }

        // Make constraints in current level hard
        for(j=start; j<end;j++) SET_HARD(j);
        
        // find first active constraint in current level 
        deactivate_constraints(work);
        reset_daqp_workspace(work);
        //for(j=0;j < work->n_active; j++)
        //    if(work->WS[j]>=start) break;
        //printf("first active in H-level: %d\n",j);
        //// reactive constraint in level current (to addresss soft->hard)
        //// TODO: can factorization be directly reused?
        //int n_active_old = (work->n_active < work->n) ? work->n_active : work->n;
        //work->n_active =j;
        //work->reuse_ind=j;
        //work->sing_ind = EMPTY_IND;
        //for(; j<n_active_old ;j++){
        //    add_constraint(work,work->WS[j],work->lam_star[j]);
        //    // Abort reactivation if WS becomes overdtermined
        //    if(work->sing_ind != EMPTY_IND){
        //        remove_constraint(work,j);
        //        printf("Singular!\n");
        //        work->sing_ind = EMPTY_IND;
        //        for(; j<n_active_old ;j++) SET_INACTIVE(work->WS[j]);
        //        break;
        //    }
        //}
        //printf("WS: {");
        //for(int ii = 0; ii < work->n_active; ii++)
        //    printf(" %d",work->WS[ii]);
        //printf(" }\n");
        // Move up hierarchy
        start = end;
    }
    work->iterations = iterations; // Append total number of iterations
    return exitflag;
}
