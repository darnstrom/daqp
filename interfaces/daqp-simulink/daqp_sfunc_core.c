#include "daqp_sfunc_core.h"
#include <stdlib.h>
#include <string.h>

DAQPSfunc* daqp_sfunc_new(int n, int m, int ms, const int* break_points, int nh,
        int warm_start, int iter_limit){
    int i;
    DAQPSfunc* s = calloc(1,sizeof(DAQPSfunc));
    if(s == NULL) return NULL;
    DAQPProblem* qp = &s->qp;
    qp->n = n;
    qp->m = m;
    qp->ms = ms;
    qp->H = calloc((size_t)n*n,sizeof(c_float));
    qp->f = calloc(n,sizeof(c_float));
    qp->A = calloc((size_t)(m-ms)*n+1,sizeof(c_float));
    qp->bupper = calloc(m+1,sizeof(c_float));
    qp->blower = calloc(m+1,sizeof(c_float));
    qp->sense = calloc(m+1,sizeof(int));
    s->sense_prev = calloc(m+1,sizeof(int));
    s->res.x = calloc(n,sizeof(c_float));
    s->res.lam = calloc(m+1,sizeof(c_float));
    if(nh > 1){
        qp->nh = nh;
        qp->break_points = malloc(nh*sizeof(int));
        if(qp->break_points != NULL)
            for(i = 0; i < nh; i++) qp->break_points[i] = break_points[i];
    }
    if(qp->H == NULL || qp->f == NULL || qp->A == NULL || qp->bupper == NULL ||
            qp->blower == NULL || qp->sense == NULL || s->sense_prev == NULL ||
            s->res.x == NULL || s->res.lam == NULL || (nh > 1 && qp->break_points == NULL)){
        daqp_sfunc_free(s);
        return NULL;
    }
    daqp_default_settings(&s->settings);
    s->settings.iter_limit = iter_limit;
    s->warm_start = warm_start;
    return s;
}

static void daqp_sfunc_teardown(DAQPSfunc* s){
    if(s->is_setup){
        s->work.settings = NULL; // Owned by s
        free_daqp_workspace(&s->work);
        free_daqp_ldp(&s->work);
    }
    memset(&s->work,0,sizeof(DAQPWorkspace));
    s->is_setup = 0;
    s->has_solved = 0;
}

int daqp_sfunc_solve(DAQPSfunc* s){
    int i, flag, mask, sense_changed = 0, structure_changed = 0;
    DAQPProblem* qp = &s->qp;
    int* sense = qp->sense;
    // A changed sense is used as is (rather than the previous working set), and soft
    // and binary constraints change what is allocated at setup => setup again
    for(i = 0; i < qp->m; i++){
        if(sense[i] != s->sense_prev[i]){
            sense_changed = 1;
            if((sense[i]^s->sense_prev[i]) & (DAQP_SOFT|DAQP_BINARY)) structure_changed = 1;
        }
        s->sense_prev[i] = sense[i];
    }
    if(structure_changed) daqp_sfunc_teardown(s);

    // Without warm start, every solve starts from scratch
    const int cold = s->warm_start == DAQP_SFUNC_COLD || !s->has_solved || sense_changed;
    const int scratch_mask = DAQP_UPDATE_unconstrained+DAQP_UPDATE_eliminate;

    if(!s->is_setup){
        memset(&s->work,0,sizeof(DAQPWorkspace));
        s->work.settings = &s->settings;
        flag = setup_daqp_main(qp,&s->work,NULL,
                s->warm_start == DAQP_SFUNC_COLD ? scratch_mask : 0);
        if(flag < 0){
            memset(&s->work,0,sizeof(DAQPWorkspace)); // (freed by setup_daqp_main)
            s->res.exitflag = flag;
            s->res.iter = 0;
            return flag;
        }
        s->is_setup = 1;
    }
    else{
        mask = DAQP_UPDATE_v+DAQP_UPDATE_d+DAQP_UPDATE_sense;
        if(s->warm_start != DAQP_SFUNC_WARM_FIXED) mask |= DAQP_UPDATE_Rinv+DAQP_UPDATE_M;
        if(s->warm_start == DAQP_SFUNC_COLD) mask |= scratch_mask;
        // Reuse the working set of the previous solve
        if(!cold) qp->sense = s->work.sense;
        flag = daqp_update_ldp(mask,&s->work,qp);
        qp->sense = sense;
        if(flag < 0){
            s->has_solved = 0;
            s->res.exitflag = flag;
            s->res.iter = 0;
            return flag;
        }
    }
    daqp_solve(&s->res,&s->work);
    s->has_solved = s->res.exitflag > 0;
    return s->res.exitflag;
}

void daqp_sfunc_free(DAQPSfunc* s){
    if(s == NULL) return;
    daqp_sfunc_teardown(s);
    free(s->qp.H);
    free(s->qp.f);
    free(s->qp.A);
    free(s->qp.bupper);
    free(s->qp.blower);
    free(s->qp.sense);
    free(s->qp.break_points);
    free(s->sense_prev);
    free(s->res.x);
    free(s->res.lam);
    free(s);
}
