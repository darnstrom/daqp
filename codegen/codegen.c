#include "codegen.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "types.h"
#include "api.h"


void render_daqp_workspace(DAQPWorkspace* work, const char *fname, const char *dir){
    char *hfname= malloc(strlen(dir)+strlen(fname) + 3); // two chars for .h and 1 for terminator
    char *hguard= malloc(strlen(fname) + 3); // two chars for _H and 1 for terminator
    char *cfname= malloc(strlen(dir)+strlen(fname) + 3); // two chars for .c and 1 for terminator

    // Header file
    strcpy(hfname, dir);
    strcat(hfname, fname);
    strcat(hfname, ".h");
    FILE* fh= fopen(hfname, "w");

    // Source file
    strcpy(cfname, dir);
    strcat(cfname, fname);
    strcat(cfname, ".c");
    FILE* fsrc= fopen(cfname, "w");

    // Header guard
    strcpy(hguard, fname);
    strcat(hguard, "_H");
    char *s = hguard;
    while(*s){
        *s = toupper((unsigned char) *s);
        s++;
    }
    fprintf(fh, "#ifndef %s\n",   hguard);
    fprintf(fh, "#define %s\n\n", hguard);

    // Include types and constants
    fprintf(fh, "#include \"types.h\"\n");
    fprintf(fh, "#include \"constants.h\"\n");

    fprintf(fsrc, "#include \"types.h\"\n");
    fprintf(fsrc, "#include \"constants.h\"\n");

    //Write settings
    fprintf(fh, "// Settings prototype\n");
    fprintf(fh, "extern DAQPSettings daqp_settings;\n\n");
    write_daqp_settings_src(fsrc,work->settings);

    // Write BnB struct
    if(work->bnb != NULL){
        write_daqp_bnb_h(fh,work->bnb,work->n);
        write_daqp_bnb_src(fsrc,work->bnb,work->n);
    }

    // Write Hierarchical
    if(work->nh > 1 ){
        fprintf(fh, "#define DAQP_HIERARCHICAL\n");
        fprintf(fh, "extern int daqp_break_points[%d];\n", work->nh);
        write_int_array(fsrc,work->break_points, work->nh,"daqp_break_points");
    }

    // TODO Check soft constraints 

    //Write workspace
    write_daqp_workspace_h(fh,work);
    write_daqp_workspace_src(fsrc,work);


    // Close header guard 
    fprintf(fh, "#endif // ifndef %s\n", hguard);

    fclose(fh);
    fclose(fsrc);

    free(hfname);
    free(cfname);
    free(hguard);
}

void write_daqp_workspace_h(FILE *f, DAQPWorkspace* work){
    int i;
    const int n = work->n;
    const int m = work->m;
    const int ms = work->ms;
    int ntot = n;
    // Account for soft constraints
    if(work->nh > 1){
        int ns = 0, start=0;
        for(i = 0; i < work->nh; i++){
            ns = (ns  > work->break_points[i]-start) ? ns : work->break_points[i]-start;
            start = work->break_points[i];
        }
        ntot+=ns;
    }
    else{// Simply count soft constraints if not hierarchical
        for(i = 0; i < m ; i++)
            if(work->sense[i] & DAQP_SOFT) ntot++;
    }

    // Refdefine NX, N_CONSTR and N_SIMPLE to static
    fprintf(f, "#undef NX\n");
    fprintf(f, "#define NX %d\n",n);
    fprintf(f, "#undef N_CONSTR\n");
    fprintf(f, "#define N_CONSTR %d\n",m);
    fprintf(f, "#undef N_SIMPLE\n");
    fprintf(f, "#define N_SIMPLE %d \n",ms);

    fprintf(f, "// Workspace prototypes\n");

    fprintf(f, "extern c_float daqp_M[%d];\n", (m-ms)*n);
    fprintf(f, "extern c_float daqp_dupper[%d];\n", m);
    fprintf(f, "extern c_float daqp_dlower[%d];\n", m);
    fprintf(f, "extern c_float daqp_Rinv[%d];\n", n*(n+1)/2);
    //fprintf(f, "extern c_float daqp_v[%d];\n", n);
    fprintf(f, "extern int daqp_sense[%d];\n\n", m);
    //fprintf(f, "extern c_float daqp_scaling[%d];\n\n", m);

    fprintf(f, "extern c_float daqp_x[%d];\n", n+1);
    fprintf(f, "extern c_float daqp_xold[%d];\n\n", n+1);

    fprintf(f, "extern c_float daqp_lam[%d];\n", ntot+1);
    fprintf(f, "extern c_float daqp_lam_star[%d];\n", ntot+1);
    fprintf(f, "extern c_float daqp_u[%d];\n\n", n+1);

    fprintf(f, "extern c_float daqp_L[%d];\n", (ntot+1)*(ntot+2)/2);
    fprintf(f, "extern c_float daqp_D[%d];\n", ntot+1);
    fprintf(f, "extern c_float daqp_xldl[%d];\n", ntot+1);
    fprintf(f, "extern c_float daqp_zldl[%d];\n\n", ntot+1);

    fprintf(f, "extern int daqp_WS[%d];\n\n", ntot+1);

    fprintf(f, "extern DAQPWorkspace daqp_work;\n\n");
}

void write_daqp_workspace_src(FILE* f, DAQPWorkspace* work){
    int i;
    int n = work->n;
    int m = work->m;
    int ms = work->ms;
    int ntot = n;
    for(i = 0; i < m ; i++) 
        if(work->sense[i] & DAQP_SOFT) ntot++;

    fprintf(f, "// Workspace\n");
    // LDP data
    write_float_array(f,work->M,(m-ms)*n,"daqp_M");
    //write_float_array(f,work->dupper,m,"daqp_dupper");
    //write_float_array(f,work->dlower,m,"daqp_dlower");
    fprintf(f, "c_float daqp_dupper[%d];\n", m);
    fprintf(f, "c_float daqp_dlower[%d];\n", m);
    write_float_array(f,work->Rinv,n*(n+1)/2,"daqp_Rinv");
    //write_float_array(f,work->v,n, "daqp_v");
    write_int_array(f,work->sense, m,"daqp_sense");
    //write_float_array(f,work->scaling, m,"daqp_scaling");

    // Iteratates
    fprintf(f, "c_float daqp_x[%d];\n", n+1);
    fprintf(f, "c_float daqp_xold[%d];\n\n", n+1);

    fprintf(f, "c_float daqp_lam[%d];\n", ntot+1);
    fprintf(f, "c_float daqp_lam_star[%d];\n", ntot+1);
    fprintf(f, "c_float daqp_u[%d];\n\n", n+1);

    fprintf(f, "c_float daqp_L[%d];\n", (ntot+1)*(ntot+2)/2);
    fprintf(f, "c_float daqp_D[%d];\n", ntot+1);
    fprintf(f, "c_float daqp_xldl[%d];\n", ntot+1);
    fprintf(f, "c_float daqp_zldl[%d];\n\n", ntot+1);

    fprintf(f, "int daqp_WS[%d];\n\n", ntot+1);

    //Workspace struct
    fprintf(f, "DAQPWorkspace daqp_work= {\n");
    fprintf(f, "NULL,\n"); // DAQPProblem
    fprintf(f, "%d, %d, %d,\n",n,m,ms); // dimensions 
    fprintf(f, "daqp_M, daqp_dupper, daqp_dlower, daqp_Rinv, NULL, daqp_sense,\n"); //LDP
    fprintf(f, "NULL,\n"); // scaling
    fprintf(f, "NULL,\n"); // RinvD
    fprintf(f, "daqp_x, daqp_xold,\n");
    fprintf(f, "daqp_lam, daqp_lam_star, daqp_u, %d,\n",-1); // fval
    fprintf(f, "daqp_L, daqp_D, daqp_xldl,daqp_zldl,%d,\n",0); // reuse_ind
    fprintf(f, "daqp_WS, %d,\n",0); //n_active
    fprintf(f, "%d,%d,\n",0,-1); //iterations + sing_id
    fprintf(f, "%f,\n",0.0); // Soft slack
    fprintf(f, "&daqp_settings, \n");
    // BnB
    if(work->bnb == NULL)
        fprintf(f, "NULL,\n");
    else
        fprintf(f, "&daqp_bnb_work,\n");
    // Hierarhical
    if(work->nh > 1)
        fprintf(f, "%d,daqp_break_points,\n",work->nh);
    else
        fprintf(f, "0, NULL,\n");
    // AVI
    fprintf(f, "NULL,\n"); // TODO: Generate for avi (also requires problem to be generated)
    // Timer
    fprintf(f, "NULL};\n\n");
}

void write_daqp_settings_src(FILE*  f, DAQPSettings* settings){

    fprintf(f, "// Settings\n");
    fprintf(f, "DAQPSettings daqp_settings = {");
    fprintf(f, "(c_float)%.20f, ", settings->primal_tol);
    fprintf(f, "(c_float)%.20f, ", settings->dual_tol);
    fprintf(f, "(c_float)%.20f, ", settings->zero_tol);
    fprintf(f, "(c_float)%.20f, ", settings->pivot_tol);
    fprintf(f, "(c_float)%.20f, ", settings->progress_tol);

    fprintf(f, "%d, ",             settings->cycle_tol);
    fprintf(f, "%d, ",             settings->iter_limit);
    fprintf(f, "(c_float)%.20f, ", settings->fval_bound);

    fprintf(f, "(c_float)%.20f, ", settings->eps_prox);
    fprintf(f, "(c_float)%.20f, ", settings->eta_prox);

    fprintf(f, "(c_float)%.20f,", settings->rho_soft);

    fprintf(f, "(c_float)%.20f,", settings->rel_subopt);
    fprintf(f, "(c_float)%.20f,", settings->abs_subopt);

    fprintf(f, "(c_float)%.20f,",  settings->sing_tol);
    fprintf(f, "(c_float)%.20f,",  settings->refactor_tol);
    fprintf(f, "(c_float)%.20f",  settings->time_limit);
    fprintf(f, "};\n\n");
}

void write_daqp_bnb_h(FILE*  f, DAQPBnB* bnb, const int n){
    fprintf(f, "#define DAQP_BNB\n");
    fprintf(f, "extern int daqp_bin_ids[%d];\n", bnb->nb);
    fprintf(f, "extern DAQPNode daqp_tree[%d];\n", bnb->nb+1);
    fprintf(f, "extern int daqp_tree_WS[%d];\n", (n+1)*(bnb->nb+1));
    fprintf(f, "extern int daqp_fixed_ids[%d];\n", bnb->nb+1);
    fprintf(f, "extern DAQPBnB daqp_bnb_work;\n\n");
}
void write_daqp_bnb_src(FILE*  f, DAQPBnB* bnb, const int n){
    if(bnb==NULL) return;
    fprintf(f, "// BnB \n");

    write_int_array(f,bnb->bin_ids, bnb->nb,"daqp_bin_ids");
    fprintf(f, "DAQPNode daqp_tree[%d];\n", bnb->nb+1);
    fprintf(f, "int daqp_tree_WS[%d];\n", (n+1)*(bnb->nb+1));
    fprintf(f, "int daqp_fixed_ids[%d];\n", bnb->nb+1);

    fprintf(f, "DAQPBnB daqp_bnb_work= {");
    fprintf(f, "daqp_bin_ids, ");
    fprintf(f, "(int)%d, ", bnb->nb);
    fprintf(f, "(int)%d, ", bnb->neq);

    fprintf(f, "daqp_tree, ");
    fprintf(f, "(int)%d, ", 0); // n_nodes

    fprintf(f, "daqp_tree_WS, ");
    fprintf(f, "(int)%d, ", 0); // nWS
    fprintf(f, "(int)%d, ", 0); // n_clean
    fprintf(f, "daqp_fixed_ids, "); // n_clean

    fprintf(f, "(int)%d, ", 0); // nodecount
    fprintf(f, "(int)%d, ", 0); // itercount

    fprintf(f, "};\n\n");
}

void write_float_array(FILE *f, c_float* a, const int N, const char *name){
    if(a == NULL)
        fprintf(f, "c_float* const %s = NULL;\n",name);
        //fprintf(f, "c_float %s[%d];\n", name, N);
    else{
        int i;
        fprintf(f, "c_float %s[%d] = {\n", name, N);
        for(i = 0; i < N; i++)
            fprintf(f, "(c_float)%.20f,\n", a[i]);
        fprintf(f, "};\n");
    }
}

void write_int_array(FILE *f, int* a, const int N, const char *name){
    if(a == NULL)
        fprintf(f, "int* const %s = NULL;\n",name);
        //fprintf(f, "int %s[%d];\n", name, N);
    else{
        int i;
        fprintf(f, "int %s[%d] = {\n", name, N);
        for(i = 0; i < N; i++)
            fprintf(f, "(int)%i,\n", a[i]);
        fprintf(f, "};\n");
    }
}
