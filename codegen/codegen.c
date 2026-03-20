#include "codegen.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "types.h"
#include "api.h"


void render_daqp_workspace(DAQPWorkspace* work, const char *fname, const char *dir, const char *prefix){
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
    fprintf(fh, "extern DAQPSettings %ssettings;\n\n", prefix);
    write_daqp_settings_src(fsrc,work->settings,prefix);

    // Write BnB struct
    if(work->bnb != NULL){
        write_daqp_bnb_h(fh,work->bnb,work->n,prefix);
        write_daqp_bnb_src(fsrc,work->bnb,work->n,prefix);
    }

    // Write Hierarchical
    if(work->nh > 1 ){
        fprintf(fh, "#define DAQP_HIERARCHICAL\n");
        fprintf(fh, "extern int %sbreak_points[%d];\n", prefix, work->nh);
        char varname[256];
        snprintf(varname, sizeof(varname), "%sbreak_points", prefix);
        write_int_array(fsrc,work->break_points, work->nh,varname);
    }

    // TODO Check soft constraints 

    //Write workspace
    write_daqp_workspace_h(fh,work,prefix);
    write_daqp_workspace_src(fsrc,work,prefix);


    // Close header guard 
    fprintf(fh, "#endif // ifndef %s\n", hguard);

    fclose(fh);
    fclose(fsrc);

    free(hfname);
    free(cfname);
    free(hguard);
}

void write_daqp_workspace_h(FILE *f, DAQPWorkspace* work, const char* prefix){
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

    fprintf(f, "extern c_float %sM[%d];\n", prefix, (m-ms)*n);
    fprintf(f, "extern c_float %sdupper[%d];\n", prefix, m);
    fprintf(f, "extern c_float %sdlower[%d];\n", prefix, m);
    fprintf(f, "extern c_float %sRinv[%d];\n", prefix, n*(n+1)/2);
    //fprintf(f, "extern c_float %sv[%d];\n", prefix, n);
    fprintf(f, "extern int %ssense[%d];\n\n", prefix, m);
    //fprintf(f, "extern c_float %sscaling[%d];\n\n", prefix, m);

    fprintf(f, "extern c_float %sx[%d];\n", prefix, n+1);
    fprintf(f, "extern c_float %sxold[%d];\n\n", prefix, n+1);

    fprintf(f, "extern c_float %slam[%d];\n", prefix, ntot+1);
    fprintf(f, "extern c_float %slam_star[%d];\n", prefix, ntot+1);
    fprintf(f, "extern c_float %su[%d];\n\n", prefix, n+1);

    fprintf(f, "extern c_float %sL[%d];\n", prefix, (ntot+1)*(ntot+2)/2);
    fprintf(f, "extern c_float %sD[%d];\n", prefix, ntot+1);
    fprintf(f, "extern c_float %sxldl[%d];\n", prefix, ntot+1);
    fprintf(f, "extern c_float %szldl[%d];\n\n", prefix, ntot+1);

    fprintf(f, "extern int %sWS[%d];\n\n", prefix, ntot+1);

    fprintf(f, "extern DAQPWorkspace %swork;\n\n", prefix);
}

void write_daqp_workspace_src(FILE* f, DAQPWorkspace* work, const char* prefix){
    int i;
    int n = work->n;
    int m = work->m;
    int ms = work->ms;
    int ntot = n;
    for(i = 0; i < m ; i++) 
        if(work->sense[i] & DAQP_SOFT) ntot++;

    char varname[256];
    fprintf(f, "// Workspace\n");
    // LDP data
    snprintf(varname, sizeof(varname), "%sM", prefix);
    write_float_array(f,work->M,(m-ms)*n,varname);
    //snprintf(varname, sizeof(varname), "%sdupper", prefix);
    //write_float_array(f,work->dupper,m,varname);
    //snprintf(varname, sizeof(varname), "%sdlower", prefix);
    //write_float_array(f,work->dlower,m,varname);
    fprintf(f, "c_float %sdupper[%d];\n", prefix, m);
    fprintf(f, "c_float %sdlower[%d];\n", prefix, m);
    snprintf(varname, sizeof(varname), "%sRinv", prefix);
    write_float_array(f,work->Rinv,n*(n+1)/2,varname);
    //snprintf(varname, sizeof(varname), "%sv", prefix);
    //write_float_array(f,work->v,n,varname);
    snprintf(varname, sizeof(varname), "%ssense", prefix);
    write_int_array(f,work->sense, m,varname);
    //snprintf(varname, sizeof(varname), "%sscaling", prefix);
    //write_float_array(f,work->scaling, m,varname);

    // Iteratates
    fprintf(f, "c_float %sx[%d];\n", prefix, n+1);
    fprintf(f, "c_float %sxold[%d];\n\n", prefix, n+1);

    fprintf(f, "c_float %slam[%d];\n", prefix, ntot+1);
    fprintf(f, "c_float %slam_star[%d];\n", prefix, ntot+1);
    fprintf(f, "c_float %su[%d];\n\n", prefix, n+1);

    fprintf(f, "c_float %sL[%d];\n", prefix, (ntot+1)*(ntot+2)/2);
    fprintf(f, "c_float %sD[%d];\n", prefix, ntot+1);
    fprintf(f, "c_float %sxldl[%d];\n", prefix, ntot+1);
    fprintf(f, "c_float %szldl[%d];\n\n", prefix, ntot+1);

    fprintf(f, "int %sWS[%d];\n\n", prefix, ntot+1);

    //Workspace struct
    fprintf(f, "DAQPWorkspace %swork= {\n", prefix);
    fprintf(f, "NULL,\n"); // DAQPProblem
    fprintf(f, "%d, %d, %d,\n",n,m,ms); // dimensions 
    fprintf(f, "%sM, %sdupper, %sdlower, %sRinv, NULL, %ssense,\n",
            prefix,prefix,prefix,prefix,prefix); //LDP
    fprintf(f, "NULL,\n"); // scaling
    fprintf(f, "NULL,\n"); // RinvD
    fprintf(f, "%sx, %sxold,\n", prefix, prefix);
    fprintf(f, "%slam, %slam_star, %su, %d,\n", prefix, prefix, prefix, -1); // fval
    fprintf(f, "%sL, %sD, %sxldl,%szldl,%d,\n",
            prefix,prefix,prefix,prefix, 0); // reuse_ind
    fprintf(f, "%sWS, %d,\n", prefix, 0); //n_active
    fprintf(f, "%d,%d,\n",0,-1); //iterations + sing_id
    fprintf(f, "%f,\n",0.0); // Soft slack
    fprintf(f, "&%ssettings, \n", prefix);
    // BnB
    if(work->bnb == NULL)
        fprintf(f, "NULL,\n");
    else
        fprintf(f, "&%sbnb_work,\n", prefix);
    // Hierarhical
    if(work->nh > 1)
        fprintf(f, "%d,%sbreak_points,\n", work->nh, prefix);
    else
        fprintf(f, "0, NULL,\n");
    // AVI
    fprintf(f, "NULL,\n"); // TODO: Generate for avi (also requires problem to be generated)
    // Timer
    fprintf(f, "NULL};\n\n");
}

void write_daqp_settings_src(FILE*  f, DAQPSettings* settings, const char* prefix){

    fprintf(f, "// Settings\n");
    fprintf(f, "DAQPSettings %ssettings = {", prefix);
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

void write_daqp_bnb_h(FILE*  f, DAQPBnB* bnb, const int n, const char* prefix){
    fprintf(f, "#define DAQP_BNB\n");
    fprintf(f, "extern int %sbin_ids[%d];\n", prefix, bnb->nb);
    fprintf(f, "extern DAQPNode %stree[%d];\n", prefix, bnb->nb+1);
    fprintf(f, "extern int %stree_WS[%d];\n", prefix, (n+1)*(bnb->nb+1));
    fprintf(f, "extern int %sfixed_ids[%d];\n", prefix, bnb->nb+1);
    fprintf(f, "extern DAQPBnB %sbnb_work;\n\n", prefix);
}
void write_daqp_bnb_src(FILE*  f, DAQPBnB* bnb, const int n, const char* prefix){
    if(bnb==NULL) return;
    char varname[256];
    fprintf(f, "// BnB \n");

    snprintf(varname, sizeof(varname), "%sbin_ids", prefix);
    write_int_array(f,bnb->bin_ids, bnb->nb,varname);
    fprintf(f, "DAQPNode %stree[%d];\n", prefix, bnb->nb+1);
    fprintf(f, "int %stree_WS[%d];\n", prefix, (n+1)*(bnb->nb+1));
    fprintf(f, "int %sfixed_ids[%d];\n", prefix, bnb->nb+1);

    fprintf(f, "DAQPBnB %sbnb_work= {", prefix);
    fprintf(f, "%sbin_ids, ", prefix);
    fprintf(f, "(int)%d, ", bnb->nb);
    fprintf(f, "(int)%d, ", bnb->neq);

    fprintf(f, "%stree, ", prefix);
    fprintf(f, "(int)%d, ", 0); // n_nodes

    fprintf(f, "%stree_WS, ", prefix);
    fprintf(f, "(int)%d, ", 0); // nWS
    fprintf(f, "(int)%d, ", 0); // n_clean
    fprintf(f, "%sfixed_ids, ", prefix); // fixed_ids

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
