#ifndef DAQP_AVI_H 
# define DAQP_AVI_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"

typedef struct{
    DAQPWorkspace* work;
    DAQPProblem problem;

    c_float* Hsym;
    c_float* H1pI;
    c_float* H2pI;
    int* P_H2;

    c_float* LU_H;
    int* P_H;

    c_float* kkt_buffer;
    int* P_S;

    c_float* xtemp;
    c_float* Hx;
    c_float* x;
    c_float* y;
}DAQPAVI;

int _daqp_avi(DAQPAVI* avi);

int daqp_lu(c_float* A, int* P, int n);
void daqp_lu_solve(c_float* LU, int* P, c_float* b, c_float* x, int n);

void daqp_solve_avi_kkt(DAQPAVI* avi);
int daqp_check_optimal_avi(DAQPAVI* avi);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_AVI_H
