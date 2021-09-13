#ifndef DAQP_DEBUG_H
# define DAQP_DEBUG_H
#include "daqp.h"
void printmatrix(c_float* matrix, int row ,int col);
void printvector(c_float* vec, int row);
void printWS(int* WS, int row);
void print_L(c_float* L, int n);

int daqp_debug(Workspace *work);
#endif //ifndef DAQP_DEBUG_H

