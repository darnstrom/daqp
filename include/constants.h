#ifndef DAQP_CONSTANTS_H
#define DAQP_CONSTANTS_H

#include <stddef.h>
typedef double c_float;

#define EMPTY_IND -1 
#define NX work->n 
#define N_CONSTR work->m 
#define N_SIMPLE work->ms 
// #define c_float double 
#define DAQP_INF ((c_float)1e30)

// DEFAULT SETTINGS 
#define DEFAULT_PRIM_TOL 1e-6
#define DEFAULT_DUAL_TOL 1e-12 
#define DEFAULT_ZERO_TOL 1e-14
#define DEFAULT_PROG_TOL 1e-14 
#define DEFAULT_PIVOT_TOL 1e-4
#define DEFAULT_CYCLE_TOL 10
#define DEFAULT_ETA 1e-6
#define DEFAULT_ITER_LIMIT 1000 
#define DEFAULT_RHO_SOFT 1e-3 

// MACROS
#define SQUARE(x) ((x)*(x))
#define ARSUM(x) ((x)*(x+1)/2)
#define R_OFFSET(X,Y) (((2*Y-X-1)*X)/2)

// EXIT FLAGS
#define EXIT_SOFT_OPTIMAL 2 
#define EXIT_OPTIMAL 1
#define EXIT_INFEASIBLE -1
#define EXIT_CYCLE -2
#define EXIT_UNBOUNDED -3
#define EXIT_ITERLIMIT -4
#define EXIT_NONCONVEX -5
#define EXIT_OVERDETERMINED_INITIAL -6

// UPDATE LDP MASKS 
#define UPDATE_Rinv 1
#define UPDATE_M 2
#define UPDATE_v 4
#define UPDATE_d 8
#define UPDATE_sense 16 

// CONSTRAINT MASKS 
#define ACTIVE 1
#define IS_ACTIVE(x) (work->sense[x]&1)
#define SET_ACTIVE(x) (work->sense[x]|=1)
#define SET_INACTIVE(x) (work->sense[x]&=~1)

#define LOWER 2 
#define IS_LOWER(x) (work->sense[x]&2)
#define SET_LOWER(x) (work->sense[x]|=2)
#define SET_UPPER(x) (work->sense[x]&=~2)

#define IMMUTABLE 4 
#define IS_IMMUTABLE(x) (work->sense[x]&4)
#define SET_IMMUTABLE(x) (work->sense[x]|=4)
#define SET_MUTABLE(x) (work->sense[x]&=~4)

#define SOFT 8 
#define IS_SOFT(x) (work->sense[x]&8)
#define SET_SOFT(x) (work->sense[x]|=8)
#define SET_HARD(x) (work->sense[x]&=~8)

#define BINARY 16 
#define IS_BINARY(x) (work->sense[x]&16)

#define IS_SIMPLE(x) (x < work->ms)


#endif //ifndef DAQP_CONSTANTS_H
