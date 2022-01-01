#ifndef DAQP_CONSTANTS_H
#define DAQP_CONSTANTS_H

#include <stddef.h>

#define EMPTY_IND -1 
#define NX work->n 
#define N_CONSTR work->m 
#define N_SIMPLE work->ms 
#define c_float double 

#define DEFAULT_PRIM_TOL 1e-6
#define DEFAULT_DUAL_TOL 1e-12 
#define DEFAULT_ZERO_TOL 1e-14
#define DEFAULT_PROG_TOL 1e-6
#define DEFAULT_FARKAS_TOL 1e-7
#define DEFAULT_PIVOT_TOL 1e-2
#define DEFAULT_CYCLE_TOL 10
#define DEFAULT_ETA 1e-6

#define DEFAULT_ITER_LIMIT 1000 
#define DEFAULT_PROX_ITER_LIMIT 1000 

#define SQUARE(x) ((x)*(x))
#define ARSUM(x) ((x)*(x+1)/2)
#define INF ((c_float)1e30)

// EXIT FLAGS
#define EXIT_OPTIMAL 1
#define EXIT_INFEASIBLE -1
#define EXIT_CYCLE -2
#define EXIT_UNBOUNDED -3
#define EXIT_ITERLIMIT -4

// CONSTRAINT MASKS 
#define INACTIVE 1
#define ACTIVE 1
#define IS_ACTIVE(x) (work->sense[x]&1)
#define SET_ACTIVE(x) (work->sense[x]|=1)
#define SET_INACTIVE(x) (work->sense[x]&=~1)

#define UPPER 0
#define LOWER 2 
#define IS_LOWER(x) (work->sense[x]&2)
#define SET_LOWER(x) (work->sense[x]|=2)
#define SET_UPPER(x) (work->sense[x]&=~2)

#define MUTABLE 0
#define IMMUTABLE 4 
#define IS_IMMUTABLE(x) (work->sense[x]&4)
#define SET_IMMUTABLE(x) (work->sense[x]|=4)
#define SET_MUTABLE(x) (work->sense[x]&=~4)

#define IS_SIMPLE(x) (x < work->ms)


#endif //ifndef DAQP_CONSTANTS_H
