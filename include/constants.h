#ifndef DAQP_CONSTANTS_H
#define DAQP_CONSTANTS_H

#define EMPTY_IND -1 
#define NX work->n 
#define N_CONSTR work->m 
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

#define EXIT_OPTIMAL 1
#define EXIT_INFEASIBLE -1
#define EXIT_CYCLE -2
#define EXIT_UNBOUNDED -3
#define EXIT_ITERLIMIT -4

#define INACTIVE_INEQUALITY 0
#define EQUALITY 1
#define ACTIVE_INEQUALITY 2
#define FREE_CONSTRAINT 3

#define INF ((c_float)1e30)

#endif //ifndef DAQP_CONSTANTS_H
