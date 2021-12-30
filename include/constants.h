#ifndef DAQP_CONSTANTS_H
#define DAQP_CONSTANTS_H

#define EMPTY_IND -1 
#define NX work->n 
#define N_CONSTR work->m 
#define c_float double 

#define INFEAS_TOL 1e-6
#define DUAL_TOL 1e-12 
#define ZERO_TOL 1e-14
#define OBJ_PROG_TOL 1e-6
#define FARKAS_TOL 1e-7
#define PIVOT_TOL 1e-2
#define CYCLE_TOL 10
#define ETA 1e-6

#define MAX_ITER 1000 
#define PROX_ITER_LIMIT 1000 
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
