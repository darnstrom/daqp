#ifndef DAQP_AUX_H 
# define DAQP_AUX_H

#include "types.h"
#include "constants.h"

void remove_constraint(Workspace* work, const int rm_ind);
void add_constraint(Workspace *work, const int add_ind, const int isupper);
//void find_constraint_to_add(Workspace *work);
void compute_primal_and_fval(Workspace *work);
int add_infeasible(Workspace *work);
int remove_blocking(Workspace *work);
void compute_CSP(Workspace *work);
void compute_singular_direction(Workspace *work);
void reorder_LDL(Workspace *work);
void pivot_last(Workspace *work);
int add_equality_constraints(Workspace *work);
#endif //ifndef DAQP_AUX_H
