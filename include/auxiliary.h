#ifndef DAQP_AUX_H 
# define DAQP_AUX_H

#include "types.h"
#include "constants.h"

void remove_constraint(Workspace* work);
void add_constraint(Workspace *work);
void find_constraint_to_add(Workspace *work);
void find_blocking_constraints(Workspace *work);
void compute_alpha_and_rm_blocking(Workspace *work);
void compute_CSP(Workspace *work);
void compute_singular_direction(Workspace *work);
void reorder_LDL(Workspace *work);
void pivot_last(Workspace *work);
void add_equality_constraints(Workspace *work);
void daqp_default_settings(DAQPSettings *settings);
#endif //ifndef DAQP_AUX_H
