#ifndef DAQP_FACTORIZATION_H
# define DAQP_FACTORIZATION_H

#include "types.h"
#include "constants.h"

void update_LDL_add(DAQPWorkspace *work, const int add_ind);
void update_LDL_remove(DAQPWorkspace *work, const int rm_ind);

#endif //ifndef DAQP_FACTORIZATION_H
