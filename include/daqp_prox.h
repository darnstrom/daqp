#ifndef DAQP_PROX_H 
# define DAQP_PROX_H 

#include "types.h"
#include "constants.h"
#include "daqp.h"

int daqp_prox(Workspace *work);
int gradient_step(Workspace* work);

#endif //ifndef DAQP_PROX_H 
