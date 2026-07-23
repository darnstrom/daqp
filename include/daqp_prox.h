#ifndef DAQP_PROX_H
# define DAQP_PROX_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"
#include "constants.h"
#include "daqp.h"

int daqp_prox(DAQPWorkspace *work);
c_float daqp_proximal_regularization(const DAQPWorkspace *work);
c_float daqp_proximal_regularization_scaled(
        const DAQPWorkspace *work, c_float hessian_scale);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_PROX_H
