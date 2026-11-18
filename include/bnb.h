#ifndef DAQP_BNB_H
# define DAQP_BNB_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"
#include "constants.h"
#include "daqp.h"

int daqp_bnb(DAQPWorkspace* work);
// The offset c between the objective J and the scale 0.5*fval of daqp_ldp (0.5*fval = J + c)
c_float daqp_bnb_offset(const DAQPWorkspace* work);
// lb is a lower bound on the objective, or -DAQP_INF (see daqp_bnb_prox)
int daqp_bnb_lb(DAQPWorkspace* work, const c_float lb);
// As daqp_bnb, with the integer-feasible solution u of the LDP (objective fval in the scale of work->fval, soft slack
// slack) as the incumbent: the solution if the tree has nothing better, or ends early (DAQP_EXIT_FEASIBLE)
int daqp_bnb_from(DAQPWorkspace* work, const c_float* u, const c_float fval, const c_float slack);
int daqp_process_node(DAQPNode* node, DAQPWorkspace* work);
int daqp_get_branch_id(DAQPWorkspace* work);
void daqp_spawn_children(DAQPNode* node, const int branch_id, DAQPWorkspace* work);

void daqp_node_cleanup_workspace(int n_clean, DAQPWorkspace* work);
int daqp_warmstart_node(DAQPNode* node, DAQPWorkspace* work);
void daqp_save_warmstart(DAQPNode* node, DAQPWorkspace* work);
int daqp_add_upper_lower(const int add_id, DAQPWorkspace* work);
int daqp_setup_cold_bnb(DAQPNode* node, DAQPWorkspace* work);

#define DAQP_LOWER_BIT 16
#define DAQP_EXTRACT_LOWER_FLAG(x) (x>>(DAQP_LOWER_BIT-1))
#define DAQP_REMOVE_LOWER_FLAG(x) (x&~(1<<DAQP_LOWER_BIT))
#define DAQP_ADD_LOWER_FLAG(x) (x|(1<<DAQP_LOWER_BIT))
#define DAQP_TOGGLE_LOWER_FLAG(x) (x^(1<<DAQP_LOWER_BIT))

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_BNB_H
