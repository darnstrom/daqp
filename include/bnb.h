#ifndef DAQP_BNB_H 
# define DAQP_BNB_H 

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"
#include "constants.h"
#include "daqp.h"

int daqp_bnb(DAQPWorkspace* work);
int process_node(DAQPNode* node, DAQPWorkspace* work);
int get_branch_id(DAQPWorkspace* work);
void spawn_children(DAQPNode* node, const int branch_id, DAQPWorkspace* work);

void node_cleanup_workspace(int n_clean, DAQPWorkspace* work);
void warmstart_node(DAQPNode* node, DAQPWorkspace* work); 
void save_warmstart(DAQPNode* node, DAQPWorkspace* work);
int add_upper_lower(const int add_id, DAQPWorkspace* work); 
void setup_cold_bnb(DAQPNode* node, DAQPWorkspace* work);

#define DAQP_LOWER_BIT 16
#define DAQP_EXTRACT_LOWER_FLAG(x) (x>>(DAQP_LOWER_BIT-1))
#define DAQP_REMOVE_LOWER_FLAG(x) (x&~(1<<DAQP_LOWER_BIT))
#define DAQP_ADD_LOWER_FLAG(x) (x|(1<<DAQP_LOWER_BIT))
#define DAQP_TOGGLE_LOWER_FLAG(x) (x^(1<<DAQP_LOWER_BIT))

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_BNB_H 
