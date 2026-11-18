#ifndef DAQP_MULTISTAGE_H
# define DAQP_MULTISTAGE_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "api.h"

struct blasfeo_dmat;

/*
 * Structure-exploiting solver for multistage QPs (stages coupled to the next stage, as in optimal control)
 * (built with the CMake option DAQP_MULTISTAGE, requires BLASFEO). See
 * docs/docs/multistage.md for the data, the constraint conventions, and the method.
 *
 * min  sum_{t<N} 0.5 x'Q x + u'S x + 0.5 u'R u + q'x + r'u + 0.5 x_N'Q_N x_N + q_N'x_N
 * s.t. x_{t+1} = A x_t + B u_t + b,
 *      lb <= [u; x][idxb] <= ub,  lg <= C x + D u <= ug,
 *      lc <= Ex x_t + Eu u_t + Fx x_{t+1} + Fu u_{t+1} <= uc
 */
typedef struct{
    int N;
    int *nx; // N+1
    int *nu; // N
    int *nb; // N+1 (NULL: none)
    int *ng; // N+1 (NULL: none)
    int *nc; // N (NULL: none)

    c_float **A, **B, **b; // Row major (NULL A/b: zero, NULL B: identity)
    c_float **Q, **S, **R, **q, **r; // (NULL: zero)

    int **idxb; // Indices into [u_t; x_t]
    c_float **lb, **ub;
    c_float **C, **D, **lg, **ug;
    c_float **Ex, **Eu, **Fx, **Fu, **lc, **uc;

    int *sense; // DAQP_SOFT or DAQP_BINARY per constraint (NULL: all hard)
    int binary_horizon; // DAQP_BINARY only at stages t < binary_horizon, relaxed at the others (0: at all stages)

    // BLASFEO matrices that the solver reads (referenced; daqp_multistage_pack forms them):
    // BAbt = [B'; A'], RSQrq = [R S; S' Q] (lower), DCt = [D'; C'], Et = [Eu'; Ex'], Ft = [Fu'; Fx']
    struct blasfeo_dmat *BAbt, *RSQrq, *DCt, *Et, *Ft;
}DAQPMultistage;

// Flags for daqp_multistage_update (changes of b, q, r and the bounds are detected: mask 0)
#define DAQP_MULTISTAGE_UPDATE_MATRICES 1 // Matrices (refactorizes)
#define DAQP_MULTISTAGE_UPDATE_SHIFT 8 // The problem moved by one stage, as in MPC: the warm start may be
                                      // the working set and the multipliers of the latest solve one stage earlier

// Set up a workspace for prob (referenced, not copied); frees it on failure
int setup_daqp_multistage(DAQPMultistage *prob, DAQPWorkspace *work);
// Bytes of memory for setup_daqp_multistage_mem, from the dimensions of prob alone (no data; ns: the largest
// number of soft constraints). Sized for the worst case (every constraint a candidate for elimination, x_0 free)
size_t daqp_multistage_memsize(const DAQPMultistage *prob, const int ns);
// As setup_daqp_multistage, in mem (size bytes): nothing is allocated. work->settings is given, and the arrays
// of the core workspace for n = daqp_multistage_nvmax(prob) (or daqp_multistage_nv(prob)) and ns soft
// constraints (lam, lam_star, WS, L, D, xldl, zldl, u, and the soft weights if ns > 0). With binary constraints,
// also xold (n), and branch and bound is allocated
int setup_daqp_multistage_mem(DAQPMultistage *prob, DAQPWorkspace *work, void *mem, size_t size);
// Largest number of variables of the LDP (x_0 free)
int daqp_multistage_nvmax(const DAQPMultistage *prob);
// Number of variables of the LDP for the data of prob (x_0 free only if it is not fixed)
int daqp_multistage_nv(const DAQPMultistage *prob);
// Refresh after the data has changed (mask: 0 or DAQP_MULTISTAGE_UPDATE_*), keeping the working set
int daqp_multistage_update(DAQPWorkspace *work, const int mask);
// res->x = [x_0; u_0; x_1; ...; x_N]; res->lam (in the order of the constraints): the multipliers of the
// working set (those of eliminated equalities and of the bounds that fix x_0 are zero, see daqp_multistage_duals)
void daqp_multistage_solve(DAQPResult *res, DAQPWorkspace *work);
// Complete multipliers lam (m) and multipliers of the dynamics pi = [pi_0; ...; pi_{N-1}] (either NULL) of the
// latest solve, valid until the next update or solve. Returns 0 if there are none (after branch and bound), 2 if the
// elimination of equalities was found inaccurate (see docs/docs/multistage.md)
int daqp_multistage_duals(DAQPWorkspace *work, c_float *lam, c_float *pi);
// Whether DAQP cannot change the value of constraint id (a zero row of M, or a bound on a fixed x_0)
int daqp_multistage_fixed_row(const DAQPWorkspace *work, const int id);
void free_daqp_multistage(DAQPWorkspace *work);

// Form the BLASFEO matrices from the row-major ones (allocated at the first call, repacked in place afterwards).
// Consecutive stages with the same row-major matrices (pointers) share one BLASFEO matrix; if that sharing has
// changed since the matrices were allocated, they are allocated again (pass DAQP_MULTISTAGE_UPDATE_MATRICES).
// Returns 0, or -1 if out of memory
int daqp_multistage_pack(DAQPMultistage *prob);
void daqp_multistage_pack_free(DAQPMultistage *prob);

// Lengths of the stacked primal solution and of the constraints
int daqp_multistage_nz(const DAQPMultistage *prob);
int daqp_multistage_m(const DAQPMultistage *prob);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef DAQP_MULTISTAGE_H
