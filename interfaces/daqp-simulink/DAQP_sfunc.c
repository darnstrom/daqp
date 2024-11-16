/*
 * Author: Christopher Schulte, RWTH Aachen University, IRT
 *         christopher.schulte@rwth-aachen.de
 *        https://irt.rwth-aachen.de/
 *
 * This file is part of the DAQP project.
 *
 * mg = number of matrix constraints
 * m = number of constraints
 * n = number of decision variables
 *
 * INPUTS:
 * 1. H: Quadratic term of the objective function [n x n]
 * 2. f: Linear term of the objective function [n x 1]
 * 3. A: Matrix of the linear constraints [mg x n]
 * 4. lb/b_lower: Lower bounds [m x 1]
 * 5. ub/b_upper: Upper bounds [m x 1]
 * 6. sense: Sense of the constraints [m x 1]
 *
 * PARAMETERS:
 * 1. maxIter: Maximum number of iterations
 *
 * OUTPUTS:
 * 1. x_opt: Optimal solution of the decision variables [n x 1]
 * 2. lambda: Lagrange multipliers [m x 1]
 * 3. fval: Value of the objective function at the optimal solution [1 x 1]
 * 4. flag: Exit flag [1 x 1]
 * 5. iteration: Number of iterations [1 x 1]
 */

#define S_FUNCTION_NAME DAQP_sfunc
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"
#include "mex.h"

/*
 * Need to include the header file for the DAQP API
 */
#include "types.h"
#include "daqp.h"
#include "api.h"
#include "utils.h"

#include <stdlib.h> // For malloc und free
#include <math.h>   // For isnan, sqrt

typedef struct
{
    int_T m;
    int_T n;
    int_T mg;
    int_T ms;
    boolean_T is_valid;
} Dimensions;

#define MIN(a, b) ((a) > (b) ? (b) : (a))

#define INFTY 1.0e20
static void removeNaN(real_T *x, int_T n)
{
    // Replace NaN values with near infinity
    int_T i;
    for (i = 0; i < n; i++)
    {
        if (isnan(x[i]))
        {
            x[i] = INFTY;
        }
    }
}
static void transposeMatrix(InputRealPtrsType x_in, c_float *x, int_T n, int_T m)
{
    // Convert Matrix column-first-order to row-first-order
    // n = number of columns
    // m = number of rows

    // check if x could be successfully allocated
    if (x == NULL)
    {
        return;
    }

    int i, j;

    for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
            x[i * m + j] = (c_float)(*x_in)[j * n + i];
}

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
#ifdef __DEBUG__
    mexPrintf("mdlInitializeSizes started.\n");
#endif

    // ----------------------------- SIMULINK -----------------------------
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);

    // --------------------------- PARAMETERS -----------------------------
    ssSetNumSFcnParams(S, 2); /* Number of expected parameters */

    // ----------------------------- INPUTS -----------------------------
    if (!ssSetNumInputPorts(S, 6))
        return;
    ssSetInputPortMatrixDimensions(S, 0, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    ssSetInputPortWidth(S, 1, DYNAMICALLY_SIZED); /* g */
    ssSetInputPortDirectFeedThrough(S, 1, 1);

    ssSetInputPortMatrixDimensions(S, 2, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
    ssSetInputPortDirectFeedThrough(S, 2, 1);

    ssSetInputPortWidth(S, 3, DYNAMICALLY_SIZED); /* lb */
    ssSetInputPortDirectFeedThrough(S, 3, 1);

    ssSetInputPortWidth(S, 4, DYNAMICALLY_SIZED); /* ub */
    ssSetInputPortDirectFeedThrough(S, 4, 1);

    ssSetInputPortWidth(S, 5, DYNAMICALLY_SIZED); /* sense */
    ssSetInputPortDirectFeedThrough(S, 5, 1);

    // ----------------------------- WORK -------------------------------
    ssSetNumPWork(S, 2); // work, result

    // ----------------------------- OUTPUTS -----------------------------
    if (!ssSetNumOutputPorts(S, 5))
        return;
    ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED); /* x_opt */
    ssSetOutputPortWidth(S, 1, DYNAMICALLY_SIZED); /* lambda */
    ssSetOutputPortWidth(S, 2, 1);                 /* fval */
    ssSetOutputPortWidth(S, 3, 1);                 /* flag */
    ssSetOutputPortWidth(S, 4, 1);                 /* iteration */

#ifdef __DEBUG__
    mexPrintf("mdlInitializeSizes finished.\n");
#endif
}

/* Function: mdlStart =========================================================
 * Abstract:
 *    This function is called once at start of the model execution. If you
 *    have states that should be initialized once, this is the place to do it.
 */
#define MDL_START
Dimensions checkDimensions(SimStruct *S)
{
    Dimensions dim;

    dim.is_valid = false;

    // Check for dimensions
    int_T dim_H = ssGetInputPortNumDimensions(S, 0);
    int_T dim_g = ssGetInputPortNumDimensions(S, 1);
    int_T dim_A = ssGetInputPortNumDimensions(S, 2);
    int_T dim_lb = ssGetInputPortNumDimensions(S, 3);
    int_T dim_ub = ssGetInputPortNumDimensions(S, 4);
    int_T dim_sense = ssGetInputPortNumDimensions(S, 5);

    // Check for dimensions
    if (dim_H != 2)
    {
        mexErrMsgTxt("ERROR (DAQP): Quadratic term of the objective function H must be a 2D matrix.");
    }
    if (dim_g != 1)
    {
        mexErrMsgTxt("ERROR (DAQP): Linear term of the objective function g must be a vector.");
    }
    if (dim_A != 2)
    {
        mexErrMsgTxt("ERROR (DAQP): Matrix of the linear constraints A must be a 2D matrix.");
    }
    if (dim_lb != 1)
    {
        mexErrMsgTxt("ERROR (DAQP): Lower bounds b_lower must be a vector.");
    }
    if (dim_ub != 1)
    {
        mexErrMsgTxt("ERROR (DAQP): Upper bounds b_upper must be a vector.");
    }
    if (dim_sense != 1)
    {
        mexErrMsgTxt("ERROR (DAQP): Sense of the constraints must be a vector.");
    }

    // Get sizes
    int_T *size_H = ssGetInputPortDimensions(S, 0);
    int_T size_g = ssGetInputPortWidth(S, 1);
    int_T *size_A = ssGetInputPortDimensions(S, 2);
    int_T size_lb = ssGetInputPortWidth(S, 3);
    int_T size_ub = ssGetInputPortWidth(S, 4);
    int_T size_sense = ssGetInputPortWidth(S, 5);

    // Check for connected ports
    if (!(ssGetInputPortConnected(S, 0)))
        mexErrMsgTxt("ERROR (DAQP): Quadratic term of the objective function H is not connected (Port 1).");
    if (!(ssGetInputPortConnected(S, 1)))
        mexErrMsgTxt("ERROR (DAQP): Linear term of the objective function g is not connected (Port 2).");

    boolean_T is_A_connected = ssGetInputPortConnected(S, 2);
    boolean_T is_lb_connected = ssGetInputPortConnected(S, 3);
    boolean_T is_ub_connected = ssGetInputPortConnected(S, 4);
    boolean_T is_sense_connected = ssGetInputPortConnected(S, 5);

    if (!(is_A_connected || is_lb_connected || is_ub_connected || is_sense_connected))
    {
        // Manually set the dimension to [1 x size_g]
        size_A[0] = 1;
        size_A[1] = size_g;
        size_lb = 1;
        size_ub = 1;
        size_sense = 1;

#ifdef __DEBUGInfo__
        mexPrintf("DAQP: A, lb, ub, sense are not connected. Automatically set A to zeros, lb to -Infty, ub to Infty, and sense to 0.\n");
#endif
    }
    if (is_A_connected && !(is_lb_connected || is_ub_connected || is_sense_connected))
    {
        int_T mg = size_A[0];
        size_lb = mg;
        size_ub = mg;
        size_sense = mg;
#ifdef __DEBUGInfo__
        mexPrintf("DAQP: lb, ub and sense are not connected. Automatically set lb to -Infty, ub to Infty, and sense to 0.\n");
#endif
    }
    else if (is_A_connected && is_ub_connected && !(is_lb_connected || is_sense_connected))
    {
        size_lb = size_ub;
        size_sense = size_ub;
#ifdef __DEBUGInfo__
        mexPrintf("DAQP: lb and sense are not connected. Automatically set lb to -Infty and sense to 0.\n");
#endif
    }
    else if (is_A_connected && !(is_lb_connected) && is_ub_connected && is_sense_connected)
    {
        size_lb = size_ub;
#ifdef __DEBUGInfo__
        mexPrintf("DAQP: lb is not connected. Automatically set lb to -Infty.\n");
#endif
    }
    else if (is_A_connected && is_lb_connected && !(is_ub_connected) && !(is_sense_connected))
    {
        size_ub = size_lb;
        size_sense = size_lb;
#ifdef __DEBUGInfo__
        mexPrintf("DAQP: ub and sense are not connected. Automatically set ub to Infty and sense to 0.\n");
#endif
    }
    else if (is_A_connected && is_lb_connected && !(is_ub_connected) && is_sense_connected)
    {
        size_ub = size_lb;
#ifdef __DEBUGInfo__
        mexPrintf("DAQP: ub is not connected. Automatically set ub to Infty.\n");
#endif
    }
    else if (is_A_connected && is_lb_connected && is_ub_connected && !(is_sense_connected))
    {
        size_sense = size_lb;
#ifdef __DEBUGInfo__
        mexPrintf("DAQP: sense is not connected. Automatically set sense to 0.\n");
#endif
    }

    // Calculate dimensions
    dim.m = size_sense;
    dim.n = size_H[0];
    dim.mg = size_A[0];
    dim.ms = dim.m - dim.mg;

    // Check for empty matrices
    if (size_H[0] == 0 || size_H[1] == 0)
        mexErrMsgTxt("ERROR (DAQP): Quadratic term of the objective function H is empty.");
    if (size_A[0] == 0 || size_A[1] == 0)
        mexErrMsgTxt("ERROR (DAQP): Matrix of the linear constraints A is empty.");
    if (size_g == 0)
        mexErrMsgTxt("ERROR (DAQP): Linear term of the objective function g is empty.");
    if (size_lb == 0)
        mexErrMsgTxt("ERROR (DAQP): Lower bounds b_lower is empty.");
    if (size_ub == 0)
        mexErrMsgTxt("ERROR (DAQP): Upper bounds b_upper is empty.");
    if (size_sense == 0)
        mexErrMsgTxt("ERROR (DAQP): Sense of the constraints is empty.");

    // Check for dimension mismatch
    if (size_H[0] != size_H[1])
        mexErrMsgTxt("ERROR (DAQP): Quadratic term H must be a square matrix.");
    if (size_g != dim.n)
        mexErrMsgTxt("ERROR (DAQP): Dimension mismatch between the linear term g of the objective function and the number of decision variables.");
    if (size_lb != dim.m)
        mexErrMsgTxt("ERROR (DAQP): Dimension mismatch in the lower bounds b_lower or sense.");
    if (size_ub != dim.m)
        mexErrMsgTxt("ERROR (DAQP): Dimension mismatch in the upper bounds b_upper or sense.");
    if (size_A[1] != dim.n)
        mexErrMsgTxt("ERROR (DAQP): Dimension 2 of the matrix of the linear constraints A must be equal to the number of decision variables.");

    // Check for invalid sizes
    if (dim.mg < 1)
        mexErrMsgTxt("ERROR (DAQP): Number of general constraints must be greater than zero.");
    if (dim.ms > dim.n)
        mexErrMsgTxt("ERROR (DAQP): Number of simple constraints must be less than or equal to the number of decision variables. Therefore the input lb, ub, and sense are too long or A has too less rows.");

    dim.is_valid = true;
    return dim;
}
DAQPWorkspace *daqp_sfun_Start_wrapper(Dimensions dim)
{
    int n = dim.n;
    int m = dim.m;
    int ms = dim.ms;

    // Allocate work struct
    DAQPWorkspace *work = (DAQPWorkspace *)calloc(1, sizeof(DAQPWorkspace));
    // Allocate settings
    allocate_daqp_settings(work);
    daqp_default_settings(work->settings);

    allocate_daqp_workspace(work, n, 0);
    work->m = m;
    work->ms = ms;

    // Allocate qp
    DAQPProblem *qp = calloc(1, sizeof(DAQPProblem));
    work->qp = qp;
    work->qp->n = work->n;
    work->qp->m = work->m;
    work->qp->ms = work->ms;

    work->qp->H = malloc(n * n * sizeof(c_float));
    work->qp->A = malloc(n * (m - ms) * sizeof(c_float));
    work->qp->f = malloc(n * sizeof(c_float));
    work->qp->bupper = malloc(m * sizeof(c_float));
    work->qp->blower = malloc(m * sizeof(c_float));

    // Setup daqp workspace
    work->Rinv = malloc(((n + 1) * n / 2) * sizeof(c_float));
    work->M = malloc(n * (m - ms) * sizeof(c_float));
    work->scaling = malloc(m * sizeof(c_float));
    work->dupper = malloc(m * sizeof(c_float));
    work->dlower = malloc(m * sizeof(c_float));
    work->v = malloc(n * sizeof(c_float));
    work->sense = calloc(1, m * sizeof(int));
    work->qp->sense = work->sense;

    // output the pointer to the workspace
    return work;
}
static void mdlStart(SimStruct *S)
{
// ----------------------------- WORK -----------------------------
// Dimensions as work
#ifdef __DEBUG__
#define __DEBUGInfo__
#endif

    Dimensions dim = checkDimensions(S);
    DAQPWorkspace *work = daqp_sfun_Start_wrapper(dim);

    DAQPResult *result = (DAQPResult *)calloc(1, sizeof(DAQPResult));
    result->lam = (c_float *)calloc(dim.m, sizeof(c_float));
    result->x = (c_float *)calloc(dim.n, sizeof(c_float));

    ssSetPWorkValue(S, 0, work);
    ssSetPWorkValue(S, 1, result);

#ifdef __DEBUG__
    mexPrintf("DAQP: Dimensions:\n");
    mexPrintf("  Signal   [n x m], Connected [Y/N]\n");
    mexPrintf(" - x_opt:  [%d x 1], %c\n", dim.n, ssGetOutputPortConnected(S, 0) ? 'Y' : 'N');
    mexPrintf(" - lambda: [%d x 1], %c\n", dim.m, ssGetOutputPortConnected(S, 1) ? 'Y' : 'N');
    mexPrintf(" - H:      [%d x %d], %c\n", dim.n, dim.n, ssGetInputPortConnected(S, 0) ? 'Y' : 'N');
    mexPrintf(" - g:      [%d x 1], %c\n", dim.n, ssGetInputPortConnected(S, 1) ? 'Y' : 'N');
    mexPrintf(" - A:      [%d x %d], %c\n", dim.mg, dim.n, ssGetInputPortConnected(S, 2) ? 'Y' : 'N');
    mexPrintf(" - lb:     [%d x 1], %c\n", dim.m, ssGetInputPortConnected(S, 3) ? 'Y' : 'N');
    mexPrintf(" - ub:     [%d x 1], %c\n", dim.m, ssGetInputPortConnected(S, 4) ? 'Y' : 'N');
    mexPrintf(" - sense:  [%d x 1], %c\n", dim.m, ssGetInputPortConnected(S, 5) ? 'Y' : 'N');
#endif
}

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block.
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    // ----------------------------- PARAMS -----------------------------
    // Get number of maximum working set iterations
    int_T WARM_START, maxIter;
    if (ssGetSFcnParamsCount(S) < 2)
    {
        WARM_START = 0; // Cold start
    }
    else
    {
        DTypeId dataType;
        ssGetSFcnParamDataType(S, 1, &dataType);
        if (dataType != SS_DOUBLE)
        {
            mexErrMsgTxt("ERROR (DAQP): WARM_START must be of type double. %d\n");
        }
        const mxArray *in_WARM_START = ssGetSFcnParam(S, 1);

        // Check if the parameter is a scalar
        if (mxGetNumberOfElements(in_WARM_START) != 1)
        {
            mexErrMsgTxt("ERROR (DAQP): WARM_START must be a scalar.");
        }

        WARM_START = (int_T)((mxGetPr(in_WARM_START))[0]); // Warm start
    }

    if (ssGetSFcnParamsCount(S) < 1)
    {
        maxIter = 100; // Default maximum number of iterations
    }
    else
    {
        DTypeId dataType;
        ssGetSFcnParamDataType(S, 0, &dataType);
        if (dataType != SS_DOUBLE)
        {
            mexErrMsgTxt("ERROR (DAQP): maxIter must be of type double.");
        }
        const mxArray *in_maxIter = ssGetSFcnParam(S, 0);

        // Check if the parameter is a scalar
        if (mxGetNumberOfElements(in_maxIter) != 1)
        {
            mexErrMsgTxt("ERROR (DAQP): maxIter must be a scalar.");
        }

        maxIter = (int_T)((mxGetPr(in_maxIter))[0]); // Maximum number of iterations

        if (maxIter <= 0)
        {
            mexErrMsgTxt("ERROR (DAQP): Maximum number of iterations must be greater than zero.");
        }
    }

    // ----------------------------- INPUTS -----------------------------
    InputRealPtrsType in_H, in_f, in_A, in_lb, in_ub, in_sense;
    in_H = ssGetInputPortRealSignalPtrs(S, 0);     // Quadratic term of the objective function
    in_f = ssGetInputPortRealSignalPtrs(S, 1);     // Linear term of the objective function
    in_A = ssGetInputPortRealSignalPtrs(S, 2);     // Matrix of the linear constraints
    in_lb = ssGetInputPortRealSignalPtrs(S, 3);    // Lower bounds
    in_ub = ssGetInputPortRealSignalPtrs(S, 4);    // Upper bounds
    in_sense = ssGetInputPortRealSignalPtrs(S, 5); // Sense of the constraints

    // ----------------------------- Get Work -----------------------------
    DAQPWorkspace *work = (DAQPWorkspace *)ssGetPWorkValue(S, 0);
    work->settings->iter_limit = maxIter;

    // ----------------------------- OUTPUTS -----------------------------
    real_T *out_x_opt = ssGetOutputPortRealSignal(S, 0);
    real_T *out_lambda = ssGetOutputPortRealSignal(S, 1);
    real_T *out_fval = ssGetOutputPortRealSignal(S, 2);
    real_T *out_flag = ssGetOutputPortRealSignal(S, 3);
    real_T *out_iteration = ssGetOutputPortRealSignal(S, 4);

    // ------------------------- DAQP -------------------------
    int ii, jj;
    int n = work->n;
    int m = work->m;
    int mg = m - work->ms;
    int update_mask = 0;

    // Update H and A
    if (WARM_START < 2 || work->iterations == 0)
    {
        // Convert input to row-first-order
        transposeMatrix(in_H, work->qp->H, work->n, work->n);
        transposeMatrix(in_A, work->qp->A, mg, work->n);

        // Check if A is connected
        if (!(ssGetInputPortConnected(S, 2)))
        {
            for (ii = mg; ii < work->m * work->n; ii++)
            {
                work->qp->A[ii] = 0.0;
            }
        }
        update_mask += UPDATE_Rinv + UPDATE_M;
    }

    // Update f
    update_mask += UPDATE_v + UPDATE_d;
    for (ii = 0; ii < work->n; ii++)
        work->qp->f[ii] = (c_float)(*in_f)[ii];

    // Update blower, bupper and sense
    boolean_T Port3Connected = ssGetInputPortConnected(S, 3);
    boolean_T Port4Connected = ssGetInputPortConnected(S, 4);
    boolean_T Port5Connected = ssGetInputPortConnected(S, 5);
    for (ii = 0; ii < work->m; ii++)
    {
        if (Port3Connected)
        {
            work->qp->blower[ii] = (c_float)(*in_lb)[ii];
        }
        else
        {
            work->qp->blower[ii] = -INFTY; // input not connected
        }

        if (Port4Connected)
        {
            work->qp->bupper[ii] = (c_float)(*in_ub)[ii];
        }
        else
        {
            work->qp->bupper[ii] = INFTY; // input not connected
        }

        if (Port5Connected)
        {
            work->qp->sense[ii] = (int_T)(*in_sense)[ii];
        }
        else
        {
            work->qp->sense[ii] = 0; // input not connected
        }
    }
    // Cleanup sense
    if (WARM_START == 0)
    {
        deactivate_constraints(work);
    }
    // Update internal problem
    update_ldp(update_mask, work);

    // Reuse previous working set
    if (WARM_START == 1)
        activate_constraints(work);

    // Solve problem
    out_flag[0] = daqp_ldp(work);
    ldp2qp_solution(work);

    DAQPResult *result = (DAQPResult *)ssGetPWorkValue(S, 1);
    daqp_extract_result(result, work);

// DEBUG
#ifdef __DEBUG__
    mexPrintf("DAQP: \n");
    mexPrintf(" - H: ");
    for (ii = 0; ii < MIN(n * n, 50); ii++)
    {
        mexPrintf("%0.2f ", work->qp->H[ii]);
    }
    mexPrintf("\n");
    mexPrintf(" - f: ");
    for (ii = 0; ii < MIN(n, 50); ii++)
    {
        mexPrintf("%0.2f ", work->qp->f[ii]);
    }
    mexPrintf("\n");
    mexPrintf(" - A: ");
    for (ii = 0; ii < MIN(mg * n, 50); ii++)
    {
        mexPrintf("%0.2f ", work->qp->A[ii]);
    }
    mexPrintf("\n");
    mexPrintf(" - lb: ");
    for (ii = 0; ii < MIN(m, 50); ii++)
    {
        mexPrintf("%0.2f ", work->qp->blower[ii]);
    }
    mexPrintf("\n");
    mexPrintf(" - ub: ");
    for (ii = 0; ii < MIN(m, 50); ii++)
    {
        mexPrintf("%0.2f ", work->qp->bupper[ii]);
    }
    mexPrintf("\n");
    mexPrintf(" - sense: ");
    for (ii = 0; ii < MIN(m, 50); ii++)
    {
        mexPrintf("%d ", work->qp->sense[ii]);
    }
    mexPrintf("\n");
    mexPrintf(" - x_opt: ");
    for (ii = 0; ii < MIN(n, 50); ii++)
    {
        mexPrintf("%0.2f ", work->x[ii]);
    }
    mexPrintf("\n");
#endif

    // ------------------------- WRITE OUTPUTS -------------------------
    if (work->x == NULL || work->lam == NULL)
    {
        mexErrMsgTxt("ERROR (DAQP): Memory allocation failed.");
    }

    for (ii = 0; ii < work->n; ++ii)
        out_x_opt[ii] = (real_T)(result->x[ii]);

    if (work->lam != NULL)
        for (ii = 0; ii < work->m; ++ii)
            out_lambda[ii] = (real_T)(result->lam[ii]);

    out_fval[0] = (real_T)(result->fval);
    out_iteration[0] = (real_T)(result->iter);

    // -------------------------- REMOVE NaN -------------------------
    removeNaN(out_x_opt, work->n);
    removeNaN(out_lambda, m - work->ms);
    removeNaN(out_fval, 1);
    removeNaN(out_flag, 1);
    removeNaN(out_iteration, 1);
}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
    // Free Work
    DAQPWorkspace *work = (DAQPWorkspace *)ssGetPWorkValue(S, 0);
    free_daqp_workspace(work);
    free_daqp_ldp(work);

    free(work->qp->H);
    free(work->qp->f);
    free(work->qp->A);
    free(work->qp->bupper);
    free(work->qp->blower);
    free(work->qp);

    free(work);

    // Free Result
    DAQPResult *result = (DAQPResult *)ssGetPWorkValue(S, 1);
    free(result->lam);
    free(result->x);
    free(result);

#ifdef __DEBUG__
    mexPrintf("mdlTerminate finished.\n");
#endif
}

#if defined(MATLAB_MEX_FILE)

/* Function: mdlSetInputPortDimensionInfo + mdlSetOutputPortDimensionInfo
 * Abstract:
 *    Needed for dynamic dimension
 */
#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
#define MDL_SET_DEFAULT_PORT_DIMENSION_INFO

static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if (!ssSetInputPortDimensionInfo(S, port, dimsInfo))
        return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if (!ssSetOutputPortDimensionInfo(S, port, dimsInfo))
        return;
}

/* Function: mdlSetDefaultPortDimensionInfo
 * Abstract:
 *    Needed for missing inputs
 */
static void setOutputPortDimensions(SimStruct *S, int_T port, int_T dim1, int_T dim2)
{
    DimsInfo_T dimsInfo;
    dimsInfo.dims = (int *)calloc(2, sizeof(int));
    if (dim2 == 1)
    {
        dimsInfo.numDims = 1;
        dimsInfo.dims[0] = dim1;
        dimsInfo.dims[1] = dim2;
        dimsInfo.width = dim1 * dim2;
    }
    else
    {
        dimsInfo.numDims = 2;
        dimsInfo.dims[0] = dim1;
        dimsInfo.dims[1] = dim2;
        dimsInfo.width = dim1 * dim2;
    }
    mdlSetOutputPortDimensionInfo(S, port, &dimsInfo);

    // Free memory
    free(dimsInfo.dims);
}
static void setInputPortDimensions(SimStruct *S, int_T port, int_T dim1, int_T dim2)
{
    DimsInfo_T dimsInfo;
    dimsInfo.dims = (int *)calloc(2, sizeof(int));
    if (dim2 == 1)
    {
        dimsInfo.numDims = 1;
        dimsInfo.dims[0] = dim1;
        dimsInfo.dims[1] = dim2;
        dimsInfo.width = dim1 * dim2;
    }
    else
    {
        dimsInfo.numDims = 2;
        dimsInfo.dims[0] = dim1;
        dimsInfo.dims[1] = dim2;
        dimsInfo.width = dim1 * dim2;
    }
    mdlSetInputPortDimensionInfo(S, port, &dimsInfo);

    // Free memory
    free(dimsInfo.dims);
}
static void mdlSetDefaultPortDimensionInfo(SimStruct *S)
{
#ifdef __DEBUG__
    mexPrintf("mdlSetDefaultPortDimensionInfo started.\n");
#endif

    // Get Dimensions
    Dimensions dim = checkDimensions(S);

    // Set Outputs
    setOutputPortDimensions(S, 0, dim.n, 1);
    setOutputPortDimensions(S, 1, dim.m, 1);

    // Set Inputs
    setInputPortDimensions(S, 0, dim.n, dim.n);
    setInputPortDimensions(S, 1, dim.n, 1);
    setInputPortDimensions(S, 2, dim.mg, dim.n);
    setInputPortDimensions(S, 3, dim.m, 1);
    setInputPortDimensions(S, 4, dim.m, 1);
    setInputPortDimensions(S, 5, dim.m, 1);

#ifdef __DEBUG__
    mexPrintf("mdlSetDefaultPortDimensionInfo finished.\n");
#endif
}

#endif

/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */
#include "simulink.c"  /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif