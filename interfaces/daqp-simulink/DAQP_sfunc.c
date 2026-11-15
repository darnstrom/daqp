/*
 * Author: Christopher Schulte, RWTH Aachen University, IRT
 *         christopher.schulte@rwth-aachen.de
 *        https://irt.rwth-aachen.de/
 *
 * This file is part of the DAQP project.
 *
 * Solves the QP
 *   min 0.5 x'*H*x + f'*x
 *   s.t. lb <= [I 0; A]*x <= ub
 * at every time step, where the first ms = m-mg elements of lb/ub are simple
 * bounds on the first ms elements of x.
 *
 * n  = number of decision variables
 * m  = number of constraints
 * mg = number of general constraints (rows of A)
 *
 * INPUTS:
 * 1. H: Quadratic term of the objective function [n x n]
 * 2. f: Linear term of the objective function [n x 1]
 * 3. A: Matrix of the general constraints [mg x n] (optional)
 * 4. lb: Lower bounds [m x 1] (optional, -inf if not connected)
 * 5. ub: Upper bounds [m x 1] (optional, +inf if not connected)
 * 6. sense: Sense of the constraints [m x 1] (optional, 0 if not connected)
 *
 * PARAMETERS:
 * 1. maxIter: Maximum number of iterations
 * 2. warm_start (optional, default 0):
 *      0 - cold start: every solve starts from the sense input
 *      1 - warm start: start from the working set of the previous solve
 *      2 - as 1, but H and A are assumed constant (only read at the first solve)
 * 3. break_points (optional, default []): hierarchical QP, where
 *    break_points(i) is the (1-based) index of the last constraint in
 *    level i, and break_points(end) = m.
 *
 * OUTPUTS:
 * 1. x: Optimal solution [n x 1]
 * 2. lambda: Lagrange multipliers [m x 1] ([1 x 1] if there are no constraints)
 * 3. fval: Objective function value at x [1 x 1]
 * 4. flag: Exit flag [1 x 1]
 * 5. iter: Number of iterations [1 x 1]
 */

#define S_FUNCTION_NAME DAQP_sfunc
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "daqp_sfunc_core.h"
#include <stdlib.h>

#define P_MAXITER 0
#define P_WARMSTART 1
#define P_BREAKPOINTS 2
#define N_PARAMS_MAX 3

enum {IN_H, IN_F, IN_A, IN_LB, IN_UB, IN_SENSE};
enum {OUT_X, OUT_LAMBDA, OUT_FVAL, OUT_FLAG, OUT_ITER};

typedef struct{
    int n, m, mg, ms;
} Dimensions;

/* ========================== Parameters ========================== */
// (The types of the parameters are validated by mdlCheckParameters)
static real_T getScalarParam(SimStruct *S, int ind, real_T default_value)
{
    if (ssGetSFcnParamsCount(S) <= ind)
        return default_value;
    return mxGetPr(ssGetSFcnParam(S, ind))[0];
}

// Returns the number of break points, which are allocated in *break_points
static int getBreakPoints(SimStruct *S, int **break_points)
{
    const mxArray *p;
    int i, nh;
    *break_points = NULL;
    if (ssGetSFcnParamsCount(S) <= P_BREAKPOINTS)
        return 0;
    p = ssGetSFcnParam(S, P_BREAKPOINTS);
    nh = (int)mxGetNumberOfElements(p);
    if (nh == 0)
        return 0;
    *break_points = (int *)malloc(nh * sizeof(int));
    if (*break_points == NULL)
        return 0;
    for (i = 0; i < nh; i++)
        (*break_points)[i] = (int)mxGetPr(p)[i];
    return nh;
}

#define MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS) && defined(MATLAB_MEX_FILE)
static int isRealDouble(const mxArray *p)
{
    return mxIsDouble(p) && !mxIsComplex(p) && !mxIsSparse(p);
}

static void mdlCheckParameters(SimStruct *S)
{
    const int np = ssGetSFcnParamsCount(S);
    const mxArray *p;
    real_T value;
    int i, nh;
    if (np < 1 || np > N_PARAMS_MAX)
    {
        ssSetErrorStatus(S, "DAQP: expected the parameters maxIter[, warm_start[, break_points]].");
        return;
    }
    p = ssGetSFcnParam(S, P_MAXITER);
    if (!isRealDouble(p) || mxGetNumberOfElements(p) != 1 || mxGetPr(p)[0] < 1)
    {
        ssSetErrorStatus(S, "DAQP: maxIter must be a positive scalar.");
        return;
    }
    if (np > P_WARMSTART)
    {
        p = ssGetSFcnParam(S, P_WARMSTART);
        value = (isRealDouble(p) && mxGetNumberOfElements(p) == 1) ? mxGetPr(p)[0] : -1;
        if (!(value == 0 || value == 1 || value == 2))
        {
            ssSetErrorStatus(S, "DAQP: warm_start must be 0 (cold start), 1 (warm start), or 2 (warm start with constant H and A).");
            return;
        }
    }
    if (np > P_BREAKPOINTS)
    {
        p = ssGetSFcnParam(S, P_BREAKPOINTS);
        if (!isRealDouble(p))
        {
            ssSetErrorStatus(S, "DAQP: break_points must be a real vector.");
            return;
        }
        nh = (int)mxGetNumberOfElements(p);
        for (i = 0; i < nh; i++)
        {
            value = mxGetPr(p)[i];
            if (value != (int)value || value < 0 || (i > 0 && value <= mxGetPr(p)[i - 1]))
            {
                ssSetErrorStatus(S, "DAQP: break_points must be nonnegative, integer, and strictly increasing.");
                return;
            }
        }
    }
}
#endif

/* ========================== Dimensions ========================== */
static void getMatrixDims(SimStruct *S, int port, int *rows, int *cols)
{
    if (ssGetInputPortNumDimensions(S, port) >= 2)
    {
        *rows = ssGetInputPortDimensions(S, port)[0];
        *cols = ssGetInputPortDimensions(S, port)[1];
    }
    else // A vector is interpreted as a row
    {
        *rows = 1;
        *cols = ssGetInputPortWidth(S, port);
    }
}

// Get the problem dimensions from the connected input ports
static int getDimensions(SimStruct *S, Dimensions *dim)
{
    int rows, cols, port;
    if (!ssGetInputPortConnected(S, IN_H))
    {
        ssSetErrorStatus(S, "DAQP: H (port 1) must be connected.");
        return 0;
    }
    if (!ssGetInputPortConnected(S, IN_F))
    {
        ssSetErrorStatus(S, "DAQP: f (port 2) must be connected.");
        return 0;
    }
    getMatrixDims(S, IN_H, &rows, &cols);
    if (rows != cols || rows < 1)
    {
        ssSetErrorStatus(S, "DAQP: H must be a nonempty square matrix.");
        return 0;
    }
    dim->n = rows;
    if (ssGetInputPortWidth(S, IN_F) != dim->n)
    {
        ssSetErrorStatus(S, "DAQP: f must have as many elements as H has rows.");
        return 0;
    }
    dim->mg = 0;
    if (ssGetInputPortConnected(S, IN_A))
    {
        getMatrixDims(S, IN_A, &rows, &cols);
        if (cols != dim->n)
        {
            ssSetErrorStatus(S, "DAQP: A must have as many columns as H.");
            return 0;
        }
        dim->mg = rows;
    }
    dim->m = -1;
    for (port = IN_LB; port <= IN_SENSE; port++)
    {
        if (!ssGetInputPortConnected(S, port))
            continue;
        if (dim->m < 0)
            dim->m = ssGetInputPortWidth(S, port);
        else if (ssGetInputPortWidth(S, port) != dim->m)
        {
            ssSetErrorStatus(S, "DAQP: lb, ub, and sense must have the same number of elements.");
            return 0;
        }
    }
    if (dim->m < 0)
        dim->m = dim->mg;
    dim->ms = dim->m - dim->mg;
    if (dim->ms < 0 || dim->ms > dim->n)
    {
        ssSetErrorStatus(S, "DAQP: lb/ub/sense must have between mg and mg+n elements (mg = number of rows of A, n = number of variables).");
        return 0;
    }
    return 1;
}

/* ========================== Initialization ========================== */
static void mdlInitializeSizes(SimStruct *S)
{
    int i;
    ssSetNumSFcnParams(S, -1); // 1-3 parameters (checked below)
#if defined(MATLAB_MEX_FILE)
    mdlCheckParameters(S);
    if (ssGetErrorStatus(S) != NULL)
        return;
#endif
    for (i = 0; i < ssGetSFcnParamsCount(S); i++)
        ssSetSFcnParamTunable(S, i, SS_PRM_NOT_TUNABLE);

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 6))
        return;
    ssSetInputPortMatrixDimensions(S, IN_H, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, IN_F, DYNAMICALLY_SIZED);
    ssSetInputPortMatrixDimensions(S, IN_A, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, IN_LB, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, IN_UB, DYNAMICALLY_SIZED);
    ssSetInputPortWidth(S, IN_SENSE, DYNAMICALLY_SIZED);
    for (i = 0; i < 6; i++)
    {
        ssSetInputPortDirectFeedThrough(S, i, 1);
        ssSetInputPortRequiredContiguous(S, i, 1);
    }

    if (!ssSetNumOutputPorts(S, 5))
        return;
    ssSetOutputPortWidth(S, OUT_X, DYNAMICALLY_SIZED);
    ssSetOutputPortWidth(S, OUT_LAMBDA, DYNAMICALLY_SIZED);
    ssSetOutputPortWidth(S, OUT_FVAL, 1);
    ssSetOutputPortWidth(S, OUT_FLAG, 1);
    ssSetOutputPortWidth(S, OUT_ITER, 1);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 1); // DAQPSfunc
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE | SS_OPTION_CALL_TERMINATE_ON_EXIT);
}

#if defined(MATLAB_MEX_FILE)
#define MDL_SET_INPUT_PORT_DIMENSION_INFO
static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    ssSetInputPortDimensionInfo(S, port, dimsInfo);
}

#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    ssSetOutputPortDimensionInfo(S, port, dimsInfo);
}

// Set the dimensions of the outputs, and of the inputs that are not connected
#define MDL_SET_DEFAULT_PORT_DIMENSION_INFO
static void mdlSetDefaultPortDimensionInfo(SimStruct *S)
{
    Dimensions dim;
    int port, m_out;
    if (!getDimensions(S, &dim))
        return;
    m_out = dim.m > 0 ? dim.m : 1;
    if (ssGetInputPortWidth(S, IN_A) == DYNAMICALLY_SIZED)
        ssSetInputPortMatrixDimensions(S, IN_A, 1, dim.n);
    for (port = IN_LB; port <= IN_SENSE; port++)
        if (ssGetInputPortWidth(S, port) == DYNAMICALLY_SIZED)
            ssSetInputPortVectorDimension(S, port, m_out);
    if (ssGetOutputPortWidth(S, OUT_X) == DYNAMICALLY_SIZED)
        ssSetOutputPortVectorDimension(S, OUT_X, dim.n);
    if (ssGetOutputPortWidth(S, OUT_LAMBDA) == DYNAMICALLY_SIZED)
        ssSetOutputPortVectorDimension(S, OUT_LAMBDA, m_out);
}
#endif

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
    Dimensions dim;
    DAQPSfunc *s;
    int *bp, nh, m, ms, i;

    ssSetPWorkValue(S, 0, NULL);
    if (!getDimensions(S, &dim))
        return;
    nh = getBreakPoints(S, &bp);

    // (A break point is the 1-based index of the last constraint in a level,
    // which is the 0-based index where the next level starts, as in DAQP)
    if (nh > 0 && bp[nh - 1] != dim.m)
    {
        free(bp);
        ssSetErrorStatus(S, "DAQP: the last element of break_points must equal the number of constraints.");
        return;
    }

    // Without constraints, an infinite bound on x(1) is used
    m = dim.m > 0 ? dim.m : 1;
    ms = dim.m > 0 ? dim.ms : 1;
    s = daqp_sfunc_new(dim.n, m, ms, bp, nh,
                       (int)getScalarParam(S, P_WARMSTART, DAQP_SFUNC_COLD),
                       (int)getScalarParam(S, P_MAXITER, 1000));
    free(bp);
    if (s == NULL)
    {
        ssSetErrorStatus(S, "DAQP: memory allocation failed.");
        return;
    }
    for (i = 0; i < m; i++)
    {
        s->qp.bupper[i] = DAQP_INF;
        s->qp.blower[i] = -DAQP_INF;
    }
    ssSetPWorkValue(S, 0, s);
}

/* ========================== Outputs ========================== */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    DAQPSfunc *s = (DAQPSfunc *)ssGetPWorkValue(S, 0);
    DAQPProblem *qp;
    int i, j, n, m, ms, mg;
    (void)tid;
    if (s == NULL)
        return;
    qp = &s->qp;
    n = qp->n;
    m = qp->m;
    ms = qp->ms;
    mg = m - ms;

    // Simulink matrices are column-major, DAQP expects row-major
    // (H and A are only read at the first solve if they are constant)
    if (s->warm_start != DAQP_SFUNC_WARM_FIXED || !s->is_setup)
    {
        const real_T *H = ssGetInputPortRealSignal(S, IN_H);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                qp->H[i * n + j] = (c_float)H[j * n + i];
        if (mg > 0)
        {
            const real_T *A = ssGetInputPortRealSignal(S, IN_A);
            for (i = 0; i < mg; i++)
                for (j = 0; j < n; j++)
                    qp->A[i * n + j] = (c_float)A[j * mg + i];
        }
    }
    {
        const real_T *f = ssGetInputPortRealSignal(S, IN_F);
        for (i = 0; i < n; i++)
            qp->f[i] = (c_float)f[i];
    }
    // (Unconnected bounds keep their infinite values from mdlStart)
    if (ssGetInputPortConnected(S, IN_LB))
    {
        const real_T *lb = ssGetInputPortRealSignal(S, IN_LB);
        for (i = 0; i < m; i++)
            qp->blower[i] = (c_float)lb[i];
    }
    if (ssGetInputPortConnected(S, IN_UB))
    {
        const real_T *ub = ssGetInputPortRealSignal(S, IN_UB);
        for (i = 0; i < m; i++)
            qp->bupper[i] = (c_float)ub[i];
    }
    if (ssGetInputPortConnected(S, IN_SENSE))
    {
        const real_T *sense = ssGetInputPortRealSignal(S, IN_SENSE);
        for (i = 0; i < m; i++)
            qp->sense[i] = (int)sense[i];
    }

    daqp_sfunc_solve(s);

    {
        real_T *x = ssGetOutputPortRealSignal(S, OUT_X);
        real_T *lambda = ssGetOutputPortRealSignal(S, OUT_LAMBDA);
        for (i = 0; i < n; i++)
            x[i] = (real_T)s->res.x[i];
        for (i = 0; i < m; i++)
            lambda[i] = (real_T)s->res.lam[i];
        ssGetOutputPortRealSignal(S, OUT_FVAL)[0] = (real_T)s->res.fval;
        ssGetOutputPortRealSignal(S, OUT_FLAG)[0] = (real_T)s->res.exitflag;
        ssGetOutputPortRealSignal(S, OUT_ITER)[0] = (real_T)s->res.iter;
    }
}

static void mdlTerminate(SimStruct *S)
{
    daqp_sfunc_free((DAQPSfunc *)ssGetPWorkValue(S, 0));
    ssSetPWorkValue(S, 0, NULL);
}

/*=============================*
 * Required S-function trailer *
 *=============================*/
#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */
#include "simulink.c"  /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif
