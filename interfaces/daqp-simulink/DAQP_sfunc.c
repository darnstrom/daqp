/*
 * Author: Christopher Schulte, RWTH Aachen University, IRT
 *         christopher.schulte@rwth-aachen.de
 *        https://www.irt.rwth-aachen.de/
 * 
 * This file is part of the DAQP project.
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
 * 1. n: Number of decision variables
 * 2. mg: Number of general constraints (mg = m - ms)
 * 3. maxIter: Maximum number of iterations (currently not used)
 * 
 * OUTPUTS:
 * 1. x_opt: Optimal solution of the decision variables [n x 1]
 * 2. lambda: Lagrange multipliers [mg x 1]
 * 3. fval: Value of the objective function at the optimal solution [1 x 1]
 * 4. flag: Exit flag [1 x 1]
 * 5. iteration: Number of iterations [1 x 1]
 */

#define S_FUNCTION_NAME  DAQP_sfunc
#define S_FUNCTION_LEVEL 2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"

#include "types.h"
#include "daqp.h"
#include "api.h"
#include <stdlib.h> // For malloc und free


/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
//     ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);
//     ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
    
    // --------------------------- PARAMETERS -----------------------------
    ssSetNumSFcnParams(S, 3);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }
    const mxArray* in_dimU = ssGetSFcnParam(S, 0);
	int_T dimU = (int_T)((mxGetPr(in_dimU))[0]);	
    
    const mxArray* in_dimA = ssGetSFcnParam(S, 1);
	int_T dimA = (int_T)((mxGetPr(in_dimA))[0]);
    
   
    // ----------------------------- INPUTS -----------------------------
    if (!ssSetNumInputPorts(S, 6)) return;
    ssSetInputPortWidth(S, 0, dimU*dimU ); /* H */
    ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    ssSetInputPortWidth(S, 1, dimU ); /* g */
    ssSetInputPortRequiredContiguous(S, 1, true); /*direct input signal access*/
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    
    ssSetInputPortWidth(S, 2, dimU*dimA ); /* A */
    ssSetInputPortRequiredContiguous(S, 2, true); /*direct input signal access*/
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    
    ssSetInputPortWidth(S, 3, DYNAMICALLY_SIZED ); /* lb */
    ssSetInputPortRequiredContiguous(S, 3, true); /*direct input signal access*/
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    
    ssSetInputPortWidth(S, 4, DYNAMICALLY_SIZED ); /* ub */
    ssSetInputPortRequiredContiguous(S, 4, true); /*direct input signal access*/
    ssSetInputPortDirectFeedThrough(S, 4, 1);
    
    ssSetInputPortWidth(S, 5, DYNAMICALLY_SIZED ); /* sense */
    ssSetInputPortRequiredContiguous(S, 5, true); /*direct input signal access*/
    ssSetInputPortDirectFeedThrough(S, 5, 1);
    	
    
    // ----------------------------- WORK -------------------------------
    ssSetNumPWork(S, 2);
    
    
    // ----------------------------- OUTPUTS -----------------------------
    if (!ssSetNumOutputPorts(S, 5)) return;
    ssSetOutputPortWidth(S, 0, dimU); /* x_opt */
    ssSetOutputPortWidth(S, 1, dimA); /* lambda */
    ssSetOutputPortWidth(S, 2, 1); /* fval */
    ssSetOutputPortWidth(S, 3, 1); /* flag */
    ssSetOutputPortWidth(S, 4, 1); /* iteration */
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
    const mxArray* in_dimU = ssGetSFcnParam(S, 0);
	int_T dimU = (int_T)((mxGetPr(in_dimU))[0]);	
    
    const mxArray* in_dimA = ssGetSFcnParam(S, 1);
	int_T dimA = (int_T)((mxGetPr(in_dimA))[0]);
    
    ssGetPWork(S)[0] = (void *) calloc( dimU, sizeof(real_T) );
    ssGetPWork(S)[1] = (void *) calloc( dimA, sizeof(real_T) );
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
static void mdlOutputs(SimStruct *S, int_T tid) {
    // ----------------------------- PARAMS -----------------------------
    // Get number of maximum working set iterations, currently not used
	const mxArray* in_maxIter = ssGetSFcnParam(S, 2);
	int maxIter = (int)((mxGetPr(in_maxIter))[0]);	
    
    // ----------------------------- INPUTS -----------------------------
    real_T *H = (real_T*) ssGetInputPortSignal(S,0);
    real_T *g = (real_T*) ssGetInputPortSignal(S,1);
    real_T *A = (real_T*) ssGetInputPortSignal(S,2);
    real_T *lb = (real_T*) ssGetInputPortSignal(S,3);
    real_T *ub = (real_T*) ssGetInputPortSignal(S,4);
    int_T *sense = (int_T*) ssGetInputPortSignal(S,5);
    
    int_T size_H = ssGetInputPortWidth(S, 0);
    int_T size_g = ssGetInputPortWidth(S, 1);
    int_T size_A = ssGetInputPortWidth(S, 2);
    int_T size_lb = ssGetInputPortWidth(S, 3);
    int_T size_ub = ssGetInputPortWidth(S, 4);
    int_T size_sense = ssGetInputPortWidth(S, 5);
    
    // ----------------------------- OUTPUTS -----------------------------
    real_T *out_x_opt   = ssGetOutputPortRealSignal(S, 0);
	real_T *out_lambda = ssGetOutputPortRealSignal(S, 1);
	real_T *out_fval = ssGetOutputPortRealSignal(S, 2);
	real_T *out_flag = ssGetOutputPortRealSignal(S, 3);
	real_T *out_iteration = ssGetOutputPortRealSignal(S, 4);

    int_T n = size_g; // Number of decision variables
    int_T mg = (int_T)((real_T) size_A / (real_T) size_g); // Number of matrix constraints (lb < Ax < ub)
    int_T m = size_lb; // Number of constraints
    int_T ms= m - mg; // Number of simple constraints
    

    // ----------------------------- WORK -----------------------------
    real_T *x_opt, *lambda;
    x_opt = (real_T *) ssGetPWork(S)[0];
	lambda = (real_T *) ssGetPWork(S)[1];
    
    DAQPResult result;
    result.x = x_opt; // primal variable
    result.lam = lambda; // dual variable

    // ----------------------------- DAQP -----------------------------
    DAQPProblem qp = {n,m,ms,H,g,A,ub,lb,sense};
    
    daqp_quadprog(&result,&qp,NULL);

    // ------------------------- WRITE OUTPUTS -------------------------
    int_T i;
    for ( i=0; i<size_g; ++i )
		out_x_opt[i] = (real_T)(x_opt[i]);
    
    for ( i=0; i<mg; ++i )
		out_lambda[i] = (real_T)(lambda[i]);
    
    out_fval[0] = (real_T)(result.fval);
    out_flag[0] = (real_T)(result.exitflag);
    out_iteration[0] = (real_T)(result.iter);
}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
    // Free memory
    int i;
	for ( i=0; i<2; ++i )
	{
		if ( ssGetPWork(S)[i] != 0 )
			free( ssGetPWork(S)[i] );
	}
}


#if defined(MATLAB_MEX_FILE)

/* Function: mdlSetInputPortDimensionInfo + mdlSetOutputPortDimensionInfo
 * Abstract:
 *    Needed for dynamic dimension
 */
#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO

static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
	if ( !ssSetInputPortDimensionInfo(S, port, dimsInfo) )
		return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
	if ( !ssSetOutputPortDimensionInfo(S, port, dimsInfo) )
		return;
}

#endif

/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif