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
#include "mex.h"
#include <math.h>

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
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
    
    // --------------------------- PARAMETERS -----------------------------
    ssSetNumSFcnParams(S, 1);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /* Return if number of expected != number of actual parameters */
        return;
    }
   
    // ----------------------------- INPUTS -----------------------------
    if (!ssSetNumInputPorts(S, 6)) return;
    ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED ); /* H */
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    ssSetInputPortWidth(S, 1, DYNAMICALLY_SIZED ); /* g */
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    
    ssSetInputPortWidth(S, 2, DYNAMICALLY_SIZED ); /* A */
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    
    ssSetInputPortWidth(S, 3, DYNAMICALLY_SIZED ); /* lb */
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    
    ssSetInputPortWidth(S, 4, DYNAMICALLY_SIZED ); /* ub */
    ssSetInputPortDirectFeedThrough(S, 4, 1);
    
    ssSetInputPortWidth(S, 5, DYNAMICALLY_SIZED ); /* sense */
    ssSetInputPortDirectFeedThrough(S, 5, 1);
    	
    
    // ----------------------------- WORK -------------------------------
    ssSetNumPWork(S, 0);
    
    
    // ----------------------------- OUTPUTS -----------------------------
    if (!ssSetNumOutputPorts(S, 5)) return;
    ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED); /* x_opt */
    ssSetOutputPortWidth(S, 1, DYNAMICALLY_SIZED); /* lambda */
    ssSetOutputPortWidth(S, 2, 1); /* fval */
    ssSetOutputPortWidth(S, 3, 1); /* flag */
    ssSetOutputPortWidth(S, 4, 1); /* iteration */
}


/* Function: mdlStart =========================================================
 * Abstract:
 *    This function is called once at start of the model execution. If you
 *    have states that should be initialized once, this is the place to do it.
 */
#define MDL_START
static void mdlStart(SimStruct *S)
{
    const mxArray* in_maxIter = ssGetSFcnParam(S, 0);
    int maxIter = (int)((mxGetPr(in_maxIter))[0]);    // Maximum number of iterations

    if (maxIter <= 0) {
        mexErrMsgTxt( "ERROR (DAQP): Maximum number of iterations must be greater than zero." );
    }
    
    int_T size_H = ssGetInputPortWidth(S, 0);
    int_T size_g = ssGetInputPortWidth(S, 1);
    int_T size_A = ssGetInputPortWidth(S, 2);
    int_T size_lb = ssGetInputPortWidth(S, 3);
    int_T size_ub = ssGetInputPortWidth(S, 4);
    int_T size_sense = ssGetInputPortWidth(S, 5);


    // ------------------------- Check for sizes -------------------------
    if (size_H == 0) mexErrMsgTxt( "ERROR (DAQP): Quadratic term of the objective function H is empty." );
    if (size_g == 0) mexErrMsgTxt( "ERROR (DAQP): Linear term of the objective function g is empty." );
    if (size_A == 0) mexErrMsgTxt( "ERROR (DAQP): Matrix of the linear constraints A is empty." );
    if (size_lb == 0) mexErrMsgTxt( "ERROR (DAQP): Lower bounds b_lower are empty." );
    if (size_ub == 0) mexErrMsgTxt( "ERROR (DAQP): Upper bounds b_upper are empty." );
    if (size_sense == 0) mexErrMsgTxt( "ERROR (DAQP): Sense of the constraints is empty." );
    if (sqrt((real_T)size_H) != (real_T) size_g) mexErrMsgTxt( "ERROR (DAQP): Dimension mismatch between the quadratic term H and the linear term g of the objective function." );
    if (size_A % size_g != 0) mexErrMsgTxt( "ERROR (DAQP): Dimension mismatch in the matrix of the linear constraints A." );
    if (size_lb != size_sense) mexErrMsgTxt( "ERROR (DAQP): Dimension mismatch in the lower bounds b_lower or sense." );
    if (size_ub != size_sense) mexErrMsgTxt( "ERROR (DAQP): Dimension mismatch in the upper bounds b_upper or sense." );
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
    // Get number of maximum working set iterations
	const mxArray* in_maxIter = ssGetSFcnParam(S, 0);
	int maxIter = (int)((mxGetPr(in_maxIter))[0]);	// Maximum number of iterations
    

    // ----------------------------- INPUTS -----------------------------
    InputRealPtrsType in_H, in_g, in_A, in_lb, in_ub, in_sense;
    in_H     = ssGetInputPortRealSignalPtrs(S, 0); // Quadratic term of the objective function
	in_g     = ssGetInputPortRealSignalPtrs(S, 1); // Linear term of the objective function
	in_A     = ssGetInputPortRealSignalPtrs(S, 2); // Matrix of the linear constraints
	in_lb    = ssGetInputPortRealSignalPtrs(S, 3); // Lower bounds
	in_ub    = ssGetInputPortRealSignalPtrs(S, 4); // Upper bounds
	in_sense = ssGetInputPortRealSignalPtrs(S, 5); // Sense of the constraints
    
    int_T size_H = ssGetInputPortWidth(S, 0); // Size of the quadratic term of the objective function
    int_T size_g = ssGetInputPortWidth(S, 1); // Size of the linear term of the objective function
    int_T size_A = ssGetInputPortWidth(S, 2); // Size of the matrix of the linear constraints
    int_T size_lb = ssGetInputPortWidth(S, 3); // Size of the lower bounds
    int_T size_ub = ssGetInputPortWidth(S, 4); // Size of the upper bounds
    int_T size_sense = ssGetInputPortWidth(S, 5); // Size of the sense of the constraints
    

    // ------------------------- sizes ------------------------------
    int n = (int) size_g; // Number of decision variables
    int mg = (int)((real_T) size_A / (real_T) size_g); // Number of matrix constraints (lb < Ax < ub)
    int m = (int) size_lb; // Number of constraints
    int ms = m - mg; // Number of simple constraints
    

    // ----------------------- Create  Work -------------------------
    c_float *x_opt  = (c_float *) calloc( n, sizeof(c_float) );
    c_float *lambda = (c_float *) calloc( mg, sizeof(c_float) );
    c_float *H      = (c_float *) calloc( n*n, sizeof(c_float) );
	c_float *g      = (c_float *) calloc( n, sizeof(c_float) );
	c_float *A      = (c_float *) calloc( mg*n, sizeof(c_float) );
	c_float *lb     = (c_float *) calloc( size_sense, sizeof(c_float) );
	c_float *ub     = (c_float *) calloc( size_sense, sizeof(c_float) );
	int *sense      = (int *) calloc( size_sense, sizeof(int) );


    // --------------------- CHECK FOR NULL ------------------------
    if ( x_opt == NULL || lambda == NULL || H == NULL || g == NULL || A == NULL || lb == NULL || ub == NULL || sense == NULL) {
        mexErrMsgTxt( "ERROR (DAQP): Memory allocation failed." );
    }
    

    // --------------------- Write to Work ----------------------
    int_T i;
    if ( H != 0 && size_H == n*n) {
		for ( i=0; i<size_H; ++i )
			H[i] = (c_float)(*in_H)[i];
	}
    if ( g != 0 && size_g == n) {
		for ( i=0; i<size_g; ++i )
			g[i] = (c_float)(*in_g)[i];
	}
    if ( A != 0 && size_A == mg*n) {
		for ( i=0; i<size_A; ++i )
			A[i] = (c_float)(*in_A)[i];
	}
    if ( lb != 0 && size_lb == size_sense) {
		for ( i=0; i<size_lb; ++i )
			lb[i] = (c_float)(*in_lb)[i];
	}
    if ( ub != 0  && size_ub == size_sense) {
		for ( i=0; i<size_ub; ++i )
			ub[i] = (c_float)(*in_ub)[i];
	}
    if ( sense != 0 ) {
		for ( i=0; i<size_sense; ++i )
			sense[i] = (int)(*in_sense)[i];
	}	
    

    // ----------------------------- OUTPUTS -----------------------------
    real_T *out_x_opt   = ssGetOutputPortRealSignal(S, 0);
	real_T *out_lambda = ssGetOutputPortRealSignal(S, 1);
	real_T *out_fval = ssGetOutputPortRealSignal(S, 2);
	real_T *out_flag = ssGetOutputPortRealSignal(S, 3);
	real_T *out_iteration = ssGetOutputPortRealSignal(S, 4);


    // --------------------- Result and settings ----------------------
    DAQPResult result;
    result.x = x_opt; // primal variable
    result.lam = lambda; // dual variable
    
    DAQPSettings settings;
    daqp_default_settings(&settings);
    settings.iter_limit = (int) maxIter;


    // ----------------------------- DAQP -----------------------------
    DAQPProblem qp = {n,m,ms,H,g,A,ub,lb,sense};
    
    daqp_quadprog(&result,&qp,&settings);


    // ------------------------- WRITE OUTPUTS -------------------------
    for ( i=0; i<size_g; ++i )
		out_x_opt[i] = (real_T)(x_opt[i]);
    
    for ( i=0; i<mg; ++i )
		out_lambda[i] = (real_T)(lambda[i]);
    
    out_fval[0] = (real_T)(result.fval);
    out_flag[0] = (real_T)(result.exitflag);
    out_iteration[0] = (real_T)(result.iter);


    // -------------------------- FREE MEMORY --------------------------
    free( H );
    free( g );
    free( A );
    free( lb );
    free( ub );
    free( sense );
    free( x_opt );
    free( lambda );
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
    // int i;
	// for ( i=0; i<8; ++i ) {
	// 	if ( ssGetPWork(S)[i] != 0 )
	// 		free( ssGetPWork(S)[i] );
	// }
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
	if ( !ssSetInputPortDimensionInfo(S, port, dimsInfo) )
		return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
	if ( !ssSetOutputPortDimensionInfo(S, port, dimsInfo) )
		return;
}

static void mdlSetDefaultPortDimensionInfo(SimStruct *S) {
    int_T size_H = ssGetInputPortWidth(S, 0);
    int_T size_g = ssGetInputPortWidth(S, 1);
    int_T size_A = ssGetInputPortWidth(S, 2);
    int_T size_lb = ssGetInputPortWidth(S, 3);
    int_T size_ub = ssGetInputPortWidth(S, 4);
    int_T size_sense = ssGetInputPortWidth(S, 5);
    
    int_T mg = (int_T)((real_T) size_A / (real_T) size_g);
    
    // if (ssGetOutputPortWidth(S, 0) == DYNAMICALLY_SIZED) {
    //     ssSetOutputPortMatrixDimensions(S, 0, size_g, 1 );
    // }
    // if (ssGetOutputPortWidth(S, 1) == DYNAMICALLY_SIZED) {
    //     ssSetOutputPortMatrixDimensions(S, 1, mg, 1 );
    // }
    if (ssGetOutputPortWidth(S, 0) == DYNAMICALLY_SIZED) {
        ssSetOutputPortMatrixDimensions(S, 0, 2, 1 );
    }
    if (ssGetOutputPortWidth(S, 1) == DYNAMICALLY_SIZED) {
        ssSetOutputPortMatrixDimensions(S, 1, 1, 1 );
    }
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