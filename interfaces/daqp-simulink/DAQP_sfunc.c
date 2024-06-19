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

#define S_FUNCTION_NAME  DAQP_sfunc
#define S_FUNCTION_LEVEL 2

// #define __DEBUG__ // uncomment for debugging

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


#include <stdlib.h> // For malloc und free
#include <math.h> // For isnan, sqrt


#define INFINITY 1.0e20
static void removeNaN(real_T *x, int_T n) {
    // Replace NaN values with near infinity
    int_T i;
    for (i = 0; i < n; i++) {
        if (isnan(x[i])) {
            x[i] = INFINITY;
        }
    }
}


static void transposeMatrix(real_T *x_in, c_float *x, int_T n, int_T m) {
    // Convert Matrix column-first-order to row-first-order
    // n = number of columns
    // m = number of rows

    // check if x could be successfully allocated
    if (x == NULL) {
        return;
    }

    int_T i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            x[j*n + i] = (c_float) x_in[i*m + j];
        }
    }
}

static void convertFortranToC(InputRealPtrsType x_in, real_T* x, int nV, int nC )
{
	int i,j;

	for ( i=0; i<nC; ++i )
		for ( j=0; j<nV; ++j )
			x[i*nV + j] = x_in[i][j];
}

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
    ssSetInputPortMatrixDimensions(S, 0, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
    ssSetInputPortDirectFeedThrough(S, 0, 1);

    ssSetInputPortWidth(S, 1, DYNAMICALLY_SIZED ); /* g */
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    
    ssSetInputPortMatrixDimensions(S, 2, DYNAMICALLY_SIZED, DYNAMICALLY_SIZED);
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
    int_T maxIter = (int_T)((mxGetPr(in_maxIter))[0]);    // Maximum number of iterations

    if (maxIter <= 0) {
        mexErrMsgTxt( "ERROR (DAQP): Maximum number of iterations must be greater than zero." );
    }
    
    int_T dim_H = ssGetInputPortNumDimensions(S,0);  
    int_T* size_H = ssGetInputPortDimensions(S, 0);
    int_T size_g = ssGetInputPortWidth(S, 1);
    int_T dim_A = ssGetInputPortNumDimensions(S,2);
    int_T* size_A = ssGetInputPortDimensions(S, 2);
    int_T size_lb = ssGetInputPortWidth(S, 3);
    int_T size_ub = ssGetInputPortWidth(S, 4);
    int_T size_sense = ssGetInputPortWidth(S, 5);


    // ------------------------- Check for sizes -------------------------
    int_T n = size_H[0]; // Number of decision variables
    int_T mg = size_A[0] ; // Number of matrix constraints (lb < Ax < ub)
    int_T m = size_sense; // Number of constraints
    int_T ms = m - mg; // Number of simple constraints

    // Check for empty matrices
    if (dim_H != 2) mexErrMsgTxt( "ERROR (DAQP): Quadratic term of the objective function H must be a matrix." );
    if (dim_A != 2) mexErrMsgTxt( "ERROR (DAQP): Matrix of the linear constraints A must be a 2D matrix." );
    if (size_H[0] == 0 || size_H[1] == 0) mexErrMsgTxt( "ERROR (DAQP): Quadratic term of the objective function H is empty." );
    if (size_A[0] == 0 || size_A[1] == 0) mexErrMsgTxt( "ERROR (DAQP): Matrix of the linear constraints A is empty." );
    if (size_g == 0) mexErrMsgTxt( "ERROR (DAQP): Linear term of the objective function g is empty." );
    if (size_lb == 0) mexErrMsgTxt( "ERROR (DAQP): Lower bounds b_lower is empty." );
    if (size_ub == 0) mexErrMsgTxt( "ERROR (DAQP): Upper bounds b_upper is empty." );
    if (size_sense == 0) mexErrMsgTxt( "ERROR (DAQP): Sense of the constraints is empty." );

    // Check for dimension mismatch
    if (size_H[0] != size_H[1]) mexErrMsgTxt( "ERROR (DAQP): Quadratic term H must be a square matrix." );
    if (size_g != n) mexErrMsgTxt( "ERROR (DAQP): Dimension mismatch between the linear term g of the objective function and the number of decision variables." );
    if (size_lb != m) mexErrMsgTxt( "ERROR (DAQP): Dimension mismatch in the lower bounds b_lower or sense." );
    if (size_ub != m) mexErrMsgTxt( "ERROR (DAQP): Dimension mismatch in the upper bounds b_upper or sense." );
    if (size_A[1] != n) mexErrMsgTxt( "ERROR (DAQP): Dimension 2 of the matrix of the linear constraints A must be equal to the number of decision variables." );

    // Check for invalid sizes
    if (mg < 1) mexErrMsgTxt( "ERROR (DAQP): Number of general constraints must be greater than zero." );
    if (ms > n) mexErrMsgTxt( "ERROR (DAQP): Number of simple constraints must be less than or equal to the number of decision variables." );

    #if defined(__DEBUG__)
    mexPrintf("mdlStart finished.\n");
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
static void mdlOutputs(SimStruct *S, int_T tid) {
    // ----------------------------- PARAMS -----------------------------
    // Get number of maximum working set iterations
	const mxArray* in_maxIter = ssGetSFcnParam(S, 0);
	int_T maxIter = (int_T)((mxGetPr(in_maxIter))[0]);	// Maximum number of iterations
    

    // ----------------------------- INPUTS -----------------------------
    InputRealPtrsType in_H, in_g, in_A, in_lb, in_ub, in_sense;
    in_H     = ssGetInputPortRealSignalPtrs(S, 0); // Quadratic term of the objective function
	in_g     = ssGetInputPortRealSignalPtrs(S, 1); // Linear term of the objective function
	in_A     = ssGetInputPortRealSignalPtrs(S, 2); // Matrix of the linear constraints
	in_lb    = ssGetInputPortRealSignalPtrs(S, 3); // Lower bounds
	in_ub    = ssGetInputPortRealSignalPtrs(S, 4); // Upper bounds
	in_sense = ssGetInputPortRealSignalPtrs(S, 5); // Sense of the constraints
    
    int_T* size_H = ssGetInputPortDimensions(S, 0); // Size of the quadratic term of the objective function
    int_T size_g = ssGetInputPortWidth(S, 1); // Size of the linear term of the objective function
    int_T* size_A = ssGetInputPortDimensions(S, 2); // Size of the matrix of the linear constraints
    int_T size_lb = ssGetInputPortWidth(S, 3); // Size of the lower bounds
    int_T size_ub = ssGetInputPortWidth(S, 4); // Size of the upper bounds
    int_T size_sense = ssGetInputPortWidth(S, 5); // Size of the sense of the constraints
    

    // ------------------------- sizes ------------------------------
    int_T n = size_g; // Number of decision variables
    int_T mg = size_A[0]; // Number of matrix constraints (lb < Ax < ub)
    int_T m = size_sense; // Number of constraints
    int_T ms = m - mg; // Number of simple constraints
    

    // ----------------------- Create  Work -------------------------
    c_float *x_opt  = (c_float *) calloc( n,          sizeof(c_float) );
    c_float *lambda = (c_float *) calloc( m,         sizeof(c_float) );
    c_float *H      = (c_float *) calloc( n*n,        sizeof(c_float) );
	c_float *g      = (c_float *) calloc( n,          sizeof(c_float) );
	c_float *A      = (c_float *) calloc( mg*n,       sizeof(c_float) );
	c_float *lb     = (c_float *) calloc( m,          sizeof(c_float) );
	c_float *ub     = (c_float *) calloc( m,          sizeof(c_float) );
	int_T *sense    = (int_T *)   calloc( m,          sizeof(int_T) );


    // --------------------- CHECK FOR NULL ------------------------
    if ( x_opt == NULL || lambda == NULL || H == NULL || g == NULL || A == NULL || lb == NULL || ub == NULL || sense == NULL) {
        mexErrMsgTxt( "ERROR (DAQP): Memory allocation failed." );
    }
    

    // --------------------- Write to Work ----------------------
    int_T i;
    transposeMatrix((real_T *) *in_H, H, n, n); // Convert Matrix column-first-order to row-first-order

    for ( i=0; i<size_g; ++i )
        g[i] = (c_float)(*in_g)[i];

    transposeMatrix((real_T *) *in_A, A, mg, n); // Convert Matrix column-first-order to row-first-order

	for ( i=0; i<size_lb; ++i )
		lb[i] = (c_float)(*in_lb)[i];
    
	for ( i=0; i<size_ub; ++i )
		ub[i] = (c_float)(*in_ub)[i];
		
    for ( i=0; i<size_sense; ++i )
		sense[i] = (int_T)(*in_sense)[i];


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
    settings.iter_limit = (int_T) maxIter;


    // ----------------------------- DAQP -----------------------------
     DAQPProblem qp = {n,m,ms,H,g,A,ub,lb,sense}; 
     daqp_quadprog(&result,&qp,&settings);


    // ------------------------- WRITE OUTPUTS -------------------------
    for ( i=0; i<n; ++i )
		out_x_opt[i] = (real_T)(x_opt[i]);
    
    for ( i=0; i<m; ++i )
		out_lambda[i] = (real_T)(lambda[i]);
    
    out_fval[0] = (real_T)(result.fval);
    out_flag[0] = (real_T)(result.exitflag);
    out_iteration[0] = (real_T)(result.iter);
    

    // -------------------------- REMOVE NaN -------------------------
    removeNaN(out_x_opt, (int_T) size_g);
    removeNaN(out_lambda, mg);
    removeNaN(out_fval, 1);
    removeNaN(out_flag, 1);
    removeNaN(out_iteration, 1);


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
    // int_T i;
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
        
    // if (ssGetOutputPortWidth(S, 0) == DYNAMICALLY_SIZED) {
    //     ssSetOutputPortMatrixDimensions(S, 0, size_g, 1 );
    // }
    // if (ssGetOutputPortWidth(S, 1) == DYNAMICALLY_SIZED) {
    //     ssSetOutputPortMatrixDimensions(S, 1, mg, 1 );
    // }
    if (ssGetOutputPortWidth(S, 0) == DYNAMICALLY_SIZED) {
        ssSetOutputPortMatrixDimensions(S, 0, size_g, 1 );
    }
    if (ssGetOutputPortWidth(S, 1) == DYNAMICALLY_SIZED) {
        ssSetOutputPortMatrixDimensions(S, 1, size_sense, 1 );
    }

    #if defined(__DEBUG__)
    mexPrintf("size_H: %d\n", size_H);
    mexPrintf("size_g: %d\n", size_g);
    mexPrintf("size_A: %d\n", size_A);
    mexPrintf("size_lb: %d\n", size_lb);
    mexPrintf("size_ub: %d\n", size_ub);
    mexPrintf("size_sense: %d\n", size_sense);
    mexPrintf("mg: %d\n", mg);
    #endif

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