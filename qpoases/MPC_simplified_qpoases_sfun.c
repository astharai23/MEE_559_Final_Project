#define S_FUNCTION_NAME  MPC_simplified_qpoases_sfun
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "MPC_simplified_qpoases.h"


static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return;
    }

    if (!ssSetNumInputPorts(S, 7)) return;

    
    /* Input port 0 */
    ssSetInputPortWidth(S, 0, 152);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortRequiredContiguous(S, 0, 1);
    
    /* Input port 1 */
    ssSetInputPortWidth(S, 1, 2);
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortRequiredContiguous(S, 1, 1);
    
    /* Input port 2 */
    ssSetInputPortWidth(S, 2, 1);
    ssSetInputPortDataType(S, 2, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortRequiredContiguous(S, 2, 1);
    
    /* Input port 3 */
    ssSetInputPortWidth(S, 3, 1);
    ssSetInputPortDataType(S, 3, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    ssSetInputPortRequiredContiguous(S, 3, 1);
    
    /* Input port 4 */
    ssSetInputPortWidth(S, 4, 1);
    ssSetInputPortDataType(S, 4, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 4, 1);
    ssSetInputPortRequiredContiguous(S, 4, 1);
    
    /* Input port 5 */
    ssSetInputPortWidth(S, 5, 1);
    ssSetInputPortDataType(S, 5, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 5, 1);
    ssSetInputPortRequiredContiguous(S, 5, 1);
    
    /* Input port 6 */
    ssSetInputPortWidth(S, 6, 1);
    ssSetInputPortDataType(S, 6, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 6, 1);
    ssSetInputPortRequiredContiguous(S, 6, 1);
    

    if (!ssSetNumOutputPorts(S, 3)) return;

    
    /* Output port 0 */
    ssSetOutputPortWidth(S, 0, 152);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    
    /* Output port 1 */
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortDataType(S, 1, SS_INT32);
    
    /* Output port 2 */
    ssSetOutputPortWidth(S, 2, 1);
    ssSetOutputPortDataType(S, 2, SS_DOUBLE);
    

    /* Single sample time */
    ssSetNumSampleTimes(S, 1);

    ssSetArrayLayoutForCodeGen(S, SS_COLUMN_MAJOR);
    ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);

    ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);
   
    ssSetOptions(S,
                 SS_OPTION_WORKS_WITH_CODE_REUSE |
                 SS_OPTION_EXCEPTION_FREE_CODE |
                 SS_OPTION_USE_TLC_WITH_ACCELERATOR);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, 0.02);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S); 
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    /* get input reference */
    double *z_op = (double *) ssGetInputPortSignal(S, 0);
    double *x0 = (double *) ssGetInputPortSignal(S, 1);
    double *pd = (double *) ssGetInputPortSignal(S, 2);
    double *Q11 = (double *) ssGetInputPortSignal(S, 3);
    double *Q22 = (double *) ssGetInputPortSignal(S, 4);
    double *R = (double *) ssGetInputPortSignal(S, 5);
    double *umax = (double *) ssGetInputPortSignal(S, 6);
    

    /* get output reference */
    double *z_opt = (double *) ssGetOutputPortSignal(S, 0);
    int *exit_flag = (int *) ssGetOutputPortSignal(S, 1);
    double *u0 = (double *) ssGetOutputPortSignal(S, 2);
    

    
    /* call function */
    MPC_simplified_call(z_op, x0, pd, Q11, Q22, R, umax, z_opt, exit_flag, u0);
    
}

static void mdlTerminate(SimStruct *S)
{

}


#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif