// HaiCS (Hamiltonians in Computational Statistics)
//            Developed in the group of
//         Professor Elena Akhmatskaya by:
//     Elena Akhmatskaya, Mario Fernández-Pendás,
//    Hristo Inouzhe, Felix Müller, Lorenzo Nagar, 
//       Jorge Pérez Heredia, Tijana Radivojević, 
//           María Xosé Rodríguez-Álvarez.
//   For references to this package, please cite:
//     * Radivojevic, T. Enhancing Sampling in
//     Computational Statistics Using Modified
//     Hamiltonians. PhD Thesis, UPV/EHU (2016)
//  * Radivojevic, T., Akhmatskaya, E. Modified
// Hamiltonian Monte Carlo for Bayesian inference.
//          Stat Comput 30, 377–404 (2020).
//
// This version provides the code for running SIR/SEIR-like
// probabilistic models (for Incidence data) with a time-dependent
// transmission rate under a Bayesian approach to infer the
// epidemiological parameters using HaiCS HMC-based samplers.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>


#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */
//#include <legendre/legendre_polynomial.h>
/* Accessor macros */

#include "Globals.h"
#include "utils.h"
#include "sir_seir_common_functions.h"
#include "read_input.h"

/*
 *--------------------------------------------------------------------
 * SPLINE Functions
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y).
 */

int mindex(int i, int j, int nrow) {
    return (j - 1) * nrow + i;
}

double bs(int nknots, int nspline, int updegree, double x, double * knots) {
    double y, y1, y2, temp1, temp2;
     if (updegree == 1) {
        if ((x >= knots[nspline - 1]) && (x < knots[nspline]))
            y = 1.0;
        else
            y = 0.0;
     }
     else {
        temp1 = 0.0;
        if ((knots[nspline + updegree - 2] - knots[nspline - 1]) > 0)
            temp1 = (x - knots[nspline - 1]) / (knots[nspline + updegree - 2] - knots[nspline - 1]);
        temp2 = 0.0;
        if ((knots[nspline + updegree - 1] - knots[nspline]) > 0)
            temp2 = (knots[nspline + updegree - 1] - x) / (knots[nspline + updegree - 1] - knots[nspline]);
        y1 = bs(nknots, nspline, updegree - 1, x, knots);
        y2 = bs(nknots, nspline + 1, updegree - 1, x, knots);
        y =  temp1 * y1 + temp2 * y2;
    }
    return y;
}

void splinebasis(int *d, int *n, int *m, double * x, double * knots, double * basis) {
    int mm = *m, dd = *d, nn = *n;
    int k = mm - dd - 1, i , j, ir, jr;
    for (i = 0; i < nn; i++) {
        ir = i + 1;
        if (x[i] == knots[mm - 1]) {
           basis [mindex (ir, k, nn) - 1] =  1.0;
            for (j = 0; j  < (k - 1);  j++) {
                jr = j + 1;
                basis [mindex (ir, jr, nn) - 1] = 0.0;
            }
        } else {
            for (j = 0; j < k ; j++) {
                jr = j + 1;
                basis[mindex (ir, jr, nn) - 1] = bs(mm, jr, dd + 1, x[i], knots);
            }
        }
    }
}

double splinefunction(double * beta, int *d, int *n, int *m, double * x, double * knots, double * basis, int *nf){
    double splineValue = 0.0;
    int ii, nnf = *nf;

    splinebasis(d,  n, m, x, knots, basis);

    for (ii = 0; ii < nnf; ii++){
        splineValue += beta[ii] * basis[ii];
    }
    return splineValue;
}

/*
 *--------------------------------------------------------------------
 * General Functions
 *--------------------------------------------------------------------
 */

double underReport(int t, double u0, double u1 , int t0, int t1){
    double res = 0.0;
    if (t <= t0) {
        res = u0;
    } else {
	if (t >= t1){
	    res = u1;
	} else {
	    res = u0 + ((u1 - u0) / (t1 - t0)) * (t - t0);
	}
    }
    return res;
}

double vectors_dot_prod(double *x, double *y, int n){
    double res = 0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

void matrix_vector_mult(double **mat, double *vec, double *result, int rows, int cols){ // in matrix form: result = mat * vec;
    int i;
    for (i = 0; i < rows; i++)
    {
        result[i] = vectors_dot_prod(mat[i], vec, cols);
    }
}

void matrix_matrix_mult(double **matrix_left, double **matrix_right, double **result_matrix, int rows, int cols, int prod_dim){ // in matrix form: result = mat * vec;
    int i,j,k;
    double result;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++){
            result = 0;
            for(k = 0; k <  prod_dim; k++){
                result += matrix_left[i][k] * matrix_right[k][j];
            }
            result_matrix[i][j] = result;
        }

    }
}


/*--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Process and verify arguments to cvsfwddenx.
 */

void ProcessArgs(char *argv[], booleantype *sensi, int *sensi_meth, booleantype *err_con)
{
  *sensi = SUNFALSE;
  *sensi_meth = -1;
  *err_con = SUNFALSE;

//  if (argc < 2) WrongArgs(argv[0]);

  if (strcmp(argv[0],"-nosensi") == 0)
    *sensi = SUNFALSE;
  else if (strcmp(argv[0],"-sensi") == 0)
    *sensi = SUNTRUE;
  else
    WrongArgs("SIR_CVODES_synth");

  if (*sensi) {

    if (strcmp(argv[1],"sim") == 0)
      *sensi_meth = CV_SIMULTANEOUS;
    else if (strcmp(argv[1],"stg") == 0)
      *sensi_meth = CV_STAGGERED;
    else if (strcmp(argv[1],"stg1") == 0)
      *sensi_meth = CV_STAGGERED1;
    else
      WrongArgs("SIR_CVODES_synth");

    if (strcmp(argv[2],"t") == 0)
      *err_con = SUNTRUE;
    else if (strcmp(argv[2],"f") == 0)
      *err_con = SUNFALSE;
    else
      WrongArgs("SIR_CVODES_synth");
  }

}

void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con]\n",name);
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = t or f\n");

    exit(0);
}

int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
