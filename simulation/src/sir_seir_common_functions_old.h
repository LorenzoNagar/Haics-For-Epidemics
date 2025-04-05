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
#include "Definitions.h"

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */

//#define NOUT  100     /* number of output times */ //substituted by g_inpset->num_out_cvodes

//#define NSpB  12 //substituted by g_inpset->num_basis_spline
//#define NP    2 + NSpB           /* number of problem parameters */
//#define NS    2 + NSpB            /* number of sensitivities computed */

#ifndef SIR_SEIR_COMMON_FUNCTIONS_H
#define SIR_SEIR_COMMON_FUNCTIONS_H

/* Type : UserData */

/* typedef struct {
  realtype p[];
  int spline_degree;
  int number_extended_knots;
  int number_basis_elements;
  double *knots;
  double *basis;
} *UserData; */

/* Prototypes of private functions */

void ProcessArgs(char *argv[], booleantype *sensi, int *sensi_meth, booleantype *err_con);
void WrongArgs(char *name);
//static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
//static void PrintOutputS(N_Vector *uS);
//static void PrintFinalStats(void *cvode_mem, booleantype sensi);
 
int check_retval(void *returnvalue, const char *funcname, int opt);

int mindex(int i, int j, int nrow);

double bs(int nknots, int nspline, int updegree, double x, double * knots);

void splinebasis(int *d, int *n, int *m, double * x, double * knots, double * basis);

double splinefunction(double * beta, int *d, int *n, int *m, double * x, double * knots, double * basis, int *nf);

double vectors_dot_prod(double *x, double *y, int n);

double underReport(int t, double u0, double u1 , int t0, int t1);

void matrix_vector_mult(double **mat, double *vec, double *result, int rows, int cols);

void matrix_matrix_mult(double **matrix_left, double **matrix_right, double **result_matrix, int rows, int cols, int prod_dim);

#endif
