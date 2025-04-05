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

#include <math.h>
#include <gsl/gsl_sf_psi.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "seir_model_synthetic.h"
#include "Globals.h"
#include "utils.h"
#include "seir_functions.h"
#include "read_input.h"
#include "sir_seir_common_functions.h"

#define M_PI 3.14159265358979323846
#define X0 0.0

#define gamma0 1e-10
#define phiInv0 1e-10
#define alpha0 1e-10

#define phiInvExp 1
//#define gammaL 0.03333333333333 
//#define gammaU 1.0

double PriorSEMIKR(){

    double t_alpha, t_gamma, t_i0, t_phi_inv;

    t_alpha =  exp(state[curr]->pos[0]) + alpha0;
    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[1]));
    // t_gamma =  exp(state[curr]->pos[1]) + gamma0;
    t_i0 =     1.0 + (num_pop - 1)/(1 + exp(-state[curr]->pos[2]));
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;

    double t_tau =   exp(state[curr]->pos[4]);
    int i, j;
    printf("------------------------------\n");
    printf("HMC proposal transformed\n");
    printf("%f %f %f %f %f ",t_alpha, t_gamma, t_i0, t_phi_inv, t_tau);
    for (i = 5; i < (5 + g_inpset->num_basis_spline); i++){
        printf("%f ", state[curr]->pos[i]);
    }
    printf("\n");
    printf("------------------------------\n");

    double spline_coordinates[g_inpset->num_basis_spline];
    for (i = 0; i < g_inpset->num_basis_spline; i++){
        spline_coordinates[i] = state[curr]->pos[i + 5];
    }

    double m[] = {g_inpset->alpha_param_1, g_inpset->gamma_param_1, g_inpset->I0_param_1, 0.1};
    double s[] = {g_inpset->alpha_param_2, g_inpset->gamma_param_2, g_inpset->I0_param_2, 75e-4};


    double logp = 0.0;
    double params[4] = {t_alpha, t_gamma, t_i0, t_phi_inv};
    double inv_gamma[] = {g_inpset->tau_param_1, g_inpset->tau_param_2};


    int nf = g_inpset->num_basis_spline;
    double ** KK = malloc((nf - 2) * sizeof(double *));
    for (i = 0; i < (nf - 2); i++){
        KK[i] = malloc(nf * sizeof(double));
    }
    double ** KKT = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++){
        KKT[i] = malloc((nf - 2) * sizeof(double));
    }
    double ** result_matrix = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++){
        result_matrix[i] = malloc(nf * sizeof(double));
    }

    KK[0][0] = 1;
    KK[0][1] = -2;
    KK[0][2] = 1;
    for(j = 3; j < nf; j++){
        KK[0][j] = 0;
    }
    for (i = 1; i < nf - 2; i++){
        for(j = 0; j < i; j++){
            KK[i][j] = 0;
        }
        KK[i][i] = 1;
        KK[i][i + 1] = -2;
        KK[i][i + 2] = 1;
        for(j = i + 3; j < nf; j++){
            KK[i][j] = 0;
        }
    }

    KKT[0][0] = 1;
    KKT[1][0] = -2;
    KKT[2][0] = 1;
    for(i = 3; i < nf; i++){
        KKT[i][0] = 0;
    }
    for (j = 1; j < nf - 2; j++){
        for(i = 0; i < j; i++){
            KKT[i][j] = 0;
        }
        KKT[j][j] = 1;
        KKT[j + 1][j] = -2;
        KKT[j + 2][j] = 1;
        for(i = j + 3; i < nf; i++){
            KKT[i][j] = 0;
        }
    }

    matrix_matrix_mult(KKT, KK, result_matrix, nf, nf, nf - 2);

    double result_mv[nf];
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);

    logp += -(params[0] - m[0])*(params[0] - m[0])/(2*s[0]*s[0]) - log(sqrt(2*M_PI)*s[0]);
    logp += -(params[1] - m[1])*(params[1] - m[1])/(2*s[1]*s[1]) - log(sqrt(2*M_PI)*s[1]);
    logp += -(params[2] - m[2])*(params[2] - m[2])/(2*s[2]*s[2]) - log(sqrt(2*M_PI)*s[2]);
    if(phiInvExp == 1){
        logp += log(g_inpset->phi_inv_param_1) - g_inpset->phi_inv_param_1 * 	params[3];
    } else {
        logp += -(params[3] - m[3])*(params[3] - m[3])/(2*s[3]*s[3]) - log(sqrt(2*M_PI)*s[3]);
    }

    logp += inv_gamma[0]*log(inv_gamma[1]) - lgamma(inv_gamma[0]) - (inv_gamma[0] + 1) * log(t_tau) - inv_gamma[1]/t_tau;
    logp += -(1/(2*t_tau*t_tau)) * vectors_dot_prod(spline_coordinates, result_mv, nf);

    /*printf("------------------------------\n");
    printf("Log-prior: %.8f \n", logp);
    printf("------------------------------\n");*/

    return logp;
}

double PriorSEMIKR_2(double *H){

    double t_alpha, t_gamma, t_i0, t_phi_inv;

    t_alpha = exp(state[curr]->pos[0] + H[0]) + alpha0;
    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-(state[curr]->pos[1] + H[1])));
    // t_gamma =  exp(state[curr]->pos[1] + H[1]) + gamma0;
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-(state[curr]->pos[2] + H[2])));
    t_phi_inv = exp(state[curr]->pos[3] + H[3]) + phiInv0;

    double t_tau =   exp(state[curr]->pos[4] + H[4]);
    int i, j;
    // /printf("------------------------------\n");
    // printf("HMC proposal transformed\n");
    // printf("%f %f %f %f %f ",t_alpha, t_gamma, t_i0, t_phi_inv, t_tau);
    // for (i = 5; i < (5 + g_inpset->num_basis_spline); i++){
    //     printf("%f ", state[curr]->pos[i]);
    // }
    // printf("\n");
    // printf("------------------------------\n");

    double spline_coordinates[g_inpset->num_basis_spline];
    for (i = 0; i < g_inpset->num_basis_spline; i++){
        spline_coordinates[i] = state[curr]->pos[i + 5] + H[i + 5];
    }

    double m[] = {g_inpset->alpha_param_1, g_inpset->gamma_param_1, g_inpset->I0_param_1, 0.1};
    double s[] = {g_inpset->alpha_param_2, g_inpset->gamma_param_2, g_inpset->I0_param_2, 75e-4};


    double logp = 0.0;
    double params[4] = {t_alpha, t_gamma, t_i0, t_phi_inv};
    double inv_gamma[] = {g_inpset->tau_param_1, g_inpset->tau_param_2};


    int nf = g_inpset->num_basis_spline;
    double ** KK = malloc((nf - 2) * sizeof(double *));
    for (i = 0; i < (nf - 2); i++){
        KK[i] = malloc(nf * sizeof(double));
    }
    double ** KKT = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++){
        KKT[i] = malloc((nf - 2) * sizeof(double));
    }
    double ** result_matrix = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++){
        result_matrix[i] = malloc(nf * sizeof(double));
    }

    KK[0][0] = 1;
    KK[0][1] = -2;
    KK[0][2] = 1;
    for(j = 3; j < nf; j++){
        KK[0][j] = 0;
    }
    for (i = 1; i < nf - 2; i++){
        for(j = 0; j < i; j++){
            KK[i][j] = 0;
        }
        KK[i][i] = 1;
        KK[i][i + 1] = -2;
        KK[i][i + 2] = 1;
        for(j = i + 3; j < nf; j++){
            KK[i][j] = 0;
        }
    }

    KKT[0][0] = 1;
    KKT[1][0] = -2;
    KKT[2][0] = 1;
    for(i = 3; i < nf; i++){
        KKT[i][0] = 0;
    }
    for (j = 1; j < nf - 2; j++){
        for(i = 0; i < j; i++){
            KKT[i][j] = 0;
        }
        KKT[j][j] = 1;
        KKT[j + 1][j] = -2;
        KKT[j + 2][j] = 1;
        for(i = j + 3; i < nf; i++){
            KKT[i][j] = 0;
        }
    }

    matrix_matrix_mult(KKT, KK, result_matrix, nf, nf, nf - 2);

    double result_mv[nf];
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);

    logp += -(params[0] - m[0])*(params[0] - m[0])/(2*s[0]*s[0]) - log(sqrt(2*M_PI)*s[0]);
    logp += -(params[1] - m[1])*(params[1] - m[1])/(2*s[1]*s[1]) - log(sqrt(2*M_PI)*s[1]);
    logp += -(params[2] - m[2])*(params[2] - m[2])/(2*s[2]*s[2]) - log(sqrt(2*M_PI)*s[2]);
    if(phiInvExp == 1){
        logp += log(g_inpset->phi_inv_param_1) - g_inpset->phi_inv_param_1 * 	params[3];
    } else {
        logp += -(params[3] - m[3])*(params[3] - m[3])/(2*s[3]*s[3]) - log(sqrt(2*M_PI)*s[3]);
    }

    logp += inv_gamma[0]*log(inv_gamma[1]) - lgamma(inv_gamma[0]) - (inv_gamma[0] + 1) * log(t_tau) - inv_gamma[1]/t_tau;
    logp += -(1/(2*t_tau*t_tau)) * vectors_dot_prod(spline_coordinates, result_mv, nf);

    /*printf("------------------------------\n");
    printf("Log-prior: %.8f \n", logp);
    printf("------------------------------\n");*/

    return logp;
}

double LoglikelihoodSEMIKR(){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp_E + g_inpset->num_comp;

    double t_alpha, t_gamma, t_i0, t_phi_inv;

    t_alpha = exp(state[curr]->pos[0]) + alpha0;
    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[1]));
    // t_gamma =  exp(state[curr]->pos[1]) + gamma0;
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-state[curr]->pos[2]));
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;

    double t_phi = 1 / t_phi_inv;

    int i;

    char *sensitivities[] = {"-nosensi"};
    double parameters_odes[3 + g_inpset->num_basis_spline];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][3 + g_inpset->num_basis_spline];


    parameters_odes[0] = t_alpha;
    parameters_odes[1] = t_gamma;
    parameters_odes[2] = t_i0;
    for(i = 3; i < (3 + g_inpset->num_basis_spline); i++){
        parameters_odes[i] = state[curr]->pos[i + 2];
    }

    int m = (g_inpset->num_basis_spline - g_inpset->spline_pol_degree - 1) + 2*(g_inpset->spline_pol_degree + 1);
    int nf = g_inpset->num_basis_spline;
    double dx, knots[m], basis[g_inpset->num_basis_spline];
    dx = (g_inpset->tfin_spl_bas - X0)/(g_inpset->num_basis_spline - g_inpset->spline_pol_degree);
    for (i = 0; i < m; i++){
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }

    SEIR_CVODES(results_y, results_s, parameters_odes, g_inpset->spline_pol_degree, m, nf, knots, basis, sensitivities);

    double C_I[kNumObs];
    C_I[0] = results_y[NNEQ - 1][0] - t_i0;
    for (i = 1; i < kNumObs; i++){
        C_I[i] = results_y[NNEQ - 1][i] - results_y[NNEQ - 1][i - 1];
    }

    double u_rep;
    int T0_under_rep, T1_under_rep;
    double U0_under_rep, U1_under_rep;
    if (g_inpset->do_under_report == 0){
        T0_under_rep = 0;
        T1_under_rep = g_inpset-> tfin_spl_bas;
        U0_under_rep = U1_under_rep = 1.0;
    } else {
        T0_under_rep = g_inpset-> T0_under;
        T1_under_rep = g_inpset-> T1_under;
        U0_under_rep = g_inpset-> U0_under;
        U1_under_rep = g_inpset-> U1_under;
    }
    double logL = 0.0;
    for (i = 0; i < kNumObs; i++){
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        logL += lgamma(g_obs[0][i] / u_rep + t_phi) - lgamma(t_phi) - lgamma(g_obs[0][i] / u_rep + 1.0);
        logL += (g_obs[0][i] / u_rep) * (log(C_I[i]) - log(C_I[i] + t_phi));
        logL += t_phi*(log(t_phi) - log(C_I[i] + t_phi));
    }

    /*printf("------------------------------\n");
    printf("Loglikelihood: %.8f \n", logL);
    printf("------------------------------\n");*/

    return logL;
}

double LoglikelihoodSEMIKR_2(double *H){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp_E + g_inpset->num_comp;

    double t_alpha, t_gamma, t_i0, t_phi_inv;

    t_alpha = exp(state[curr]->pos[0] + H[0]) + alpha0;
    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-(state[curr]->pos[1] + H[1])));
    // t_gamma =  exp(state[curr]->pos[1] + H[1]) + gamma0;
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-(state[curr]->pos[2] + H[2])));
    t_phi_inv = exp(state[curr]->pos[3] + H[3]) + phiInv0;

    double t_phi = 1 / t_phi_inv;

    int i;

    char *sensitivities[] = {"-nosensi"};
    double parameters_odes[3 + g_inpset->num_basis_spline];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][3 + g_inpset->num_basis_spline];


    parameters_odes[0] = t_alpha;
    parameters_odes[1] = t_gamma;
    parameters_odes[2] = t_i0;
    for(i = 3; i < (3 + g_inpset->num_basis_spline); i++){
        parameters_odes[i] = state[curr]->pos[i + 2] + H[i + 2];
    }

    int m = (g_inpset->num_basis_spline - g_inpset->spline_pol_degree - 1) + 2*(g_inpset->spline_pol_degree + 1);
    int nf = g_inpset->num_basis_spline;
    double dx, knots[m], basis[g_inpset->num_basis_spline];
    dx = (g_inpset->tfin_spl_bas - X0)/(g_inpset->num_basis_spline - g_inpset->spline_pol_degree);
    for (i = 0; i < m; i++){
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }

    SEIR_CVODES(results_y, results_s, parameters_odes, g_inpset->spline_pol_degree, m, nf, knots, basis, sensitivities);

    double C_I[kNumObs];
    C_I[0] = results_y[NNEQ - 1][0] - t_i0;
    for (i = 1; i < kNumObs; i++){
        C_I[i] = results_y[NNEQ - 1][i] - results_y[NNEQ - 1][i - 1];
    }

    double u_rep;
    int T0_under_rep, T1_under_rep;
    double U0_under_rep, U1_under_rep;
    if (g_inpset->do_under_report == 0){
        T0_under_rep = 0;
        T1_under_rep = g_inpset-> tfin_spl_bas;
        U0_under_rep = U1_under_rep = 1.0;
    } else {
        T0_under_rep = g_inpset-> T0_under;
        T1_under_rep = g_inpset-> T1_under;
        U0_under_rep = g_inpset-> U0_under;
        U1_under_rep = g_inpset-> U1_under;
    }
    double logL = 0.0;
    for (i = 0; i < kNumObs; i++){
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        logL += lgamma(g_obs[0][i] / u_rep + t_phi) - lgamma(t_phi) - lgamma(g_obs[0][i] / u_rep + 1.0);
        logL += (g_obs[0][i] / u_rep) * (log(C_I[i]) - log(C_I[i] + t_phi));
        logL += t_phi*(log(t_phi) - log(C_I[i] + t_phi));
    }

    /*printf("------------------------------\n");
    printf("Loglikelihood: %.8f \n", logL);
    printf("------------------------------\n");*/

    return logL;
}

void GradientSEMIKR(TypeState *st){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp_E + g_inpset->num_comp;

    int i,j;
    double spline_coordinates[g_inpset->num_basis_spline];
    for (i = 0; i < g_inpset->num_basis_spline; i++){
        spline_coordinates[i] = st->pos[i + 5];
    }

    double t_alpha, t_gamma, t_i0, t_phi_inv;

    t_alpha = exp(st->pos[0]) + alpha0;
    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-st->pos[1]));
    // t_gamma =  exp(st->pos[1]) + gamma0;
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-st->pos[2]));
    t_phi_inv = exp(st->pos[3]) + phiInv0;

    double t_phi = 1 / t_phi_inv;
    double t_tau =   exp(st->pos[4]);

    double J1;
    if(isinf(st->pos[1])){
        J1 = 0;
    } else {
        J1 = ((g_inpset->gammaUpper - g_inpset->gammaLower) * exp(-st->pos[1])) / ((1 + exp(-st->pos[1])) * (1 + exp(-st->pos[1])));
        if(isnan(J1)){
            J1 = 0;
        }
    }
    double J2;
    if(isinf(st->pos[2])){
        J2 = 0;
    } else {
        J2 = ((num_pop - 1)*exp(-st->pos[2]))/((1 + exp(-st->pos[2]))*(1 + exp(-st->pos[2])));
        if(isnan(J2)){
            J2 = 0;
        }
    }

    double Jacobian[5 + g_inpset->num_basis_spline];
    Jacobian[0] = t_alpha - alpha0;
    // Jacobian[1] = t_gamma - gamma0;
    Jacobian[1] = J1;
    Jacobian[2] = J2;
    Jacobian[3] = t_phi_inv - phiInv0;
    Jacobian[4] = t_tau;
    for (i = 5; i < (5 + g_inpset->num_basis_spline); i++){
        Jacobian[i] = 1;
    }
    double DJacobian[5 + g_inpset->num_basis_spline];
    DJacobian[0] = 1;
    // DJacobian[1] = 1;
    DJacobian[1] = (1 - exp(st->pos[1])) / (1 + exp(st->pos[1]));
    DJacobian[2] = (1 - exp(st->pos[2])) / (1 + exp(st->pos[2]));
    DJacobian[3] = 1;
    DJacobian[4] = 1;
    for (i = 5; i < (5 + g_inpset->num_basis_spline); i++){
        DJacobian[i] = 0;
    }

    char *sensitivities[] = {"-sensi", "stg", "t"};
    double parameters_odes[3 + g_inpset->num_basis_spline];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][3 + g_inpset->num_basis_spline];

    parameters_odes[0] = t_alpha;
    parameters_odes[1] = t_gamma;
    parameters_odes[2] = t_i0;
    for(i = 3; i < (3 + g_inpset->num_basis_spline); i++){
        parameters_odes[i] = st->pos[i + 2];
    }

    int mm = (g_inpset->num_basis_spline - g_inpset->spline_pol_degree - 1) + 2*(g_inpset->spline_pol_degree + 1);
    int nf = g_inpset->num_basis_spline;
    double dx, knots[mm], basis[g_inpset->num_basis_spline];
    dx = (g_inpset->tfin_spl_bas - X0)/(g_inpset->num_basis_spline - g_inpset->spline_pol_degree);
    for (i = 0; i < mm; i++){
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }

    SEIR_CVODES(results_y, results_s, parameters_odes, g_inpset->spline_pol_degree, mm, nf, knots, basis, sensitivities);

    double C_I[kNumObs];
    C_I[0] = results_y[NNEQ - 1][0] - t_i0;
    for (i = 1; i < kNumObs; i++){
        C_I[i] = results_y[NNEQ - 1][i] - results_y[NNEQ - 1][i - 1];
    }

    double m[] = {g_inpset->alpha_param_1, g_inpset->gamma_param_1, g_inpset->I0_param_1, 0.1};
    double s[] = {g_inpset->alpha_param_2, g_inpset->gamma_param_2, g_inpset->I0_param_2, 75e-4};

    double inv_gamma[] = {g_inpset->tau_param_1, g_inpset->tau_param_2};
    double u_rep;
    int T0_under_rep, T1_under_rep;
    double U0_under_rep, U1_under_rep;
    if (g_inpset->do_under_report == 0){
        T0_under_rep = 0;
        T1_under_rep = g_inpset-> tfin_spl_bas;
        U0_under_rep = U1_under_rep = 1.0;
    } else {
        T0_under_rep = g_inpset-> T0_under;
        T1_under_rep = g_inpset-> T1_under;
        U0_under_rep = g_inpset-> U0_under;
        U1_under_rep = g_inpset-> U1_under;
    }

    for (j = 0; j < 1; j++){
        u_rep = underReport(1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        st->grad[j] = results_s[NNEQ - 1][0][j]*((g_obs[0][0] / u_rep) / C_I[0] - (g_obs[0][0] / u_rep + t_phi)/(C_I[0] + t_phi));
        for (i = 1; i < kNumObs; i++){
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            st->grad[j] += (results_s[NNEQ - 1][i][j] - results_s[NNEQ - 1][i - 1][j])*((g_obs[0][i] / u_rep) / C_I[i] - (g_obs[0][i] / u_rep + t_phi)/(C_I[i] + t_phi));
       }
       st->grad[j] += -(parameters_odes[j] - m[j])/(s[j]*s[j]);
       st->grad[j] *= Jacobian[j];
       st->grad[j] += DJacobian[j];
    }

    for (j = 1;j < 2; j++){
        u_rep = underReport(1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        st->grad[j] = results_s[NNEQ - 1][0][j] * ((g_obs[0][0] / u_rep) / C_I[0] - (g_obs[0][0] / u_rep + t_phi)/(C_I[0] + t_phi));
        for (i = 1; i < kNumObs; i++){
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            st->grad[j] += (results_s[NNEQ - 1][i][j] - results_s[NNEQ - 1][i - 1][j])*((g_obs[0][i] / u_rep) / C_I[i] - (g_obs[0][i] / u_rep + t_phi)/(C_I[i] + t_phi));
        }
        st->grad[j] += -(parameters_odes[j] - m[j])/(s[j]*s[j]);
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }

    for (j = 2;j < 3; j++){
        u_rep = underReport(1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        st->grad[j] = (results_s[NNEQ - 1][0][j] - 1)*((g_obs[0][0] / u_rep) / C_I[0] - (g_obs[0][0] / u_rep + t_phi)/(C_I[0] + t_phi));
        for (i = 1; i < kNumObs; i++){
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            st->grad[j] += (results_s[NNEQ - 1][i][j] - results_s[NNEQ - 1][i - 1][j])*((g_obs[0][i] / u_rep) / C_I[i] - (g_obs[0][i] / u_rep + t_phi)/(C_I[i] + t_phi));
        }
        st->grad[j] += -(parameters_odes[j] - m[j])/(s[j]*s[j]);
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }

    st->grad[3] = 0;
    for (i = 0; i < kNumObs; i++){
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        st->grad[3] += gsl_sf_psi(g_obs[0][i] / u_rep + t_phi) - gsl_sf_psi(t_phi);
        st->grad[3] += (C_I[i] - (g_obs[0][i] / u_rep))/(C_I[i] + t_phi);
        st->grad[3] += log(t_phi) - log(C_I[i] + t_phi);
    }
    st->grad[3] *= -1 / (t_phi_inv * t_phi_inv);
    if(phiInvExp == 1){
        st->grad[3] += -g_inpset->phi_inv_param_1;
    } else {
        st->grad[3] += -(t_phi_inv - m[3])/(s[3]*s[3]);
    }
    st->grad[3] *= Jacobian[3];
    st->grad[3] += DJacobian[3];

    double ** KK = malloc((nf - 2) * sizeof(double *));
    for (i = 0; i < (nf - 2); i++){
        KK[i] = malloc(nf * sizeof(double));
    }
    double ** KKT = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++){
        KKT[i] = malloc((nf - 2) * sizeof(double));
    }
    double ** result_matrix = malloc(nf * sizeof(double *));
    for (i = 0; i < nf; i++){
        result_matrix[i] = malloc(nf * sizeof(double));
    }

    KK[0][0] = 1;
    KK[0][1] = -2;
    KK[0][2] = 1;
    for(j = 3; j < nf; j++){
        KK[0][j] = 0;
    }
    for (i = 1; i < nf - 2; i++){
        for(j = 0; j < i; j++){
            KK[i][j] = 0;
        }
        KK[i][i] = 1;
        KK[i][i + 1] = -2;
        KK[i][i + 2] = 1;
        for(j = i + 3; j < nf; j++){
            KK[i][j] = 0;
        }
    }

    KKT[0][0] = 1;
    KKT[1][0] = -2;
    KKT[2][0] = 1;
    for(i = 3; i < nf; i++){
        KKT[i][0] = 0;
    }
    for (j = 1; j < nf - 2; j++){
        for(i = 0; i < j; i++){
            KKT[i][j] = 0;
        }
        KKT[j][j] = 1;
        KKT[j + 1][j] = -2;
        KKT[j + 2][j] = 1;
        for(i = j + 3; i < nf; i++){
            KKT[i][j] = 0;
        }
    }

    matrix_matrix_mult(KKT, KK, result_matrix, nf, nf, nf - 2);

    double result_mv[nf];
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);

    st->grad[4] = 0;
    st->grad[4] += -(inv_gamma[0] + 1)/t_tau + inv_gamma[1]/(t_tau * t_tau) + vectors_dot_prod(spline_coordinates, result_mv, nf)/(t_tau * t_tau * t_tau);
    st->grad[4] *= Jacobian[4];
    st->grad[4] += DJacobian[4];

    for (j = 5; j < (5 + g_inpset->num_basis_spline); j++){
        u_rep = underReport(1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        st->grad[j] = results_s[NNEQ - 1][0][j - 2]*((g_obs[0][0] / u_rep) / C_I[0] - (g_obs[0][0] / u_rep + t_phi)/(C_I[0] + t_phi));
        for (i = 1; i < kNumObs; i++){
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            st->grad[j] += (results_s[NNEQ - 1][i][j - 2] - results_s[NNEQ - 1][i - 1][j - 2])*((g_obs[0][i] / u_rep)/C_I[i] - (g_obs[0][i] / u_rep + t_phi)/(C_I[i] + t_phi));
        }
        st->grad[j] += -(1/(t_tau * t_tau)) * result_mv[j - 5];
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }

//     printf("------------------------------\n");
//     printf("gradient: \n");
//     for (i = 0; i < (5 + g_inpset->num_basis_spline); i++){
//         printf("%f ", st->grad[i]);
//     }
//     printf("\n");
//     printf("------------------------------\n");
//     //return 0;
//     fprintf(stderr, "Salimos Gradient\n");
}

double LogJacobianSEMIKR(){
    double Jacobian[5 + g_inpset->num_basis_spline];
    double t_alpha, t_phi_inv;

    t_alpha = exp(state[curr]->pos[0]) + alpha0;
    Jacobian[0] = t_alpha - alpha0;

    double J1;
    if(isinf(state[curr]->pos[1])){
        J1 = 0;
    } else {
        J1 = ((g_inpset->gammaUpper - g_inpset->gammaLower) * exp(-state[curr]->pos[1])) / ((1 + exp(-state[curr]->pos[1])) * (1 + exp(-state[curr]->pos[1])));
        if(isnan(J1)){
            J1 = 0;
        }
    }
    Jacobian[1] = J1;
    double J2;
    if(isinf(state[curr]->pos[2])){
        J2 = 0;
    } else {
        J2 = ((num_pop - 1)*exp(-state[curr]->pos[2]))/((1 + exp(-state[curr]->pos[2]))*(1 + exp(-state[curr]->pos[2])));
        if(isnan(J2)){
            J2 = 0;
        }
    }
    Jacobian[2] = J2;
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;
    Jacobian[3] = t_phi_inv - phiInv0;

    double t_tau =   exp(state[curr]->pos[4]);
    int i;

    Jacobian[4] = t_tau;
    for (i = 5; i < (5 + g_inpset->num_basis_spline); i++){
        Jacobian[i] = 1;
    }

    double logJ = 1.0;
    for(i = 0; i < 5; i++){
        logJ *= Jacobian[i];
    }
    logJ = log(logJ);
    /*printf("------------------------------\n");
    printf("Log-Jacobian: %.8f \n", logJ);
    printf("------------------------------\n");*/
    return logJ;
}

double LogJacobianSEMIKR_2(double *H){
    double Jacobian[5 + g_inpset->num_basis_spline];
    double t_alpha, t_phi_inv;

    t_alpha = exp(state[curr]->pos[0] + H[0]) + alpha0;
    Jacobian[0] = t_alpha - alpha0;

    double J1;
    if(isinf(state[curr]->pos[1])){
        J1 = 0;
    } else {
        J1 = ((g_inpset->gammaUpper - g_inpset->gammaLower) * exp(-(state[curr]->pos[1] + H[1]))) / ((1 + exp(-(state[curr]->pos[1] + H[1]))) * (1 + exp(-(state[curr]->pos[1] + H[1]))));
        if(isnan(J1)){
            J1 = 0;
        }
    }
    Jacobian[1] = J1;
    double J2;
    if(isinf(state[curr]->pos[2])){
        J2 = 0;
    } else {
        J2 = ((num_pop - 1)*exp(-(state[curr]->pos[2] + H[2])))/((1 + exp(-(state[curr]->pos[2] + H[2])))*(1 + exp(-(state[curr]->pos[2] + H[2]))));
        if(isnan(J2)){
            J2 = 0;
        }
    }
    Jacobian[2] = J2;
    t_phi_inv = exp(state[curr]->pos[3] + H[3]) + phiInv0;
    Jacobian[3] = t_phi_inv - phiInv0;

    double t_tau =   exp(state[curr]->pos[4] + H[4]);
    int i;

    Jacobian[4] = t_tau;
    for (i = 5; i < (5 + g_inpset->num_basis_spline); i++){
        Jacobian[i] = 1;
    }

    double logJ = 1.0;
    for(i = 0; i < 5; i++){
        logJ *= Jacobian[i];
    }
    logJ = log(logJ);
    /*printf("------------------------------\n");
    printf("Log-Jacobian: %.8f \n", logJ);
    printf("------------------------------\n");*/
    return logJ;
}

void JacobianMatSEMIKR(double *jac_mat) {

    double t_alpha, t_phi_inv;

    t_alpha = exp(state[curr]->pos[0]) + alpha0;
    jac_mat[0] = 1 / (t_alpha - alpha0);

    double J1;
    if(isinf(state[curr]->pos[1])){
        J1 = 0;
    } else {
        J1 = ((g_inpset->gammaUpper - g_inpset->gammaLower) * exp(-state[curr]->pos[1])) / ((1 + exp(-state[curr]->pos[1])) * (1 + exp(-state[curr]->pos[1])));
        if(isnan(J1)){
            J1 = 0;
        }
    }
    jac_mat[1] = 1 / J1;
    double J2;
    if(isinf(state[curr]->pos[2])){
        J2 = 0;
    } else {
        J2 = ((num_pop - 1)*exp(-state[curr]->pos[2]))/((1 + exp(-state[curr]->pos[2]))*(1 + exp(-state[curr]->pos[2])));
        if(isnan(J2)){
            J2 = 0;
        }
    }
    jac_mat[2] = 1 / J2;
    t_phi_inv = exp(state[curr]->pos[3]) + phiInv0;
    jac_mat[3] = 1 / (t_phi_inv - phiInv0);

    double t_tau = exp(state[curr]->pos[4]);
    int i;

    jac_mat[4] = 1 / t_tau;
    for (i = 5; i < (5 + g_inpset->num_basis_spline); i++){
        jac_mat[i] = 1;
    }
}

void TransfSEMIKR() {
    state[curr]->pos[0] = log(state[curr]->pos[0] - alpha0);
    // state[curr]->pos[1] = log(state[curr]->pos[1] - gamma0);
    state[curr]->pos[1] = log((state[curr]->pos[1] - g_inpset->gammaLower) / (g_inpset->gammaUpper - state[curr]->pos[1]));
    state[curr]->pos[2] = log((state[curr]->pos[2] - 1.0)/(num_pop - state[curr]->pos[2]));
    state[curr]->pos[3] = log(state[curr]->pos[3] - phiInv0);
    state[curr]->pos[4] = log(state[curr]->pos[4]);
}

void InvTransfSEMIKR() {
    state[curr]->pos[0] = exp(state[curr]->pos[0]) + alpha0;
    // state[curr]->pos[1] = exp(state[curr]->pos[1]) + gamma0;
    state[curr]->pos[1] = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower)/(1.0 + exp(-state[curr]->pos[1]));
    state[curr]->pos[2] = 1.0 + (num_pop - 1.0)/(1.0 + exp(-state[curr]->pos[2]));
    state[curr]->pos[3] = exp(state[curr]->pos[3]) + phiInv0;
    state[curr]->pos[4] = exp(state[curr]->pos[4]);
}

void HessianSEMIKR(){
    int iii, kk;
    double h = 1e-8;

    //printf("Entering Hessian: \n");
    for (kk = 0; kk < kNumParam; kk++){
        for (iii = 0; iii < kNumParam; iii++){
            if (iii != kk){
                temporalStatePlusH[curr]->pos[iii] = state[curr]->pos[iii];
                temporalStateMinusH[curr]->pos[iii] = state[curr]->pos[iii];
            } else {
                temporalStatePlusH[curr]->pos[iii] = (state[curr]->pos[iii] + h);
                temporalStateMinusH[curr]->pos[iii] = (state[curr]->pos[iii] - h);
            }
            //printf("%f ", temporalStatePlusH[curr]->pos[iii]);
        }
        //printf("\n");
        GradientSEMIKR(temporalStatePlusH[curr]);
        GradientSEMIKR(temporalStateMinusH[curr]);
        /*printf("gradient %d: \n", k);
        for (i = 0; i < (5 + g_inpset->num_basis_spline); i++){
            printf("%f ", temporalStatePlusH[curr]->grad[i]);
        }
        printf("\n");
        for (i = 0; i < (5 + g_inpset->num_basis_spline); i++){
            printf("%f ", temporalStateMinusH[curr]->grad[i]);
        }
        printf("\n");*/
        //printf("------------------------------\n");
        CSM_STORE(state[curr]->hes, kNumParam, kk, kk, temporalStatePlusH[curr]->grad[kk] / (2 * h) - temporalStateMinusH[curr]->grad[kk] / (2 * h));
        for (iii = kk + 1; iii < kNumParam; iii++){
            CSM_STORE(state[curr]->hes, kNumParam, kk, iii, temporalStatePlusH[curr]->grad[iii] / (2 * h) - temporalStateMinusH[curr]->grad[iii] / (2 * h));
            //state[curr]->hes[kk][iii] = temporalStatePlusH[curr]->grad[iii] / (2 * h) - temporalStateMinusH[curr]->grad[iii] / (2 * h);
            //hessian[k * (4 + g_inpset->num_basis_spline)  + i] = state[curr]->hes[k][i];
        }
    }
    /*printf("------------------------------\n");
    printf("Hessian: \n");
    for (kk = 0; kk < 4 + g_inpset->num_basis_spline; kk++){
        for (iii = 0; iii < (4 + g_inpset->num_basis_spline); iii++){
            printf("%f ", state[curr]->hes[kk][iii]);
        }
        printf("\n");
    }
    printf("------------------------------\n");*/
}

void GradientDescentSEMIKR(){
    printf("Starting Gradient Descent\n");
    printf("------------------------------\n");
    int i,j;
    double logPosterior;

    logPosterior = PriorSEMIKR() + LoglikelihoodSEMIKR() + LogJacobianSEMIKR();
    printf("%d%%\n", 100 * 0 / g_inpset->max_iter);
    printf("------------------------------\n");
    printf("Log-Posterior: %.8f \n", logPosterior);
    printf("------------------------------\n");

    // double h = 1e-7;
    // double H_values[1 + 5 + g_inpset->num_basis_spline][5 + g_inpset->num_basis_spline];
    // double LogPosterior_values[1 + 5 + g_inpset->num_basis_spline];
    // double Num_Derivative[5 + g_inpset->num_basis_spline];

    // for (i = 0; i < (1 + 5 + g_inpset->num_basis_spline); i++){
    //     if (i == 0){
    //         for (j = 0; j < (5 + g_inpset->num_basis_spline); j++){
    //             H_values[i][j] = 0;
    //         }
    //     } else{
    //         for (j = 0; j < (i - 1); j++){
    //             H_values[i][j] = 0;
    //         }
    //         H_values[i][i - 1] = h;
    //         for (j = i; j < (5 + g_inpset->num_basis_spline); j++){
    //             H_values[i][j] = 0;
    //         }
    //     }
    // }

    for(i = 0; i < g_inpset->max_iter; i++){
        GradientSEMIKR(state[curr]);

        // for(j = 0; j < (1 + 5 + g_inpset->num_basis_spline); j++){
        //     LogPosterior_values[j] = PriorSEMIKR_2(H_values[j]) + LoglikelihoodSEMIKR_2(H_values[j]) + LogJacobianSEMIKR_2(H_values[j]);
        // }


        // for(j = 1; j < (1 + 5 + g_inpset->num_basis_spline); j++){
        //     Num_Derivative[j - 1] = (LogPosterior_values[j] - LogPosterior_values[0]) / h;
        // }
        // printf("Numerical Derivative\n");
        // for (j = 0; j < (5 + g_inpset->num_basis_spline); j++){
        //     printf("%f ", Num_Derivative[j]);
        // }
        // printf("\n");
        printf("------------------------------\n");
        printf("gradient: \n");
        for (j = 0; j < (5 + g_inpset->num_basis_spline); j++){
            printf("%f ", state[curr]->grad[j]);
        }
        printf("\n");
        printf("------------------------------\n");

        for(j = 0; j < kNumParam; j++){
            state[curr]->pos[j] = state[curr]->pos[j] + g_inpset->learning_rate * state[curr]->grad[j];
        }

        if (i % (g_inpset->max_iter / 10) == 0 && i > 0){
            logPosterior = PriorSEMIKR() + LoglikelihoodSEMIKR() + LogJacobianSEMIKR();
            printf("%d%%\n", 100 * i / g_inpset->max_iter);
            printf("------------------------------\n");
            printf("Log-Posterior: %.8f \n", logPosterior);
            printf("------------------------------\n");
        }
    }
    printf("%d%%\n", 100);
    printf("------------------------------\n");
}

int ReadModelSEMIKR_Incidence() {

    //printf("Entramos: ReadModelSEMIKR_synth\n");

    char filename[100];
    FILE *fptr_obs;
    FILE *fptr_initialpoint;
    FILE *fptr_splinebasis;
    int i = 0, ii = 0;
    int symbol;
    double num;

    kNumObs = 0;

    //check if num_basis_spline < 3
    if(g_inpset->num_basis_spline < 3){
        fprintf(stderr, "seir_model_synthetic.c - num_basis_spline has to at least 3, but it is '%d'", g_inpset->num_basis_spline);
        exit(EXIT_FAILURE);
    }

    //open file with dataset
    sprintf(filename, "./benchmarks/%s_data.txt", g_inpset->data);

    fptr_obs = fopen(filename, "r");
    if((fptr_obs) == NULL){
        fprintf(stderr, "seir_model_synthetic.c - Failed to open the dataset file '%s'", filename);
        exit(EXIT_FAILURE);
    }

    //count newlines in the whole file (that is the number of observations kNumObs)
    while(!feof(fptr_obs)){
        symbol = fgetc(fptr_obs);
        if(symbol == '\n'){ kNumObs++; }
    }

    rewind(fptr_obs); //go to the beginning of the file

    //allocate memory for the data, only one column (new infectious)
    g_obs = MatrixAlloc(1, kNumObs);

    //read in the file and write data in g_obs
    for(ii = 0; ii < kNumObs; ii++){
        if(fscanf(fptr_obs, "%lf", &num) != EOF){
            g_obs[0][ii] = num;
        }
    }

    num_pop = g_obs[0][kNumObs - 1];

    kNumObs = kNumObs - 1;

    fprintf(red_lfp, "Number of observations %d\n", kNumObs);

    fprintf(red_lfp, "Population size %f\n", num_pop);

    rewind(fptr_obs); //go to the beginning of the file

    /* assign number of model parameters and functions */
    kNumParam = 5 + g_inpset->num_basis_spline;
    Prior = PriorSEMIKR;
    Loglikelihood = LoglikelihoodSEMIKR;
    Gradient = GradientSEMIKR;
    Transf = TransfSEMIKR;
    InvTransf = InvTransfSEMIKR;
    LogJacobian = LogJacobianSEMIKR;
    JacobianMat = JacobianMatSEMIKR;
    Hessian = HessianSEMIKR;

    /** initialize state of model parameters ***/
    for (i = back2; i <= forw2; i++) {
        SAlloc(&(state[i]), kNumParam);
        SAlloc(&(temporalStatePlusH[i]), kNumParam);
        SAlloc(&(temporalStateMinusH[i]), kNumParam);
    }

    ///////////////////////////////
    //initialize position variables

    //check if file exists. if yes, read in the initial point. if no, generate a random one.
    int bool_initialpointfile_exists = -1;
    sprintf(filename, "./output/%s/initialpoint_%s.txt", id_global, id_global);
    fptr_initialpoint = fopen(filename, "r");
    if ((fptr_initialpoint) == NULL) {
        bool_initialpointfile_exists = 0;
    } else{
        bool_initialpointfile_exists = 1;
    }

    //check if spline basis file exists
    int bool_splinebasis_exists = -1;
    sprintf(filename, "./output/%s/splinebasis_%s.txt", id_global, id_global);
    fptr_splinebasis = fopen(filename, "r");
    if ((fptr_splinebasis) == NULL) {
        bool_splinebasis_exists = 0;
    } else{
        bool_splinebasis_exists = 1;
    }

    fprintf(red_lfp, "bool_initialpointfile_exists is %i\n", bool_initialpointfile_exists);

    fprintf(red_lfp, "bool_splinebasis_exists is %i\n", bool_splinebasis_exists);

    //depending on inputfile generate a random starting point or read in from file
    if((bool_initialpointfile_exists == 0) || (bool_splinebasis_exists == 0)){
        fprintf(stderr, "Error in SEIR model input: initialpoints and/or splinebasis files not found. \n");
        exit(EXIT_FAILURE);
    } else if((bool_initialpointfile_exists == 1) && (bool_splinebasis_exists == 1)){ //read in from file

        // count number of spaces in the first line
        int numSizeInitialPoint = 0;
        symbol = 1; //initialise symbol
        while(!feof(fptr_initialpoint) && symbol != '\n'){
            symbol = fgetc(fptr_initialpoint);
            if(isspace(symbol)){ numSizeInitialPoint++; }
        }

        fprintf(red_lfp, "numSizeInitialPoint is %i\n", numSizeInitialPoint);

        int numSizeSplineBasis = 0;
        symbol = 1; //initialise symbol
        while(!feof(fptr_splinebasis) && symbol != '\n'){
            symbol = fgetc(fptr_splinebasis);
            if(isspace(symbol)){ numSizeSplineBasis++; }
        }

        fprintf(red_lfp, "numSizeSplineBasis is %i\n", numSizeSplineBasis);

        if(g_inpset->num_basis_spline != numSizeSplineBasis){
            fprintf(stderr, "seir_model_synthetic.c - num_basis_spline (%d) is not equal to numSizeSplineBasis (%i) in seir_model_synthetic.c. Exiting...\n", g_inpset->num_basis_spline, numSizeSplineBasis);
            exit(EXIT_FAILURE);
        }

        if(kNumParam != numSizeInitialPoint + numSizeSplineBasis){
            fprintf(stderr, "seir_model_synthetic.c - kNumParam (%d) is not equal to numSizeInitialPoint + numSizeSplineBasis (%i) in seir_model_synthetic.c. Exiting...\n", kNumParam, numSizeInitialPoint + numSizeSplineBasis);
            exit(EXIT_FAILURE);
        }

        //go to beginning of file
        rewind(fptr_initialpoint);
        rewind(fptr_splinebasis);

        //read initial points in one by one
        for (ii = 0; ii < numSizeInitialPoint; ii++) {
            if (fscanf(fptr_initialpoint, "%lf", &num) != EOF) {
            state[curr]->pos[ii] = num;
            }
        }

        fclose(fptr_initialpoint);

        //read spline basis in one by one
        for (ii = 0; ii < numSizeSplineBasis; ii++) {
            if (fscanf(fptr_splinebasis, "%lf", &num) != EOF) {
                state[curr]->pos[numSizeInitialPoint + ii] = num;
            }
        }

        fprintf(red_lfp, "initialpoints are");
        for (i=0; i < numSizeInitialPoint + numSizeSplineBasis; i++) {
            fprintf(red_lfp, "%lf ", state[curr]->pos[i]);
        }

        fclose(fptr_splinebasis);

    } else{
        fprintf(stderr, "seir_model_synthetic.c - If-else not matched for initial position\n");
        exit(EXIT_FAILURE);
    }
    // end initialize position variables

    //always allocate memory for current hessian
    global_misc->dim1_hes = kNumParam;
    global_misc->dim2_hes = kNumParam;
    state[curr]->hes = MatrixAlloc(global_misc->dim1_hes, global_misc->dim2_hes); //allocate and init with zeroes

    if (g_inpset->do_grad_desc == 1){
        TransfSEMIKR();
        GradientDescentSEMIKR();
        InvTransfSEMIKR();
        for(ii=0 ; ii < 5 + g_inpset->num_basis_spline; ii++ ){
            fprintf(fp_startingpoint, "%lf\n", state[curr]->pos[ii]);
        }
    }

    fclose(fp_startingpoint);

    /*** initialize momenta vectors from a gaussian distr ****/
    if (state[curr]->mom == NULL) {
        printf("seir_model_synthetic.c - Error in allocating memory (state mom)");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < kNumParam; i++) {
        state[curr]->mom[i] = gsl_ran_gaussian_ziggurat(g_rng, 1.);
    }

    /***** END initialization ******/
    Ham = (TypeHamiltonian *)malloc(sizeof(TypeHamiltonian));

    return 0;
}
