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

#include "sir_stand_model_synthetic.h"
#include "sir_seir_common_functions.h"
#include "Globals.h"
#include "utils.h"
#include "sir_stand_functions.h"
#include "read_input.h"

#define M_PI 3.14159265358979323846
#define X0 0.0

#define gamma0 1e-10
#define phiInv0 1e-10
#define beta0 1e-10

#define phiInvExp 1
//#define gammaL 0.03333333333333 
//#define gammaU 1.0

double PriorSIKR_stand_Incidence(){

    double t_gamma, t_i0, t_phi_inv, t_beta;

    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[0]));
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-state[curr]->pos[1]));
    t_phi_inv = exp(state[curr]->pos[2]) + phiInv0;
    t_beta = exp(state[curr]->pos[3]) + beta0;

    //double t_tau =   exp(state[curr]->pos[3]);
    //int i, j;
    printf("------------------------------\n");
    printf("HMC proposal transformed\n");
    printf("%f %f %f %f ", t_gamma, t_i0, t_phi_inv, t_beta);
    /*for (i = 4; i < (kNumParam); i++){
        printf("%f ", state[curr]->pos[i]);
    }*/
    printf("\n");
    printf("------------------------------\n");

    /*double spline_coordinates[g_inpset->num_basis_spline];
    for (i = 0; i < g_inpset->num_basis_spline; i++){
        spline_coordinates[i] = state[curr]->pos[i + 4];
    }*/

    double m[] = {g_inpset->gamma_param_1, g_inpset->I0_param_1, 0.1, g_inpset->beta_mean};
    double s[] = {g_inpset->gamma_param_2, g_inpset->I0_param_2, 75e-4, g_inpset->beta_sd};

    double logp = 0.0;
    double params[4] = {t_gamma, t_i0, t_phi_inv, t_beta};
    //double inv_gamma[] = {g_inpset->tau_param_1, g_inpset->tau_param_2};

    /*int nf = g_inpset->num_basis_spline;
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
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);*/
    
    logp += -(params[0] - m[0])*(params[0] - m[0])/(2*s[0]*s[0]) - log(sqrt(2*M_PI)*s[0]);
    logp += -(params[1] - m[1])*(params[1] - m[1])/(2*s[1]*s[1]) - log(sqrt(2*M_PI)*s[1]);
    if(phiInvExp == 1){
	logp += log(g_inpset->phi_inv_param_1) - g_inpset->phi_inv_param_1 * 	params[2];
    } else {
	logp += -(params[2] - m[2])*(params[2] - m[2])/(2*s[2]*s[2]) - log(sqrt(2*M_PI)*s[2]);
    }
    //logp += inv_gamma[0]*log(inv_gamma[1]) - lgamma(inv_gamma[0]) - (inv_gamma[0] + 1) * log(t_tau) - inv_gamma[1]/t_tau;
    //logp += -(1/(2*t_tau*t_tau)) * vectors_dot_prod(spline_coordinates, result_mv, nf);
    logp += -(params[3] - m[3])*(params[3] - m[3])/(2*s[3]*s[3]) - log(sqrt(2*M_PI)*s[3]);

    return logp;
}

double PriorSIKR_stand_Incidence_2(double *H){

    double t_gamma, t_i0, t_phi_inv, t_beta;

    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-(state[curr]->pos[0] + H[0])));
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-(state[curr]->pos[1] + H[1])));
    t_phi_inv = exp(state[curr]->pos[2] + H[2]) + phiInv0;
    t_beta = exp(state[curr]->pos[3] + H[3]) + beta0;

    //double t_tau = exp(state[curr]->pos[3] + H[3]);
    //int i, j;
    /*printf("------------------------------\n");
    printf("HMC proposal transformed\n");
    printf("%f %f %f %f ", t_gamma, t_i0, t_phi_inv, t_tau);
    for (i = 4; i < (kNumParam); i++){
        printf("%f ", state[curr]->pos[i]);
    }
    printf("\n");
    printf("------------------------------\n");*/

    /*double spline_coordinates[g_inpset->num_basis_spline];
    for (i = 0; i < g_inpset->num_basis_spline; i++){
        spline_coordinates[i] = state[curr]->pos[i + 4] + H[i + 4];
    }*/

    double m[] = {g_inpset->gamma_param_1, g_inpset->I0_param_1, 0.1, g_inpset->beta_mean};
    double s[] = {g_inpset->gamma_param_2, g_inpset->I0_param_2, 75e-4, g_inpset->beta_sd};

    double logp = 0.0;
    double params[4] = {t_gamma, t_i0, t_phi_inv, t_beta};
    //double inv_gamma[] = {g_inpset->tau_param_1, g_inpset->tau_param_2};

    /*int nf = g_inpset->num_basis_spline;
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
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);*/

    logp += -(params[0] - m[0])*(params[0] - m[0])/(2*s[0]*s[0]) - log(sqrt(2*M_PI)*s[0]);
    logp += -(params[1] - m[1])*(params[1] - m[1])/(2*s[1]*s[1]) - log(sqrt(2*M_PI)*s[1]);
    if(phiInvExp == 1){
	logp += log(g_inpset->phi_inv_param_1) - g_inpset->phi_inv_param_1 * 	params[2];
    } else {
	logp += -(params[2] - m[2])*(params[2] - m[2])/(2*s[2]*s[2]) - log(sqrt(2*M_PI)*s[2]);
    }
    //logp += inv_gamma[0]*log(inv_gamma[1]) - lgamma(inv_gamma[0]) - (inv_gamma[0] + 1) * log(t_tau) - inv_gamma[1]/t_tau;
    //logp += -(1/(2*t_tau*t_tau)) * vectors_dot_prod(spline_coordinates, result_mv, nf);
    logp += -(params[3] - m[3])*(params[3] - m[3])/(2*s[3]*s[3]) - log(sqrt(2*M_PI)*s[3]);

    return logp;
}

double LoglikelihoodSIKR_stand_Incidence(){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp;

    double t_gamma, t_i0, t_phi_inv, t_beta;

    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-state[curr]->pos[0]));
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-state[curr]->pos[1]));
    t_phi_inv = exp(state[curr]->pos[2]) + phiInv0;
    t_beta = exp(state[curr]->pos[3]) + beta0;


    double t_phi = 1 / t_phi_inv;

    int i;

    char *sensitivities[] = {"-nosensi"};
    double parameters_odes[2 + 1];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][2 + 1];


    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_i0;
    parameters_odes[2] = t_beta;

    /*int m = (g_inpset->num_basis_spline - g_inpset->spline_pol_degree - 1) + 2*(g_inpset->spline_pol_degree + 1);
    int nf = g_inpset->num_basis_spline;
    double dx, knots[m], basis[g_inpset->num_basis_spline];
    dx = (g_inpset->tfin_spl_bas - X0)/(g_inpset->num_basis_spline - g_inpset->spline_pol_degree);
    for (i = 0; i < m; i++){
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }*/

    SIR_CVODES_stand_Incidence(results_y, results_s, parameters_odes, sensitivities);

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

    return logL;
}

double LoglikelihoodSIKR_stand_Incidence_2(double *H){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp;

    double t_gamma, t_i0, t_phi_inv, t_beta;

    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-(state[curr]->pos[0] + H[0])));
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-(state[curr]->pos[1] + H[1])));
    t_phi_inv = exp(state[curr]->pos[2]+ H[2]) + phiInv0;
    t_beta = exp(state[curr]->pos[3]+ H[3]) + beta0;


    double t_phi = 1 / t_phi_inv;

    int i;

    char *sensitivities[] = {"-nosensi"};
    double parameters_odes[2 + 1];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][2 + 1];


    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_i0;
    parameters_odes[2] = t_beta;

    /*int m = (g_inpset->num_basis_spline - g_inpset->spline_pol_degree - 1) + 2*(g_inpset->spline_pol_degree + 1);
    int nf = g_inpset->num_basis_spline;
    double dx, knots[m], basis[g_inpset->num_basis_spline];
    dx = (g_inpset->tfin_spl_bas - X0)/(g_inpset->num_basis_spline - g_inpset->spline_pol_degree);
    for (i = 0; i < m; i++){
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }*/

    SIR_CVODES_stand_Incidence(results_y, results_s, parameters_odes, sensitivities);

    double C_I[kNumObs];
    C_I[0] = results_y[NNEQ - 1][0] - t_i0;
    //C_I[0] = results_y[NNEQ - 1][0];
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

    return logL;
}

void GradientSIKR_stand_Incidence(TypeState *st){
    //choose the number of components with respect to the model
    int NNEQ = 3 + g_inpset->num_comp;

    int i,j;
    /*double spline_coordinates[g_inpset->num_basis_spline];
    for (i = 0; i < g_inpset->num_basis_spline; i++){
        spline_coordinates[i] = st->pos[i + 4];
    }*/

    double t_gamma, t_i0, t_phi_inv, t_beta;

    t_gamma = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower) / (1 + exp(-st->pos[0]));
    t_i0 = 1.0 + (num_pop - 1)/(1 + exp(-st->pos[1]));
    t_phi_inv = exp(st->pos[2]) + phiInv0;
    t_beta = exp(st->pos[3]) + beta0;

    double t_phi = 1 / t_phi_inv;
    //double t_tau =   exp(st->pos[3]);
    double Jacobian[kNumParam];

    double J1;
    if(isinf(st->pos[0])){
        J1 = 0;
    } else {
        J1 = ((g_inpset->gammaUpper - g_inpset->gammaLower) * exp(-st->pos[0])) / ((1 + exp(-st->pos[0])) * (1 + exp(-st->pos[0])));
        if(isnan(J1)){
            J1 = 0;
        }
    }
    
    double J2;
    if(isinf(st->pos[1])){
        J2 = 0;
    } else {
        J2 = ((num_pop - 1)*exp(-st->pos[1]))/((1 + exp(-st->pos[1]))*(1 + exp(-st->pos[1])));
        if(isnan(J2)){
            J2 = 0;
        }
    }  
    Jacobian[0] = J1;
    Jacobian[1] = J2;
    Jacobian[2] = t_phi_inv - phiInv0;
    Jacobian[3] = t_beta - beta0;
    for (i = 4; i < (kNumParam); i++){
        Jacobian[i] = 1;
    }
    double DJacobian[kNumParam];
    DJacobian[0] = (1 - exp(st->pos[0])) / (1 + exp(st->pos[0]));
    DJacobian[1] = (1 - exp(st->pos[1])) / (1 + exp(st->pos[1]));
    DJacobian[2] = 1;
    DJacobian[3] = 1;
    for (i = 4; i < (kNumParam); i++){
        DJacobian[i] = 0;
    }

    char *sensitivities[] = {"-sensi", "stg", "t"};
    double parameters_odes[2 + 1];
    double results_y[NNEQ][kNumObs];
    double results_s[NNEQ][kNumObs][2 + 1];


    parameters_odes[0] = t_gamma;
    parameters_odes[1] = t_i0;
    parameters_odes[2] = t_beta;

    /*int mm = (g_inpset->num_basis_spline - g_inpset->spline_pol_degree - 1) + 2*(g_inpset->spline_pol_degree + 1);
    int nf = g_inpset->num_basis_spline;
    double dx, knots[mm], basis[g_inpset->num_basis_spline];
    dx = (g_inpset->tfin_spl_bas - X0)/(g_inpset->num_basis_spline - g_inpset->spline_pol_degree);
    for (i = 0; i < mm; i++){
        knots[i] = X0 + (i - g_inpset->spline_pol_degree) * dx;
    }*/

    SIR_CVODES_stand_Incidence(results_y, results_s, parameters_odes, sensitivities);

    double C_I[kNumObs];
    C_I[0] = results_y[NNEQ - 1][0] - t_i0;
    for (i = 1; i < kNumObs; i++){
        C_I[i] = results_y[NNEQ - 1][i] - results_y[NNEQ - 1][i - 1];
    }

    double m[] = {g_inpset->gamma_param_1, g_inpset->I0_param_1, 0.1, g_inpset->beta_mean};
    double s[] = {g_inpset->gamma_param_2, g_inpset->I0_param_2, 75e-4, g_inpset->beta_sd};

    //double inv_gamma[] = {g_inpset->tau_param_1, g_inpset->tau_param_2};
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

    for (j = 0;j < 1; j++){
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
        st->grad[j] = (results_s[NNEQ - 1][0][j] - 1)*((g_obs[0][0] / u_rep) / C_I[0] - (g_obs[0][0] / u_rep + t_phi)/(C_I[0] + t_phi));
        for (i = 1; i < kNumObs; i++){
	    u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
	    st->grad[j] += (results_s[NNEQ - 1][i][j] - results_s[NNEQ - 1][i - 1][j])*((g_obs[0][i] / u_rep) / C_I[i] - (g_obs[0][i] / u_rep + t_phi)/(C_I[i] + t_phi));
        }
        st->grad[j] += -(parameters_odes[j] - m[j])/(s[j]*s[j]);
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }

    st->grad[2] = 0;
    for (i = 0; i < kNumObs; i++){
        u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        st->grad[2] += gsl_sf_psi(g_obs[0][i] / u_rep + t_phi) - gsl_sf_psi(t_phi);
        st->grad[2] += (C_I[i] - (g_obs[0][i] / u_rep))/(C_I[i] + t_phi);
        st->grad[2] += log(t_phi) - log(C_I[i] + t_phi);
    }
    st->grad[2] *= -1 / (t_phi_inv * t_phi_inv);
    if(phiInvExp == 1){
        st->grad[2] += -g_inpset->phi_inv_param_1;
    } else {
        st->grad[2] += -(t_phi_inv - m[2])/(s[2]*s[2]);
    }
    st->grad[2] *= Jacobian[2];
    st->grad[2] += DJacobian[2];

    /*double ** KK = malloc((nf - 2) * sizeof(double *));
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
    matrix_vector_mult(result_matrix, spline_coordinates, result_mv, nf, nf);*/

    /*st->grad[3] = 0;
    st->grad[3] += -(inv_gamma[0] + 1)/t_tau + inv_gamma[1]/(t_tau * t_tau) + vectors_dot_prod(spline_coordinates, result_mv, nf)/(t_tau * t_tau * t_tau);
    st->grad[3] *= Jacobian[3];
    st->grad[3] += DJacobian[3];*/

    for (j = 3; j < (kNumParam); j++){
        u_rep = underReport(1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
        st->grad[j] = results_s[NNEQ - 1][0][j - 1]*((g_obs[0][0] / u_rep)/C_I[0] - (g_obs[0][0] / u_rep + t_phi)/(C_I[0] + t_phi));
        for (i = 1; i < kNumObs; i++){
            u_rep = underReport(i + 1, U0_under_rep, U1_under_rep, T0_under_rep, T1_under_rep);
            st->grad[j] += (results_s[NNEQ - 1][i][j - 1] - results_s[NNEQ - 1][i - 1][j - 1])*((g_obs[0][i] / u_rep)/C_I[i] - (g_obs[0][i] / u_rep + t_phi)/(C_I[i] + t_phi));
        }
        st->grad[j] += -(parameters_odes[j - 1] - m[j])/(s[j]*s[j]);
        st->grad[j] *= Jacobian[j];
        st->grad[j] += DJacobian[j];
    }
}

double LogJacobianSIKR_stand_Incidence(){
    double Jacobian[kNumParam];
    double t_phi_inv;
    
    double J1;
    if(isinf(state[curr]->pos[0])){
        J1 = 0;
    } else {
        J1 = ((g_inpset->gammaUpper - g_inpset->gammaLower) * exp(-state[curr]->pos[0])) / ((1 + exp(-state[curr]->pos[0])) * (1 + exp(-state[curr]->pos[0])));
        if(isnan(J1)){
            J1 = 0;
        }
    }
    Jacobian[0] = J1;

    double J2;
    if(isinf(state[curr]->pos[1])){
        J2 = 0;
    } else {
        J2 = ((num_pop - 1)*exp(-state[curr]->pos[1]))/((1 + exp(-state[curr]->pos[1]))*(1 + exp(-state[curr]->pos[1])));
        if(isnan(J2)){
            J2 = 0;
        }
    }
    Jacobian[1] = J2;

    t_phi_inv = exp(state[curr]->pos[2]) + phiInv0;
    Jacobian[2] = t_phi_inv - phiInv0;

    double t_beta =   exp(state[curr]->pos[3]) + beta0;
    int i;

    Jacobian[3] = t_beta - beta0;
    for (i = 4; i < (kNumParam); i++){
        Jacobian[i] = 1;
    }

    double logJ = 1.0;
    for(i = 0; i < 4; i++){
        logJ *= Jacobian[i];
    }
    logJ = log(logJ);

    return logJ;
}

double LogJacobianSIKR_stand_Incidence_2(double *H){
    double Jacobian[kNumParam];
    double t_phi_inv;

    double J1;
    if(isinf(state[curr]->pos[0])){
        J1 = 0;
    } else {
        J1 = ((g_inpset->gammaUpper - g_inpset->gammaLower) * exp(-(state[curr]->pos[0] + H[0]))) / ((1 + exp(-(state[curr]->pos[0] + H[0]))) * (1 + exp(-(state[curr]->pos[0] + H[0]))));
        if(isnan(J1)){
            J1 = 0;
        }
    }
    Jacobian[0] = J1;

    double J2;
    if(isinf(state[curr]->pos[1])){
        J2 = 0;
    } else {
        J2 = ((num_pop - 1)*exp(-(state[curr]->pos[1] + H[1])))/((1 + exp(-(state[curr]->pos[1] + H[1])))*(1 + exp(-(state[curr]->pos[1] + H[1]))));
        if(isnan(J2)){
            J2 = 0;
        }
    }
    Jacobian[1] = J2;

    t_phi_inv = exp(state[curr]->pos[2] + H[2]) + phiInv0;
    Jacobian[2] = t_phi_inv - phiInv0;

    double t_beta = exp(state[curr]->pos[3] + H[3]) + beta0;
    int i;

    Jacobian[3] = t_beta - beta0;
    for (i = 4; i < (kNumParam); i++){
        Jacobian[i] = 1;
    }

    double logJ = 1.0;
    for(i = 0; i < 4; i++){
        logJ *= Jacobian[i];
    }
    logJ = log(logJ);

    return logJ;
}

void JacobianMatSIKR_stand_Incidence(double *jac_mat) {

    double t_phi_inv, t_beta;

    double J1;
    if(isinf(state[curr]->pos[0])){
        J1 = 0;
    } else {
        J1 = ((g_inpset->gammaUpper - g_inpset->gammaLower) * exp(-state[curr]->pos[0])) / ((1 + exp(-state[curr]->pos[0])) * (1 + exp(-state[curr]->pos[0])));
        if(isnan(J1)){
            J1 = 0;
        }
    }
    jac_mat[0] = 1 / J1;

    double J2;
    if(isinf(state[curr]->pos[1])){
        J2 = 0;
    } else {
        J2 = ((num_pop - 1)*exp(-state[curr]->pos[1]))/((1 + exp(-state[curr]->pos[1]))*(1 + exp(-state[curr]->pos[1])));
        if(isnan(J2)){
            J2 = 0;
        }
    }
    jac_mat[1] = 1 / J2;

    t_phi_inv = exp(state[curr]->pos[2]) + phiInv0;
    jac_mat[2] = 1 / (t_phi_inv - phiInv0);

    t_beta = exp(state[curr]->pos[3]) + beta0;
    int i;

    jac_mat[3] = 1 / (t_beta - beta0);
    for (i = 4; i < (kNumParam); i++){
        jac_mat[i] = 1;
    }
}

void TransfSIKR_stand_Incidence() {
    state[curr]->pos[0] = log((state[curr]->pos[0] - g_inpset->gammaLower) / (g_inpset->gammaUpper - state[curr]->pos[0]));
    state[curr]->pos[1] = log((state[curr]->pos[1] - 1.0)/(num_pop - state[curr]->pos[1]));
    state[curr]->pos[2] = log(state[curr]->pos[2] - phiInv0);
    state[curr]->pos[3] = log(state[curr]->pos[3] - beta0);
}

void InvTransfSIKR_stand_Incidence() {
    state[curr]->pos[0] = g_inpset->gammaLower + (g_inpset->gammaUpper - g_inpset->gammaLower)/(1.0 + exp(-state[curr]->pos[0]));
    state[curr]->pos[1] = 1.0 + (num_pop - 1.0)/(1.0 + exp(-state[curr]->pos[1]));
    state[curr]->pos[2] = exp(state[curr]->pos[2]) + phiInv0;
    state[curr]->pos[3] = exp(state[curr]->pos[3]) + beta0;
}

void HessianSIKR_stand_Incidence(){
    int iii, kk;
    double h = 1e-8;

    for (kk = 0; kk < kNumParam; kk++){
        for (iii = 0; iii < kNumParam; iii++){
            if (iii != kk){
                temporalStatePlusH[curr]->pos[iii] = state[curr]->pos[iii];
                temporalStateMinusH[curr]->pos[iii] = state[curr]->pos[iii];
            } else {
                temporalStatePlusH[curr]->pos[iii] = (state[curr]->pos[iii] + h);
                temporalStateMinusH[curr]->pos[iii] = (state[curr]->pos[iii] - h);
            }
        }

        GradientSIKR_stand_Incidence(temporalStatePlusH[curr]);
        GradientSIKR_stand_Incidence(temporalStateMinusH[curr]);

        CSM_STORE(state[curr]->hes, kNumParam, kk, kk, temporalStatePlusH[curr]->grad[kk] / (2 * h) - temporalStateMinusH[curr]->grad[kk] / (2 * h));
        for (iii = kk + 1; iii < kNumParam; iii++){
            CSM_STORE(state[curr]->hes, kNumParam, kk, iii, temporalStatePlusH[curr]->grad[iii] / (2 * h) - temporalStateMinusH[curr]->grad[iii] / (2 * h));
        }
    }
}

void GradientDescentSIKR_stand_Incidence(){
    printf("Starting Gradient Descent\n");
    printf("------------------------------\n");
    int i,j;
    double logPosterior;
    logPosterior = PriorSIKR_stand_Incidence() + LoglikelihoodSIKR_stand_Incidence() + LogJacobianSIKR_stand_Incidence();
    printf("%d%%\n", 100 * 0 / g_inpset->max_iter);
    printf("------------------------------\n");
    printf("Log-Posterior: %.8f \n", logPosterior);
    printf("------------------------------\n");

    double h = 1e-6;
    double H_values[1 + kNumParam][kNumParam];
    double LogPosterior_values[1 + kNumParam];
    double Num_Derivative[kNumParam];

    for (i = 0; i < (1 + kNumParam); i++){
        if (i == 0){
            for (j = 0; j < (kNumParam); j++){
                H_values[i][j] = 0;
            }
        } else{
            for (j = 0; j < (i - 1); j++){
                H_values[i][j] = 0;
            }
            H_values[i][i - 1] = h;
            for (j = i; j < (kNumParam); j++){
                H_values[i][j] = 0;
            }
        }
    }

    for(i = 0; i < g_inpset->max_iter; i++){
        GradientSIKR_stand_Incidence(state[curr]);

        for(j = 0; j < (1 + kNumParam); j++){
            LogPosterior_values[j] = PriorSIKR_stand_Incidence_2(H_values[j]) + LoglikelihoodSIKR_stand_Incidence_2(H_values[j]) + LogJacobianSIKR_stand_Incidence_2(H_values[j]);
        }


        for(j = 1; j < (1 + kNumParam); j++){
            Num_Derivative[j - 1] = (LogPosterior_values[j] - LogPosterior_values[0]) / h;
        }
        printf("Numerical Derivative\n");
        for (j = 0; j < (kNumParam); j++){
            printf("%f ", Num_Derivative[j]);
        }
        printf("\n");
        printf("------------------------------\n");
        printf("gradient: \n");
        for (j = 0; j < (kNumParam); j++){
            printf("%f ", state[curr]->grad[j]);
        }
        printf("\n");
        printf("------------------------------\n");
        
        for(j = 0; j < (kNumParam); j++){
            state[curr]->pos[j] = state[curr]->pos[j] + g_inpset->learning_rate * state[curr]->grad[j];
        }

        if (i % (g_inpset->max_iter / 10) == 0 && i > 0){
            logPosterior = PriorSIKR_stand_Incidence() + LoglikelihoodSIKR_stand_Incidence() + LogJacobianSIKR_stand_Incidence();
            printf("%d%%\n", 100 * i / g_inpset->max_iter);
            printf("------------------------------\n");
            printf("Log-Posterior: %.8f \n", logPosterior);
            printf("------------------------------\n");
        }
    }
    printf("%d%%\n", 100);
    printf("------------------------------\n");
}

int ReadModelSIKR_stand_Incidence() {

    char filename[100];
    FILE *fptr_obs;
    FILE *fptr_initialpoint;
    int i = 0, ii = 0;
    int symbol;
    double num;

    kNumObs = 0;

    //open file with dataset
    sprintf(filename, "./benchmarks/%s_data.txt", g_inpset->data);

    fptr_obs = fopen(filename, "r");
    if((fptr_obs) == NULL){
        fprintf(stderr, "sir_model_synthetic.c - Failed to open the dataset file '%s'", filename);
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
    kNumParam = 3 + 1;
    Prior = PriorSIKR_stand_Incidence;
    Loglikelihood = LoglikelihoodSIKR_stand_Incidence;
    Gradient = GradientSIKR_stand_Incidence;
    Transf = TransfSIKR_stand_Incidence;
    InvTransf = InvTransfSIKR_stand_Incidence;
    LogJacobian = LogJacobianSIKR_stand_Incidence;
    JacobianMat = JacobianMatSIKR_stand_Incidence;
    Hessian = HessianSIKR_stand_Incidence;

    /** initialize state of model parameters ***/
    for (i = back2; i <= forw2; i++) {
        SAlloc(&(state[i]), kNumParam);
        SAlloc(&(temporalStatePlusH[i]), kNumParam);
        SAlloc(&(temporalStateMinusH[i]), kNumParam);
    }

    ///////////////////////////////
    //initialize position variables

    //check if initialpoint file exists
    int bool_initialpointfile_exists = -1;
    sprintf(filename, "./output/%s/initialpoint_%s.txt", id_global, id_global);
    fptr_initialpoint = fopen(filename, "r");
    if ((fptr_initialpoint) == NULL) {
        bool_initialpointfile_exists = 0;
    } else{
        bool_initialpointfile_exists = 1;
    }

    fprintf(red_lfp, "bool_initialpointfile_exists is %i\n", bool_initialpointfile_exists);

    //depending on inputfile generate a random starting point or read in from file
    if(bool_initialpointfile_exists == 0){
        fprintf(stderr, "Error in SIR model input: initialpoints file not found. \n");
        exit(EXIT_FAILURE);
    } else if(bool_initialpointfile_exists == 1){ //read in from file

        // count number of spaces in the first line
        int numSizeInitialPoint = 0;
        symbol = 1; //initialise symbol
        while(!feof(fptr_initialpoint) && symbol != '\n'){
            symbol = fgetc(fptr_initialpoint);
            if(isspace(symbol)){ numSizeInitialPoint++; }
        }

        fprintf(red_lfp, "numSizeInitialPoint is %i\n", numSizeInitialPoint);

        if(kNumParam != numSizeInitialPoint){
            fprintf(stderr, "sir_model_synthetic.c - kNumParam (%d) is not equal to numSizeInitialPoint (%i) in sir_model_synthetic.c. Exiting...\n", kNumParam, numSizeInitialPoint);
            exit(EXIT_FAILURE);
        }

        //go to beginning of file
        rewind(fptr_initialpoint);

        //read initial points in one by one
        for (ii = 0; ii < numSizeInitialPoint; ii++) {
            if (fscanf(fptr_initialpoint, "%lf", &num) != EOF) {
                state[curr]->pos[ii] = num;
            }
        }

        fclose(fptr_initialpoint);

        fprintf(red_lfp, "initialpoints are ");
        for (i = 0; i < numSizeInitialPoint; i++) {
            fprintf(red_lfp, "%lf ", state[curr]->pos[i]);
        }
        fprintf(red_lfp, "\n");

    } else{
        fprintf(stderr, "sir_stand_model_synthetic.c - If-else not matched for initial position\n");
        exit(EXIT_FAILURE);
    }
    // end initialize position variables

    // Perform gradient descent if enabled
    if (g_inpset->do_grad_desc == 1) {
        // Transform parameters to the real line before gradient descent
        TransfSIKR_stand_Incidence();

        // Perform gradient descent and get positions
        GradientDescentSIKR_stand_Incidence(); // positions

        // Transform parameters back from the real line
        InvTransfSIKR_stand_Incidence();

        // Save the final position from gradient descent to the starting point file
        sprintf(filename, "./output/%s/startingpoint_%s.txt", id_global, id_global);
        FILE *fp_startingpoint = fopen(filename, "w");
        if (fp_startingpoint == NULL) {
            fprintf(stderr, "ReadModelSIKR_stand_Incidence - Failed to open starting point file '%s'\n", filename);
            exit(EXIT_FAILURE);
        }
        for (ii = 0; ii < kNumParam; ii++) {
            fprintf(fp_startingpoint, "%lf\n", state[curr]->pos[ii]);
        }
        fclose(fp_startingpoint);
    }

    /*** initialize momenta vectors from a gaussian distr ****/
    if (state[curr]->mom == NULL) {
        printf("sir_stand_model_synthetic.c - Error in allocating memory (state mom)");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < kNumParam; i++) {
        state[curr]->mom[i] = gsl_ran_gaussian_ziggurat(g_rng, 1.);
    }

    /***** END initialization ******/
    Ham = (TypeHamiltonian *)malloc(sizeof(TypeHamiltonian));

    return 0;
}
