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

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>

#include "hmc.h"
#include "Globals.h"
#include "integrators.h"
#include "Definitions.h"
#include "utils.h"
#include "sir_seir_common_functions.h"
#include "sir_model_functions.h" // Cote: check whether this can be done in a different way (need function printODE() from sir_model.c)
#include "sir_functions.h"
#include "seir_functions.h"
#include "sir_stand_functions.h"

void Hamiltonian(double logpd, TypeState **st, int dimension, double stepsize, TypeHamiltonian *H, int iter);

//florian--wayne: metropolis random walk
void HamiltonianMetropolis(double logpd, TypeState **st, int dimension, double stepsize, TypeHamiltonian *H, int iter);

// pointer to the Hamiltonian function (shadow or true)
void (*HamiltonianF)(double /*logpd*/, TypeState **/*state*/, int /*dimension*/, double /*stepsize*/, TypeHamiltonian *H, int iter);

//florian--wayne: pointer to the Hamiltonian function to be called in the PMU step. This is different for Metropolis as the kinetic energy in the Metropolis RW is set to 0
void (*HamiltonianPMU)(double /*logpd*/, TypeState **/*state*/, int /*dimension*/, double /*stepsize*/, TypeHamiltonian *H, int iter);

static void MomFlip(double *mom, int dimension);
static void MomFlipNo(double *mom, int dimension);
/* momentum flip function */
static void (*MomentumFlipF)(double */*mom*/, int /*dimension*/);

void MetropolisTest(double dH, bool *acc, bool *acc_ed, bool *prev_acc, double *acc_prob_prev_inv, int *numRF);

void AlwaysAcceptTest(double dH, bool *acc, bool *acc_ed, bool *PREV_ACC, double *acc_prob_prev_inv, int *numRF); //added by Felix, 12.07.2018

static void (*AcceptanceTest)(double /*dH*/, bool */*acc*/, bool */*acc_ed*/, bool */*prev_acc*/, double */*accProbPrevInv*/, int */*numRF*/);

static void PMU(double phi, TypeState **st, double *mom_old, int dimension, double logpd, double ss, bool *Accepted, int *numAM, TypeHamiltonian *H, int iter);

static void (*MomentumUpdateF)(double /*phi*/, TypeState **st, double */*mom_old*/, int /*dimension*/, double /*logpd*/, double /*ss*/, bool */*mom_acc*/, int */*numAM*/, TypeHamiltonian *H, int iter);

static void NoWeights(TypeHamiltonian Ham);
static void (*CalculateWeights)(TypeHamiltonian Ham);

/* book keeping */
static int num_proposed, num_accepted, num_accepted_ed, num_prop_mom, num_acc_mom, num_red_flips;
static int num_acc_tune, iter_tuning = 0; // for tuning the time step
static double ar_tuning; //for tuning the time step
static int flag_dt_tuned = 0; // flag for checking if a stepsize has been tuned
static double dt_tuned; //tuned time step for the scaling burn-in
static bool ACCEPTED, ACCEPTED_ed, MOM_ACCEPTED, PREV_ACCEPTED;
static double acc_prob_prev_inv;
static bool flag = true; //for normalization time
static int iter;

//variables for normalization time
double t0 = 0.;
double t1 = 0.;
double time_iter = 0.;
double NormalizationTime = 0.;

//variables for TUNE method
double ar_last = 0.;

TypeState *state_old[5] = {NULL,NULL,NULL,NULL,NULL};

// Hamiltonian()
// mom[]: momenta
// dimension: number of elements in the momenta vector
void Hamiltonian(double logpd, TypeState **st, int dimension, double stepsize, TypeHamiltonian *H, int iter) {

  H->pot = -logpd;
  
  H->kin = 0.5*cblas_ddot(dimension, st[curr]->mom, 1, st[curr]->mom, 1);

  H->ham = H->pot + H->kin;

}

// florian--wayne:
/*************************************************************
 *     hamiltonian for metropolis RW
 *     kinetic energy is 0
 *     mom[]: momenta
 *     dimension: number of elements in the momenta vector
 ************************************************************/

void HamiltonianMetropolis(double logpd, TypeState **st, int dimension, double stepsize, TypeHamiltonian *H, int iter) {

  H->pot = -logpd;

  H->kin = 0;

  H->ham = H->pot;

}

//momentum flip - automatic
static void MomFlip(double *mom, int dimension) {

  double alpha = -1.0;
  cblas_dscal(dimension, alpha, mom, 1);

}

//momentum flip - no flip
static void MomFlipNo(double *mom, int dimension) {
}

//metropolis test
void MetropolisTest(double dH, bool *acc, bool *acc_ed, bool *prev_acc, double *acc_prob_prev_inv, int *numRF) {

  double exp_dH,
  r = -1;

  if (dH <= 0) {
    *acc = TRUE;
    *acc_ed = TRUE;
  } else {
    *acc_ed = FALSE;
    r = gsl_rng_uniform(g_rng);
    exp_dH = exp(-dH);
    if (r < exp_dH) {
      *acc = TRUE;
    } else {
      *acc = FALSE;
    }
  }
}

//AlwaysAcceptTest
void AlwaysAcceptTest(double dH, bool *acc, bool *acc_ed, bool *prev_acc, double *acc_prob_prev_inv, int *numRF) {

  *acc = TRUE;

} //end of function AlwaysAcceptTest

//partial momentum update step
static void PMU(double phi, TypeState **st, double *mom_old, int dimension, double logpd, double ss, bool *Accepted, int *numAM, TypeHamiltonian *H, int iter) {

  int i;
  double aa = 0;
  double r = -1;
  double rphi;
  double rphicomp;
  double u[dimension];
  int INCX = 1,
  INCY = 1;

  for (i = 0; i < dimension; i++){
    u[i] = gsl_ran_gaussian_ziggurat(g_rng, 1.);
  }

  //check if phi is adaptive
  if ((g_inpset->t_varphi == 4) || (g_inpset->t_varphi == 5) || (g_inpset->t_varphi == 6)) {
    aa = cblas_ddot(dimension, mom_old, INCX, mom_old, INCY) - cblas_ddot(dimension, u, INCX, u, INCY);
  }
  
  //asign value if phi is adaptive
  if (g_inpset->t_varphi == 4) { //adaptive maximal MD AR
    if (aa < 0) {
      phi = 1 - phi;
    } else if (aa == 0) {
      r = gsl_rng_uniform(g_rng);
      if (r < 0.5) {
        phi = 1 - phi;
      }
    }
  } else if (g_inpset->t_varphi == 5) { //adaptive minimal MD AR
    if (aa > 0) {
      phi = 1 - phi;
    } else if (aa == 0) {
      r = gsl_rng_uniform(g_rng);
      if (r < 0.5) {
        phi = 1 - phi;
      }
    }
  } else if (g_inpset->t_varphi == 6) {}//adaptive optimal MD AR
    
  /// perform the momentum update
  if ( DOUBLE_EQ(phi,1) ) {
    cblas_dcopy(dimension, u, INCX, st[curr]->mom, INCY);//mom=u
  } else {
    rphi = sqrt(phi);
    rphicomp = sqrt(1-phi);
    for (i = 0; i < dimension; i++){ st[curr]->mom[i] =  rphicomp * mom_old[i] + rphi * u[i]; }
  }
  
  *numAM += 1;
  *Accepted = TRUE;
  
  //florian--wayne: 
  HamiltonianPMU(logpd, st, dimension, ss, H, iter);

}

//no weights needed
static void NoWeights(TypeHamiltonian Ham) {

}

//hmc_update()
//Update a state of the Markov chain
static void HMCUpdate(int iter) {

  TypeHamiltonian Ham_curr;

  fflush(stdout);

  double logpd = 0.,
  logpd_old = 0.,
  deltaH = 0.;
  int i;

  ////////////////////////////
  //MODEL PARAMETERS

  logpd_old = logpd = deltaH = 0.0;
  
  // keep the old position, momenta, gradient
  CopyState(state_old[curr], state[curr], kNumParam);

  logpd_old = Prior() + Loglikelihood() + LogJacobian();

  num_prop_mom++;


  // printf("Position and Momentum for gamma prior momentum update: %f, %f\n", state[curr]->pos[1], state[curr]->mom[1]);
  //update momenta and calculate current Hamiltonian; Hamiltonian here does not inlcude Jacobian, as it is constant for the PMMC step
  MomentumUpdateF(varphi[iter], state, state_old[curr]->mom, kNumParam, logpd_old, stepsize[iter], &MOM_ACCEPTED, &num_acc_mom, Ham, iter);

  //store the current Ham into the global Ham; it used some values from the previous state of the Ham
  Ham_curr = *Ham;

  // keep the current momenta, needed in case of rejection
  memcpy(state_old[curr]->mom, state[curr]->mom, kNumParam*sizeof(double));

  fprintf(tsfp, "%f\n", stepsize[iter]);
  fprintf(tlfp, "%d\n", traj_length[iter]);

  // compute normalization time multipliying one integration time duration for L
  if(flag){
    t0 = clock();
  }

  // For fixed gamma in SIR-like models we need to supress the momentum from the Hamiltonian calculation.
  // For the logposterior it is handled internaly.
  if(g_inpset->gammaFixed == 1){
    if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0 || strcmp(g_inpset->model, "SIKR_Incidence") == 0){
      state[curr]->mom[0] = 0;
    }
    if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0){
      state[curr]->mom[1] = 0;
    }
  }
  // printf("Position and Momentum for gamma prior MDMove: %f, %f\n", state[curr]->pos[1], state[curr]->mom[1]);

  // run the trajectory / make the proposal
  MDMove(Gradient, kNumParam, state, stepsize[iter], &(g_coeff[iter][0]), traj_length[iter]);
  
  // Check if the proposed point fulfils that all states sum to N (the size of the poulation.)
  InvTransf();
  if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0){
    double R_compartment = num_pop - state[curr]-> pos[1] - state[curr]-> pos[2];
    if (R_compartment < 0){
      state[curr]-> pos[1] = num_pop - state[curr]-> pos[2];
    }
    if(R_compartment > num_pop){
      if(state[curr]-> pos[1] < 0){
        state[curr]-> pos[1] = 0;
      } else if (state[curr]-> pos[2] < 0){
        state[curr]-> pos[2] = 0;
      }    
    }
    printf("Proposed step R0: %e\n", num_pop - state[curr]-> pos[1] - state[curr]-> pos[2]);
  } else if (strcmp(g_inpset->model, "SIKR_Incidence") == 0){
    double R_compartment = num_pop - state[curr]-> pos[1] - state[curr]-> pos[2];
    if (R_compartment < 0){
      state[curr]-> pos[1] = num_pop - state[curr]-> pos[2] ;
    }
    if(R_compartment > num_pop){
      if(state[curr]-> pos[1] < 0){
        state[curr]-> pos[1] = 0;
      } else if (state[curr]-> pos[2] < 0){
        state[curr]-> pos[2] = 0;
      }    
    }
    printf("Proposed step R0: %e\n", num_pop - state[curr]-> pos[1] - state[curr]-> pos[2]);
  } else if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0){
    double R_compartment = num_pop - state[curr]-> pos[2] - state[curr]-> pos[3] - state[curr]-> pos[4];
    if (R_compartment < 0){
      state[curr]-> pos[2] = num_pop - state[curr]-> pos[3] - state[curr]-> pos[4];
    }
    if(R_compartment > num_pop){
      if(state[curr]-> pos[2] < 0){
          state[curr]-> pos[2] = 0;
      } else if (state[curr]-> pos[3] < 0){
          state[curr]-> pos[3] = 0;
      } else {
          state[curr]-> pos[4] = 0;
      }        
    }
    printf("Original R0: %e\n", R_compartment);
    printf("Proposed step R0: %e\n", num_pop - state[curr]-> pos[2] - state[curr]-> pos[3] - state[curr]-> pos[4]);
  }
  Transf();
  // printf("Position and Momentum for gamma after MDMove: %f, %f\n", state[curr]->pos[1], state[curr]->mom[1]);

  // compute normalization time part2
  if(flag){
    t1 = clock();
    time_iter = ((t1 - t0)/CLOCKS_PER_SEC/traj_length[iter])*(g_inpset->iter_sampling);
    if(time_iter != 0){  //in order to avoid NormalizationTime equal to 0
      flag = false;
    }
  }

  num_proposed++;

  logpd = Prior() + Loglikelihood() + LogJacobian();

  HamiltonianF(logpd, state, kNumParam, stepsize[iter], Ham, iter);//proposed Ham is stored in Ham global

  //deltaH = (Ham->ham + logJac) - (Ham_curr.ham + logJac_old);
  deltaH = (Ham->ham) - (Ham_curr.ham);
  if(deltaH < 0){
    num_proposals_dHneg_arr[iter] = 1; //store the number of proposals with negative energy error
  }else{
    num_proposals_dHneg_arr[iter] = 0;
  }

  // acceptance test
  AcceptanceTest(deltaH, &ACCEPTED, &ACCEPTED_ed, &PREV_ACCEPTED, &acc_prob_prev_inv, &num_red_flips);

  //HessianSIKR_Incidence();

  if (ACCEPTED) {
    num_accepted++;
    accept_arr[iter] = 1;
    if((flag_dt_tuned == 0) && (iter < g_inpset->iter_burn_in)){ //during the burn-in, AR is used both for tuning a time step and for calculating fitting factors
      num_acc_tune++;
    }
    if (ACCEPTED_ed){
      num_accepted_ed++;
    }
    CalculateWeights(*Ham);
    fprintf(hfp, "%f\n", Ham->ham);
  } else {
    accept_arr[iter] = 0;
    CopyState(state[curr], state_old[curr], kNumParam);
    logpd = logpd_old;
    *Ham = Ham_curr;
    MomentumFlipF(state[curr]->mom, kNumParam);
    CalculateWeights(*Ham);
    fprintf(hfp, "%f\n", Ham->ham);
  }

  // tuning the time step during the burn-in stage
  if(flag_dt_tuned == 0){
    if((iter+1)%(g_inpset->iter_tune) == 0){
      ar_tuning = (double)(num_acc_tune)/(g_inpset->iter_tune); // corrected acceptance rate
      printf("AR_TUNING: %f\n", ar_tuning);
      if(ar_tuning < g_inpset->AR_target - g_inpset->delta_AR_target){ //if the AR is smaller than the target, lower a stepsize
        for (i = iter + 1; i < g_inpset->iter_burn_in + g_inpset->iter_sampling; i = i + 1){
          stepsize[i] = stepsize[i]-stepsize[i]/10 * (g_inpset->iter_burn_in + g_inpset->iter_sampling)/(g_inpset->iter_burn_in + g_inpset->iter_sampling + iter + 1 - iter_tuning);
        }
      // the acceptance rate is restarted after a change in the time steps
      num_acc_tune = 0;
      iter_tuning = iter + 1;
      }else if(ar_tuning > g_inpset->AR_target + g_inpset->delta_AR_target){ //if the AR is bigger than the target, increase a stepsize
        for (i = iter + 1; i < g_inpset->iter_burn_in + g_inpset->iter_sampling; i = i + 1){
          stepsize[i] = stepsize[i] + stepsize[i]/10 * (g_inpset->iter_burn_in + g_inpset->iter_sampling)/(g_inpset->iter_burn_in + g_inpset->iter_sampling + iter + 1 - iter_tuning);
        }
      // the acceptance rate is restarted after a change in the time steps
      num_acc_tune = 0;
      iter_tuning = iter + 1;  
      }else{
        iter_tuning = iter + 1;
        flag_dt_tuned = 1;
        for (i = 0; i < kNumParam; i = i + 1){
          InvTransf();
          fprintf(final_point_tune_fp, "%f ", state[curr]->pos[i]);
          Transf();
        }
        //fprintf(final_point_tune_fp, "\n");
        dt_tuned = stepsize[iter + 1];
        for(i = 0; i < g_inpset->iter_burn_in + g_inpset->iter_sampling; i = i + 1){
          stepsize[i] = dt_tuned;
        }
        //print the last AR tuned
        printf("Tuned AR: %lf\n", ar_tuning);
        fprintf(red_lfp, "Tuned AR: %lf\n", ar_tuning);
        fprintf(artfp, "AR_tuned %lf\n", ar_tuning);

        //print the tuned dt
        printf("Tuned dt: %lf\n", dt_tuned);
        fprintf(red_lfp, "Tuned dt: %lf\n", dt_tuned);
        fprintf(artfp, "dt_tuned %lf\n", dt_tuned);

        //number of iterations needed for tuning
        printf("Iterations needed for tuning: %i\n", iter_tuning);
        fprintf(red_lfp, "Iterations needed for tuning: %i\n", iter_tuning);
        fprintf(artfp, "iter_tuning %i\n", iter_tuning);
      }
    }
  }

  fprintf(lpdfp, "%f\n", logpd);

} //end of function "HMCUpdate"

//hmc()
double HMC()
{

  int i;
  int ii;
  int tot_iter = g_inpset->iter_burn_in + g_inpset->iter_sampling;
  double AR;
  clock_t start0, end0, startT, endT;
  double elapsed0, elapsedT;

  ACCEPTED = ACCEPTED_ed = MOM_ACCEPTED = FALSE;
  PREV_ACCEPTED = FALSE;

  acc_prob_prev_inv = 0;

  num_proposed = 0;
  num_accepted = 0;
  num_accepted_ed = 0;
  num_prop_mom = 0;
  num_acc_mom = 0;
  num_red_flips = 0;

  // allocate memory for the old state keeping
  for (i = back2; i <= forw2; i++) {
    SAlloc(&(state_old[i]), kNumParam);
  }

  // setup of functions
  if (strcmp(g_inpset->method, "HMC") == 0){

    AcceptanceTest = MetropolisTest;
    MomentumFlipF = MomFlipNo;

  } else if (strcmp(g_inpset->method, "GHMC") == 0){
    
    AcceptanceTest = MetropolisTest;
    MomentumFlipF = MomFlip;

  }else{

    fprintf(stderr, "hmc.c - no match found in if-elseif for method. Exiting.");
    exit(EXIT_FAILURE);

  } //end of if-elseif for method

  CalculateWeights = NoWeights;

  //allocating memory for stepsize, traj_length, g_coeff, g_SH_coeff, varphi
  //stepsize = (double *)malloc(tot_iter*sizeof(double)); //allocate and init with zeroes
  //traj_length = (int *)malloc(tot_iter*sizeof(int)); //allocate and init with zeroes
  //g_coeff = MatrixAlloc(tot_iter, 8); //allocate and init with zeroes
  //varphi = (double *)malloc(tot_iter*sizeof(double)); //allocate and init with zeroes

  //allocation memory for acceptance/rejection arrays
  accept_arr = (int *)malloc(tot_iter*sizeof(int)); //allocate and init with zeroes
  num_proposals_dHneg_arr = (int *)malloc(tot_iter*sizeof(int)); //allocate and init with zeroes

  Transf();//transformation of the model parameters at the beginning
  start0 = clock();

  ////////////////////////////////////////////////////////////
  // burn-in phase

  HamiltonianF = Hamiltonian;
  HamiltonianPMU = Hamiltonian;
  MomentumUpdateF = PMU;
  
  ///////////////////////////////////////////////////////
  //assign quantities for burn-in phase
  //1. stepsize = 1/D, 2. number of integration steps L = 1, 3. integrator coefficients of 1-stage VV, 4. Shadow Hamiltonian coefficients, 5. phi chosen from the input

  if(g_inpset->iter_burn_in == 0){
    flag_dt_tuned = 1;
    printf("The number of burn-in iterations chosen is 0, only production stage.\n");
  }else{

    printf("------ Starting burn-in phase\n");
    fprintf(red_lfp, "------ Starting burn-in phase\n");

    double delta_t_start = fmin(1. / kNumParam, 0.002); // dt for starting the burn-in

    for (i = 0; i < g_inpset->iter_burn_in; i++) {
      stepsize[i] = delta_t_start;
      traj_length[i] = 1;
      g_coeff[i][0] = 0.5;
      g_coeff[i][1] = 1.0;
    }
  
    //phi assignation in GHMC
    fnAssignPhi(varphi, g_inpset->varphi, g_inpset->t_varphi, 0, g_inpset->iter_burn_in - 1); //choice of phi given in the inputfile, if HMC is set to be constant equal to 1 in the fnAssignPhi function in integrators.c

    printf("The first time step for the tuning is %lf\n", stepsize[0]);
    fprintf(red_lfp, "The first time step for the tuning is %lf\n", stepsize[0]);

    printf("The trajectory used for the tuning is %d\n", traj_length[0]);
    fprintf(red_lfp, "The trajectory used for the tuning is %d\n", traj_length[0]);

    //assign MDMove for the burn-in, i.e. 1-stage Velocity Verlet
    MDMove = VV;

    //tuning the time step
    while ((flag_dt_tuned == 0) && (iter < g_inpset->iter_burn_in)) {

      HMCUpdate(iter);

      iter++;
    }

    if(flag_dt_tuned == 0){
      fprintf(stderr, "hmc.c - Stepsize cannot be tuned with the burn-in iterations chosen in input. Exiting. Try both to increase the number of iterations and/or iter_tune.\n");
      fprintf(stderr, "The last AR obtained is: %lf\n", ar_tuning);
      fprintf(stderr, "The last dt used is: %lf\n", stepsize[iter-1]);
      exit(EXIT_FAILURE);
    }

    for (iter = 0; iter < g_inpset->iter_burn_in; iter++) {
      HMCUpdate(iter);

      //write trajectory file
      if ((iter+1)%(g_inpset->thinning) == 0) {
        //main trajectory
        for (i = 0; i < kNumParam; i++) { fprintf(tfp, "%f ", state[curr]->pos[i]); }
        fprintf(tfp, "\n");

        //SIKR Incidence
        if (strcmp(g_inpset->model, "SIKR_Incidence") == 0) {
          printODE_SIR(odefp);
        }
        //SEMIKR synthetic
        if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0) {
          printODE_SEIR(odefp);
        }
        //SIR standard Incidence
        if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0) {
          printODE_SIR_stand(odefp);
        }
      }
    }
    // end of the loop for the real burn-in phase

    end0 = clock();
    elapsed0 = ((double)(end0-start0))/CLOCKS_PER_SEC;

    //store the final point of the burn-in stage
    InvTransf();
    for (i = 0; i < kNumParam; i++) { fprintf(final_point_burnin_fp, "%f ", state[curr]->pos[i]); }
    Transf();
    fprintf(final_point_burnin_fp, "\n");

    //fitting factor calculation
    int accept_burnin = 0;
    for(ii = 0; ii < g_inpset->iter_burn_in; ii++){
      accept_burnin += accept_arr[ii];
    }

    int accept_dHneg_burnin = 0;
    for(ii = 0; ii < g_inpset->iter_burn_in; ii++){
      accept_dHneg_burnin += num_proposals_dHneg_arr[ii];
    }

    int accept_dHpos_burnin = accept_burnin - accept_dHneg_burnin;
    double AR_filt_new_burnin = (double)accept_dHpos_burnin/g_inpset->iter_burn_in;
    double AR_burnin = (double)accept_burnin/g_inpset->iter_burn_in;

    double exp_dH_filtAR = - log(AR_filt_new_burnin);
    double exp_dH_Gupta = 4 * PI * (1 - AR_burnin) * (1 - AR_burnin);
    double fitting_filtAR = fmax(1, 1 / dt_tuned * pow(32*exp_dH_filtAR/kNumParam, 1.0/6.0));
    double fitting_Gupta = fmax(1, 1 / dt_tuned * pow(32*exp_dH_Gupta/kNumParam, 1.0/6.0));
    double fitting_factor = 2 * fitting_Gupta * fitting_filtAR / (fitting_Gupta + fitting_filtAR);

    //print acceptance rate burn-in
    printf("Acceptance rate burn-in: %lf\n", AR_burnin);
    fprintf(red_lfp, "Acceptance rate burn-in: %lf\n", AR_burnin);
    fprintf(artfp, "AR_burnin %lf\n", AR_burnin);

    //print fitting factor to art.txt
    fprintf(red_lfp, "Fitting factor: %lf\n", fitting_factor);
    fprintf(artfp, "fitting_factor %lf\n", fitting_factor);

    //print HSL in 1-stage units to art.txt
    double hsl = 1./fitting_factor;
    fprintf(red_lfp, "Estimated HSL: %lf\n", hsl);
    fprintf(artfp, "hsl %lf\n", hsl);

    //print HSL filtAR in 1-stage units to art.txt
    double hsl_filtAR = 1./fitting_filtAR;
    fprintf(red_lfp, "Estimated HSL_filtAR: %lf\n", hsl_filtAR);
    fprintf(artfp, "hsl_filtAR %lf\n", hsl_filtAR);

    //print HSL Gupta in 1-stage units to art.txt
    double hsl_Gupta = 1./fitting_Gupta;
    fprintf(red_lfp, "Estimated HSL_Gupta: %lf\n", hsl_Gupta);
    fprintf(artfp, "hsl_Gupta %lf\n", hsl_Gupta);

    //time elapsed in burn-in
    printf("Elapsed time burnin: %lf sec\n", elapsed0);
    fprintf(red_lfp, "Elapsed time burnin: %lf\n", elapsed0);
    fprintf(artfp, "TimeBurnin %lf\n", elapsed0);
  }

  ///////////////////////////
  //start production phase loop

  if(g_inpset->iter_sampling == 0){
    printf("The number of production iterations chosen is 0, only burnin stage.\n");
  }else{

    printf("------ Starting Production phase (collecting posterior samples...)\n");
    fprintf(red_lfp, "------ Starting Production phase (collecting posterior samples...)\n");

    flag_dt_tuned = 1; //skip the tuning part in HMCUpdate
  
    //assign HMC parameters for sampling
    
    fnAssignStepsizes(stepsize, delta_t, g_inpset->t_stepsize, g_inpset->iter_burn_in, tot_iter - 1);
    fnAssignTrajectoryLengths(traj_length, &g_inpset->L, g_inpset->t_L, g_inpset->iter_burn_in, tot_iter - 1);
    fnAssignIntegratorCoeffs(g_coeff, g_inpset->iter_burn_in, tot_iter - 1, 0);
    fnAssignPhi(varphi, g_inpset->varphi, g_inpset->t_varphi, g_inpset->iter_burn_in, tot_iter - 1);

    //assign MDMove for the production stage

    //find number of stages as first chracter of g_inpset->integrator for assignment later
    char strNumStages[256];
    strncpy(strNumStages, g_inpset->integrator, 1);
    strNumStages[1] = 0; // null terminate destination

    if(strcmp(strNumStages, "1") == 0){ //1-stage integrators
      MDMove = VV;
    }else if(strcmp(strNumStages, "2") == 0){ //2-stage integrators
      MDMove = V2S;
    }else if(strcmp(strNumStages, "3") == 0){ //3-stage integrators
      MDMove = V3S;
    }else{
      fprintf(stderr, "HMC.c - No match found in cases for integrator for assigning MDMove. Exiting.");
      exit(EXIT_FAILURE);
    }

    startT = clock();

    num_proposed  = 0;
    num_accepted  = 0;
    num_prop_mom = 0;
    num_acc_mom = 0;
    num_red_flips = 0;

    for (iter = g_inpset->iter_burn_in; iter < tot_iter; iter++) {

      HMCUpdate(iter);

      //write trajectory file
      if ((iter+1)%(g_inpset->thinning) == 0) {
        //main trajectory
        for (i = 0; i < kNumParam; i++) { fprintf(tfp, "%f ", state[curr]->pos[i]); }
        fprintf(tfp, "\n");

        //SIKR Incidence
        if (strcmp(g_inpset->model, "SIKR_Incidence") == 0) {
          printODE_SIR(odefp);
        }
        //SEMIKR synthetic
        if (strcmp(g_inpset->model, "SEMIKR_Incidence") == 0) {
          printODE_SEIR(odefp);
        }
        //SIR standard Incidence
        if (strcmp(g_inpset->model, "SIR_standard_Incidence") == 0) {
          printODE_SIR_stand(odefp);
        }
      }
    }
    
    endT = clock();
    elapsedT = ((double)(endT-startT))/CLOCKS_PER_SEC;

    /////////////////////////
    //print results to files

    printf("------ Printing statistics.\n");

    fprintf(red_lfp, "------ Printing statistics to art file.\n");

    //average acceptance rate
    AR = (double)num_accepted/(double)num_proposed;
    printf("Acceptance rate: %lf\n", AR);
    fprintf(red_lfp, "Acceptance rate: %lf\n", AR);
    fprintf(artfp, "AR %lf\n", AR);

    //print total number of proposals
    printf("Num proposals: %i\n", num_proposed);
    fprintf(red_lfp, "Num proposals: %i\n", num_proposed);
    fprintf(artfp, "num_proposed: %i\n", num_proposed);

    //time elapsed total
    printf("Elapsed time: %lf sec\n", elapsedT);
    fprintf(red_lfp, "Elapsed time: %lf\n", elapsedT);
    fprintf(artfp, "TimeTotal %lf\n", elapsedT);

    //mean of integrator coefficient b1
    double mean_g_coeff_i_0 = 0;
    for (ii = 0; ii < tot_iter; ii++){ mean_g_coeff_i_0 += g_coeff[ii][0]; }
    mean_g_coeff_i_0 = mean_g_coeff_i_0 / tot_iter;
    printf("Average b: %lf\n", mean_g_coeff_i_0);
    fprintf(red_lfp, "Average b: %lf\n", mean_g_coeff_i_0);
    fprintf(artfp, "mean_b %lf\n", mean_g_coeff_i_0);

    //standard deviation of integrator coefficient b1
    double sd_g_coeff_i_0 = 0;
    for (ii = 0; ii < tot_iter; ii++){ sd_g_coeff_i_0 += pow(g_coeff[ii][0] - mean_g_coeff_i_0, 2.0); }
    sd_g_coeff_i_0 = sqrt(sd_g_coeff_i_0 / tot_iter);
    fprintf(red_lfp, "standard deviation b: %lf\n", sd_g_coeff_i_0);
    fprintf(artfp, "sd_b %lf\n", sd_g_coeff_i_0);

    //mean of angle phi
    double mean_phi = 0;
    for (ii = 0; ii < tot_iter; ii++){ mean_phi += varphi[ii]; }
    mean_phi = mean_phi / tot_iter;
    printf("Average angle phi: %lf\n", mean_phi);
    fprintf(red_lfp, "Average angle phi: %lf\n", mean_phi);
    fprintf(artfp, "mean_phi %lf\n", mean_phi);

    //standard deviation of angle phi
    double sd_phi = 0;
    for (ii = 0; ii < tot_iter; ii++){ sd_phi += pow(varphi[ii] - mean_phi, 2.0); }
    sd_phi = sqrt(sd_phi / tot_iter);
    fprintf(red_lfp, "standard deviation phi: %lf\n", sd_phi);
    fprintf(artfp, "sd_phi %lf\n", sd_phi);
  }
  //end of main sampling loop

  ///////////////////////////////////////////////
  //close opened file pointers
  fclose(tfp);
  fclose(lpdfp);
  fclose(hfp);
  fclose(tsfp);
  fclose(artfp);
  if(g_inpset->iter_burn_in > 0){
    fclose(final_point_tune_fp);
    fclose(final_point_burnin_fp);
  }
  fclose(odefp);

  ///////////////////////////////////////////////
  //free allocated memory
  for (i = back2; i <= forw2; i++) {
    SFree(&state[i]);
    SFree(&state_old[i]);
  }

  MatrixFree(g_obs);

  free(stepsize);
  free(traj_length);
  free(varphi);

  gsl_rng_free(g_rng);

  free(g_inpset);

  return 0;
}
