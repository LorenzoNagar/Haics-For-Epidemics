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

#include <string.h>
#include <math.h>

#include "Globals.h"
#include "Definitions.h"

// scale user given time step
// scale user given time step
double scale_timestep(double ss) {
  double ss_scaled = 0.;

  ss_scaled = g_inpset->scaling_value * ss;

  return ss_scaled;
}

double FindIntegrator2stage(double *b_coeff, double ss) {
  double ss_scaled, b;
  int ind;

  ss_scaled = scale_timestep(ss);
  
  if (ss_scaled >= 4){
    b = 0.25;
    fprintf(red_lfp, "WARNING: Scaled timestep (%f) exceeds stability limit. b set to 0.25 (VV).\n", ss_scaled);
  } else if (ss_scaled > 0){

    ind = round(ss_scaled * 10000);
    
    //safety checks
    if(ind > hbar_to_b_num_disc_steps){
      ind = hbar_to_b_num_disc_steps;
    }else if(ind == 0){
      ind = 1;
    }

    //assign b
    b = b_coeff[ind-1];
    
  } else{
    perror("adapt_integrator.c - There was a problem in the else-if clauses while assigning the coefficient b in the FindIntegrator2stage function.\n");
    exit(EXIT_FAILURE);
  }

  return b;
}

double FindIntegrator3stage(double *b_coeff, double ss) {
  double ss_scaled, b;
  int ind;

  ss_scaled = scale_timestep(ss);
  
  if (ss_scaled >= 6){
    b = 1.0f/6.0f;
    fprintf(red_lfp, "WARNING: Scaled timestep (%f) exceeds stability limit. b set to 1/6 (Strang).\n", ss_scaled);
  } else if (ss_scaled > 0){
    
    ind = round(ss_scaled * 10000);

    //safety checks
    if(ind > hbar_to_b_num_disc_steps_3s){
      ind = hbar_to_b_num_disc_steps_3s;
    }else if(ind == 0){
      ind = 1;
    }

    //assign b
    b = b_coeff[ind-1];
    
  } else{
    perror("adapt_integrator.c - There was a problem in the else-if clauses while assigning the coefficient b in the FindIntegrator3stage function.\n");
    exit(EXIT_FAILURE);
  }

  return b;
}

double FindIntegratorAIA(double *b_coeff, double ss, int num_integrator_stages) {
  double ss_scaled, b;
  int ind;

  ss_scaled = sqrt(2) * ss;
  
  if (ss_scaled >= 2 * num_integrator_stages){
    b = 0.25;
    fprintf(red_lfp, "WARNING: Scaled timestep (%f) exceeds stability limit. b set to 0.25 (VV).\n", ss_scaled);
  } else if (ss_scaled > 0){

      ind = round(ss_scaled * 10000);

    //safety checks
    if(ind > hbar_to_b_num_disc_steps){
      ind = hbar_to_b_num_disc_steps;
    }else if(ind == 0){
      ind = 1;
    }

    //assign b
    b = b_coeff[ind-1];
    
  } else{
    perror("adapt_integrator.c - There was a problem in the else-if clauses while assigning the coefficient b in the FindIntegratorAIA function.\n");
    exit(EXIT_FAILURE);
  }

  return b;
}
