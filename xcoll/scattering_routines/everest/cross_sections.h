// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CROSS_SECTIONS_H
#define XCOLL_EVEREST_CROSS_SECTIONS_H

#include <stddef.h>
#include <math.h>

// =======================================================
// ====== Coulomb Cross Sections =====================
// =======================================================
/*gpufun*/
void theta_cut_from_step(double* theta_cut, double length_step, double X0, double pc) {
    // Lynch-Dahl cutoff
    double theta_ld = 2.326 * (13.6e-3/pc * sqrt(length_step/X0) * (1 + 0.038*log(length_step/X0)));

    if (theta_ld <= 1e-6) {
        *theta_cut = 5e-6;
    } else {
        *theta_cut = theta_ld;
    }
}

/*gpufun*/
void get_coulomb_cross_section(double Z, double length_step,
                               double theta_cut, double KE, double* cs_coulomb){
    double hbar_c         = 0.1973269804 ;  // [GeV*fm] 
    double constant = (1./137. * hbar_c * 1*Z)/(KE); // fm  // change 1 so it matches the charge number of the particle.

    *cs_coulomb = 10*(M_PI/4. * constant * constant * (1 + cos(theta_cut))/(1 - cos(theta_cut))); //in fm2, so convert to mb * 10s
}


// =======================================================
// ====== Nuclear interaction length =====================
// =======================================================

/*gpufun*/
void get_interaction_length(MaterialData restrict material, double interaction_lengths[4], 
                            double N, double sqrt_s, double particle_id) {

    double cs_prod_hA    = MaterialData_evaluate_glauber_spline(material, sqrt_s, 1, particle_id); // production
    double cs_el_hA      = MaterialData_evaluate_glauber_spline(material, sqrt_s, 2, particle_id); // elastic nucleus
    double cs_el_nucleon = MaterialData_evaluate_glauber_spline(material, sqrt_s, 3, particle_id); // el nucleon
    double cs_sd_hA      = MaterialData_evaluate_glauber_spline(material, sqrt_s, 4, particle_id); // single diffractive

    // Production
    interaction_lengths[0] = (1)/(N*cs_prod_hA*1.0e-31);   // [m]

    // Elastic Nucleus: Total - Inelastic
    if (cs_el_hA < 1e-15){
        cs_el_hA = 1e-12; // In case. Makes Lambda large
        interaction_lengths[1] = (1)/(N*cs_el_hA*1.0e-31);
    } else {
        interaction_lengths[1] = (1)/(N*cs_el_hA*1.0e-31); // [m]
    }

    // Elastic nucleon
    interaction_lengths[2] = 1/(N*cs_el_nucleon*1.0e-31);   // [m]
    // Single diffractive
    interaction_lengths[3] = (1)/(N*cs_sd_hA*1.0e-31);      // [m]
}

#endif // XCOLL_EVEREST_CROSS_SECTIONS_H