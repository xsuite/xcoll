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
void get_coulomb_cross_section(double Z, double pc, double N,
                               double theta_rms, double nuclear_slope, double KE, double* cs_coulomb){
    double hbar_c         = 0.197;  // [GeV*fm] 
    double constant = (1./137. * hbar_c * 1*Z)/(KE); // fm
    *cs_coulomb = 10*(M_PI/4. * constant * constant * (1 + cos(theta_rms))/(1 - cos(theta_rms))); //in fm2, so convert to mb * 10

}


// =======================================================
// ====== Nuclear interaction length =====================
// =======================================================

/*gpufun*/
void get_interaction_length(MaterialData restrict material, double interaction_lengths[4], 
                            double cs_tot, double Z, double N, double sqrt_s, double pc, double particle_id) {

    double cs_prod_hA    = MaterialData_evaluate_glauber_spline(material, sqrt_s, 1, particle_id); // production
    double cs_el_hA      = MaterialData_evaluate_glauber_spline(material, sqrt_s, 2, particle_id); // elastic nucleus
    double cs_el_nucleon = MaterialData_evaluate_glauber_spline(material, sqrt_s, 3, particle_id); // el nucleon
    double cs_sd_hA      = MaterialData_evaluate_glauber_spline(material, sqrt_s, 4, particle_id); // single diffractive
    printf("cross sections = %f, %f, %f, %f\n", cs_prod_hA, cs_el_hA, cs_el_nucleon, cs_sd_hA); // --- IGNORE ---

     // Production
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