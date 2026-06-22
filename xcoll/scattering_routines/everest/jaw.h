// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_JAW_H
#define XCOLL_EVEREST_JAW_H
#include <math.h>
#include <stdio.h>

// This function needs to be rewritten to work with the geometry. 
// Note that with this version (to be compliant with the old code),
// we sample length step from the interaction length. In the new code that 
// should be removed such that we compare path length of MCS with sampled total interaction length to decide.
// The nuclear interaction function has this logic commented out; so it shall ALWAYS be called in the new framework.
// See thesis Simone.


/*gpufun*/
double jaw(EverestData restrict everest, MaterialData restrict material, LocalParticle* part,
           double pc, double length, int edge_check) {
    if (LocalParticle_get_state(part) < 1){
        // Do nothing if already absorbed
        return pc;
    }

    pc /= 1.e9; // [GeV]
    double rlen = length;
    double s0 = LocalParticle_get_s(part);
    double steps = 0;
    double theta_cut;
    while (1) {
        // Length of the step until nuclear interaction
        double length_step;

            double cs_coulomb;
            double nuclear_slope;
            double N           = MaterialData_get__atoms_per_volume(material);
            double Z           = sqrt(MaterialData_get__Z2_eff(material));
            double X0          = MaterialData_get__radiation_length(material);
            double M           = MaterialData_get__molar_mass(material) * 0.931494103;

            double KE = sqrt((M*M + 1) * (1e-3*XC_PROTON_MASS)*(1e-3*XC_PROTON_MASS) + 2*M*1e-3*XC_PROTON_MASS*pc) - (M+1)*(1e-3*XC_PROTON_MASS); // [GeV]
            double particle_id = 0;
            everest->ecmsq     = 2*XC_PROTON_MASS*1.0e-3*pc;
            double sqrt_s      = sqrt(everest->ecmsq);
            theta_cut_from_step(&theta_cut, rlen, X0, pc);
            get_coulomb_cross_section(Z, rlen, theta_cut, KE, &cs_coulomb); // [mb]
            double cross_section_tot  = MaterialData_evaluate_glauber_spline(material, sqrt_s, 0, particle_id); //+ cs_coulomb; // [mb]
            length_step = RandomExponential_generate(part) *(1.)/(N*cross_section_tot*1.0e-31); // [m]

        if (length_step > rlen) {
            // Length to nuclear interaction is longer than remaining: MCS to end and exit collimator
            mcs(everest, material, part, rlen, pc, edge_check);
            break;
        }

        mcs(everest, material, part, length_step, pc, edge_check);
        if (LocalParticle_get_state(part) < 1 || (edge_check && LocalParticle_get_x(part) <= 0)){
            // Particle lost all energy due to ionisation, or left the collimator
            break;
        }
        pc = do_nuclear_interaction_and_ionisation_loss(everest, part, length_step, material, pc, theta_cut);
        if (LocalParticle_get_state(part) < 1){
            // Particle was absorbed
            break;
        }
        // Calculate the remaining interaction length and close the iteration loop.
        rlen = rlen - length_step;
        steps++;
    }
    return pc*1e9;  // Back to eV
}

#endif /* XCOLL_EVEREST_JAW_H */
