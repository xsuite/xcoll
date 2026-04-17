// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H
#define XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H
#include <math.h>

/*gpufun*/
double do_nuclear_interaction_and_ionisation_loss(EverestData restrict everest, LocalParticle* part, double length,// FindRoot finder, 
                                                  MaterialData restrict material, double pc){
                                                //   LocalTrajectory traj, MaterialData restrict material, double pc) {

    double interaction_length_tot;
    double cs_coulomb;
    double cross_section_tot;
    double sqrt_t_p;
    int64_t i_slot = -1;

    double N          = MaterialData_get__atoms_per_volume(material);
    double Z          = sqrt(MaterialData_get__Z2_eff(material));
    double X0         = MaterialData_get__radiation_length(material);
    everest->ecmsq    = 2*XC_PROTON_MASS*1.0e-3*pc;
    double sqrt_s      = sqrt(everest->ecmsq);
    double nuclear_slope = MaterialData_get__nuclear_slope(material);
    double theta_init = (13.6e-3 / pc) * sqrt(length / X0) * (1.0 + 0.038 * log(length / X0));
    // double theta_init = atan2(MultipleCoulombTrajectory_get_tan_t0( (MultipleCoulombTrajectory) LocalTrajectory_member(traj) ), 1);

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t scatter               = everest->coll->record_scatterings;

    get_coulomb_cross_section(Z, pc, N, theta_init, nuclear_slope, &cs_coulomb); // [mb]
    cross_section_tot  = MaterialData_evaluate_glauber_spline(material, sqrt_s, 0) + cs_coulomb; // [mb]
    interaction_length_tot = RandomExponential_generate(part) *(1.)/(N*cross_section_tot*1.0e-27); // [m]
    //double mcs_path_length = FindRoot_get_path_length(finder);

    if ((length - (interaction_length_tot)) < 1e-12) {
        // MCS to exit
        calculate_ionisation_properties(everest, material, pc);
        pc = calcionloss(everest, material, part, length, pc, 1);
        return pc;
    } else {
        double interaction_lengths[5];
        get_interaction_length(material, interaction_lengths, cross_section_tot, Z, N, theta_init, sqrt_s, pc);
        // // 1 = inel, 2 = el, 3 = prod, 4 = sd, (5 = pp/pn), 5(6) = coulomb
        double min_length = (RandomExponential_generate(part) * interaction_lengths[0]); // Inelastic
        int chosen = 1;
        double elastic_length = (RandomExponential_generate(part) * interaction_lengths[1]);
        double SD_length = (RandomExponential_generate(part) * interaction_lengths[3]);
        double coulomb_length = (RandomExponential_generate(part) * 1./(N*cs_coulomb*1.0e-27));
    
        if (elastic_length < min_length) { // Elastic
            min_length = elastic_length;
            chosen = 2;
        }
        if ((RandomExponential_generate(part) * interaction_lengths[2]) < min_length) { // Production
            min_length = (RandomExponential_generate(part) * interaction_lengths[2]);
            chosen = 3;
        }
        if (SD_length < min_length) { // Single diffractive
            min_length = SD_length;
            chosen = 4;
        }
        if (coulomb_length < min_length) { // Coulomb
            min_length = coulomb_length;
            chosen = 5;
        }

        calculate_ionisation_properties(everest, material, pc);
        pc = calcionloss(everest, material, part, min_length, pc, 1);
        // pc = calcionloss(everest, material, part, 1./(fractions[chosen]*N), pc, 1);

        if (chosen == 1 || chosen == 3){
            // Inelastic or Production: particle is absorbed
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
            LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
            pc = 1.e-9; 
            sqrt_t_p = 0;

        } else if (chosen == 2){
            // Elastic
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_PN_ELASTIC);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/nuclear_slope)/pc;

        } else if (chosen == 4){
            // Single diffractive
            double pc_in = pc;
            double xm2 = exp(RandomUniform_generate(part)*everest->xln15s);
            pc = pc*(1 - xm2/everest->ecmsq);
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_SINGLE_DIFFRACTIVE);

            if (pc <= 1.e-9 || pc != pc) {
                // Very small (<1eV) or NaN
                if (scatter) InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
                LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
                pc = 1.e-9; 
                sqrt_t_p = 0;
            } else {
                printf("Nuclear slope is: %e GeV^-2\n", nuclear_slope);
                double bsd;
                if (xm2 < 2.) {
                    bsd = 2*everest->bpp;
                } else if (xm2 >= 2. && xm2 <= 5.) {
                    bsd = ((106.0 - 17.0*xm2)*everest->bpp)/36.0;
                } else {
                    bsd = (7*everest->bpp)/12.0;
                } // THIS IS THE REASON FOR THE TAILS
                sqrt_t_p = sqrt(RandomExponential_generate(part)/bsd)/sqrt(pc_in*pc);
            }

        } else if (chosen == 5){
            // Coulomb
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_COULOMB);
            sqrt_t_p = sqrt(RandomRutherford_generate(everest->coll->rng, part))/pc;

        } else {
            printf("Error in nuclear interaction choice.\n");
            return pc; // No interaction, return original momentum
        }

        double tan_theta = sqrt_t_p * sqrt(1 - sqrt_t_p*sqrt_t_p/4)/(1 - sqrt_t_p*sqrt_t_p/2);
        double alpha = 2*M_PI*RandomUniform_generate(part);
        double tan_theta_x = tan_theta*cos(alpha);
        double tan_theta_y = tan_theta*sin(alpha);

        // Change the angles
        #ifdef XCOLL_USE_EXACT
                LocalParticle_add_to_exact_xp_yp(part, tan_theta_x, tan_theta_y);
        #else
                LocalParticle_add_to_xp_yp(part, tan_theta_x, tan_theta_y);
        #endif

        if (scatter) InteractionRecordData_log_child(record, i_slot, part);
    }
    return pc;
}
#endif /* XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H */