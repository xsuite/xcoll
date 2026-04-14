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
    // 0 get lengths
    // 1 compare lengths
    // 2 either return for MCS or do nuclear interaction
    // 3 if nuclear interaction, ionisation loss, get remanining lengths and do nucl int
    double interaction_length_tot;
    double cross_section_tot;
    double sqrt_t_p;
    int64_t i_slot = -1;
    // double A          = MaterialData_get__A(material);
    double N          = MaterialData_get__atoms_per_volume(material);
    double Z          = sqrt(MaterialData_get__Z2_eff(material));
    double X0         = MaterialData_get__radiation_length(material);
    everest->ecmsq    = 2*XC_PROTON_MASS*1.0e-3*pc;
    double sqrt_s      = sqrt(everest->ecmsq);

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t scatter               = everest->coll->record_scatterings;

    // I think we handled A inside the CS, so this is only CS
    cross_section_tot  = MaterialData_evaluate_glauber_spline(material, sqrt_s, 0); // [mb]
    interaction_length_tot = (1.)/(N*cross_section_tot*1.0e-27); // [m]
    //double mcs_path_length = FindRoot_get_path_length(finder);
    printf("sqrt(s) = %e GeV, total cross section = %e mb, interaction length = %e m\n", sqrt_s, cross_section_tot, interaction_length_tot);
    double mcs_path_length = length;
    double P_int = 1.0 - exp(-mcs_path_length / interaction_length_tot);
    if (RandomUniform_generate(part) > P_int) {
        // Ionisation loss to interaction point
        calculate_ionisation_properties(everest, material, pc);
        pc = calcionloss(everest, material, part, mcs_path_length, pc, 1);
        return pc; // false for nucl int.
    // if ( (mcs_path_length - interaction_length_tot) < 1e-12) {
    //     // MCS to exit
    //     // Ionisation loss to interaction point
    //     calculate_ionisation_properties(everest, material, pc);
    //     pc = calcionloss(everest, material, part, mcs_path_length, pc, 1);
    //     return pc; // false for nucl int.
    } else {
        double interaction_lengths[6];
        double nuclear_slope = MaterialData_get__nuclear_slope(material);
        double Neff = MaterialData_get__num_nucleons_eff(material);
        // double theta_init = atan2(MultipleCoulombTrajectory_get_tan_t0( (MultipleCoulombTrajectory) LocalTrajectory_member(traj) ), 1);
        double theta_init = (13.6e-3 / pc) * sqrt(length / X0) * (1.0 + 0.038 * log(length / X0));

        // TESTING> THIS IS OLD 
        // get_interaction_length(material, interaction_lengths, cross_section_tot, Z, N, Neff, theta_init, sqrt_s, pc);
        // // 1 = inel, 2 = el, 3 = prod, 4 = sd, (5 = pp/pn), 5(6) = coulomb
        // printf("Interaction lengths (m): Inel: %e, El: %e, Prod: %e, SD: %e, pp/pn: %e, Coulomb: %e\n", 
        //         interaction_lengths[0], interaction_lengths[1], interaction_lengths[2], interaction_lengths[3], interaction_lengths[4], interaction_lengths[5]);
        // // Finding the smallest length
        // double min_length = interaction_lengths[0]; // Inelastic
        // int min_index = 1;

        // if (interaction_lengths[1] < min_length) { // Elastic
        //     min_length = interaction_lengths[1];
        //     min_index = 2;
        // }
        // if (interaction_lengths[2] < min_length) { // Production
        //     min_length = interaction_lengths[2];
        //     min_index = 3;
        // }
        // if (interaction_lengths[3] < min_length) { // Single diffractive
        //     min_length = interaction_lengths[3];
        //     min_index = 4;
        // }
        // // if (interaction_lengths[4] < min_length) { // Proton-proton / proton-neutron
        // //     min_length = interaction_lengths[4];
        // //     min_index = 5;
        // // }
        // if (interaction_lengths[4] < min_length) { // Coulomb
        //     min_length = interaction_lengths[4];
        //     min_index = 5;
        // }
        // printf("Chosen interaction index: %d with length %e m\n", min_index, min_length);
// TESTING OLD ENDS HERE

        // Draw a uniform random number and select process by weight
        double fractions[5];
        for (int i = 1; i < 5; i++) {
            fractions[i] = MaterialData_evaluate_glauber_spline(material, sqrt_s, i); // [mb]
        }

        double r = RandomUniform_generate(part) * cross_section_tot;
        int chosen = 0;
        double cumulative_sum = 0.0;
        for (int i = 1; i < 5; i++) {
            cumulative_sum += fractions[i];
            if (r < cumulative_sum) { 
                chosen = i;    
                break; 
            }
        }
        printf("Chosen interaction index: %d\n", chosen);
        double min_index = chosen; // 1 = inel, 2 = el, 3 = prod, 4 = sd, (5 = pp/pn), 5(6) = coulomb
        // Ionisation loss to interaction point
        calculate_ionisation_properties(everest, material, pc);
        // pc = calcionloss(everest, material, part, min_length, pc, 1);
        pc = calcionloss(everest, material, part, 1./(fractions[chosen]*N), pc, 1);
        // Doing the nuclear interaction
        if (min_index == 1 || min_index == 3){
            // Inelastic or Production: particle is absorbed
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
            LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
            pc = 1.e-9; 
            sqrt_t_p = 0;

        } else if (min_index == 2){
            // Elastic
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_PN_ELASTIC);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/nuclear_slope)/pc;

        } else if (min_index == 4){
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
                sqrt_t_p = sqrt(RandomExponential_generate(part)/nuclear_slope)/sqrt(pc_in*pc);
            }

        // } else if (min_index == 5){
        //     // Proton-proton / proton-neutron
        //     if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_PP_ELASTIC);
        //     sqrt_t_p = sqrt(RandomExponential_generate(part)/nuclear_slope)/pc;

        } else if (min_index == 5){
            // Coulomb
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_COULOMB);
            sqrt_t_p = sqrt(RandomRutherford_generate(everest->coll->rng, part))/pc;

        } else {
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