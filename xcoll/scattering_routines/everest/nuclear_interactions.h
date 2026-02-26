// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H
#define XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H

/*gpufun*/
double do_nuclear_interaction(EverestData restrict everest, LocalParticle* part, FindRoot finder, 
                              LocalTrajectory traj, MaterialData restrict material, double pc) {
    // 0 get lengths
    // 1 compare lengths
    // 2 either return for MCS or do nuclear interaction
    // 3 if nuclear interaction, get remanining lengths and do nucl int
    double interaction_length_tot;
    double A          = MaterialData_get__atomic_mass(material);
    double molar_mass = MaterialData_get__molar_mass(material);
    double Z          = MaterialData_get__atomic_number(material);
    double rho        = MaterialData_get__density(material);
    everest->ecmsq    = 2*XC_PROTON_MASS*1.0e-3*pc;
    double sqrt_s      = sqrt(everest->ecmsq);

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t scatter               = everest->coll->record_scatterings;

    total_cross_section(material, &interaction_length_tot, &A, &molar_mass, &rho, &sqrt_s);
    FindRoot_find_path_length(finder, traj);
    double mcs_path_length = FindRoot_get_path_length(finder);

    if ( (mcs_path_length - interaction_length_tot) < 1e-12) {
        // MCS to exit
        return 0; // false for nucl int.

    } else {
        // Nuclear interaction dominates: return path length for nuclear interaction and do the interaction
        double interaction_lengths[6];
        double theta_init = atan2(MultipleCoulombTrajectory_get_tan_t0( (MultipleCoulombTrajectory) LocalTrajectory_member(traj) ), 1);

        // Get interaction lengths for all types of interactions, to find the dominant one
        get_interaction_length(everest, material, pc, &theta_init, interaction_lengths, sqrt_s, interaction_length_tot);
        // 1 = inel, 2 = el, 3 = prod, 4 = sd, 5 = pp/pn, 6 = coulomb

        // Finding the smallest length
        double min_length = interaction_lengths[0];
        int min_index = 1;

        if (interaction_lengths[1] < min_length) {
            min_length = interaction_lengths[1];
            min_index = 2;
        }
        if (interaction_lengths[2] < min_length) {
            min_length = interaction_lengths[2];
            min_index = 3;
        }
        if (interaction_lengths[3] < min_length) {
            min_length = interaction_lengths[3];
            min_index = 4;
        }
        if (interaction_lengths[4] < min_length) {
            min_length = interaction_lengths[4];
            min_index = 5;
        }
        if (interaction_lengths[5] < min_length) {
            min_length = interaction_lengths[5];
            min_index = 6;
        }

        // Doing the nuclear interaction
        if (min_index == 1 || min_index == 3){
            // Inelastic or Production: particle is absorbed
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_LOST_ON_EVEREST_COLL);
            LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
            pc = 1.e-9; 
            sqrt_t_p = 0;

        } else if (min_index == 2){
            // Elastic
            double b_nuclear_elastic;
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_PN_ELASTIC);
            get_slope_hadron_nucleus(&A, &b_nuclear_elastic);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/b_nuclear_elastic)/pc;

        } else if (min_index == 4){
            // Single diffractive
            double pc_in = pc;
            pc = pc*(1 - xm2/everest->ecmsq);
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_SINGLE_DIFFRACTIVE);

            if (pc <= 1.e-9 || pc != pc) {
                // Very small (<1eV) or NaN
                if (sc) InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
                LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
                pc = 1.e-9; 
                sqrt_t_p = 0;
            } else {
                get_slope_single_diffraction(everest->ecmsq, &b_sd);
                sqrt_t_p = sqrt(RandomExponential_generate(part)/b_sd)/sqrt(pc_in*pc);
            }

        } else if (min_index == 5){
            // Proton-proton / proton-neutron
            double b_pp_pn;
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_PP_ELASTIC);
            get_slope_proton_proton(everest->ecmsq, &b_pp_pn);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/b_pp_pn)/pc;

        } else if (min_index == 6){
            // Coulomb
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_COULOMB);
            sqrt_t_p = sqrt(RandomRutherford_generate(everest->coll->rng, part))/pc;

        } else {
            // Unsupported interaction type
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

            if (sc) InteractionRecordData_log_child(record, i_slot, part);
    }
    return pc;
}