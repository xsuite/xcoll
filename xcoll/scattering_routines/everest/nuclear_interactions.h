// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H
#define XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H
#include <math.h>
// /*gpufun*/
// double get_particle_id(double pdg_id) {
//     if (pdg_id == 2212) {
//         return 0; // proton
//     } else if (pdg_id == 211) {
//         return 4; // pion plus
//     } else if (pdg_id == 321) {
//         return 2; // kaon plus
//     } else if (pdg_id == -211) {
//         return 3; // pion minus
//     } else if (pdg_id == -321) {
//         return 1; // kaon minus
//     } else {
//         return -1; // unknown particle
//     }
// }
/*gpufun*/
double do_nuclear_interaction_and_ionisation_loss(EverestData restrict everest, LocalParticle* part, double length,// FindRoot finder, 
                                                  MaterialData restrict material, double pc){
                                                //   LocalTrajectory traj, MaterialData restrict material, double pc) {

    double interaction_length_tot;
    double cs_coulomb;
    double nuclear_slope;
    double cross_section_tot;
    double sqrt_t_p;
    int64_t i_slot = -1;

    double N           = MaterialData_get__atoms_per_volume(material);
    double Z           = sqrt(MaterialData_get__Z2_eff(material));
    double X0          = MaterialData_get__radiation_length(material);
    // double particle_id = get_particle_id(LocalParticle_get_pdg_id(part));
    double particle_id = 0;
    everest->ecmsq     = 2*XC_PROTON_MASS*1.0e-3*pc;
    double sqrt_s      = sqrt(everest->ecmsq);
    double theta_init = (13.6e-3 / pc) * sqrt(length / X0) * (1.0 + 0.038 * log(length / X0)); // add random part 
    // double theta_init = atan2(MultipleCoulombTrajectory_get_tan_t0( (MultipleCoulombTrajectory) LocalTrajectory_member(traj) ), 1);
 
    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t scatter               = everest->coll->record_scatterings;

    if (particle_id == 3 || particle_id == 4){
        nuclear_slope = MaterialData_get__nuclear_slope_pion(material);
    } else {
        nuclear_slope = MaterialData_get__nuclear_slope(material);
    }
    get_coulomb_cross_section(Z, pc, N, theta_init, nuclear_slope, &cs_coulomb); // [mb]
    cross_section_tot  = MaterialData_evaluate_glauber_spline(material, sqrt_s, 0, particle_id) + cs_coulomb; // [mb]
    interaction_length_tot = RandomExponential_generate(part) *(1.)/(N*cross_section_tot*1.0e-31); // [m]
    //double mcs_path_length = FindRoot_get_path_length(finder);

    if ((length - (interaction_length_tot)) < 1e-12) {
        // MCS to exit
        calculate_ionisation_properties(everest, material, pc);
        pc = calcionloss(everest, material, part, length, pc, 1);
        return pc;
    } else {
        double interaction_lengths[5];
        get_interaction_length(material, interaction_lengths, cross_section_tot, Z, N, theta_init, sqrt_s, pc, particle_id);
        // 1: Prod, 2: elastic nucleus, 3: elastic nucleon, 4: single diffractive, 5: Coulomb
        int chosen                    = 1;
        double min_length             = (RandomExponential_generate(part) * interaction_lengths[0]);
        double elastic_length         = (RandomExponential_generate(part) * interaction_lengths[1]);
        double elastic_nucleon_length = (RandomExponential_generate(part) * interaction_lengths[2]);
        double SD_length              = (RandomExponential_generate(part) * interaction_lengths[3]);
        double coulomb_length         = (RandomExponential_generate(part) * 1./(N*cs_coulomb*1.0e-31));
        if (elastic_length < min_length) {         // Elastic
            min_length = elastic_length;
            chosen = 2;
        }
        if (elastic_nucleon_length < min_length) { // Elastic nucleon 
            min_length = elastic_nucleon_length;
            chosen = 3;
        }
        if (SD_length < min_length) {              // Single diffractive
            min_length = SD_length;
            chosen = 4;
        }
        if (coulomb_length < min_length) {         // Coulomb
            min_length = coulomb_length;
            chosen = 5;
        }

        calculate_ionisation_properties(everest, material, pc);
        pc = calcionloss(everest, material, part, min_length, pc, 1); // should be along mcs traj, how

        if (chosen == 1){
            // Production: particle is absorbed
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
            LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
            pc = 1.e-9; 
            sqrt_t_p = 0;

        } else if (chosen == 2){
            // Elastic Nucleus
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_PN_ELASTIC);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/nuclear_slope)/pc;

        } else if (chosen == 3){
            // Elastic Nucleon
            double pp_new = 9.3 + 0.22 * log(everest->ecmsq) + 0.03*(log(everest->ecmsq))*(log(everest->ecmsq)); //TODO: EVEREST->bpp
            if (scatter) i_slot = InteractionRecordData_log(record, record_index, part, XC_PP_ELASTIC);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/pp_new)/pc;
 
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
                double bsd;
                double pp_new = 9.3 + 0.22 * log(everest->ecmsq) + 0.03*(log(everest->ecmsq))*(log(everest->ecmsq));
                if (xm2 < 2.) {
                    bsd = 2* pp_new;//everest->bpp;
                } else if (xm2 >= 2. && xm2 <= 5.) {
                    bsd = ((106.0 - 17.0*xm2)*pp_new)/36.0;
                } else {
                    bsd = (7*pp_new)/12.0;
                } // THIS IS THE REASON FOR THE TAILS say ok here
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