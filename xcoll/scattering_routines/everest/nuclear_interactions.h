// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H
#define XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H

double do_nuclear_interaction(EverestData restrict everest, FindRoot finder, LocalTrajectory traj, MaterialData restrict material, double pc) {
    // 0 get lengths
    // 1 compare lengths
    // 2 either return for MCS or do nuclear interaction
    // 3 if nuclear interaction, get remanining lengths and do nucl int

    double interaction_length_tot;
    double A = MaterialData_get__atomic_mass(material);
    double molar_mass = MaterialData_get__molar_mass(material);
    double rho = MaterialData_get__density(material);
    everest->ecmsq = 2*XC_PROTON_MASS*1.0e-3*pc;
    double ecmsq = everest->ecmsq;
    total_cross_section(material, &interaction_length_tot, &A, &molar_mass, &rho, &ecmsq);
    FindRoot_find_path_length(finder, traj);
    double mcs_path_length = FindRoot_get_path_length(finder);

    if ( (mcs_path_length - interaction_length_tot) < 1e-12) {
        // MCS to exit
        return 0; // false for nucl int.

    } else {
        // Nuclear interaction dominates: return path length for nuclear interaction and do the interaction
        double interaction_length_inel, interaction_length_el, interaction_length_prod, interaction_length_sd;
        double interaction_length_pp_pn;
        get_interaction_length(everest, &interaction_length_inel, &A, &molar_mass, &rho, &ecmsq, &interaction_length_tot, 2);
        get_interaction_length(everest, &interaction_length_el, &A, &molar_mass, &rho, &ecmsq, &interaction_length_tot, 3);
        get_interaction_length(everest, &interaction_length_prod, &A, &molar_mass, &rho, &ecmsq, &interaction_length_tot, 4);
        get_interaction_length(everest, &interaction_length_sd, &A, &molar_mass, &rho, &ecmsq, &interaction_length_tot, 6);
        get_interaction_length(everest, &interaction_length_pp_pn, &A, &molar_mass, &rho, &ecmsq, &interaction_length_tot, 7);

        double min_length = interaction_length_inel;
        int min_index = 1;  // 1 = inel, 2 = el, 3 = prod, 4 = sd, 5 = pp/pn

        if (interaction_length_el < min_length) {
            min_length = interaction_length_el;
            min_index = 2;
        }
        if (interaction_length_prod < min_length) {
            min_length = interaction_length_prod;
            min_index = 3;
        }
        if (interaction_length_sd < min_length) {
            min_length = interaction_length_sd;
            min_index = 4;
        }
        if (interaction_length_pp_pn < min_length) {
            min_length = interaction_length_pp_pn;
            min_index = 5;
        }
        if (min_index == 1){
            // Inelastic
        } else if (min_index == 2){
            // Elastic
        } else if (min_index == 3){
            // Production
        } else if (min_index == 4){
            // Single diffractive
        } else if (min_index == 5){
            // Proton-proton / proton-neutron
        }
    }
}