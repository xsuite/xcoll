// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H
#define XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H
#include <math.h>


// Allen & Hastings Approximation for the Exponential Integral function E1(x) = int_x^inf (e^-t / t) dt
/*gpufun*/
inline void E1_approx(double* E1_approx, double x) {
    if (x <= 0.0) {
        *E1_approx = 1e21;  // E1 undefined for x <= 0
        return;
    }
    if (x <= 1.0) {

        // Small x expansion
        const double a0 = -0.57722;
        const double a1 =  0.99999;
        const double a2 = -0.24991;
        const double a3 =  0.05519;
        const double a4 = -0.00976;
        const double a5 =  0.00108;

        double x2 = x * x;
        double x3 = x2 * x;
        double x4 = x3 * x;
        double x5 = x4 * x;

        *E1_approx = -log(x) + (a0 + a1*x + a2*x2 + a3*x3 + a4*x4 + a5*x5);

    } else {

        // Large x rational approximation
        const double b0 =  0.26777;
        const double b1 =  8.63476;
        const double b2 = 18.05902;
        const double b3 =  8.57333;

        const double c0 =  3.95850;
        const double c1 = 21.09965;
        const double c2 = 25.63296;
        const double c3 =  9.57332;

        double x2 = x * x;
        double x3 = x2 * x;

        double numerator   = b0 + b1*x + b2*x2 + b3*x3;
        double denominator = c0 + c1*x + c2*x2 + c3*x3;
        *E1_approx = (exp(-(x)) / x) * (numerator / denominator);
    }
}

void get_coulomb_interaction_length(double Z, double A, double pc, double theta_init, double* lambda_coulomb){
    double b_coulomb;
    double E1, cs_coulomb;
    double R;
    double t_cut = ((pc)*2.325*(theta_init))*((pc)*2.325*(theta_init));
    double hbar_c = sqrt(0.389); // [mb*GeV^2]
    double constant = (4*M_PI*Z*Z*(1./137.)*(1./137.)*(hbar_c*hbar_c));

    get_slope_hadron_nucleus(A, &b_coulomb);
    R = 2*hbar_c*sqrt(b_coulomb);
    E1_approx(&E1, (R*R*(b_coulomb)*t_cut));
    cs_coulomb = -constant * (R*R*(b_coulomb)*(E1) - exp(-R*R*(b_coulomb)*t_cut)/t_cut);
    *lambda_coulomb = 1. / (N * cs_coulomb);
}

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
    double A          = MaterialData_get_A(material);
    double N          = MaterialData_get__atomic_per_volume(material);
    double Z          = sqrt(MaterialData_get__Z2_eff(material));
    everest->ecmsq    = 2*XC_PROTON_MASS*1.0e-3*pc;
    double sqrt_s      = sqrt(everest->ecmsq);

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t scatter               = everest->coll->record_scatterings;

    total_cross_section(&interaction_length_tot, &cross_section_tot, A, Z, N, sqrt_s);
    //FindRoot_find_path_length(finder, traj);
    //double mcs_path_length = FindRoot_get_path_length(finder);
    double mcs_path_length = length;
    if ( (mcs_path_length - interaction_length_tot) < 1e-12) {
        // MCS to exit
        // since mcs is already done.. then either we do ionsiaton loss here outside
        // Ionisation loss to interaction point
        calculate_ionisation_properties(everest, material, pc);
        pc = calcionloss(everest, material, part, mcs_path_length, pc, 1);
        return pc; // false for nucl int.
    } else {
        double interaction_lengths[6], coulomb_length;
        double Neff = MaterialData_get__num_nucleons_eff(material);
        // double theta_init = atan2(MultipleCoulombTrajectory_get_tan_t0( (MultipleCoulombTrajectory) LocalTrajectory_member(traj) ), 1);
        double theta_init = atan2(LocalParticle_get_xp(part), 1);

        // Get interaction lengths for all types of interactions, to find the dominant one
        get_interaction_length(material, interaction_lengths, cross_section_tot, A, Z, N, Neff, theta_init, sqrt_s, pc);
        // 1 = inel, 2 = el, 3 = prod, 4 = sd, 5 = pp/pn, 6 = coulomb

        // Finding the smallest length
        double min_length = interaction_lengths[0]; // Inelastic
        int min_index = 1;

        if (interaction_lengths[1] < min_length) { // Elastic
            min_length = interaction_lengths[1];
            min_index = 2;
        }
        if (interaction_lengths[2] < min_length) { // Production
            min_length = interaction_lengths[2];
            min_index = 3;
        }
        if (interaction_lengths[3] < min_length) { // Single diffractive
            min_length = interaction_lengths[3];
            min_index = 4;
        }
        if (interaction_lengths[4] < min_length) { // Proton-proton / proton-neutron
            min_length = interaction_lengths[4];
            min_index = 5;
        }
        if (interaction_lengths[5] < min_length) { // Coulomb
            min_length = interaction_lengths[5];
            min_index = 6;
        }

        // Ionisation loss to interaction point
        calculate_ionisation_properties(everest, material, pc);
        pc = calcionloss(everest, material, part, min_length, pc, 1);

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
            get_slope_hadron_nucleus(A, &b_nuclear_elastic);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/b_nuclear_elastic)/pc;

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
                double b_sd;
                get_slope_single_diffraction(everest->ecmsq, &b_sd, part);
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
#endif // XCOLL_EVEREST_NUCLEAR_INTERACTIONS_H