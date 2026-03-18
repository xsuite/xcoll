// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CROSS_SECTIONS_H
#define XCOLL_EVEREST_CROSS_SECTIONS_H

#include <stddef.h>
#include <math.h>

// ====== Splines & Helper functions =============
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

// =======================================================
// ====== Cross sections & Slopes =======================
// =======================================================
/*gpufun*/
void get_slope_hadron_nucleus(double A, double* b){
    // we can add pions, its a bit different
    // From GEANT4
    if (A <= 62){
        *b = 14.5 * pow(A, 2.0/3.0); 
    } else {
        *b = 60.0 * pow(A, 1.0/3.0);
    }
}
/*gpufun*/
void get_slope_proton_proton(double s, double* b){
    // from russian guy
    double B0 = 9.3; // +- 0.3 
    double alpha_1 = 0.11; // +- 0.06
    double alpha_2 = 0.03; // +- 0.01
    *b = B0 + 2*alpha_1*log(s) + alpha_2*pow(log(sqrt(s)), 2);
}
/*gpufun*/
void get_slope_single_diffraction(double s, double* b, LocalParticle* part){
    // from pythia
    double M_2 = exp(RandomUniform_generate(part)*(log(0.15*s)));
    *b = 2*2.3 + 2*0.25*log((s/M_2));
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

// =======================================================
// ====== Nuclear interaction length =====================
// =======================================================

/*gpufun*/
void get_interaction_length(MaterialData restrict material, double interaction_lengths[6], double cs_tot, double A, double Z, 
                            double N, double Neff, double theta_init, double sqrt_s, double pc) {
    // cs_type: Nucleus: 1 = inelastic, 2 = elastic, 3 = production, 4 = quasi-elastic, 
    //          Nucleon: 5 = single diffractive, 6 = proton-proton/proton-neutron, 7 = Coulomb

    // we dont need to take A into account because its already done.
    double cs_inel_hA = MaterialData_get__cs_inel_hA(material);
    double cs_prod_hA = MaterialData_get__cs_prod_hA(material);
    double cs_sd_hA = MaterialData_get__cs_sd_hA(material);
    double cs_el_hA = MaterialData_get__cs_el_hA(material);
    double cs_tot_hA = MaterialData_get__cs_tot_hA(material);
    double cs_tot_pp = MaterialData_get__cs_tot_pp(material);
    double cs_el_pp = MaterialData_get__cs_el_pp(material);
    double cs_inel_pp = MaterialData_get__cs_inel_pp(material);
    double cs_tot_pn = MaterialData_get__cs_tot_pn(material);

    // // Inelastic
    // cs_inel_hA = M_PI*R*R * log(1 + (cs_tot)/(M_PI*(R*R))); // Glauber-Gribov approximation
    interaction_lengths[0] = ((A))/(N*cs_inel_hA);

    // Elastic: Total - Inelastic
    if (cs_el_hA < 1e-15){
        cs_el_hA = 1e-12; // In case. Makes Lambda large
    } else {
        interaction_lengths[1] = (A)/(N*cs_el_hA);
    }

    // Production
    interaction_lengths[2] = ((A))/(N*cs_prod_hA);

    // Single diffractive
    interaction_lengths[3] = ((A))/(N*cs_sd_hA);

    // // Proton-proton / proton-neutron
    // interaction_lengths[4] = ((A))/(N*cs_pp_hN*Neff);
    interaction_lengths[4] = 1e21;

    // Coulomb
    calculate_coulomb_cross_section(Z, A, pc, theta_init, &lambda_coulomb);
    interaction_lengths[5] = lambda_coulomb;
}

#endif // XCOLL_EVEREST_CROSS_SECTIONS_H