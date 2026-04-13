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


/*gpufun*/
void E2_approx(double* E2, double x) {
    const double gamma = 0.5772156649015328606;
    // small x: use series (stable)
    if (x < 1e-3) {
        *E2 = 1.0 - x * (log(x) + gamma - 1.0);
        return;
    }
    double E1;
    E1_approx(&E1, x);
    *E2 = exp(-x) - x * E1;
}
// =======================================================
// ====== Cross sections & Slopes =======================
// =======================================================
// /*gpufun*/
// void get_slope_hadron_nucleus(double A, double* b){
//     // we can add pions, its a bit different
//     // From GEANT4
//     if (A <= 62){
//         *b = 14.5 * pow(A, 2.0/3.0); 
//     } else {
//         *b = 60.0 * pow(A, 1.0/3.0);
//     }
// }

// I dont think we need this anymore
// /*gpufun*/
// void get_slope_proton_proton(double s, double* b){
//     // from russian guy
//     double B0 = 9.3; // +- 0.3 
//     double alpha_1 = 0.11; // +- 0.06
//     double alpha_2 = 0.03; // +- 0.01
//     *b = B0 + 2*alpha_1*log(s) + alpha_2*log(sqrt(s))*log(sqrt(s));
// }
// /*gpufun*/
// void get_slope_single_diffraction(double s, double* b, LocalParticle* part){
//     // from pythia
//     double M_2 = exp(RandomUniform_generate(part)*(log(0.15*s)));
//     *b = 2*2.3 + 2*0.25*log((s/M_2));
// }

void get_coulomb_interaction_length(double Z, double pc, double N,
                                    double theta_init, double nuclear_slope, double* lambda_coulomb){
    double E2, cs_coulomb;
    double R;
    double t_cut = ((pc)*2.325*(theta_init))*((pc)*2.325*(theta_init));
    double hbar_c = sqrt(0.389); // [mb*GeV^2]
    double constant = (4*M_PI*Z*Z*(1./137.)*(1./137.)*(hbar_c*hbar_c));
    R = 2*hbar_c*sqrt(nuclear_slope);
    E2_approx(&E2, (R*R*(856.)*t_cut));
    cs_coulomb = -constant * (R*R*856.*(E2));
    *lambda_coulomb = 1. / (N * cs_coulomb);
}

// =======================================================
// ====== Nuclear interaction length =====================
// =======================================================

/*gpufun*/
void get_interaction_length(MaterialData restrict material, double interaction_lengths[6], 
                            double cs_tot, double Z, double N, double Neff,
                            double theta_init, double sqrt_s, double pc) {
    // cs_type: Nucleus: 1 = inelastic, 2 = elastic, 3 = production, 4 = quasi-elastic, 
    //          Nucleon: 5 = single diffractive, 6 = proton-proton/proton-neutron, 7 = Coulomb

    // double cs_tot_hA  = MaterialData_evaluate_glauber_spline(material, sqrt_s, 0);
    double cs_inel_hA = MaterialData_evaluate_glauber_spline(material, sqrt_s, 1);
    double cs_el_hA   = MaterialData_evaluate_glauber_spline(material, sqrt_s, 2);
    double cs_prod_hA = MaterialData_evaluate_glauber_spline(material, sqrt_s, 3);
    double cs_sd_hA   = MaterialData_evaluate_glauber_spline(material, sqrt_s, 4);
    double nuclear_slope = MaterialData_get__nuclear_slope(material);
    double lambda_coulomb;
    // // Inelastic
    // cs_inel_hA = M_PI*R*R * log(1 + (cs_tot)/(M_PI*(R*R))); // Glauber-Gribov approximation
    interaction_lengths[0] = (1)/(N*cs_inel_hA);

    // Elastic: Total - Inelastic
    if (cs_el_hA < 1e-15){
        cs_el_hA = 1e-12; // In case. Makes Lambda large
    } else {
        interaction_lengths[1] = (1)/(N*cs_el_hA);
    }

    // Production
    interaction_lengths[2] = (1)/(N*cs_prod_hA);

    // Single diffractive
    interaction_lengths[3] = (1)/(N*cs_sd_hA);

    // // Proton-proton / proton-neutron
    // interaction_lengths[4] = ((A))/(N*cs_pp_hN*Neff);
    interaction_lengths[4] = 1e21;

    // Coulomb
    get_coulomb_interaction_length(Z, pc, N, theta_init, nuclear_slope, &lambda_coulomb);
    interaction_lengths[5] = lambda_coulomb;
}

#endif // XCOLL_EVEREST_CROSS_SECTIONS_H