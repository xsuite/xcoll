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
    if (x > 10.0) {
        double invx = 1.0 / x;
        double invx2 = invx * invx;

        *E2 = invx * (1.0
            - 2.0 * invx
            + 6.0 * invx2
            - 24.0 * invx2 * invx);
        return;
    }
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
// ====== Coulomb Cross Sections =====================
// =======================================================
/*gpufun*/
void get_coulomb_interaction_length(double Z, double pc, double N,
                                    double theta_init, double nuclear_slope, double* lambda_coulomb){
    double E2, cs_coulomb;
    double R;
    double t_cut = ((pc)*2.325*(theta_init))*((pc)*2.325*(theta_init));
    double hbar_c_squared = 0.389;                                              // [mb*GeV^2]
    double hbar_c         = 0.197;                                              // [GeV*fm]
    double constant       = (4*M_PI*Z*Z*(1./137.)*(1./137.)*(hbar_c_squared));
    double R_fm           = 2.0 * hbar_c * sqrt(nuclear_slope);                 // fm
    double R_GeV          = R_fm / hbar_c;                                      // GeV^-1

    E2_approx(&E2, (R_GeV*R_GeV*(856.)*t_cut));
    cs_coulomb = constant * (R_GeV*R_GeV*856.*(E2));
    if (cs_coulomb < 1e-15){
        printf("Coulomb cross section is very small: %e mb.\n", cs_coulomb);
    }
    *lambda_coulomb = 1. / (N * cs_coulomb*1.0e-27);
}
/*gpufun*/
void get_coulomb_cross_section(double Z, double pc, double N,
                               double theta_init, double nuclear_slope, double* cs_coulomb){
    double E2;
    double R;
    double t_cut = ((pc)*2.325*(theta_init))*((pc)*2.325*(theta_init));
    double hbar_c_squared = 0.389;                                              // [mb*GeV^2]
    double hbar_c         = 0.197;                                              // [GeV*fm]
    double constant       = (4*M_PI*Z*Z*(1./137.)*(1./137.)*(hbar_c_squared));
    double R_fm           = 2.0 * hbar_c * sqrt(nuclear_slope);                 // fm
    double R_GeV          = R_fm / hbar_c;                                      // GeV^-1

    E2_approx(&E2, (R_GeV*R_GeV*(856.)*t_cut));
    *cs_coulomb = constant * (R_GeV*R_GeV*856.*(E2))/((R_GeV*R_GeV*(856.)*t_cut));
    if (*cs_coulomb < 1e-15){
        *cs_coulomb = 1e-15; // Avoid negative cross section
    }
}


// =======================================================
// ====== Nuclear interaction length =====================
// =======================================================

/*gpufun*/
void get_interaction_length(MaterialData restrict material, double interaction_lengths[5], 
                            double cs_tot, double Z, double N,
                            double theta_init, double sqrt_s, double pc) {
    double cs_inel_hA = MaterialData_evaluate_glauber_spline(material, sqrt_s, 1);
    double cs_el_hA   = MaterialData_evaluate_glauber_spline(material, sqrt_s, 2);
    double cs_el_nucleon = cs_inel_hA - MaterialData_evaluate_glauber_spline(material, sqrt_s, 3); // el nucleon
    double cs_sd_hA   = MaterialData_evaluate_glauber_spline(material, sqrt_s, 4);
    double nuclear_slope = MaterialData_get__nuclear_slope(material);
    double lambda_coulomb;
    get_coulomb_cross_section(Z, pc, N, theta_init, nuclear_slope, &lambda_coulomb);

    // Inelastic
    interaction_lengths[0] = (1)/(N*cs_inel_hA*1.0e-27);   // [m]

    // Elastic: Total - Inelastic
    if (cs_el_hA < 1e-15){
        cs_el_hA = 1e-12; // In case. Makes Lambda large
    } else {
        interaction_lengths[1] = (1)/(N*cs_el_hA*1.0e-27); // [m]
    }

    // Elastic nucleon
    interaction_lengths[2] = 1/(N*cs_el_nucleon*1.0e-27);   // [m]
    // Single diffractive
    interaction_lengths[3] = (1)/(N*cs_sd_hA*1.0e-27);     // [m]
    }

#endif // XCOLL_EVEREST_CROSS_SECTIONS_H