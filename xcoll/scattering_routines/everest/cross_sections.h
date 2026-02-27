// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CROSS_SECTIONS_H
#define XCOLL_EVEREST_CROSS_SECTIONS_H

#include <splines.h>
#include <stddef.h>
#include <math.h>

// ====== Splines & Helper functions =============
/*gpufun*/
void _get_R(double A, double* R) {
    if (A > 21) {
        *R = 1.1*pow(A, 1.0/3.0)*1e-15 * 0.9;  // [m], 0.8 < f(A) < 1.0
    } else {
        *R = 1.1*pow(A, 1.0/3.0)*1e-15 * 1.05; // [m], 1.0 < f(A) < 1.1
    }
}

// Horner's method for evaluating polynomials
static inline double eval_interval(const Spline* s, double sqrt_s) {
    double dx = sqrt_s - s->x0;
    return ((s->a * dx + s->b) * dx + s->c) * dx + s->d;
}
/*gpufun*/
double eval_spline(const Spline* spline, size_t n, double sqrt_s) {
    // Outside range
    if (sqrt_s <= spline[0].x0)
        return eval_interval(&spline[0], sqrt_s);
    if (sqrt_s >= spline[n-1].x1)
        return eval_interval(&spline[n-1], sqrt_s);

    // Binary search for interval
    size_t left = 0;
    size_t right = n - 1;

    while (left <= right) {
        size_t mid = (left + right) / 2;

        if (sqrt_s < spline[mid].x0)
            right = mid - 1;
        else if (sqrt_s > spline[mid].x1)
            left = mid + 1;
        else
            return eval_interval(&spline[mid], sqrt_s);
    }
    return 1e21; // should never reach here
}


// Allen & Hastings Approximation for the Exponential Integral function E1(x) = int_x^inf (e^-t / t) dt
/*gpufun*/
static inline void E1_approx(double* E1_approx, double x) {
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

// ================ Hadron - Nucleus cross sections ================ 
// Proton - Nucleus
/*gpufun*/
void total_cross_section(double* lambda, double* cs_tot, double A, double Z, double molar_mass, 
                         double rho, double sqrt_s){
    double cs_tot_pp  = eval_spline(spline_tot_pp, N_spline_tot_pp, sqrt_s);
    double cs_tot_pn, R;

    if (sqrt_s < 30) {
        cs_tot_pn  = eval_spline(spline_tot_pn, N_spline_tot_pn, sqrt_s); // There is only world data to 30 GeV
    } else {
        cs_tot_pn = cs_tot_pp;
    }
    if (cs_tot_pp > 1e20){
        cs_tot_pp = 0;
    }
    if (cs_tot_pn > 1e20){  // pn is proton - neutron
        cs_tot_pn = 0;
    }

    double cs_tot_hN = Z*cs_tot_pp + (A - Z)*cs_tot_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn

    if (A < 4) {
        *lambda = (molar_mass*A)/(XC_AVOGADRO*rho*cs_tot_hN);
        *cs_tot = cs_tot_hN;
        return;
    } else {
        _get_R(A, &R);
        double cs_tot_hA = 2*M_PI*pow(R,2) * log(1 + (cs_tot_hN)/(2*M_PI*pow(R,2))); // Glauber-Gribov approximation
        *lambda = (molar_mass)/(XC_AVOGADRO*rho*cs_tot_hA);
        *cs_tot = cs_tot_hA;
    }
}

/*gpufun*/
void calculate_elastic_cross_section(double Z, double* cs_el_hN, double A, double sqrt_s){
    double cs_el_pp  = eval_spline(spline_el_pp, N_spline_el_pp, sqrt_s);
    double cs_el_pn = 0;
    if (sqrt_s < 30) {
        cs_el_pn  = eval_spline(spline_el_pn, N_spline_el_pn, sqrt_s); // There is only world data to 30 GeV
    } else {
        cs_el_pn = cs_el_pp;
    }
    if (cs_el_pp > 1e20){
        cs_el_pp = 0;
    }
    if (cs_el_pn > 1e20){  // pn is proton`- neutron
        cs_el_pn = 0;
    }
    *cs_el_hN = Z*cs_el_pp + (A - Z)*cs_el_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn
}
/*gpufun*/
void calculate_inelastic_cross_section(double Z, double* cs_inel_hN, double A, double sqrt_s){
    double cs_inel_pp  = eval_spline(spline_inel_pp, N_spline_inel_pp, sqrt_s);
    double cs_inel_pn;
    if (sqrt_s < 30) {
        cs_inel_pn  = eval_spline(spline_inel_pn, N_spline_inel_pn, sqrt_s); // There is only world data to 30 GeV
    } else {
        cs_inel_pn = cs_inel_pp;
    }
    if (cs_inel_pp > 1e20){
        cs_inel_pp = 0;
    }
    if (cs_inel_pn > 1e20){  // pn is proton - neutron
        cs_inel_pn = 0;
    }
    *cs_inel_hN = Z*cs_inel_pp + (A - Z)*cs_inel_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn
}

void calculate_coulomb_cross_section(double Z, double A, double pc, double theta_init, double* cs_coulomb){
    double* b_coulomb;
    double* E1;
    double R;
    double t_cut = ((pc)*2.325*(theta_init))*((pc)*2.325*(theta_init));
    double hbar_c = sqrt(0.389); // [mb*GeV^2]
    double constant = (4*M_PI*Z*Z*(1./137.)*(1./137.)*(hbar_c*hbar_c));

    get_slope_hadron_nucleus(A, b_coulomb);
    R = 2*hbar_c*sqrt(*b_coulomb);
    E1_approx(E1, (R*R*(*b_coulomb)*t_cut));
    *cs_coulomb = -constant * (R*R*(*b_coulomb)*(*E1) - exp(-R*R*(*b_coulomb)*t_cut)/t_cut);
}

// ======== Slopes ============
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
    *b = 2*2.3 + 2*0.25*log((*s/M_2));
}

// =======================================================
// ====== Nuclear interaction length =====================
// =======================================================

/*gpufun*/
void get_interaction_length(double interaction_lengths[6], double cs_tot, double A, double Z, 
                            double molar_mass, double rho, double theta_init, double sqrt_s, double pc) {
    // cs_type: Nucleus: 1 = inelastic, 2 = elastic, 3 = production, 4 = quasi-elastic, 
    //          Nucleon: 5 = single diffractive, 6 = proton-proton/proton-neutron, 7 = Coulomb
    double R;

    if (A < 4) {
        // Semi-supported material
        double cs_inel_hN, cs_el_hN;
        calculate_inelastic_cross_section(&Z, &cs_inel_hN, &A, sqrt_s);
        calculate_elastic_cross_section(&Z, &cs_el_hN, &A, sqrt_s); // TODO: decide either this or tot - inel

        // Inelastic
        interaction_lengths[0] = (molar_mass*(A))/(XC_AVOGADRO*(rho)*cs_inel_hN);

        // Elastic
        interaction_lengths[1] = (molar_mass*(A))/(XC_AVOGADRO*(rho)*cs_el_hN);

        // Production
        interaction_lengths[2] = (molar_mass*(A))/(XC_AVOGADRO*(rho)*cs_inel_hN);

        // Single diffractive
        interaction_lengths[3] = 1e20; // single diffractive not supported for A < 4

        // Proton-proton / proton-neutron
        double Neff = MaterialData_get__num_nucleons_effective(material);
        double cs_pp_hN = eval_spline(spline_tot_pp, N_spline_tot_pp, sqrt_s);
        interaction_lengths[4] = (molar_mass)/(XC_AVOGADRO*(rho)*cs_pp_hN*Neff);

        // Coulomb
        double cs_coulomb;
        calculate_coulomb_cross_section(&Z, &A, pc, theta_init, &cs_coulomb);
        interaction_lengths[5] = (molar_mass*(A))/(XC_AVOGADRO*(rho)*cs_coulomb);
    
    } else {
        double cs_inel_hA, cs_el_hA;// A after GG, N before
        double cs_prod_hA, cs_inel_hN;
        double cs_coulomb;
        double cs_sd_hA, alpha;
        calculate_inelastic_cross_section(Z, &cs_inel_hN, A, sqrt_s);
        _get_R(A, &R);
        
        // Inelastic
        cs_inel_hA = M_PI*R*R * log(1 + (cs_tot)/(M_PI*(R*R))); // Glauber-Gribov approximation
        interaction_lengths[0] = (molar_mass)/(XC_AVOGADRO*(rho)*cs_inel_hA);

        // Elastic: Total - Inelastic
        cs_el_hA = cs_tot - cs_inel_hA;
        if (cs_el_hA < 0){
            cs_el_hA = 1e-10; // In case. Makes Lambda large
        } else {
            interaction_lengths[1] = (molar_mass)/(XC_AVOGADRO*(rho)*cs_el_hA);
        }

        // Production
        cs_prod_hA = M_PI*R*R * log(1 + (cs_inel_hN)/(M_PI*R*R)); // Glauber-Gribov approximation
        interaction_lengths[2] = (molar_mass)/(XC_AVOGADRO*(rho)*cs_prod_hA);

        // Single diffractive
        alpha = cs_tot/(2*M_PI*R*R + cs_tot);
        cs_sd_hA = M_PI*R*R * (alpha - log(1 + alpha)); // Glauber-Gribov approximation
        interaction_lengths[3] = (molar_mass)/(XC_AVOGADRO*(rho)*cs_sd_hA);

        // Proton-proton / proton-neutron
        double Neff = MaterialData_get__num_nucleons_effective(material);
        double cs_pp_hN = eval_spline(spline_tot_pp, N_spline_tot_pp, sqrt_s);
        interaction_lengths[4] = (molar_mass)/(XC_AVOGADRO*(rho)*cs_pp_hN*Neff);

        // Coulomb
        calculate_coulomb_cross_section(Z, A, pc, theta_init, &cs_coulomb);
        interaction_lengths[5] = (molar_mass)/(XC_AVOGADRO*(rho)*cs_coulomb);
        // Quasi-Elastic: Inelastic - Production
        // double cs_inel_hA, cs_tot_hN;
        // double cs_prod_hA, cs_inel_hN;
        // calculate_inelastic_cross_section(&Z, &cs_inel_hN, &A, sqrt_s);
        // _get_R(&A, &R);
        // cs_inel_hA = M_PI*R*R * log(1 + (*cs_tot)/(M_PI*R*R)); // Glauber-Gribov approximation
        // cs_prod_hA = M_PI*R*R * log(1 + (cs_inel_hN)/(M_PI*R*R)); // Glauber-Gribov approximation
        // cs_qel_hA = cs_inel_hA - cs_prod_hA; // quasi-elastic is inelastic - production
        // if (cs_qel_hA < 0){
        //     cs_qel_hA = 1e-10; // In case. Makes Lambda large
        // } else { 
        //     *lambda = (molar_mass)/(N_A*(rho)*cs_qel_hA);
        // }
    }
}
#endif // XCOLL_EVEREST_CROSS_SECTIONS_H