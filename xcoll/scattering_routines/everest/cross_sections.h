// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CROSS_SECTIONS_H
#define XCOLL_EVEREST_CROSS_SECTIONS_H

#include <splines.h>
#include <stddef.h>
// ====== Splines =============

// Horner's method for evaluating polynomials
static inline double eval_interval(const Spline* s, double* sqrt_s) {
    double dx = *sqrt_s - s->x0;
    return ((s->a * dx + s->b) * dx + s->c) * dx + s->d;
}

double eval_spline(const Spline* spline, size_t n, double* sqrt_s) {

    // Outside range
    if (*sqrt_s <= spline[0].x0)
        return eval_interval(&spline[0], sqrt_s);
    if (*sqrt_s >= spline[n-1].x1)
        return eval_interval(&spline[n-1], sqrt_s);

    // Binary search for interval
    size_t left = 0;
    size_t right = n - 1;

    while (left <= right) {
        size_t mid = (left + right) / 2;

        if (*sqrt_s < spline[mid].x0)
            right = mid - 1;
        else if (*sqrt_s > spline[mid].x1)
            left = mid + 1;
        else
            return eval_interval(&spline[mid], sqrt_s);
    }
    return 1e21; // should never reach here
}


// =======================================================
// ====== Cross sections ======
// =======================================================

// ================ Hadron - Nucleus cross sections ================ 
// Proton - Nucleus
void total_cross_section(MaterialData restrict material, double* lambda, double* A, double* molar_mass, double* rho, double* sqrt_s){
    double Z = MaterialData_get__atomic_number(material);
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
    double cs_tot_hN = Z*cs_tot_pp + (*A - Z)*cs_tot_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn
    if (*A < 4) {
        *lambda = (*molar_mass*(*A))/(N_A*(*rho)*cs_tot_hN);
        return;
    } else {
    _get_R(A, &R);
    double cs_tot_hA = 2*M_PI*pow(R,2) * log(1 + (cs_tot_hN)/(2*M_PI*pow(R,2))); // Glauber-Gribov approximation
    *lambda = (*molar_mass)/(N_A*(*rho)*cs_tot_hA);
    }
}


void calculate_elastic_cross_section(MaterialData restrict material, double* cs_el_hN, double* A, double* sqrt_s){
    double Z = MaterialData_get__atomic_number(material);
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
    *cs_el_hN = Z*cs_el_pp + (*A - Z)*cs_el_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn
}

void calculate_inelastic_cross_section(MaterialData restrict material, double* cs_inel_hN, double* A, double* sqrt_s){
    double Z = MaterialData_get__atomic_number(material);
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
    *cs_inel_hN = Z*cs_inel_pp + (*A - Z)*cs_inel_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn
}

// Pion - Nucleus
// Kaon - Nucleus
// etc.
// Proton - Proton

// ======== slopes ============
// angles, and slopes part 2 :) 
void get_slope_hadron_nucleus(MaterialData restrict material, double* sqrt_s, double* b){
    // we can add pions, its a bit different
    // From GEANT4
    double A = MaterialData_get__atomic_mass(material);
    if (A <= 62){
        *b = 14.5 * pow(A, 2.0/3.0); 
    } else {
        *b = 60.0 * pow(A, 1.0/3.0);
    }
}

void get_slope_proton_proton(double* sqrt_s, double* b){
    // from russian guy
    double B0 = 9.3; // +- 0.3 
    double alpha_1 = 0.11; // +- 0.06
    double alpha_2 = 0.03; // +- 0.01
    *b = B0 + 2*alpha_1*log(*sqrt_s) + alpha_2*pow(log(*sqrt_s), 2);
}

void get_slope_single_diffraction(double* sqrt_s, double* b){
    // unsure
}

// =======================================================
// ====== Nuclear interaction length =====================
// =======================================================

void _get_R(double* A, double* R) {
    if (*A > 21) {
        *R = 1.1*pow(*A, 1.0/3.0)*1e-15 * 0.9;  // [m], 0.8 < f(A) < 1.0
    } else {
        *R = 1.1*pow(*A, 1.0/3.0)*1e-15 * 1.05; // [m], 1.0 < f(A) < 1.1
    }
}

/*gpufun*/
void get_interaction_length(EverestData restrict everest, double* lambda, double* A, double* molar_mass, double* rho, double* sqrt_s, double* cs_tot, double cs_type){
    // cs_type: Nucleus: 1 = total, 2 = inelastic, 3 = elastic, 4 = production, 5 = quasi-elastic, 
    //          Nucleon: 6 = single diffractive, 7 = proton-proton/proton-neutron
    double N_A = 6.02214076e23;  // Avogadro's number
    double R;

    if (A < 4) {
        // Semi-supported material
        if (cs_type == 1){ // get total cs 
            // pass
        } else if (cs_type == 2){ // get inelastic cs
            double cs_inel_hN;
            calculate_inelastic_cross_section(material, &cs_inel_hN, &A, sqrt_s);
            *lambda = (*molar_mass*(*A))/(N_A*(*rho)*cs_inel_hN);

        } else if (cs_type == 3){ // get elastic cs
            double cs_el_hN;
            calculate_elastic_cross_section(material, &cs_el_hN, &A, sqrt_s); // TODO: decide either this or tot - inel
            *lambda = (*molar_mass*(*A))/(N_A*(*rho)*cs_el_hN);

        } else if (cs_type == 4){ // get production cs
            double cs_prod_hN;
            calculate_inelastic_cross_section(material, &cs_prod_hN, &A, sqrt_s);
            *lambda = (*molar_mass*(*A))/(N_A*(*rho)*cs_prod_hN);

        } else if (cs_type == 7){ 
            // Proton-proton / proton-neutron
            double Neff = MaterialData_get__num_nucleons_effective(material);
            double cs_pp_hN = eval_spline(spline_tot_pp, N_spline_tot_pp, sqrt_s);
            *lambda = (*molar_mass)/(N_A*(*rho)*cs_pp_hN*Neff);

        } else {
            // Unsupported cs type 
            *lambda = 1e21; // effectively infinite interaction length
        }
    } else {
         if (cs_type == 1){ 
            // // Total 
            // double cs_tot_hA, cs_tot_hN;
            // calculate_total_cross_section(material, &cs_tot_hN, &A, sqrt_s);
            // _get_R(&A, &R);
            // cs_tot_hA = 2*M_PI*pow(R,2) * log(1 + (cs_tot_hN)/(2*M_PI*pow(R,2))); // Glauber-Gribov approximation
            // *lambda = (*molar_mass)/(N_A*(*rho)*cs_tot_hA);

        } else if (cs_type == 2){ 
            // Inelastic
            double cs_inel_hA, cs_tot_hN;
            _get_R(&A, &R);
            cs_inel_hA = M_PI*pow(R,2) * log(1 + (*cs_tot)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            *lambda = (*molar_mass)/(N_A*(*rho)*cs_inel_hA);

        } else if (cs_type == 3){ 
            // Elastic: Total - Inelastic   // TODO: we do have the spline tho, so we can use it directly
            double cs_inel_hA;
            double cs_tot_hA, cs_tot_hN;
            double cs_el_hA;
            _get_R(&A, &R);
            cs_inel_hA = M_PI*pow(R,2) * log(1 + (*cs_tot)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            cs_tot_hA = 2*M_PI*pow(R,2) * log(1 + (*cs_tot)/(2*M_PI*pow(R,2))); // Glauber-Gribov approximation
            cs_el_hA = cs_tot_hA - cs_inel_hA; // elastic is total - inelastic
            if (cs_el_hA < 0){
                cs_el_hA = 1e-10; // In case. Makes Lambda large
            } else{
                *lambda = (*molar_mass)/(N_A*(*rho)*cs_el_hA);
            }

        } else if (cs_type == 4){ 
            // Production
            double cs_inel_hA;
            double cs_prod_hA, cs_inel_hN;
            calculate_inelastic_cross_section(material, &cs_inel_hN, &A, sqrt_s);
            _get_R(&A, &R);
            cs_prod_hA = M_PI*pow(R,2) * log(1 + (cs_inel_hN)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            *lambda = (*molar_mass)/(N_A*(*rho)*cs_prod_hA);
            
        } else if (cs_type == 5){ 
            // Quasi-Elastic: Inelastic - Production
            double cs_inel_hA, cs_tot_hN;
            double cs_prod_hA, cs_inel_hN;
            calculate_inelastic_cross_section(material, &cs_inel_hN, &A, sqrt_s);
            _get_R(&A, &R);
            cs_inel_hA = M_PI*pow(R,2) * log(1 + (*cs_tot)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            cs_prod_hA = M_PI*pow(R,2) * log(1 + (cs_inel_hN)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            cs_qel_hA = cs_inel_hA - cs_prod_hA; // quasi-elastic is inelastic - production
            if (cs_qel_hA < 0){
                cs_qel_hA = 1e-10; // In case. Makes Lambda large
            } else{ 
                *lambda = (*molar_mass)/(N_A*(*rho)*cs_qel_hA);
            }

        } else if (cs_type == 6){ 
            // Single diffractive
            double cs_sd_hA, cs_tot_hN, alpha;
            _get_R(&A, &R);
            alpha = (*cs_tot)/(2*M_PI*pow(R,2) + *cs_tot);
            cs_sd_hA = M_PI*pow(R,2) * (alpha - log(1 + alpha)); // Glauber-Gribov approximation
            *lambda = (*molar_mass)/(N_A*(*rho)*cs_sd_hA);

        } else if (cs_type == 7){ 
            // Proton-proton / proton-neutron
            double Neff = MaterialData_get__num_nucleons_effective(material);
            double cs_pp_hN = eval_spline(spline_tot_pp, N_spline_tot_pp, sqrt_s);
            *lambda = (*molar_mass)/(N_A*(*rho)*cs_pp_hN*Neff);

        } else {
            // Unsupported cs type 
            *lambda = 1e21; // effectively infinite interaction length
        }
    }
}
