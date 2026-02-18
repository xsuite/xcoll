// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCL_H
#define XCOLL_EVEREST_NUCL_H

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
// ====== Nuclear interaction length and cross sections ======
// =======================================================
void _get_R(double* A, double* R) {
    if (*A > 21) {
        *R = 1.1*pow(*A, 1.0/3.0)*1e-15 * 0.9; // [m], 0.8 < f(A) < 1.0
    } else {
        *R = 1.1*pow(*A, 1.0/3.0)*1e-15 * 1.05; // [m], 1.0 < f(A) < 1.1
    }
}

/*gpufun*/
void get_interaction_length(MaterialData restrict material, double* sqrt_s, double* lambda, double cs_type){
    double N_A = 6.02214076e23;  // Avogadro's number
    double A = MaterialData_get__atomic_mass(material);
    double molar_mass = MaterialData_get__molar_mass(material);
    double rho = MaterialData_get__density(material);
    double R;

    if (A < 4) {
        // Semi-supported material
        if (cs_type == 1){ // get total cs 
            double cs_tot_A;
            _calculate_total_cross_section(material, &cs_tot_A, &A, sqrt_s);
            *lambda = (molar_mass*A)/(N_A*rho*cs_tot_A);

        } else if (cs_type == 2){ // get elastic cs
            double cs_el_A;
            _calculate_elastic_cross_section(material, &cs_el_A, &A, sqrt_s);
            *lambda = (molar_mass*A)/(N_A*rho*cs_el_A);

        } else if (cs_type == 3){ // get inelastic cs
            double cs_inel_A;
            _calculate_inelastic_cross_section(material, &cs_inel_A, &A, sqrt_s);
            *lambda = (molar_mass*A)/(N_A*rho*cs_inel_A);

        } else if (cs_type == 4){ // get single diffractive cs
            double cs_sd_A;
            _calculate_single_diff_cross_section(material, &cs_sd_A, &A, sqrt_s);
            *lambda = (molar_mass*A)/(N_A*rho*cs_sd_A);

        } else {
            // Unsupported cs type 
            *lambda = 1e21; // effectively infinite interaction length
        }
    } else {
         if (cs_type == 1){ 
            // Total 
            double cs_tot_hA, cs_tot_hN;
            _calculate_total_cross_section(material, &cs_tot_hN, &A, sqrt_s);
            _get_R(&A, &R);
            cs_tot_hA = 2*M_PI*pow(R,2) * log(1 + (cs_tot_hN)/(2*M_PI*pow(R,2))); // Glauber-Gribov approximation
            *lambda = (molar_mass)/(N_A*rho*cs_tot_hA);

        } else if (cs_type == 2){ 
            // Inelastic
            double cs_inel_hA, cs_tot_hN;
            _calculate_total_cross_section(material, &cs_tot_hN, &A, sqrt_s);
            _get_R(&A, &R);
            cs_inel_hA = M_PI*pow(R,2) * log(1 + (cs_tot_hN)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            *lambda = (molar_mass)/(N_A*rho*cs_inel_hA);

        } else if (cs_type == 3){ 
            // Elastic: Total - Inelastic
            double cs_inel_hA;
            double cs_tot_hA, cs_tot_hN;
            double cs_el_hA;
            _calculate_total_cross_section(material, &cs_tot_hN, &A, sqrt_s);
            _get_R(&A, &R);
            cs_inel_hA = M_PI*pow(R,2) * log(1 + (cs_tot_hN)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            cs_tot_hA = 2*M_PI*pow(R,2) * log(1 + (cs_tot_hN)/(2*M_PI*pow(R,2))); // Glauber-Gribov approximation
            cs_el_hA = cs_tot_hA - cs_inel_hA; // elastic is total - inelastic
            if (cs_el_hA < 0){
                cs_el_hA = 1e-10; // In case. Makes Lambda large
            } else{
                *lambda = (molar_mass)/(N_A*rho*cs_el_hA);
            }
        } else if (cs_type == 5){ 
            // Production
            double cs_inel_hA;
            double cs_prod_hA, cs_inel_hN;
            _calculate_inelastic_cross_section(material, &cs_inel_hN, &A, sqrt_s);
            _get_R(&A, &R);
            cs_prod_hA = M_PI*pow(R,2) * log(1 + (cs_inel_hN)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            *lambda = (molar_mass)/(N_A*rho*cs_prod_hA);
            
        } else if (cs_type == 4){ 
            // Quasi-Elastic: Inelastic - Production
            double cs_inel_hA, cs_tot_hN;
            double cs_prod_hA, cs_inel_hN;
            _calculate_total_cross_section(material, &cs_tot_hN, &A, sqrt_s);
            _calculate_inelastic_cross_section(material, &cs_inel_hN, &A, sqrt_s);
            _get_R(&A, &R);
            cs_inel_hA = M_PI*pow(R,2) * log(1 + (cs_tot_hN)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            cs_prod_hA = M_PI*pow(R,2) * log(1 + (cs_inel_hN)/(M_PI*pow(R,2))); // Glauber-Gribov approximation
            cs_qel_hA = cs_inel_hA - cs_prod_hA; // quasi-elastic is inelastic - production
            if (cs_qel_hA < 0){
                cs_qel_hA = 1e-10; // In case. Makes Lambda large
            } else{ 
                *lambda = (molar_mass)/(N_A*rho*cs_qel_hA);
            }
        } else {
            // Unsupported cs type 
            *lambda = 1e21; // effectively infinite interaction length
        }
    }
}

// ================ Hadron - Nucleus cross sections ================
// Proton - Nucleus
void _calculate_total_cross_section_pH(MaterialData restrict material, double* cs_tot_hN, double* A, double* sqrt_s){
    double Z = MaterialData_get__atomic_number(material);
    double cs_tot_pp  = eval_spline(spline_tot_pp, N_spline_tot_pp, sqrt_s);
    double cs_tot_pn  = eval_spline(spline_tot_pn, N_spline_tot_pn, sqrt_s);
    if (cs_tot_pp > 1e20){
        cs_tot_pp = 0;
    }
    if (cs_tot_pn > 1e20){  // pn is proton - neutron
        *cs_tot_pn = 0;
    }
    *cs_tot_hN = Z*cs_tot_pp + (*A - Z)*cs_tot_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn
}


void calculate_elastic_cross_section(MaterialData restrict material, double* cs_el_hN, double* A, double* sqrt_s){
    double Z = MaterialData_get__atomic_number(material);
    double cs_el_pp  = eval_spline(spline_el_pp, N_spline_el_pp, sqrt_s);
    double cs_el_pn  = eval_spline(spline_el_pn, N_spline_el_pn, sqrt_s);
    if (cs_el_pp > 1e20){
        cs_el_pp = 0;
    }
    if (cs_el_pn > 1e20){  // pn is proton - neutron
        *cs_el_pn = 0;
    }
    *cs_el_hN = Z*cs_el_pp + (*A - Z)*cs_el_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn
}

void calculate_inelastic_cross_section(MaterialData restrict material, double* cs_inel_hN, double* A, double* sqrt_s){
    double Z = MaterialData_get__atomic_number(material);
    double cs_inel_pp  = eval_spline(spline_inel_pp, N_spline_inel_pp, sqrt_s);
    double cs_inel_pn  = eval_spline(spline_inel_pn, N_spline_inel_pn, sqrt_s);
    if (cs_inel_pp > 1e20){
        cs_inel_pp = 0;
    }
    if (cs_inel_pn > 1e20){  // pn is proton - neutron
        *cs_inel_pn = 0;
    }
    *cs_inel_hN = Z*cs_inel_pp + (*A - Z)*cs_inel_pn; //  cs * A = Z*cs_pp + (A-Z)*cs_pn
}

void calculate_single_diff_cross_section(MaterialData restrict material, double* cs_sd_hN, double* A, double* sqrt_s){
    //...
}

    //....get the splines and the poly of the cs. 
    // if A get file A get polynomial for A
    // if B get file B get polynomial for B
    // if C get file C get polynomial for C
    // ...

// calc cs el ++
// compare lengths (but probs in jaw ? )

// ======== slopes ============
// angles, and slopes part 2 :) 