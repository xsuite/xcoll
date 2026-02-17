// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCL_H
#define XCOLL_EVEREST_NUCL_H

// ====== lengths =============

// Get length function from lambda
void get_interaction_length(MaterialData restrict material, double pc, double* lambda, double cs_type){
    // get the interaction length of 
    double N_A = 6.02214076e23;  // Avogadro's number
    double A = MaterialData_get__atomic_mass(material);
    double molar_mass = MaterialData_get__molar_mass(material);
    double rho = MaterialData_get__density(material);

    if (cs_type == 1){ // get total cs 
        double cs_tot_A;
        _calculate_total_cross_section(material, &cs_tot_A, &A, pc); // pointer!!
        *lambda = (molar_mass*A)/(N_A*rho*cs_tot_A); // probably make this a pointer somewhere
    } else if (cs_type == 2){ // get elastic cs
        double cs_el_A;
        _calculate_elastic_cross_section(material, &cs_el_A, &A, pc); // pointer!!
        *lambda = (molar_mass*A)/(N_A*rho*cs_el_A); // probably make this a pointer somewhere
    } else if (cs_type == 3){ // get inelastic cs
        double cs_inel_A;
        _calculate_inelastic_cross_section(material, &cs_inel_A, &A, pc); // pointer!!
        *lambda = (molar_mass*A)/(N_A*rho*cs_inel_A); // probably make this a pointer somewhere
    } else if (cs_type == 4){ // get single diffractive cs
        double cs_sd_A;
        _calculate_single_diff_cross_section(material, &cs_sd_A, &A, pc); // pointer!!
        *lambda = (molar_mass*A)/(N_A*rho*cs_sd_A); // probably make this a pointer somewhere
    } else {
        // Unsupported cs type 
        *lambda = 1e21; // effectively infinite interaction length
    }
    //....
}
// calculate cs tot -- POINTER
void _calculate_total_cross_section(MaterialData restrict material, double* cs_tot, double* A, double pc){
    double Z = MaterialData_get__atomic_number(material);
    //.... get file with polynomials for the splines.. --> add to files
    // read files --> and get the get_splines... with pc
    // A cs = Z cs pp + (A-Z) cs pn
    // 
}


void calculate_elastic_cross_section(){
    //.... should we say energy here, or when we actually use it? 
}

void calculate_inelastic_cross_section(){
    //....
}

void calculate_single_diff_cross_section(){
    //....
}

void get_splines(...what file, pc){
    //....get the splines and the poly of the cs. 
    // if A get file A get polynomial for A
    // if B get file B get polynomial for B
    // if C get file C get polynomial for C
    // ...
}

// calc cs el ++
// compare lengths (but probs in jaw ? )

// ======== slopes ============
// angles, and slopes part 2 :) 