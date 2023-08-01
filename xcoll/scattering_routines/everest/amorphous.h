// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_AMORPHOUS_H
#define XCOLL_EVEREST_AMORPHOUS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// Processes
int const proc_out         =  -1;     // Crystal not hit
int const proc_AM          =   1;     // Amorphous
int const proc_VR          =   2;     // Volume reflection
int const proc_CH          =   3;     // Channeling
int const proc_VC          =   4;     // Volume capture
int const proc_absorbed    =   5;     // Absorption
int const proc_DC          =   6;     // Dechanneling
int const proc_pne         =   7;     // Proton-neutron elastic interaction
int const proc_ppe         =   8;     // Proton-proton elastic interaction
int const proc_diff        =   9;     // Single diffractive
int const proc_ruth        =  10;     // Rutherford scattering
int const proc_ch_absorbed =  15;     // Channeling followed by absorption
int const proc_ch_pne      =  17;     // Channeling followed by proton-neutron elastic interaction
int const proc_ch_ppe      =  18;     // Channeling followed by proton-proton elastic interaction
int const proc_ch_diff     =  19;     // Channeling followed by single diffractive
int const proc_ch_ruth     =  20;     // Channeling followed by Rutherford scattering
int const proc_TRVR        = 100;     // Volume reflection in VR-AM transition region
int const proc_TRAM        = 101;     // Amorphous in VR-AM transition region

/*gpufun*/
double* nuclear_interaction(RandomRutherfordData rng, LocalParticle* part, double pc,
                           struct ScatteringParameters* scat, double iProc) {

    double* result = (double*)malloc(2 * sizeof(double));

    //Choose nuclear interaction
    double aran = RandomUniform_generate(part);
    int ichoix = 1;

    while (aran > scat->cprob[ichoix]) {
        ichoix += 1;
    }

    //Do the interaction
    double teta = 0 ; //default value to cover ichoix=1

    if (ichoix==1) {
        iProc = proc_ch_absorbed; //deep inelastic, impinging p disappeared
    } else if (ichoix==2) { //p-n elastic
        iProc = proc_ch_pne;
        teta  = sqrt(RandomExponential_generate(part)/scat->bn)/pc;
    } else if (ichoix==3) { //p-p elastic
        iProc = proc_ch_ppe;
        teta  = sqrt(RandomExponential_generate(part)/scat->bpp)/pc;
    } else if (ichoix==4) { //Single diffractive
        iProc = proc_ch_diff;
        double xm2 = exp(RandomUniform_generate(part)*scat->xln15s);

        double bsd = 0.0;
        if (xm2 < 2.) {
            bsd = 2*scat->bpp;
        } else if (xm2 >= 2. && xm2 <= 5.) {
            bsd = ((106.0 - 17.0*xm2)*scat->bpp)/36.0;
        } else if (xm2 > 5.) {
            bsd = (7*scat->bpp)/12.0;
        }
        teta = sqrt(RandomExponential_generate(part)/bsd)/pc;
        pc = pc*(1 - xm2/scat->ecmsq);
    } else { //(ichoix==5)
        iProc = proc_ch_ruth;
        teta  = sqrt(RandomRutherford_generate(rng, part))/pc;
    }

    double tx = teta*RandomNormal_generate(part);
    double tz = teta*RandomNormal_generate(part);

    //Change p angle
    double rpp = LocalParticle_get_rpp(part);
    LocalParticle_add_to_px(part, tx/rpp);
    LocalParticle_add_to_py(part, tz/rpp);

    result[0] = pc;
    result[1] = iProc;
    return result;
}

/*gpufun*/
double* moveam(RandomRutherfordData rng, LocalParticle* part, double dz, double pc,
               struct ScatteringParameters* scat, CrystalMaterialData material, double iProc) {

    double* result = (double*)malloc(3 * sizeof(double));

//     calculate_scattering(scat, pc, (GeneralMaterialData) material, 1.);

    // Material properties
    double const dlr      = CrystalMaterialData_get_crystal_radiation_length(material);
    double const collnt   = CrystalMaterialData_get_nuclear_collision_length(material);

    // Multiple Coulomb Scattering
    double dya   = (13.6/pc)*sqrt(dz/dlr); // RMS of coloumb scattering MCS (mrad)
    double kxmcs = dya*RandomNormal_generate(part)*1.0e-3;
    double kymcs = dya*RandomNormal_generate(part)*1.0e-3;

    double rpp = LocalParticle_get_rpp(part);
    LocalParticle_add_to_px(part, kxmcs/rpp);
    LocalParticle_add_to_py(part, kymcs/rpp);

    // Can nuclear interaction happen?
    double zlm = collnt*RandomExponential_generate(part);

    if (zlm < dz) {
        double* result_ni = nuclear_interaction(rng, part, pc, scat, iProc);
        pc    = result_ni[0];
        iProc = result_ni[1];
        free(result_ni);
    } else {
        zlm = -1; // TODO: this is a hack, to say that no interaction happened
    }

    result[0] = pc;
    result[1] = iProc;
    result[2] = zlm;
    return result; // Turn on/off nuclear interactions
}



// double* Drift_amorphous_4d(RandomRutherfordData rng, LocalParticle* part, double length, double pc, CrystalMaterialData material, double iProc, struct IonisationProperties properties){

//     double* result = (double*)malloc(2 * sizeof(double));

//     Drift_single_particle_4d(part, 0.5*length);
//     double energy_loss = calcionloss(part, length, properties);
//     double* result_am = moveam(rng, part, length, energy_loss, pc, material, iProc);

//     pc = result_am[0];
//     iProc = result_am[1];
//     free(result_am);

//     Drift_single_particle_4d(part, 0.5*length);

//     result[0] = pc;
//     result[1] = iProc;
//     return result;
// }

/*gpufun*/
double* Amorphous(RandomRutherfordData rng, LocalParticle* part, double length, double pc, CrystalMaterialData material,
                 struct IonisationProperties* properties, struct ScatteringParameters* scat, int do_moveam) {

    double* result = (double*)malloc(2 * sizeof(double));
    double iProc = proc_AM;

    if (do_moveam == 1) {
        Drift_single_particle_4d(part, 0.5*length);
        double energy_loss = calcionloss(part, length, properties);
        pc  = pc - energy_loss*length; // Energy lost because of ionization process[GeV]
        double* result_am = moveam(rng, part, length, pc, scat, material, iProc);
        pc = result_am[0];
        iProc = result_am[1];
        free(result_am);
        Drift_single_particle_4d(part, 0.5*length);
    } else{
        Drift_single_particle_4d(part, length);
        double energy_loss = calcionloss(part, length, properties);
        pc  = pc - energy_loss*length; // Energy lost because of ionization process[GeV]
    }

    result[0] = pc;
    result[1] = iProc;
    return result;
}

#endif /* XCOLL_EVEREST_AMORPHOUS_H */
