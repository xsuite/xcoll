// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRY_PROP_H
#define XCOLL_EVEREST_CRY_PROP_H
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

struct CrystalProperties {
    double const_dech;
    double prob_tail_c1;
    double prob_tail_c2;
    double prob_tail_c3;
    double prob_tail_c4;
    double energy_loss;
    double energy_loss_tail;
};


struct CrystalProperties calculate_crystal_properties(double pc, CrystalMaterialData material) {

    struct CrystalProperties prop;

    double const pmap = 938.271998;
    double const pmae = 0.51099890;
    double const k = 0.307075; // Constant in front bethe-bloch [mev g^-1 cm^2]

    double const aTF = 0.194e-10; // Screening function [m]
    double const dP  = 1.920e-10; // Distance between planes (110) [m]
    double const u1  = 0.075e-10; // Thermal vibrations amplitude

    double crade = 2.817940285e-15;

    // Material properties
    double const exenergy = CrystalMaterialData_get_excitation_energy(material);
    double const rho      = CrystalMaterialData_get_density(material);
    double const anuc     = CrystalMaterialData_get_A(material);
    double const zatom    = CrystalMaterialData_get_Z(material);

    // Energy variables
    double mom    = pc*1.0e3;                // [GeV]
    double enr    = sqrt(pow(mom, 2.) + pow(pmap, 2.)); // [MeV]
    double gammar = enr/pmap;
    double betar  = mom/enr;
    double bgr    = betar*gammar;
    double mep    = pmae/pmap;  // Electron/proton


    // Dechanneling length calculation
    prop.const_dech = ((256.0/(9*pow(M_PI,2.))) * (1./(log(((2.*pmae)*gammar)/(exenergy*1.0e3)) - 1.))) * ((aTF*dP)/(crade*pmae)); // [m/MeV]
    prop.const_dech = prop.const_dech*1.0e3; // [m/GeV]


    // Ionisation loss
    double tmax = (2.*pmae*pow(bgr,2.))/(1. + 2.*gammar*mep + pow(mep,2.));  // [MeV]
    double plen = sqrt((rho*zatom)/anuc)*28.816e-6; // [MeV]

    prop.energy_loss = (k*zatom/(anuc* pow(betar,2.))) * (0.5*log(((((2.*pmae)*bgr)*bgr)*tmax)/(1.0e6* pow(exenergy,2.))) -
            pow(betar,2.) - log(plen/(exenergy*1.0e3)) - log(bgr) + 0.5);
    prop.energy_loss = prop.energy_loss*rho*1.0e-1; // [GeV]
    prop.energy_loss_tail = ((k*zatom)/(anuc*pow(betar,2.))) * (0.5*log((2*pmae*bgr*bgr*tmax)/(1.0e6*pow(exenergy,2.)))
                                                           - pow(betar,2.) - 
                log(plen/(exenergy*1.0e3)) - log(bgr) + 0.5 + pow(tmax,2.)/(8.*(pow(gammar,2.))*(pow(pmap,2.))));
    prop.energy_loss_tail = prop.energy_loss_tail*rho*1.0e-1; // [GeV/m]

    double Tt = prop.energy_loss*1.0e3 + k*zatom*4.0e2*rho/(anuc*pow(betar, 2.)); // [MeV]

    // Calculate different coefficients for terms in dz (length) to get the tail probability
    double const prob_factor = rho*1.e2*k*zatom/(anuc*pow(betar, 2.));
    prop.prob_tail_c1 = prob_factor * 0.5 / Tt;
    prop.prob_tail_c2 = prob_factor * (
                        tmax/(4.*pow(gammar,2.)*pow(pmap,2.)) - 0.5/tmax
                        - log(tmax/Tt)*pow(betar,2.)/(2.*tmax)
                   );  // * dz
    prop.prob_tail_c3 = prob_factor * pow(betar,2.)/(2.*tmax);  // * dz * log(dz)
    prop.prob_tail_c4 = -prob_factor * Tt/(4.*pow(gammar,2.)*pow(pmap,2.));  // * dz * dz

    return prop;
}
    
    
    
    
#endif /* XCOLL_EVEREST_CRY_PROP_H */