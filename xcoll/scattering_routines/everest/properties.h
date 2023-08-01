// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_PROP_H
#define XCOLL_EVEREST_PROP_H
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


struct ScatteringParameters {
    double cprob[6];
    double xintl;
    double bn;
    double ecmsq;
    double xln15s;
    double bpp;
};

void calculate_scattering(struct ScatteringParameters* scat, double p, GeneralMaterialData material, double rescale) {

    // Material properties
    double const anuc  = GeneralMaterialData_get_A(material);
    double const rho   = GeneralMaterialData_get_density(material);
    double const bnref = GeneralMaterialData_get_nuclear_elastic_slope(material);
    double csref[6];
    csref[0]           = GeneralMaterialData_get_cross_section(material, 0);
    csref[1]           = GeneralMaterialData_get_cross_section(material, 1);
    csref[5]           = GeneralMaterialData_get_cross_section(material, 5);

    double csect[6];

    scat->ecmsq = 2*XC_PROTON_MASS*1.0e-3*p;
    scat->xln15s = log(0.15*scat->ecmsq);

    // Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
    double pptot = 0.041084 - 0.0023302*log(scat->ecmsq) + 0.00031514*pow(log(scat->ecmsq), 2.);

    // Claudia used the fit from TOTEM for ppel (in barn)
    double ppel = (11.7 - 1.59*log(scat->ecmsq) + 0.134*pow(log(scat->ecmsq), 2.))/1e3;   // TODO /1.e3

    // Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
    double ppsd = (4.3 + 0.3*log(scat->ecmsq))/1e3;

    // Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
    scat->bpp = 7.156 + 1.439*log(sqrt(scat->ecmsq));

    // freep: number of nucleons involved in single scattering
    double freep = XC_FREE_CO * pow(anuc,(1./3.));

    // Rescale reference cross sections (only used for crystals)
    csref[0] = csref[0]*rescale;
    csref[1] = csref[1]*rescale;

    // compute pp and pn el+single diff contributions to cross-section
    // (both added : quasi-elastic or qel later)
    csect[3] = freep * ppel;
    csect[4] = freep * ppsd;

    // correct TOT-CSec for energy dependence of qel
    // TOT CS is here without a Coulomb contribution
    csect[0] = csref[0] + freep * (pptot - XC_PP_CS_REF);
    //Also correct inel-CS
    if(csref[0] == 0) {   // TODO: is this needed?
        csect[1] = 0;   // TODO: Is this correct? It seems as if it should be huge instead
        scat->bn = 1.e10;
    } else {
        csect[1] = csref[1] * csect[0] / csref[0];
        scat->bn = bnref * csect[0] / csref[0];
    }

    // Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect[2] = csect[0] - csect[1] - csect[3] - csect[4];
    csect[5] = csref[5];

    // Now add Coulomb
    csect[0] += csect[5];

    // Interaction length in meter
    scat->xintl = 1.0e-2*anuc/(XC_AVOGADRO*rho*csect[0]*1e-24);

    // Filling CProb with cumulated normalised Cross-sections
    int i;
    scat->cprob[0] = 0;
    if(csect[0] == 0) {   // TODO: is this needed?
        for (i=1; i<5; ++i){
            scat->cprob[i] = scat->cprob[i-1];  // TODO: seems wrong
        }
    } else {
        for (i=1; i<5; ++i){
            scat->cprob[i] = scat->cprob[i-1] + csect[i]/csect[0];
        }
    }
    scat->cprob[5] = 1;
}


struct IonisationProperties {
    double prob_tail_c1;
    double prob_tail_c2;
    double prob_tail_c3;
    double prob_tail_c4;
    double energy_loss;
    double energy_loss_tail;
};

void calculate_ionisation_properties(struct IonisationProperties* prop, double pc, GeneralMaterialData material) {

    // Material properties
    double const exenergy = GeneralMaterialData_get_excitation_energy(material)*1.0e3; // MeV
    double const rho      = GeneralMaterialData_get_density(material);
    double const anuc     = GeneralMaterialData_get_A(material);
    double const zatom    = GeneralMaterialData_get_Z(material);

    // Energy variables
    double momentum = pc*1.0e3;   // [MeV]
    double energy   = sqrt(pow(momentum, 2.) + pow(XC_PROTON_MASS, 2.)); // [MeV]
    double gammar   = energy/XC_PROTON_MASS;
    double betar    = momentum/energy;
    double bgr      = betar*gammar;
    double mep      = XC_ELECTRON_MASS/XC_PROTON_MASS;  // Electron/proton

    // tmax is max energy loss from kinematics
    double tmax = 2.*XC_ELECTRON_MASS*pow(bgr, 2.)/ (1. + 2.*gammar*mep + pow(mep, 2.));  // [MeV]
    double plen = sqrt(rho*zatom/anuc)*28.816e-6; // [MeV]

    prop->energy_loss = 0.5*log(2.*XC_ELECTRON_MASS*bgr*bgr*tmax/pow(exenergy, 2.)) + 0.5;
    prop->energy_loss -= pow(betar, 2.);
    prop->energy_loss -= log(plen/exenergy);
    prop->energy_loss -= log(bgr);
    prop->energy_loss *= XC_BETHE_BLOCH*zatom/(anuc* pow(betar, 2.))*rho*1.0e-1; // [GeV/m]

    prop->energy_loss_tail = 0.5*log(2*XC_ELECTRON_MASS*bgr*bgr*tmax/pow(exenergy, 2.)) + 0.5;
    prop->energy_loss_tail -= pow(betar, 2.);
    prop->energy_loss_tail -= log(plen/exenergy);
    prop->energy_loss_tail -= log(bgr);
    prop->energy_loss_tail += pow(tmax, 2.)/(8.*pow(energy, 2.));
    prop->energy_loss_tail *= XC_BETHE_BLOCH*zatom/(anuc*pow(betar, 2.))*rho*1.0e-1; // [GeV/m]

    double Tt = prop->energy_loss*1.0e3 + XC_BETHE_BLOCH*zatom*4.0e2*rho/(anuc*pow(betar, 2.)); // [MeV/m]

    // Calculate different coefficients for terms in dz (length) to get the tail probability
    double const prob_factor = rho*1.e2*XC_BETHE_BLOCH*zatom/(anuc*pow(betar, 2.));
    prop->prob_tail_c1 = prob_factor * 0.5 / Tt;
    prop->prob_tail_c2 = prob_factor * (
                        tmax/(4.*pow(energy, 2.)) - 0.5/tmax
                        - log(tmax/Tt)*pow(betar,2.)/(2.*tmax)
                   );  // * dz
    prop->prob_tail_c3 = prob_factor * pow(betar,2.)/(2.*tmax);  // * dz * log(dz)
    prop->prob_tail_c4 = -prob_factor * Tt/(4.*pow(energy, 2.));  // * dz * dz

    return prop;
}


/*gpufun*/
double calcionloss(LocalParticle* part, double length, struct IonisationProperties* prop) {

    double prob_tail = prop->prob_tail_c1 + prop->prob_tail_c2 * length
                     + prop->prob_tail_c3 * length * log(length) + prop->prob_tail_c4 * length * length;

    if (RandomUniform_generate(part) < prob_tail) {
        return prop->energy_loss_tail;
    } else {
        return prop->energy_loss;
    }
}


double calculate_dechanneling_length(double pc, CrystalMaterialData material) {
    // Material properties
    double const exenergy = CrystalMaterialData_get_excitation_energy(material)*1.0e3; // MeV

    // Energy variables
    double momentum = pc*1.0e3;   // [MeV]
    double energy   = sqrt(pow(momentum, 2.) + pow(XC_PROTON_MASS, 2.)); // [MeV]
    double gammar   = energy/XC_PROTON_MASS;

    double const_dech = 256.0/(9.*pow(M_PI, 2.)) / (log(2.*XC_ELECTRON_MASS*gammar/exenergy) - 1.);
    const_dech       *= (XC_SCREENING*XC_PLANE_DISTANCE)/(XC_CRADE*XC_ELECTRON_MASS)*1.0e3; // [m/GeV]
    return const_dech;
}

#endif /* XCOLL_EVEREST_PROP_H */
