// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_PROP_H
#define XCOLL_EVEREST_PROP_H
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <headers/track.h>
#include <xcoll/scattering_routines/everest/everest.h>
#include <xcoll/scattering_routines/everest/constants.h>


GPUFUN
void calculate_scattering(EverestData restrict everest, double pc) {
    if (everest->coll->only_mcs){   // TODO: this should be done smarter
        return;
    }
    // Material properties
    double const anuc  = everest->coll->anuc;
    double const rho   = everest->coll->rho;
    double const bnref = everest->coll->bnref;
    double csref[6];
    csref[0] = everest->coll->csref[0];
    csref[1] = everest->coll->csref[1];
    csref[5] = everest->coll->csref[5];

    double csect[6];

    everest->ecmsq = 2*XC_PROTON_MASS*1.0e-3*pc;
    double ecmsq = everest->ecmsq;
    everest->xln15s = log(0.15*ecmsq);

    // Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
    double pptot = 0.041084 - 0.0023302*log(ecmsq) + 0.00031514*pow(log(ecmsq), 2.);

    // Claudia used the fit from TOTEM for ppel (in barn)
    double ppel = (11.7 - 1.59*log(ecmsq) + 0.134*pow(log(ecmsq), 2.))/1e3;   // TODO /1.e3

    // Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
    double ppsd = (4.3 + 0.3*log(ecmsq))/1e3;

    // Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
    everest->bpp = 7.156 + 1.439*log(sqrt(ecmsq));

    // freep: number of nucleons involved in single scattering
    double freep = XC_FREE_CO * pow(anuc,(1./3.));

    // Rescale reference cross sections (only used for crystals)
    csref[0] = csref[0]*everest->rescale_scattering;
    csref[1] = csref[1]*everest->rescale_scattering;

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
        everest->bn = 1.e10;
    } else {
        csect[1] = csref[1] * csect[0] / csref[0];
        everest->bn = bnref * csect[0] / csref[0];
    }

    // Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect[2] = csect[0] - csect[1] - csect[3] - csect[4];
    csect[5] = csref[5];

    // Now add Coulomb
    csect[0] += csect[5];

    // Interaction length in meter
    everest->xintl = 1.0e-2*anuc/(XC_AVOGADRO*rho*csect[0]*1e-24);

    // Filling CProb with cumulated normalised Cross-sections
    int i;
    everest->cprob[0] = 0;
    if(csect[0] == 0) {   // TODO: is this needed?
        for (i=1; i<5; ++i){
            everest->cprob[i] = everest->cprob[i-1];  // TODO: seems wrong
        }
    } else {
        for (i=1; i<5; ++i){
            everest->cprob[i] = everest->cprob[i-1] + csect[i]/csect[0];
        }
    }
    everest->cprob[5] = 1;
}


GPUFUN
double calculate_dechannelling_length(EverestData restrict everest, double pc) {
    // Material properties
    double const exenergy = everest->coll->exenergy;

    // Energy variables
    double momentum = pc*1.0e3;   // [MeV]
    double energy   = sqrt(pow(momentum, 2.) + pow(XC_PROTON_MASS, 2.)); // [MeV]
    double gammar   = energy/XC_PROTON_MASS;

    double const_dech = 256.0/(9.*pow(M_PI, 2.)) / (log(2.*XC_ELECTRON_MASS*gammar/exenergy/1000) - 1.);
    const_dech       *= (XC_SCREENING*XC_PLANE_DISTANCE)/(XC_CRADE*XC_ELECTRON_MASS)*1.0e3; // [m/GeV]
    return const_dech;
}


#endif /* XCOLL_EVEREST_PROP_H */
