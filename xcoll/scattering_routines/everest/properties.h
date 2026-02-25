// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_PROP_H
#define XCOLL_EVEREST_PROP_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#endif  // XO_CONTEXT_CPU

#include "xobjects/headers/common.h"
#include "xtrack/headers/constants.h"
#include "xcoll/scattering_routines/everest/everest.h"
#include "xcoll/scattering_routines/everest/constants.h"


GPUFUN
void calculate_scattering(EverestData restrict everest, MaterialData restrict material,double pc) {
    // Material properties
    double const atoms = MaterialData_get__atoms_per_volume(material);
    double const bnref = MaterialData_get__nuclear_elastic_slope(material);
    double const freep = MaterialData_get__num_nucleons_eff(material);
    double csref[6];
    csref[0] = MaterialData_get__cross_section(material, 0);
    csref[1] = MaterialData_get__cross_section(material, 1);
    csref[5] = MaterialData_get__cross_section(material, 5);

    double csect[6];

    everest->ecmsq = 2*XC_PROTON_MASS*1.0e-3*pc;
    double ecmsq = everest->ecmsq;
    everest->xln15s = log(0.15*ecmsq);

    // Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
    double pptot = 0.041084 - 0.0023302*log(ecmsq) + 0.00031514*POW2(log(ecmsq));

    // Claudia used the fit from TOTEM for ppel (in barn)
    double ppel = (11.7 - 1.59*log(ecmsq) + 0.134*POW2(log(ecmsq)))/1e3;   // TODO /1.e3

    // Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
    double ppsd = (4.3 + 0.3*log(ecmsq))/1e3;

    // Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
    everest->bpp = 7.156 + 1.439*log(sqrt(ecmsq));

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
    everest->xintl = 1./(atoms*csect[0]*1e-28);

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
double calculate_dechannelling_length(EverestData restrict everest, MaterialData restrict material, double pc) {

    // Material properties
    double exenergy = MaterialData_get__excitation_energy(material);
    if (exenergy < 0){
        // Unsupported material for ionisation loss
        return 1.e21;
    }
    exenergy *= 1.0e-6; // [MeV]

    // Energy variables
    double momentum = pc*1.0e3;   // [MeV]
    double energy   = sqrt(POW2(momentum) + POW2(XC_PROTON_MASS)); // [MeV]
    double gammar   = energy/XC_PROTON_MASS;

    double const_dech = 256.0/(9.*POW2(PI)) / (log(2.*XC_ELECTRON_MASS*gammar/exenergy/1000) - 1.);
    const_dech       *= (XC_SCREENING*XC_PLANE_DISTANCE)/(XC_CRADE*XC_ELECTRON_MASS)*1.0e3; // [m/GeV]
    return const_dech;
}

#endif /* XCOLL_EVEREST_PROP_H */
