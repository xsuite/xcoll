#ifndef XCOLL_EVEREST_SCAT_INIT_H
#define XCOLL_EVEREST_SCAT_INIT_H
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



struct ScatteringParameters calculate_scattering(double p0, MaterialData material) {

    struct ScatteringParameters scat;

    // Material properties
    double const anuc  = MaterialData_get_A(material);
    double const rho   = MaterialData_get_density(material);
    double const bnref = MaterialData_get_nuclear_elastic_slope(material);
    double csref[6];
    csref[0]           = MaterialData_get_cross_section(material, 0);
    csref[1]           = MaterialData_get_cross_section(material, 1);
    csref[5]           = MaterialData_get_cross_section(material, 5);

    double csect[6];
    
    // Constants 
    double pptref = 0.04;
    double freeco = 1.618;
    double pmap   = 938.271998;
    double fnavo  = 6.02214076e23;

    scat.ecmsq = (2*(pmap*1.0e-3)) * p0;
    scat.xln15s = log(0.15*scat.ecmsq);
    // Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
    double pptot = 0.041084 - 0.0023302*log(scat.ecmsq) + 0.00031514*pow(log(scat.ecmsq),2);
    // Claudia used the fit from TOTEM for ppel (in barn)
    double ppel = (11.7-1.59*log(scat.ecmsq)+0.134*pow(log(scat.ecmsq),2))/1e3;                      // TODO /1.e3
    // Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
    double ppsd = (4.3+0.3*log(scat.ecmsq))/1e3;

    // Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
    scat.bpp = 7.156 + 1.439*log(sqrt(scat.ecmsq));

    // freep: number of nucleons involved in single scattering
    double freep = freeco * pow(anuc,(1./3.));

    // compute pp and pn el+single diff contributions to cross-section
    // (both added : quasi-elastic or qel later)
    csect[3] = freep * ppel;
    csect[4] = freep * ppsd;

    // correct TOT-CSec for energy dependence of qel
    // TOT CS is here without a Coulomb contribution
    csect[0] = csref[0] + freep * (pptot - pptref);
    scat.bn = bnref * csect[0] / csref[0];

    // also correct inel-CS
    csect[1] = csref[1] * csect[0] / csref[0];

    // Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect[2] = csect[0] - csect[1] - csect[3] - csect[4];
    csect[5] = csref[5];

    // Now add Coulomb
    csect[0] += csect[5];

    // Interaction length in meter
    scat.xintl = (1.0e-2*anuc)/(((fnavo * rho)*csect[0])*1e-24);

    // Filling CProb with cumulated normalised Cross-sections
//     double cprob[1] = cprob0 + csect1/csect0;
//     double cprob[2] = cprob1 + csect2/csect0;
//     double cprob[3] = cprob2 + csect3/csect0;
//     double cprob[4] = cprob3 + csect4/csect0;

    int i;
    scat.cprob[0] = 0;
    for (i=1; i<5; ++i){
        scat.cprob[i] = scat.cprob[i-1] + csect[i]/csect[0];
    }
    scat.cprob[5] = 1;
    
    return scat;

}


#endif
