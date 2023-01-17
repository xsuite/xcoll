#include <stdlib.h>
#include <math.h>
#include <stdio.h>



struct ScatteringParameters {
    double cprob0;
    double cprob1;
    double cprob2;
    double cprob3;
    double cprob4;
    double cprob5;
    double xintl;
    double bn;
    double ecmsq;
    double xln15s;
    double bpp;
};



struct ScatteringParameters calculate_scattering(double p0, MaterialData material) {

    struct ScatteringParameters scat;

    // Material properties
    double const zatom    = MaterialData_get_Z(material);
    double const anuc     = MaterialData_get_A(material);
    double const rho      = MaterialData_get_density(material);
    double const emr      = MaterialData_get_nuclear_radius(material);
    double const bnref    = MaterialData_get_nuclear_elastic_slope(material);
    double const csref0   = MaterialData_get_cross_section(material, 0);
    double const csref1   = MaterialData_get_cross_section(material, 1);
    double const csref5   = MaterialData_get_cross_section(material, 5);
    
    // Constants 
    double pptref = 0.04;
    double freeco = 1.618;
    double pmap   = 938.271998;
    double fnavo  = 6.02214076e23;

    double ecmsq = (2*(pmap*1.0e-3)) * p0;
    double xln15s = log(0.15*ecmsq);
    // Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
    double pptot = 0.041084 - 0.0023302*log(ecmsq) + 0.00031514*pow(log(ecmsq),2);
    // Claudia used the fit from TOTEM for ppel (in barn)
    double ppel = (11.7-1.59*log(ecmsq)+0.134*pow(log(ecmsq),2))/1e3;
    // Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
    double ppsd = (4.3+0.3*log(ecmsq))/1e3;

    // Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
    double bpp = 7.156 + 1.439*log(sqrt(ecmsq));

    // freep: number of nucleons involved in single scattering
    double freep = freeco * pow(anuc,(1./3.));

    // compute pp and pn el+single diff contributions to cross-section
    // (both added : quasi-elastic or qel later)
    double csect3 = freep * ppel;
    double csect4 = freep * ppsd;

    // correct TOT-CSec for energy dependence of qel
    // TOT CS is here without a Coulomb contribution
    double csect0 = csref0 + freep * (pptot - pptref);
    double bn = (bnref * csect0) / csref0;

    // also correct inel-CS
    double csect1 = (csref1 * csect0) / csref0;

    // Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    double csect2 = ((csect0 - csect1) - csect3) - csect4;
    double csect5 = csref5;

    // Now add Coulomb
    csect0 = csect0 + csect5;

    // Interaction length in meter
    double xintl = (1.0e-2*anuc)/(((fnavo * rho)*csect0)*1e-24);

    // Filling CProb with cumulated normalised Cross-sections
    double cprob5 = 1;
    double cprob0 = 0;
    double cprob1 = cprob0 + csect1/csect0;
    double cprob2 = cprob1 + csect2/csect0;
    double cprob3 = cprob2 + csect3/csect0;
    double cprob4 = cprob3 + csect4/csect0;
    // for i in range(1,5,1):
    //     cprob[i] = cprob[i-1] + csect[i]/csect[0]

    scat.cprob0 = cprob0;
    scat.cprob1 = cprob1;
    scat.cprob2 = cprob2;
    scat.cprob3 = cprob3;
    scat.cprob4 = cprob4;
    scat.cprob5 = cprob5;
    scat.xintl  = xintl;
    scat.bn     = bn;
    scat.ecmsq  = ecmsq;
    scat.xln15s = xln15s;
    scat.bpp    = bpp;
    
    return scat;

}

