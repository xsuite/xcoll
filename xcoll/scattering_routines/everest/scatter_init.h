#include <stdlib.h>
#include <math.h>

double* calculate_scattering(double p0, double anuc, double rho, double zatom, double emr, double csref0, double csref1, double csref5, double bnref) {
    
    static double result[11];
    
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

    result[0]=cprob0;
    result[1]=cprob1;
    result[2]=cprob2;
    result[3]=cprob3;
    result[4]=cprob4;
    result[5]=cprob5;
    result[6]=xintl;
    result[7]=bn;
    result[8]=ecmsq;
    result[9]=xln15s;
    result[10]=bpp;
    
    return result;

}