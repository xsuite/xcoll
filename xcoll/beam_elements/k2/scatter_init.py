import numpy as np

def calculate_scattering(p0,anuc,rho,zatom,emr,csref0,csref1,csref5,bnref):
    
    # Output parameters
    cprob = np.array([0,0,0,0,0,0], dtype=np.float64)
    csect = np.array([0,0,0,0,0,0], dtype=np.float64) # Cross section
    
    # Constants 
    pptref = 0.04
    freeco = 1.618
    pmap  = 938.271998
    fnavo  = 6.02214076e23

    ecmsq = (2*(pmap*1.0e-3)) * p0
    xln15s = np.log(0.15*ecmsq)
    # Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
    pptot = 0.041084 - 0.0023302*np.log(ecmsq) + 0.00031514*np.log(ecmsq)**2
    # Claudia used the fit from TOTEM for ppel (in barn)
    ppel = (11.7-1.59*np.log(ecmsq)+0.134*np.log(ecmsq)**2)/1e3
    # Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
    ppsd = (4.3+0.3*np.log(ecmsq))/1e3

    # Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
    bpp = 7.156 + 1.439*np.log(np.sqrt(ecmsq))

    # freep: number of nucleons involved in single scattering
    freep = freeco * anuc**(1/3)

    # compute pp and pn el+single diff contributions to cross-section
    # (both added : quasi-elastic or qel later)
    csect[3] = freep * ppel
    csect[4] = freep * ppsd

    # correct TOT-CSec for energy dependence of qel
    # TOT CS is here without a Coulomb contribution
    csect[0] = csref0 + freep * (pptot - pptref)
    bn = (bnref * csect[0]) / csref0

    # also correct inel-CS
    csect[1] = (csref1 * csect[0]) / csref0

    # Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect[2] = ((csect[0] - csect[1]) - csect[3]) - csect[4]
    csect[5] = csref5

    # Now add Coulomb
    csect[0] = csect[0] + csect[5]

    # Interaction length in meter
    xintl = (1.0e-2*anuc)/(((fnavo * rho)*csect[0])*1e-24)

    # Filling CProb with cumulated normalised Cross-sections
    cprob[5] = 1
    for i in range(1,5,1):
        cprob[i] = cprob[i-1] + csect[i]/csect[0]

    return cprob, xintl, bn, ecmsq, xln15s, bpp, csect