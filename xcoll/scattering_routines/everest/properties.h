// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_PROP_H
#define XCOLL_EVEREST_PROP_H
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


/*gpufun*/
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


/*gpufun*/
void calculate_ionisation_properties(EverestData restrict everest, double pc) {
    if (everest->coll->only_mcs){   // TODO: this should be done smarter
        return;
    }

    // Material properties
    double const exenergy = everest->coll->exenergy*1.0e3; // [MeV]
    double const rho      = everest->coll->rho;
    double const anuc     = everest->coll->anuc;
    double const zatom    = everest->coll->zatom;

    // Energy variables
    double momentum = pc*1.0e3;   // [MeV]
    double energy   = sqrt(pow(momentum, 2.) + pow(XC_PROTON_MASS, 2.)); // [MeV]
    double gammar   = energy/XC_PROTON_MASS;
    double betar    = momentum/energy;
    double bgr      = betar*gammar;
    double mep      = XC_ELECTRON_MASS/XC_PROTON_MASS;  // Electron/proton

    // tmax is max energy loss from kinematics
    // mep = mass_electron / mass_beam
    double tmax = 2.*XC_ELECTRON_MASS*pow(bgr, 2.)/ (1. + 2.*gammar*mep + pow(mep, 2.));  // [MeV]
    double plen = sqrt(rho*zatom/anuc)*28.816e-6; // [MeV]

    // TODO: Bete-Bloch needs to be adapted for high energy (more loss due to delta ray radiation)
    //       and low energy (nuclear scattering)
    // TODO: should scale with z^2 (incoming charge)
    double BB_fac = XC_BETHE_BLOCH*zatom/(anuc* pow(betar, 2.));  // [MeV*cm^2/g]
    everest->energy_loss = 0.5*log(2.*XC_ELECTRON_MASS*bgr*bgr*tmax/pow(exenergy, 2.)) + 0.5;
    everest->energy_loss -= pow(betar, 2.);
    everest->energy_loss -= log(plen/exenergy);
    everest->energy_loss -= log(bgr);
    everest->energy_loss *= BB_fac*rho*1.0e-1; // [GeV/m]

    // Straggling: energy loss is Gaussian for thick absorbers (TODO: it's Landau-Vavilov for thin absorbers when eps << Tmax)
    // TODO: should be multiplied with z^2
    // TODO: does it come naturally from variations in trajectory (when using NewGeometry)?
    everest->straggling = sqrt(BB_fac/2.*tmax*(1-pow(betar, 2.)/2.))*rho*1.0e-1; // [GeV/m]

    everest->energy_loss_tail = 0.5*log(2*XC_ELECTRON_MASS*bgr*bgr*tmax/pow(exenergy, 2.)) + 0.5;
    everest->energy_loss_tail -= pow(betar, 2.);
    everest->energy_loss_tail -= log(plen/exenergy);
    everest->energy_loss_tail -= log(bgr);
    everest->energy_loss_tail += pow(tmax, 2.)/(8.*pow(energy, 2.));
    everest->energy_loss_tail *= XC_BETHE_BLOCH*zatom/(anuc*pow(betar, 2.))*rho*1.0e-1; // [GeV/m]

    double Tt = everest->energy_loss*1.0e3 + XC_BETHE_BLOCH*zatom*4.0e2*rho/(anuc*pow(betar, 2.)); // [MeV/m]

    // Calculate different coefficients for terms in dz (length) to get the tail probability
    double const prob_factor = rho*1.e2*XC_BETHE_BLOCH*zatom/(anuc*pow(betar, 2.));
    everest->prob_tail_c1 = prob_factor * 0.5 / Tt;
    everest->prob_tail_c2 = prob_factor * (
                        tmax/(4.*pow(energy, 2.)) - 0.5/tmax
                        - log(tmax/Tt)*pow(betar,2.)/(2.*tmax)
                   );  // * dz
    everest->prob_tail_c3 = prob_factor * pow(betar,2.)/(2.*tmax);   // * dz * log(dz)
    everest->prob_tail_c4 = -prob_factor * Tt/(4.*pow(energy, 2.));  // * dz * dz
}


/*gpufun*/
double calcionloss(EverestData restrict everest, LocalParticle* part, double ionisation_length, double pc, double scale_factor) {

#ifdef XCOLL_REFINE_ENERGY
    calculate_ionisation_properties(everest, pc);
#endif

    double ionisation_loss;
    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;

    double mp = LocalParticle_get_mass0(part)/1.e9;  // [GeV]   TODO: update when allowing other than protons
    double kinetic_energy = sqrt(pow(pc, 2.) + pow(mp, 2.)) - mp;
    int8_t dead = 0;
    double cutoff = 1.e-6; // Lower cutoff of 1keV
    double new_pc = 0;
    double new_pc_2 = 0;

    double prob_tail = everest->prob_tail_c1 + everest->prob_tail_c2 * ionisation_length
                     + everest->prob_tail_c3 * ionisation_length * log(ionisation_length)
                     + everest->prob_tail_c4 * ionisation_length * ionisation_length;

    if (RandomUniform_generate(part) < prob_tail) {
        ionisation_loss = scale_factor*ionisation_length*everest->energy_loss_tail;
    } else {
        ionisation_loss = scale_factor*ionisation_length*everest->energy_loss;
    }

    // Straggling: energy loss is Gaussian for thick absorbers (TODO: it's Landau-Vavilov for thin absorbers when eps << Tmax)
    double ran = RandomNormal_generate(part);
    ionisation_loss += fabs(ran)*everest->straggling;

    if (ionisation_loss > kinetic_energy - cutoff) {
        // All energy lost due to ionisation!
        dead = 1;
    } else {
        new_pc_2 = pow(kinetic_energy - ionisation_loss + mp, 2.) - pow(mp, 2.);
        if (new_pc_2 <= 1.e-12 || new_pc_2 != new_pc_2){
            // Rounding error. Kill particle to avoid NaN  (a != a is true only for NaN)
            dead = 1;
        } else {
            new_pc = sqrt(new_pc_2);
            if (new_pc != new_pc){
                // NaN
                dead = 1;
            }
        }
    }
    if (dead) {
        if (sc) InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
        LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
        return cutoff;
    } else {
        return new_pc;
    }
}


/*gpufun*/
double calculate_dechanneling_length(EverestData restrict everest, double pc) {
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
