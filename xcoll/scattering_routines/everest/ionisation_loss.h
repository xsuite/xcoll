// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_IONLOSS_H
#define XCOLL_EVEREST_IONLOSS_H
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

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
    double BB_fac = XC_BETHE_BLOCH/2*zatom/anuc/betar/betar;  // [MeV*cm^2/g]
    everest->energy_loss_xi_m = BB_fac*rho*1.0e-1; // [GeV/m]   xi per meter
    double e_ionisation = log(2.*XC_ELECTRON_MASS*bgr*bgr/exenergy);
    double density = 2*log(plen/exenergy) + 2*log(bgr) - 1;
    everest->energy_loss_m = e_ionisation + log(tmax/exenergy) - 2*pow(betar, 2.) - density;
    everest->energy_loss_m *= everest->energy_loss_xi_m; // [GeV/m]

    // Straggling: energy loss is Gaussian for thick absorbers (TODO: it's Landau-Vavilov for thin absorbers when eps << Tmax)
    // TODO: does it come naturally from variations in trajectory (when using NewGeometry)?
    everest->energy_loss_most_probable_m = e_ionisation - log(exenergy*1.e-3) + 0.194 - pow(betar, 2.) - density;
    everest->energy_loss_most_probable_m *= everest->energy_loss_xi_m; // [GeV/m]   missing factor xi ln xi

    everest->energy_loss_tail_m = e_ionisation+ log(tmax/exenergy) - 2*pow(betar, 2.) - density;  // Bethe-Bloch
    everest->energy_loss_tail_m += 2*tmax*tmax/8/energy/energy;
    everest->energy_loss_tail_m *= everest->energy_loss_xi_m; // [GeV/m]

    double Tt = everest->energy_loss_m*1.0e3 + 8*everest->energy_loss_xi_m*1e3; // [MeV/m]

    // Calculate different coefficients for terms in dz (length) to get the tail probability
    double const prob_factor = everest->energy_loss_xi_m*1.e3; // [MeV/m]
    everest->prob_tail_c1 = prob_factor / Tt;
    everest->prob_tail_c2 = prob_factor * (
                        2*tmax/(4.*pow(energy, 2.)) - 1/tmax
                        - log(tmax/Tt)*pow(betar,2.)/tmax
                   );  // * dz
    everest->prob_tail_c3 = prob_factor * pow(betar,2.)/tmax;   // * dz * log(dz)
    everest->prob_tail_c4 = -prob_factor * Tt/(2.*pow(energy, 2.));  // * dz * dz
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

    double xi = ionisation_length*everest->energy_loss_xi_m;
    if (RandomUniform_generate(part) < prob_tail) {
        ionisation_loss = ionisation_length*everest->energy_loss_tail_m;
    } else {
        ionisation_loss = ionisation_length*everest->energy_loss_most_probable_m;
        ionisation_loss += xi*log(xi);
    }

    // Straggling: energy loss is double Gaussian-like for thick absorbers (TODO: it's Landau-Vavilov for thin absorbers when xi << Tmax)
    double ran1 = RandomNormal_generate(part);
    double ran2 = RandomUniform_generate(part);
    // Sample from two Gaussians, to mimic the fat tail of the Landau distribution
    if (ran2 >= 0.8) {
        ionisation_loss += 3*xi + ran1*xi*3.39729; // 3.39729 = 4/sqrt(2ln2)
    } else {
        ionisation_loss += xi + ran1*xi*1.69864; // 1.69864 = 2/sqrt(2ln2)
    }
    if (ionisation_loss < 1.e-12) {
        ionisation_loss = 0;
    } else {
        ionisation_loss *= scale_factor;
    }

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

#endif /* XCOLL_EVEREST_IONLOSS_H */