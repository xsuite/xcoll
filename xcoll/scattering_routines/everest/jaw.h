// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_JAW_H
#define XCOLL_EVEREST_JAW_H
#include <math.h>
#include <stdio.h>

#include <headers/track.h>
#include <xcoll/scattering_routines/everest/everest.h>
#include <xcoll/scattering_routines/everest/ionisation_loss.h>
#include <xcoll/scattering_routines/everest/nuclear_interaction.h>
#include <xcoll/scattering_routines/everest/multiple_coulomb_scattering.h>


GPUFUN
double jaw(EverestData restrict everest, LocalParticle* part, double pc, double length, int edge_check) {
    if (LocalParticle_get_state(part) < 1){
        // Do nothing if already absorbed
        return pc;
    }

    pc /= 1.e9; // [GeV]

    if (everest->coll->only_mcs) {
        // TODO: ionisation loss should also be calculated when only_mcs
        mcs(everest, part, length, pc, edge_check);

    } else {
        double rlen = length;
        double s0 = LocalParticle_get_s(part);
        while (1) {
            // Length of the step until nuclear interaction
            double length_step = everest->xintl*RandomExponential_generate(part);

            if (length_step > rlen) {
                // Length to nuclear interaction is longer than remaining: MCS to end and exit collimator
                mcs(everest, part, rlen, pc, edge_check);
                break;
            }

            mcs(everest, part, length_step, pc, edge_check);
            if (LocalParticle_get_state(part) < 1 || (edge_check && LocalParticle_get_x(part) <= 0)){
                // Particle lost all energy due to ionisation, or left the collimator
                break;
            }

            pc = nuclear_interaction(everest, part, pc);
            if (LocalParticle_get_state(part) < 1){
                // Particle was absorbed
                break;
            }

            // Calculate the remaining interaction length and close the iteration loop.
            rlen = rlen - length_step;
        }
        calculate_ionisation_properties(everest, pc);
        double ionisation_length = LocalParticle_get_s(part) - s0;
        pc = calcionloss(everest, part, ionisation_length, pc, 1);
    }
    return pc*1e9;  // Back to eV
}

#endif /* XCOLL_EVEREST_JAW_H */
