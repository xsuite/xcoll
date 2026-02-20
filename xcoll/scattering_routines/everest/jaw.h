// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_JAW_H
#define XCOLL_EVEREST_JAW_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#endif  // XO_CONTEXT_CPU


/*gpufun*/
double jaw(EverestData restrict everest, MaterialData restrict material, LocalParticle* part,
           double pc, double length, int edge_check) {
    if (LocalParticle_get_state(part) < 1){
        // Do nothing if already absorbed
        return pc;
    }

    pc /= 1.e9; // [GeV]

    double rlen = length;
    double s0 = LocalParticle_get_s(part);
    while (1) {
        // Length of the step until nuclear interaction
        double length_step;
        if (MaterialData_get__cross_section(material, 0) < 0){
            // No cross-section defined: only do MCS
            length_step = 1.e21;
        } else {
            length_step = everest->xintl*RandomExponential_generate(part);
        }

        if (length_step > rlen) {
            // Length to nuclear interaction is longer than remaining: MCS to end and exit collimator
            mcs(everest, material, part, rlen, pc, edge_check);
            break;
        }

        mcs(everest, material, part, length_step, pc, edge_check);
        if (LocalParticle_get_state(part) < 1 || (edge_check && LocalParticle_get_x(part) <= 0)){
            // Particle lost all energy due to ionisation, or left the collimator
            break;
        }

        pc = nuclear_interaction(everest, material, part, pc);
        if (LocalParticle_get_state(part) < 1){
            // Particle was absorbed
            break;
        }

        // Calculate the remaining interaction length and close the iteration loop.
        rlen = rlen - length_step;
    }
    calculate_ionisation_properties(everest, material, pc);
    double ionisation_length = LocalParticle_get_s(part) - s0;
    pc = calcionloss(everest, material, part, ionisation_length, pc, 1);
    return pc*1e9;  // Back to eV
}

#endif /* XCOLL_EVEREST_JAW_H */
