// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_JAW_H
#define XCOLL_EVEREST_JAW_H
#include <math.h>
#include <stdio.h>


/*gpufun*/
double jaw(EverestData restrict everest, LocalParticle* part, double p, double length, int edge_check) {
    if (LocalParticle_get_state(part) < 1){
        // Do nothing if already absorbed
        return p;
    }

    double rlen = length;
    double s0 = LocalParticle_get_s(part);
    double x0 = LocalParticle_get_x(part);
    double y0 = LocalParticle_get_y(part);
    p /= 1e9;   // Energy (not momentum) in GeV

    if (everest->coll->only_mcs) {
        mcs(everest, part, rlen, p, edge_check);

    } else {
        // Do a step for a point-like interaction.
        // Get monte-carlo interaction length.
        while (1) {
            double Ax, Ay;
            double s1 = LocalParticle_get_s(part) - s0;
            calculate_ionisation_properties(everest, p);
            double length_step = everest->xintl*RandomExponential_generate(part);

            // If the monte-carlo interaction length is longer than the remaining
            // length, then put it to the remaining length, do mcs and return.
            if (length_step > rlen) {
                mcs(everest, part, rlen, p, edge_check);
                break;
            }
// changes from here -----------------------------------------------------  
            A(everest->coll->radl, p, &Ax, &Ay);
            mcs_trajectory = MultipleCoulombTrajectory(everest->coll->radl, Ax, Ay); // like this? Then check how if we do it or nucl? 
            // get crossing
            //curve_length = MultipleCoulombTrajectory_length(s0, x0, y0, s1, s2=crossing, Ax, Ay);
            // if we do then get crossing 
            // if not then nucl int, and start over again 
            // before we restart we update s? 

            // Otherwise do multi-coulomb scattering.
            mcs(everest, part, length_step, p, edge_check);

            if(LocalParticle_get_x(part) <= 0) {
                // PARTICLE LEFT COLLIMATOR BEFORE ITS END.
                break;
            }
// changes to here -----------------------------------------------------
            p = nuclear_interaction(everest, part, p);
            if (LocalParticle_get_state(part) < 1){
                // PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.
                break;
            }

            // Calculate the remaining interaction length and close the iteration loop.
            rlen = rlen - length_step;
        }
        // TODO: ionisation loss should also be calculated when only_mcs
        double m_dpodx = calcionloss(everest, part, rlen);  // DM routine to include tail // TODO: should not be rlen but s after updating
        double s = LocalParticle_get_s(part) - s0;
        p = p-m_dpodx*s; // TODO: This is correct: ionisation loss is only calculated and applied at end of while (break)
    }
    return p*1e9;  // Back to eV
}

#endif /* XCOLL_EVEREST_JAW_H */
