// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCL_H
#define XCOLL_EVEREST_NUCL_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU

#include <xtrack/headers/track.h>
#include <xcoll/headers/particle_states.h>
#include <xcoll/scattering_routines/everest/everest.h>
#include <xcoll/scattering_routines/everest/properties.h>


/*gpufun*/
double nuclear_interaction(EverestData restrict everest, MaterialData restrict material,
                           LocalParticle* part, double pc) {
    if (MaterialData_get__cross_section(material, 0) < 0){
        // Unsupported material for nuclear interaction
        return pc;
    }
    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;

#ifdef XCOLL_REFINE_ENERGY
    calculate_scattering(everest, material, pc);
#endif

    //Choose nuclear interaction
    double aran = RandomUniform_generate(part);
    int ichoix = 1;

    while (aran > everest->cprob[ichoix]) {
        ichoix += 1;
    }

    // Do the interaction
    // Scattered angle is cos theta = 1 + t / (2p^2) for elastic scattering
    //    from Mandelstam t = (p1-p3)^2) = 2m^2 - 2E1E3 + 2p1.p3
    //    if elastic, p1 = p3, and hence t = 2m^2 - 2(m^2 + p^2) + 2p^2 cos(theta)
    int64_t i_slot = -1;
    if (ichoix==1) {
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
        LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);

    } else {
        double sqrt_t_p;
        if (ichoix==2) {
            // p-n elastic
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_PN_ELASTIC);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/everest->bn)/pc;

        } else if (ichoix==3) {
            // p-p elastic
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_PP_ELASTIC);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/everest->bpp)/pc;

        } else if (ichoix==4) {
            // Single diffractive
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_SINGLE_DIFFRACTIVE);
            double xm2 = exp(RandomUniform_generate(part)*everest->xln15s);
            double bsd;
            if (xm2 < 2.) {
                bsd = 2*everest->bpp;
            } else if (xm2 >= 2. && xm2 <= 5.) {
                bsd = ((106.0 - 17.0*xm2)*everest->bpp)/36.0;
            } else {
                bsd = (7*everest->bpp)/12.0;
            }
            double pc_in = pc;
            pc = pc*(1 - xm2/everest->ecmsq);
            if (pc <= 1.e-9 || pc != pc) {
                // Very small (<1eV) or NaN
                if (sc) InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
                LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);
                pc = 1.e-9;
                sqrt_t_p = 0;
            } else {
                // Corrected 1/p into 1/sqrt(pp')
                sqrt_t_p = sqrt(RandomExponential_generate(part)/bsd)/sqrt(pc_in*pc);
            }

        } else {
            // Coulomb
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_COULOMB);
            sqrt_t_p = sqrt(RandomRutherford_generate(everest->coll->rng, part))/pc;
        }

        // theta = arccos(1 + t/(2p^2))  =>  tan(theta) = sqrt( -t/p^2 * (1 + t/(4p^2)) ) / (1 + t/(2p^2))
        // Note that in elastic scattering, t < 0, but we sampled t > 0 so we need to flip the sign
        double tan_theta = sqrt_t_p * sqrt(1 - sqrt_t_p*sqrt_t_p/4)/(1 - sqrt_t_p*sqrt_t_p/2);
        double alpha = 2*M_PI*RandomUniform_generate(part);
        double tan_theta_x = tan_theta*cos(alpha);
        double tan_theta_y = tan_theta*sin(alpha);

        // Change the angles
#ifdef XCOLL_USE_EXACT
        LocalParticle_add_to_exact_xp_yp(part, tan_theta_x, tan_theta_y);
#else
        LocalParticle_add_to_xp_yp(part, tan_theta_x, tan_theta_y);
#endif

        if (sc) InteractionRecordData_log_child(record, i_slot, part);
    }

    return pc;
}

#endif /* XCOLL_EVEREST_NUCL_H */
