// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCLEAR_INTERACTION_H
#define XCOLL_EVEREST_NUCLEAR_INTERACTION_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#endif /* XO_CONTEXT_CPU */

#include "xobjects/headers/common.h"
#include "xtrack/random/random_src/uniform.h"
// #include "xtrack/random/random_src/uniform_accurate.h"
#include "xtrack/random/random_src/exponential.h"
#include "xtrack/random/random_src/rutherford.h"
#include "xcoll/scattering_routines/everest/everest.h"
#include "xcoll/scattering_routines/everest/properties.h"
#include "xcoll/interaction_record/interaction_record_src/interaction_record.h"


GPUFUN
double nuclear_interaction(EverestData RESTRICT everest, MaterialData RESTRICT material,
                           LocalParticle* part, double pc) {
    if (MaterialData_get__cross_section(material, 0) < 0) {
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
    // double aran = RandomUniformAccurate_generate(part); // TODO: do we need 64bit randoms?
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
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_ABSORBED, everest->shape_id);
        if (LocalParticle_get_state(part) == XC_SECONDARY_PARTICLE) {
            LocalParticle_set_state(part, XC_LOST_ON_MATERIAL_SEC);
        } else {
            LocalParticle_set_state(part, XC_LOST_ON_MATERIAL);
        }

    } else {
        double sqrt_t_p;
        if (ichoix==2) {
            // p-n elastic
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_PN_ELASTIC, everest->shape_id);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/everest->bn)/pc; // TODO: do we need 64bit randoms?

        } else if (ichoix==3) {
            // p-p elastic
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_PP_ELASTIC, everest->shape_id);
            sqrt_t_p = sqrt(RandomExponential_generate(part)/everest->bpp)/pc; // TODO: do we need 64bit randoms?

        } else if (ichoix==4) {
            // Single diffractive
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_SINGLE_DIFFRACTIVE, everest->shape_id);
            double xm2 = exp(RandomUniform_generate(part)*everest->xln15s); // TODO: do we need 64bit randoms?
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
                if (sc) InteractionRecordData_log(record, record_index, part, XC_ABSORBED, everest->shape_id);
                if (LocalParticle_get_state(part) == XC_SECONDARY_PARTICLE) {
                    LocalParticle_set_state(part, XC_LOST_ON_MATERIAL_SEC);
                } else {
                    LocalParticle_set_state(part, XC_LOST_ON_MATERIAL);
                }
                pc = 1.e-9;
                sqrt_t_p = 0;
            } else {
                // Corrected 1/p into 1/sqrt(pp')
                sqrt_t_p = sqrt(RandomExponential_generate(part)/bsd)/sqrt(pc_in*pc); // TODO: do we need 64bit randoms?
            }

        } else {
            // Coulomb
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_COULOMB, everest->shape_id);
            sqrt_t_p = sqrt(RandomRutherford_generate(everest->coll->rng, part))/pc; // TODO: do we need 64bit randoms?
        }

        // theta = arccos(1 + t/(2p^2))  =>  tan(theta) = sqrt( -t/p^2 * (1 + t/(4p^2)) ) / (1 + t/(2p^2))
        // Note that in elastic scattering, t < 0, but we sampled t > 0 so we need to flip the sign
        double tan_theta = sqrt_t_p * sqrt(1 - sqrt_t_p*sqrt_t_p/4)/(1 - sqrt_t_p*sqrt_t_p/2);
        double alpha = 2*M_PI*RandomUniform_generate(part); // TODO: do we need 64bit randoms?
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

#endif /* XCOLL_EVEREST_NUCLEAR_INTERACTION_H */
