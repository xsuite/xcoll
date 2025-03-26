// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_NUCL_H
#define XCOLL_EVEREST_NUCL_H
#include <math.h>
#include <stdio.h>


/*gpufun*/
double nuclear_interaction(EverestData restrict everest, LocalParticle* part, double pc) {
    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    int8_t sc = everest->coll->record_scatterings;

#ifdef XCOLL_REFINE_ENERGY
    calculate_scattering(everest, pc);
#endif

    //Choose nuclear interaction
    double aran = RandomUniform_generate(part);
    int ichoix = 1;

    while (aran > everest->cprob[ichoix]) {
        ichoix += 1;
    }

    //Do the interaction
    int64_t i_slot = -1;
    if (ichoix==1) {
        if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_ABSORBED);
        LocalParticle_set_state(part, XC_LOST_ON_EVEREST_COLL);

    } else {
        double teta;
        if (ichoix==2) { // p-n elastic
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_PN_ELASTIC);
            teta  = sqrt(RandomExponential_generate(part)/everest->bn)/pc;

        } else if (ichoix==3) { // p-p elastic
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_PP_ELASTIC);
            teta  = sqrt(RandomExponential_generate(part)/everest->bpp)/pc;

        } else if (ichoix==4) { // Single diffractive
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
                teta = 0;
            } else {
                // Corrected 1/p into 1/sqrt(pp')
                teta = sqrt(RandomExponential_generate(part)/bsd)/sqrt(pc_in*pc);
            }

        } else { // Coulomb
            if (sc) i_slot = InteractionRecordData_log(record, record_index, part, XC_COULOMB);
            teta = sqrt(RandomRutherford_generate(everest->coll->rng, part))/pc;
        }

        // TODO: I am not convinced that we can just sample two independent random numbers
        // I believe it should be tan(tx) = cos(phi) * tan(teta)    and    tan(ty) = sin(phi) * tan(teta)
        // with phi uniformly sampled between 0 and 2 pi
        double tx = teta*RandomNormal_generate(part);
        double tz = teta*RandomNormal_generate(part);

        //Change the angles
#ifdef XCOLL_USE_EXACT
        LocalParticle_add_to_exact_xp_yp(part, tx, tz);
#else
        LocalParticle_add_to_xp_yp(part, tx, tz);
#endif

        if (sc) InteractionRecordData_log_child(record, i_slot, part);
    }

    return pc;
}

#endif /* XCOLL_EVEREST_NUCL_H */