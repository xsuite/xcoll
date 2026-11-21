// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_RF_SWEEP_H
#define XCOLL_RF_SWEEP_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#endif /* XO_CONTEXT_CPU */

#include "xobjects/headers/common.h"
#include "xtrack/headers/track.h"


GPUFUN
void RFSweep_track_local_particle(RFSweepData el, LocalParticle* part0) {
    double  const step = RFSweepData_get__rf_sweep_df_step(el); // L*df/(f0 + df)
    int64_t const start_turn = RFSweepData_get_start_turn(el);
    int64_t const stop_turn = RFSweepData_get_stop_turn(el);

    START_PER_PARTICLE_BLOCK(part0, part);
        int64_t at_turn = LocalParticle_get_at_turn(part);
        if (at_turn >= start_turn && at_turn <= stop_turn) {
            LocalParticle_add_to_zeta(part, -(at_turn - start_turn + 1)*step);
        }
    END_PER_PARTICLE_BLOCK;
}

#endif /* XCOLL_RF_SWEEP_H */
