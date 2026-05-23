// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_RF_SWEEP_H
#define XCOLL_RF_SWEEP_H

/*gpufun*/
void RFSweep_track_local_particle(RFSweepData el, LocalParticle* part0){
    double  const step = RFSweepData_get__rf_sweep_df_step(el); // L*df/(f0 + df)
    int64_t const start_turn = RFSweepData_get_start_turn(el);
    int64_t const stop_turn = RFSweepData_get_stop_turn(el);

    //start_per_particle_block (part0->part)
        int64_t at_turn = LocalParticle_get_at_turn(part);
        if (at_turn >= start_turn && at_turn <= stop_turn){
            LocalParticle_add_to_zeta(part, -(at_turn - start_turn + 1)*step);
        }
    //end_per_particle_block
}

#endif /* XCOLL_RF_SWEEP_H */
