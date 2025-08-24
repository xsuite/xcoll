// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_BLOWUP_H
#define XCOLL_BLOWUP_H

#include <headers/track.h>
#include <xcoll/headers/particle_states.h>


GPUFUN
void BlowUp_track_local_particle(BlowUpData el, LocalParticle* part0){

    int8_t plane          = BlowUpData_get__plane(el);
    double max_kick       = BlowUpData_get__max_kick(el);
    int8_t active         = BlowUpData_get__active(el);
    int8_t individual     = BlowUpData_get_use_individual_kicks(el);
    int64_t start_at_turn = BlowUpData_get_start_at_turn(el);
    int64_t stop_at_turn  = BlowUpData_get_stop_at_turn(el);

    START_PER_PARTICLE_BLOCK(part0, part);
        if (active){
            int64_t at_turn = LocalParticle_get_at_turn(part);
            if (at_turn >= start_at_turn && at_turn < stop_at_turn){
                double ran;
                if (individual){
                    ran = (2*RandomUniform_generate(part) - 1);
                } else {
                    ran = BlowUpData_get__rans(el, at_turn - start_at_turn);
                }
                double kick = max_kick * ran;
                if (plane == 1){
                    LocalParticle_add_to_px(part, kick);
                } else if (plane == -1){
                    LocalParticle_add_to_py(part, kick);
                } else {
                    LocalParticle_kill_particle(part, XC_ERR_INVALID_XOFIELD);
                }
            }
        }
    END_PER_PARTICLE_BLOCK;
}

#endif /* XCOLL_BLOWUP_H */
