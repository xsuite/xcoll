// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_ADT_H
#define XCOLL_ADT_H

/*gpufun*/
void ADT_track_local_particle(ADTData el, LocalParticle* part0){

    int8_t plane     = ADTData_get__plane(el);
    double amplitude = ADTData_get__amplitude(el);
    int8_t active    = ADTData_get__active(el);

    //start_per_particle_block (part0->part)
        if (active){
            double kick = amplitude * (2.*RandomUniform_generate(part) - 1.);
            if (plane == 1){
                LocalParticle_add_to_xp(part, kick);
            } else if (plane == -1){
                LocalParticle_add_to_yp(part, kick);
            } else {
                LocalParticle_kill_particle(part, XC_ERR_INVALID_XOFIELD);
            }
        }
    //end_per_particle_block
}

#endif /* XCOLL_ADT_H */
