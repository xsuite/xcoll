// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_INVALID_H
#define XCOLL_INVALID_H


/*gpufun*/
void InvalidCollimator_track_local_particle(InvalidCollimatorData el, LocalParticle* part0) {
    kill_all_particles(part0, xcoll_state_invalid_tracking); // xcoll lost state error
}

#endif
