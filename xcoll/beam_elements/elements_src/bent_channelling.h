// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef BENT_CHANNELLING_H
#define BENT_CHANNELLING_H

// Xsuite core includes
#include <headers/particle_states.h>
#include <headers/checks.h>



/*gpufun*/
void BentChannellingDev_track_local_particle(
    BentChannellingDevData el,
    LocalParticle* part0
) {
    track_bent_channelling_particles(el, part0);
}

#endif

