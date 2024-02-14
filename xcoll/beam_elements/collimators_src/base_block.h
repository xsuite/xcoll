// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_BASE_BLOCK_H
#define XCOLL_BASE_BLOCK_H

/*gpufun*/
void BaseBlock_track_local_particle(BaseBlockData el, LocalParticle* part0) {
    kill_all_particles(part0, XC_ERR_INVALID_TRACK);
}

#endif /* XCOLL_BASE_BLOCK_H */
