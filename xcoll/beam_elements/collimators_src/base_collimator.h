// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_BASE_COLL_H
#define XCOLL_BASE_COLL_H

/*gpufun*/
void BaseCollimator_track_local_particle(BaseCollimatorData el, LocalParticle* part0) {
    // UNUSED(el);
    kill_all_particles(part0, XC_ERR_INVALID_TRACK);
}

#endif /* XCOLL_BASE_COLL_H */
