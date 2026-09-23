// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_BEAMGAS_H
#define XCOLL_BEAMGAS_H

// A BeamGasScattering element is a passive marker during tracking: the
// beam-residual-gas interactions are generated up front by
// BeamGasScattering.scatter(), which produces a weighted macro-particle
// sample that is then tracked through the line. Keeping the element passive
// (instead of collective) means that installing beam-gas scattering centres
// in a line costs nothing in tracking speed and does not split the line into
// collective chunks.

/*gpufun*/
void BeamGasScattering_track_local_particle(BeamGasScatteringData el, LocalParticle* part0){
    (void) el;
    (void) part0;
    return;
}

#endif /* XCOLL_BEAMGAS_H */
