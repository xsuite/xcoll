// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EMITTANCE_MONITOR_H
#define XCOLL_EMITTANCE_MONITOR_H

#include <xtrack/headers/track.h>


/*gpufun*/
void EmittanceMonitor_track_local_particle(EmittanceMonitorData el, LocalParticle* part0){
    ParticleStatsMonitor_track_local_particle((ParticleStatsMonitorData) el, part0);
}

#endif /* XCOLL_EMITTANCE_MONITOR_H */
