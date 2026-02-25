// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_EMITTANCE_MONITOR_H
#define XCOLL_EMITTANCE_MONITOR_H

#ifdef XO_CONTEXT_CPU
#include <math.h>    // for round
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU

#include "xobjects/headers/common.h"
#include "xtrack/headers/track.h"
#include "xcoll/beam_elements/elements_src/monitor.h"


GPUFUN
void EmittanceMonitor_track_local_particle(
    EmittanceMonitorData el,
    LocalParticle* part0
){
    // Use the stats monitor to track
    ParticleStatsMonitor_track_local_particle(
        (ParticleStatsMonitorData) el, part0
    );
}

#endif /* XCOLL_EMITTANCE_MONITOR_H */
