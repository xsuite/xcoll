#ifndef XCOLL_BASE_H
#define XCOLL_BASE_H
    
/*gpufun*/
void BaseCollimator_track_local_particle(BaseCollimatorData el, LocalParticle* part0) {
    xcoll_kill_all_particles(part0);
}

#endif