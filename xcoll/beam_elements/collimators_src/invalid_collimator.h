#ifndef XCOLL_INVALID_H
#define XCOLL_INVALID_H

/*gpufun*/
void xcoll_kill_particle(LocalParticle* part) {
    LocalParticle_set_x(part, 1e30);
    LocalParticle_set_px(part, 1e30);
    LocalParticle_set_y(part, 1e30);
    LocalParticle_set_py(part, 1e30);
    LocalParticle_set_zeta(part, 1e30);
    LocalParticle_update_delta(part, -1);  // zero energy
    LocalParticle_set_state(part, -399);   // xcoll lost state error
}

/*gpufun*/
void xcoll_kill_all_particles(LocalParticle* part0) {
    //start_per_particle_block (part0->part)
        xcoll_kill_particle(part)
    //end_per_particle_block
}

/*gpufun*/
void InvalidCollimator_track_local_particle(InvalidCollimatorData el, LocalParticle* part0) {
    xcoll_kill_all_particles(part0);
}

#endif
