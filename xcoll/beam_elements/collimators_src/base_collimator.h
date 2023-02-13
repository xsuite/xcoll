#ifndef XCOLL_BASE_H
#define XCOLL_BASE_H

/*gpufun*/
int8_t xcoll_assert_tracking(LocalParticle* part){
    // Whenever we are not tracking, e.g. in a twiss, the particle will be at_turn < 0.
    // We test this to distinguish genuine tracking from twiss, and kill the particle in the latter case.
    if (LocalParticle_get_at_turn(part) < 0){
        xcoll_kill_particle(part);
        return 0;
    }
    return 1;
}

/*gpufun*/
int8_t xcoll_assert_rng_set(LocalParticle* part){
    int64_t rng_s1 = LocalParticle_get__rng_s1(part);
    int64_t rng_s2 = LocalParticle_get__rng_s2(part);
    int64_t rng_s3 = LocalParticle_get__rng_s3(part);
    int64_t rng_s4 = LocalParticle_get__rng_s4(part);
    if (rng_s1==0 && rng_s2==0 && rng_s3==0 && rng_s4==0) {
        xcoll_kill_particle(part);
        return 0;
    }
    return 1;
}

/*gpufun*/
void BaseCollimator_track_local_particle(BaseCollimatorData el, LocalParticle* part0) {
    xcoll_kill_all_particles(part0);
}

#endif
