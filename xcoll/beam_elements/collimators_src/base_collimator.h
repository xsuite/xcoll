#ifndef XCOLL_BASE_H
#define XCOLL_BASE_H


// TODO:
//    Do not split 4d and zeta in drifts
//    Use drift function from xtrack Drift element (call its C function)
//    Use rotation function from xtrack XYRotation element (call its C function)

/*gpufun*/
void xcoll_drift_6d(LocalParticle* part0, double length) {
    //start_per_particle_block (part0->part)
        double const rpp    = LocalParticle_get_rpp(part);
        double const rv0v   = 1./LocalParticle_get_rvv(part);
        double const xp     = LocalParticle_get_px(part) * rpp;
        double const yp     = LocalParticle_get_py(part) * rpp;
        double const dzeta  = 1 - rv0v * ( 1. + ( xp*xp + yp*yp ) / 2. );

        LocalParticle_add_to_x(part, xp * length );
        LocalParticle_add_to_y(part, yp * length );
        LocalParticle_add_to_s(part, length);
        LocalParticle_add_to_zeta(part, length * dzeta );
    //end_per_particle_block
}

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
int8_t xcoll_assert_rutherford_set(RandomGeneratorData rng, LocalParticle* part){
    double A = RandomGeneratorData_get_rutherford_A(rng);
    double B = RandomGeneratorData_get_rutherford_B(rng);
    if (A==0. && B==0.) {
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
