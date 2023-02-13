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
int8_t xcoll_assert_rutherford_set(RandomGeneratorData rng, LocalParticle* part){
    double A = RandomGeneratorData_get_rutherford_A(rng);
    double B = RandomGeneratorData_get_rutherford_B(rng);
    if (A==0. && B==0.) {
        xcoll_kill_particle(part);
        return 0;
    }
    return 1;
}

/*gpukern*/
void RandomGeneratorData_set_rutherford_by_xcoll_material(RandomGeneratorData ran, GeneralMaterialData material){
    double const zatom    = GeneralMaterialData_get_Z(material);
    double const emr      = GeneralMaterialData_get_nuclear_radius(material);
    double const hcut     = GeneralMaterialData_get_hcut(material);
    double const lcut     = 0.0009982;
    double const c = 0.8561e3; // TODO: Where tha fuck does this come from??
    double A = pow(zatom,2);
    double B = c*pow(emr,2);
    RandomGeneratorData_set_rutherford(ran, A, B, lcut, hcut);
}

/*gpufun*/
void BaseCollimator_track_local_particle(BaseCollimatorData el, LocalParticle* part0) {
    xcoll_kill_all_particles(part0);
}

#endif /* XCOLL_BASE_H */
