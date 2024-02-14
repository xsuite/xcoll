// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_CHECKS_H
#define XCOLL_CHECKS_H

/*gpufun*/
int8_t xcoll_check_particle_init(RandomRutherfordData rng, LocalParticle* part) {
    int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
    if (!is_tracking){
        printf("Collimator tracking code is called, but we are not supposed to be tracking!");
        fflush(stdout);
    }
    int8_t rng_is_set  = assert_rng_set(part, RNG_ERR_SEEDS_NOT_SET);
    if (!rng_is_set){
        printf("Random generator seeds in particles object are not set!");
        fflush(stdout);
    }
    int8_t ruth_is_set = assert_rutherford_set(rng, part, RNG_ERR_RUTH_NOT_SET);
    if (!ruth_is_set){
        printf("Rutherford random generator not initialised!");
        fflush(stdout);
    }
    return is_tracking*rng_is_set*ruth_is_set;
}
    
#endif /* XCOLL_CHECKS_H */