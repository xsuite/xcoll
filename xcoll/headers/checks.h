// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_CHECKS_H
#define XCOLL_CHECKS_H

#include <headers/track.h>
#include <headers/checks.h>
#include <headers/particle_states.h>
#include <random/random_src/rutherford.h>
#include <xcoll/headers/particle_states.h>

// This is a quick macro to use inside a function body on a parameter that is not
// used inside the function (this avoids throwing warnings at compilation time).
#ifndef UNUSED
#define UNUSED(expr) (void)(expr)
#endif

GPUFUN
int8_t xcoll_check_particle_init(RandomRutherfordData rng, LocalParticle* part) {
    int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
    if (!is_tracking){
        printf("Collimator tracking code is called, but we are not supposed to be tracking!"); //only_for_context cpu_serial
        fflush(stdout);                                                                        //only_for_context cpu_serial
    }
    int8_t rng_is_set  = assert_rng_set(part, RNG_ERR_SEEDS_NOT_SET);
    if (!rng_is_set){
        printf("Random generator seeds in particles object are not set!"); //only_for_context cpu_serial
        fflush(stdout);                                                    //only_for_context cpu_serial
    }
    int8_t ruth_is_set = assert_rutherford_set(rng, part, RNG_ERR_RUTH_NOT_SET);
    if (!ruth_is_set){
        printf("Rutherford random generator not initialised!"); //only_for_context cpu_serial
        fflush(stdout);                                         //only_for_context cpu_serial
    }
    return is_tracking*rng_is_set*ruth_is_set;
}

#endif /* XCOLL_CHECKS_H */
