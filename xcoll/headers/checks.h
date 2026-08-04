// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_CHECKS_H
#define XCOLL_CHECKS_H

#ifdef XO_CONTEXT_CPU
#include <stdio.h>
#include <stdint.h>  // for int64_t etc
#endif /* XO_CONTEXT_CPU */

#include "xobjects/headers/common.h"
#include "xtrack/headers/checks.h"


// This is a quick macro to use inside a function body on a parameter that is not
// used inside the function (this avoids throwing warnings at compilation time).
#ifndef UNUSED
#define UNUSED(expr) (void)(expr)
#endif


GPUFUN
void xcoll_kill_with_message(LocalParticle* part0, int64_t kill_state,
                             GPUGLMEM const char *message) {
    kill_all_particles(part0, kill_state);
#ifdef XO_CONTEXT_CPU
        printf("Error: %s", message);
        fflush(stdout);
#endif /* XO_CONTEXT_CPU */
}


GPUFUN
int8_t xcoll_check_particle_init(RandomRutherfordData rng, LocalParticle* part) {
    int8_t is_tracking = assert_tracking(part, XC_ERR_INVALID_TRACK);
#ifdef XO_CONTEXT_CPU
    if (!is_tracking) {
        printf("Collimator tracking code is called, but we are not supposed to be tracking!");
        fflush(stdout);
    }
#endif /* XO_CONTEXT_CPU */
    int8_t rng_is_set  = assert_rng_set(part, RNG_ERR_SEEDS_NOT_SET);
#ifdef XO_CONTEXT_CPU
    if (!rng_is_set) {
        printf("Random generator seeds in particles object are not set!");
        fflush(stdout);
    }
#endif /* XO_CONTEXT_CPU */
    int8_t ruth_is_set = assert_rutherford_set(rng, part, RNG_ERR_RUTH_NOT_SET);
#ifdef XO_CONTEXT_CPU
    if (!ruth_is_set) {
        printf("Rutherford random generator not initialised!");
        fflush(stdout);
    }
#endif /* XO_CONTEXT_CPU */
    return is_tracking*rng_is_set*ruth_is_set;
}

#endif /* XCOLL_CHECKS_H */
