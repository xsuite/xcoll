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


#endif /* XCOLL_CHECKS_H */
