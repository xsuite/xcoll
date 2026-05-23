// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_ROT_H
#define XCOLL_GEOM_ROT_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#endif  // XO_CONTEXT_CPU


/*gpufun*/
double YRotation_single_particle_rotate_only(LocalParticle* part, double s, double angle){
    double x   = LocalParticle_get_x(part);
    double rpp = LocalParticle_get_rpp(part);
    double sin_y = sin(angle);
    double cos_y = cos(angle);
    LocalParticle_set_x(part, x*cos_y - s*sin_y);
    LocalParticle_add_to_px(part,-angle/rpp);
    return x*sin_y + s*cos_y;  // new s
}

#endif /* XCOLL_GEOM_ROT_H */
