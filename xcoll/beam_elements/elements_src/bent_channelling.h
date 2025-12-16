// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef BENT_CHANNELLING_H
#define BENT_CHANNELLING_H

// Xsuite core includes
#include <headers/particle_states.h>
#include <headers/checks.h>



/*gpufun*/
void BentChannellingDev_track_local_particle(
    BentChannellingDevData el,
    LocalParticle* part0
) {
    // pre-calculations
    const double length   = BentChannellingDevData_get_length(el);

    double U0      = BentChannellingDevData_get_U0(el);
    double Umax    = BentChannellingDevData_get_Umax(el);
    double R       = BentChannellingDevData_get_R(el);

    double dp      = BentChannellingDevData_get_dp(el);
    double aTF     = BentChannellingDevData_get_aTF(el);
    double uT      = BentChannellingDevData_get_uT(el);

    double alpha_i = BentChannellingDevData_get_alpha_i(el);
    double beta_i  = BentChannellingDevData_get_beta_i(el);

    START_PER_PARTICLE_BLOCK(part0, part);
        track_bent_channelling_body_single_particle(
                part, el
        );
    END_PER_PARTICLE_BLOCK;
}

#endif

