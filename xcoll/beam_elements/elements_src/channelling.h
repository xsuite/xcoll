// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_CHANNELLINGDEV_H
#define XCOLL_CHANNELLINGDEV_H

/*gpufun*/
void ChannellingDev_track_local_particle(ChannellingDevData el, LocalParticle* part0){
    double length = ChannellingDevData_get_length(el);

    // Some constants
    double beta_i = 0.573481;
    double a_TF = 0.194E-10;
    double U_N = 3.526347;

    //start_per_particle_block (part0->part)
        double x0 = LocalParticle_get_x(part);
        double theta0 = LocalParticle_get_xp(part);
        double x, theta;

        // Here should come the actual channelling code, that takes x0 and theta0 as
        // input, calculates Et and phi from them, and then applies the Jacobi functions
        // to calculate the x and theta after a certain length of channelling.

        LocalParticle_set_x(part, x);
        LocalParticle_set_xp(part, theta);
        // TODO: should drift other coordinates.
    //end_per_particle_block
}

#endif /* XCOLL_CHANNELLINGDEV_H */
