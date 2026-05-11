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
     // -----------------------------
    // Read element parameters
    // -----------------------------
    const double length   = BentChannellingDevData_get_length(el);

    double U0      = BentChannellingDevData_get_U0(el);
    double Umax    = BentChannellingDevData_get_Umax(el);
    double R       = BentChannellingDevData_get_bending_radius(el);

    double dp      = BentChannellingDevData_get_dp(el);
    double aTF     = BentChannellingDevData_get_aTF(el);
    double uT      = BentChannellingDevData_get_uT(el);

    double alpha_i = BentChannellingDevData_get_alpha_i(el);
    double beta_i  = BentChannellingDevData_get_beta_i(el);
    // caclulate constants 
    double aTF_over_beta = aTF/beta_i;
    double two_aTF_over_beta = 2.0*aTF_over_beta;
    double beta_over_aTF = 1/aTF_over_beta;
    double U_N = U_N_(Umax, dp, alpha_i, beta_i, beta_over_aTF);

    START_PER_PARTICLE_BLOCK(part0, part);
        // -----------------------------
        //  Beam & particle physics
        // -----------------------------
        double p0c          = LocalParticle_get_p0c(part);
        double delta        = LocalParticle_get_delta(part);
        double rvv          = LocalParticle_get_rvv(part);
        double beta0        = LocalParticle_get_beta0(part);
        double charge_ratio = LocalParticle_get_charge_ratio(part);
        double chi          = LocalParticle_get_chi(part);

        // Effective beam momentum factor (bpc)
        double bpc = beta0 * rvv * (1.0 + delta) * p0c* charge_ratio / chi;
        //constants for each particle
        double U  = (U_N / bpc)*beta_over_aTF*beta_over_aTF;  
        double sqrt_U = sqrt(U);
        track_bent_channelling_body_single_particle(
                part, el,
                length, Umax, dp, aTF,
            alpha_i, beta_i,
            aTF_over_beta, two_aTF_over_beta, beta_over_aTF,
            U_N, R, bpc, U, sqrt_U
        );
    END_PER_PARTICLE_BLOCK;
}

#endif

