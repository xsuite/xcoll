// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef TRACK_BENT_CHANNELLING_H
#define TRACK_BENT_CHANNELLING_H

#include <stdint.h>


// ================================================================
//  Single-particle tracking (BODY)
// ================================================================

GPUFUN
void track_bent_channelling_body_single_particle(
        LocalParticle* part,
         BentChannellingDevData el)
{
    // Particle dead? skip.
    if (LocalParticle_get_state(part) <= 0)
        return;
// -----------------------------
    // Read element parameters
    // -----------------------------
    const double length   = BentChannellingDevData_get_length(el);

    double U0      = BentChannellingDevData_get_U0(el);
    double Umax    = BentChannellingDevData_get_Umax(el);
    double R       = BentChannellingDevData_get_R(el);

    double dp      = BentChannellingDevData_get_dp(el);
    double aTF     = BentChannellingDevData_get_aTF(el);
    double uT      = BentChannellingDevData_get_uT(el);

    double alpha_i = BentChannellingDevData_get_alpha_i(el);
    double beta_i  = BentChannellingDevData_get_beta_i(el);

    int8_t method  = BentChannellingDevData_get_method(el);
    int8_t order   = BentChannellingDevData_get_order(el);
    int8_t variant = BentChannellingDevData_get_variant(el);

    int n_steps = BentChannellingDevData_get_n_steps(el);
    // Apply defaults (0 -> use default)
    if (method  == 0) method  = BENTCH_DEFAULT_METHOD;
    if (order   == 0) order   = BENTCH_DEFAULT_ORDER;
    if (variant == 0) variant = BENTCH_DEFAULT_VARIANT;
    if (n_steps == 0) n_steps = BENTCH_DEFAULT_NSTEPS;

// APPLY DEFAULTS 
    if (U0      == 0.0) U0     = BENTCH_DEFAULT_U0;
    if (Umax    == 0.0) Umax   = BENTCH_DEFAULT_UMAX;

    if (aTF     == 0.0) aTF     = BENTCH_DEFAULT_ATF;
    if (uT      == 0.0) uT      = BENTCH_DEFAULT_UT;
    if (dp      == 0.0) dp     = BENTCH_DEFAULT_DP;
    if (alpha_i == 0.0) alpha_i = BENTCH_DEFAULT_ALPHA;
    if (beta_i  == 0.0) beta_i  = BENTCH_DEFAULT_BETA;
    // Maybe -1, because R can accually be 0.
    if (R       == -1.0) R       = BENTCH_DEFAULT_R;

    // Initial phase-space coordinates
    const double x0  = LocalParticle_get_x(part);
    const double px0 = LocalParticle_get_xp(part);

    // -----------------------------
    //  Beam & particle physics
    // -----------------------------
    const double p0c          = LocalParticle_get_p0c(part);
    const double delta        = LocalParticle_get_delta(part);
    const double rvv          = LocalParticle_get_rvv(part);
    const double beta0        = LocalParticle_get_beta0(part);
    const double charge_ratio = LocalParticle_get_charge_ratio(part);
    const double chi          = LocalParticle_get_chi(part);

    // Effective beam momentum factor (bpc) - so far seems to not work
    double bpc = beta0 * rvv * (1.0 + delta) * p0c* charge_ratio / chi;
    //double bpc = beta0*rvv*p0c;
    //const double bpc = 150e9;


    // Working variables
    double x  = x0;
    double px = px0;

    // ============================================================
    //   STEP SUBDIVISION
    // ============================================================
    if (n_steps <= 0)
        n_steps = 1;

    const double ds = length / (double)n_steps;

    for (int i = 0; i < n_steps; ++i) {

        if (method == 2) {

            fM2_apply_yoshida(
                ds, Umax, aTF, alpha_i, beta_i, dp,
                x, px, bpc, R, order, variant, &x, &px
            );
        }
        else if (method == 3) {

            fM3_apply_yoshida(
                ds, Umax, aTF, alpha_i, beta_i, dp,
                x, px, bpc, R, order, variant,
                &x, &px
            );
        }
        else if (method == 4) {

            fM4_apply_yoshida(
                ds, Umax, aTF, alpha_i, beta_i, dp,
                x, px, bpc, R, order, variant,
                &x, &px
            );
        }
    }

    // ------------------------------------------------------------
    //  should i do something for s as well?
    // yes! I am doing it exactly like the magnets do, updating s and zeta.
   
    // ------------------------------------------------------------
    LocalParticle_set_x(part,  x);
    LocalParticle_set_xp(part, px);
    LocalParticle_add_to_s(part, length);
    LocalParticle_add_to_zeta(part, length);
}



// ================================================================
//  Multi-particle version 
// ================================================================
GPUFUN
void track_bent_channelling_particles(
        BentChannellingDevData el,
        LocalParticle* part0)
{
    
    START_PER_PARTICLE_BLOCK(part0, part);

        track_bent_channelling_body_single_particle(
            part, el
        );

    END_PER_PARTICLE_BLOCK;
}

#endif // TRACK_BENT_CHANNELLING_H



