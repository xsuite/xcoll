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
        const double length,
        const double U0,
        const double eta,
        const double R,
        const double dp,
        const double aTF,
        const double uT,
        const double alpha_i,
        const double beta_i,
        int8_t method,
        int8_t order,
        int8_t variant)
{
    // Particle dead? skip.
    if (LocalParticle_get_state(part) <= 0)
        return;

    // Apply defaults (0 -> use default)
    if (method  == 0) method  = BENTCH_DEFAULT_METHOD;
    if (order   == 0) order   = BENTCH_DEFAULT_ORDER;
    if (variant == 0) variant = BENTCH_DEFAULT_VARIANT;
    // COPY PARAMETERS INTO MUTABLE LOCALS
double U0_      = U0;
double aTF_     = aTF;
double uT_      = uT;
double dp_      = dp;
double alpha_i_ = alpha_i;
double beta_i_  = beta_i;
double R_       = R;

// APPLY DEFAULTS (modify the locals, not the parameters!)
if (U0_      == 0.0) U0_      = BENTCH_DEFAULT_U0;
if (aTF_     == 0.0) aTF_     = BENTCH_DEFAULT_ATF;
if (uT_      == 0.0) uT_      = BENTCH_DEFAULT_UT;
if (dp_      == 0.0) dp_      = BENTCH_DEFAULT_DP;
if (alpha_i_ == 0.0) alpha_i_ = BENTCH_DEFAULT_ALPHA;
if (beta_i_  == 0.0) beta_i_  = BENTCH_DEFAULT_BETA;
if (R_       == 0.0) R_       = BENTCH_DEFAULT_R;

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

    // Effective beam momentum factor (bpc)
    //double bpc = beta0 * rvv * (1.0 + delta) * p0c* charge_ratio / chi;
    //double bpc = beta0*rvv*p0c;
    const double bpc = 150e9;


    // Working variables
    double x  = x0;
    double px = px0;

    // ============================================================
    //   STEP SUBDIVISION
    // ============================================================
    int n_steps = BENTCH_DEFAULT_NSTEPS;
    if (n_steps <= 0)
        n_steps = 1;

    const double ds = length / (double)n_steps;

    for (int i = 0; i < n_steps; ++i) {

        if (method == 2) {

            fM2_apply_yoshida(
                ds, x, px, bpc, order, variant, &x, &px
            );
        }
        else if (method == 3) {

            fM3_apply_yoshida(
                ds, x, px, bpc, order, variant,
                &x, &px
            );
        }
        else if (method == 4) {

            fM4_apply_yoshida(
                ds, x, px, bpc, order, variant,
                &x, &px
            );
        }
    }

    // ------------------------------------------------------------
    //  Write back to particle
    // ------------------------------------------------------------
    LocalParticle_set_x(part,  x);
    LocalParticle_set_xp(part, px);
}



// ================================================================
//  Multi-particle version USING XSUiTE MACROS (CORRECT!!!)
// ================================================================
GPUFUN
void track_bent_channelling_particles(
        BentChannellingDevData el,
        LocalParticle* part0)
{
    double  length   = BentChannellingDevData_get_length(el);
    double  U0       = BentChannellingDevData_get_U0(el);
    double  eta      = BentChannellingDevData_get_eta(el);
    double  R        = BentChannellingDevData_get_R(el);
    int8_t  method   = BentChannellingDevData_get_method(el);
    int8_t  order    = BentChannellingDevData_get_order(el);
    int8_t  variant  = BentChannellingDevData_get_variant(el);

    // parameters from the element 
    // In the future, i should include the material option,
    // so I can change them accordingly.
    double  aTF      = BentChannellingDevData_get_aTF(el);
    double  uT       = BentChannellingDevData_get_uT(el);      // maybe unused for now
    double  alpha_i  = BentChannellingDevData_get_alpha_i(el);
    double  beta_i   = BentChannellingDevData_get_beta_i(el);

    START_PER_PARTICLE_BLOCK(part0, part);

        track_bent_channelling_body_single_particle(
            part,
            length,
            U0, eta,
            R, 
            dp, aTF, uT,
             alpha_i, beta_i,
            method, order, variant
        );

    END_PER_PARTICLE_BLOCK;
}

#endif // TRACK_BENT_CHANNELLING_H



