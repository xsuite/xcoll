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

    int8_t method  = BentChannellingDevData_get_method(el);
    int8_t order   = BentChannellingDevData_get_order(el);
    int8_t variant = BentChannellingDevData_get_variant(el);

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

    int n_steps = BentChannellingDevData_get_n_steps(el);
    if (n_steps < 0){
        double _n_steps_auto = BentChannellingDevData_get__n_steps_auto(el);
        n_steps = fmax(3, ceil(_n_steps_auto/sqrt(bpc)));
    }

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

#endif // TRACK_BENT_CHANNELLING_H



