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
        BentChannellingDevData el,
        // element-level
        double length,
        double Umax,
        double dp,
        double aTF,
        double alpha_i,
        double beta_i,
        double aTF_over_beta,
        double two_aTF_over_beta,
        double beta_over_aTF,
        double U_N,
        double R,
        // particle-level
        double bpc,
        double U,
        double sqrt_U)

{
    // Particle dead? skip.
    if (LocalParticle_get_state(part) <= 0)
        return;
    
    
    
 
    
    int8_t method  = BentChannellingDevData_get_method(el);
    int8_t order   = BentChannellingDevData_get_order(el);
    int8_t variant = BentChannellingDevData_get_variant(el);

    // Initial phase-space coordinates
    const double x0  = LocalParticle_get_x(part);
    const double px0 = LocalParticle_get_xp(part);


//-----to be deleted from here-------start
        
    
    

    //-----to be deleted from here-------end

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


    const double ds = length / (double)n_steps;

    for (int i = 0; i < n_steps; ++i) {

        if (method == 2) {

            fM2_apply_yoshida(
                ds, x, px, bpc, R, 
                aTF_over_beta, beta_over_aTF, U_N, U, sqrt_U, 
                order, variant, &x, &px
            );
        }
        else if (method == 3) {

            fM3_apply_yoshida(
                ds, x, px, bpc, R, 
                aTF_over_beta, beta_over_aTF, U_N, U, sqrt_U, 
                order, variant, &x, &px
            );
        }
        else if (method == 4) {

            fM4_apply_yoshida(
                ds, x, px, bpc, R, 
                beta_over_aTF, two_aTF_over_beta,U_N,
                order, variant, &x, &px
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



