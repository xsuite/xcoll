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
    
    // Geometry and Position collimator-like parameters
    double width  = BentChannellingDevData_get_width(el);
    double height = BentChannellingDevData_get_height(el);
    double angle  = BentChannellingDevData_get_angle(el);
    double jaw_U  = BentChannellingDevData_get_jaw_U(el);
    double tilt   = BentChannellingDevData_get_tilt(el);
    
    
    // Initial phase-space coordinates
    const double x0  = LocalParticle_get_x(part);
    const double px0 = LocalParticle_get_xp(part);
    
    // We need y coordinates as well, because sometimes the channelling can happen in these planes instead.
    const double y0  = LocalParticle_get_y(part);
    const double py0 = LocalParticle_get_yp(part);


    // ------------------------------------------------------------
    // Transform from machine frame to collimator frame
    // ------------------------------------------------------------
    // The collimator angle is given in degrees.
    // In the collimator frame, x_col is the active jaw coordinate.
    // For angle = 90 deg, x_col is approximately y.

    double alpha = angle * 3.14159265358979323846 / 180.0;
    double ca = cos(alpha);
    double sa = sin(alpha);

    double x_col  =  x0 * ca + y0 * sa;
    double y_col  = -x0 * sa + y0 * ca;

    double px_col =  px0 * ca + py0 * sa;
    double py_col = -px0 * sa + py0 * ca;

    // ------------------------------------------------------------
    // Jaw-local coordinates
    // ------------------------------------------------------------
    // For a left jaw, the crystal material starts at x_col = jaw_U.
    // Therefore q_in > 0 means that the particle is already inside
    // the crystal material at the element entrance.

    double q_in = x_col - jaw_U;
    double p0   = px_col - tilt;

    // Passive coordinate: used only for checking the finite crystal height.
    double q_passive = y_col;

    // ------------------------------------------------------------
    // First-impact approximation
    // ------------------------------------------------------------

    double q0 = q_in;          // local depth inside the crystal
    double s_impact = 0.0;     // distance travelled before impact
    double length_eff = length;
    int hits_crystal = 0;

    // Case 1: already inside the crystal at the element entrance
    if ((q_in >= 0.0) && (q_in <= width)) {

        hits_crystal = 1;
        q0 = q_in;
    }

    // Case 2: outside at entrance, but moving towards the left jaw
    else if ((q_in < 0.0) && (p0 > 0.0)) {

        s_impact = -q_in / p0;

        if ((s_impact >= 0.0) && (s_impact <= length)) {

            hits_crystal = 1;

            // At first impact, the particle is at the crystal surface.
            q0 = 0.0;

            // Only the remaining crystal length contributes to the interaction.
            length_eff = length - s_impact;

            // Drift to the impact point.
            LocalParticle_add_to_s(part, s_impact);
            LocalParticle_add_to_zeta(part, s_impact);
        }
    }

    // If the particle does not hit the crystal, or is outside the passive size,
    // it passes through without interaction.
    if ((!hits_crystal) || (fabs(q_passive) > 0.5 * height)) {

        LocalParticle_add_to_s(part, length_eff);
        LocalParticle_add_to_zeta(part, length_eff);
        return;
    }


    int n_steps = BentChannellingDevData_get_n_steps(el);
    if (n_steps < 0){
        double _n_steps_auto = BentChannellingDevData_get__n_steps_auto(el);
        n_steps = 1000; //temporary
    }

    // Finding all the local channel areas inside the crystal
    double channel_index_f = floor(q0 / dp);
    double channel_center  = (channel_index_f + 0.5) * dp;
    double x_local         = q0 - channel_center;
    
    
    // Checking if particle position allows channelling after all
    double x_c = sqrt(0.9) * 0.5 * dp;

    if (fabs(x_local) > x_c) {
    // Too close to an atomic plane: not channelled.
    // Temporary behaviour: pass through without channelling.
        LocalParticle_add_to_s(part, length_eff);
        LocalParticle_add_to_zeta(part, length_eff);
    return;
    }

// // Checking if particle energy allows channelling after all

    double arg = 0.5 * beta_over_aTF * x_local;
    double sinh_arg = sinh(arg);
    double U_local = 2.0 * U_N * sinh_arg * sinh_arg;
    
    double E_T = 0.5 * bpc * p0 * p0 + U_local;

    if (E_T > Umax) {
        LocalParticle_add_to_s(part, length_eff);
        LocalParticle_add_to_zeta(part, length_eff);
    return;
    }
 


    // Working variables
    double x  = x_local;
    double px = p0;

    // ============================================================
    //   STEP SUBDIVISION
    // ============================================================


    const double ds = length_eff / (double)n_steps;

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
    double q_final = channel_center + x;

    // Final coordinates in the collimator frame
    double x_col_out  = q_final + jaw_U;
    double bend_angle = length_eff / R;
    //it worked, after trial and error
    double px_col_out = px + bend_angle + tilt;
    
    // Passive coordinate is unchanged in this first approximation
    double y_col_out  = y_col;
    double py_col_out = py_col;

    // Transform back from collimator frame to machine frame
    double x_out  = x_col_out * ca - y_col_out * sa;
    double y_out  = x_col_out * sa + y_col_out * ca;

    double px_out = px_col_out * ca - py_col_out * sa;
    double py_out = px_col_out * sa + py_col_out * ca;

    LocalParticle_set_x(part,  x_out);
    LocalParticle_set_y(part,  y_out);
    LocalParticle_set_xp(part, px_out);
    LocalParticle_set_yp(part, py_out);
        
    
    LocalParticle_add_to_s(part, length_eff);
    LocalParticle_add_to_zeta(part, length_eff);
    }

#endif // TRACK_BENT_CHANNELLING_H



