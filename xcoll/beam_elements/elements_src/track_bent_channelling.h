// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef TRACK_BENT_CHANNELLING_H
#define TRACK_BENT_CHANNELLING_H

#include <stdint.h>


// Helpers for recording interactions

/*gpufun*/
int8_t BentChannellingDevData_get_record_impacts(BentChannellingDevData el){
    return BentChannellingDevData_get__record_interactions(el) % 2;
}

/*gpufun*/
int8_t BentChannellingDevData_get_record_exits(BentChannellingDevData el){
    return (BentChannellingDevData_get__record_interactions(el) >> 1) % 2;
}

/*gpufun*/
int8_t BentChannellingDevData_get_record_scatterings(BentChannellingDevData el){
    return (BentChannellingDevData_get__record_interactions(el) >> 2) % 2;
}

/*gpufun*/
void BentChannellingDev_log_in_jaw_frame(
        InteractionRecordData record,
        RecordIndex record_index,
        LocalParticle* part,
        int64_t interaction_type,
        double q_active,
        double p_active,
        double q_passive,
        double p_passive
){
    // Save lab-frame coordinates
    double x_lab  = LocalParticle_get_x(part);
    double px_lab = LocalParticle_get_xp(part);
    double y_lab  = LocalParticle_get_y(part);
    double py_lab = LocalParticle_get_yp(part);

    // Temporarily write jaw-frame coordinates:
    // active plane -> x/px, passive plane -> y/py.
    LocalParticle_set_x(part,  q_active);
    LocalParticle_set_xp(part, p_active);
    LocalParticle_set_y(part,  q_passive);
    LocalParticle_set_yp(part, p_passive);

    // Log in jaw frame, Everest-compatible convention.
    InteractionRecordData_log(record, record_index, part, interaction_type);

    // Restore lab-frame coordinates
    LocalParticle_set_x(part,  x_lab);
    LocalParticle_set_xp(part, px_lab);
    LocalParticle_set_y(part,  y_lab);
    LocalParticle_set_yp(part, py_lab);
}

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
    
    // ------------------------------------------------------------
    // Interaction record
    // ------------------------------------------------------------
    InteractionRecordData record = BentChannellingDevData_getp_internal_record(el, part);
    RecordIndex record_index = NULL;

    int8_t record_impacts = 0;
    int8_t record_scatterings = 0;

    if (record) {
        record_index = InteractionRecordData_getp__index(record);
        record_impacts = BentChannellingDevData_get_record_impacts(el);
        record_scatterings = BentChannellingDevData_get_record_scatterings(el);
    }
 
    
    int8_t method  = BentChannellingDevData_get_method(el);
    int8_t order   = BentChannellingDevData_get_order(el);
    int8_t variant = BentChannellingDevData_get_variant(el);
    
    // Geometry and Position collimator-like parameters
    double width  = BentChannellingDevData_get__width(el);
    double height = BentChannellingDevData_get__height(el);
    // angle
    double sin_z = BentChannellingDevData_get__sin_z(el);
    double cos_z = BentChannellingDevData_get__cos_z(el);
    int is_vertical = (fabs(sin_z) > fabs(cos_z));
    
    double jaw_U  = BentChannellingDevData_get__jaw_U(el);
    //double tilt   = BentChannellingDevData_get__tan_y(el);
    double tilt = atan(BentChannellingDevData_get__tan_y(el));
    
    // Initial phase-space coordinates
    const double x0  = LocalParticle_get_x(part);
    const double px0 = LocalParticle_get_xp(part);
    
    // We need y coordinates as well, because sometimes the channelling can happen in these planes instead.
    const double y0  = LocalParticle_get_y(part);
    const double py0 = LocalParticle_get_yp(part);



    // ------------------------------------------------------------
    // Select active plane and define jaw-local coordinates
    // ------------------------------------------------------------
    // For the TWOCRYST vertical TCCP case:
    //   active coordinate = y
    //   active angle      = yp

    // For a horizontal crystal:
    //   active coordinate = x
    //   active angle      = xp

    // Here, angle/orientation is used only to choose the active plane.
    // We do not rotate the coordinates explicitly.

    double q_in;
    double p0;
    double q_passive;
    double p_passive;

    if (is_vertical) {

        // Vertical crystal:
        // channelling happens in y/yp.
        q_in = y0 - jaw_U;
        p0   = py0 - tilt;

        // x is only used to check the finite transverse size.
        q_passive = x0;
        p_passive = px0;
    }
    else {

        // Horizontal crystal:
        // channelling happens in x/xp.
        q_in = x0 - jaw_U;
        p0   = px0 - tilt;

        // y is only used to check the finite transverse size.
        q_passive = y0;
        p_passive = py0;
    }

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



    // particle has entered the crystal jaw
    if (record && record_impacts) {
    BentChannellingDev_log_in_jaw_frame(
        record,
        record_index,
        part,
        XC_ENTER_JAW_L,
        q0,
        p0,
        q_passive,
        p_passive
    );
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




// // Checking if particle energy allows channelling after all plus Bending angular acceptance condition

    double arg = 0.5 * beta_over_aTF * x_local;
    double sinh_arg = sinh(arg);
    double U_local = 2.0 * U_N * sinh_arg * sinh_arg;

    // Critical bending radius estimate.
    // If R <= Rcrit, the effective potential well disappears.
    double Rcrit = bpc * (0.5 * dp) / (2.0 * Umax);

    if (R <= Rcrit) {
        LocalParticle_add_to_s(part, length_eff);
        LocalParticle_add_to_zeta(part, length_eff);
        return;
    }


    
    double U_bend = bpc * x_local / R;
    double U_eff = U_local + U_bend;


    // This avoids using the straight-crystal Umax when the well is tilted.
    double Umax_eff = Umax * (1.0 - Rcrit / R);

    if (Umax_eff <= 0.0) {
        LocalParticle_add_to_s(part, length_eff);
        LocalParticle_add_to_zeta(part, length_eff);
        return;
    }

    double E_T_eff = 0.5 * bpc * p0 * p0 + U_eff;

    if (E_T_eff > Umax_eff) {
        LocalParticle_add_to_s(part, length_eff);
        LocalParticle_add_to_zeta(part, length_eff);
        return;
    }
    
    // particle passed the channelling acceptance conditions
    if (record && record_scatterings) {
    BentChannellingDev_log_in_jaw_frame(
        record,
        record_index,
        part,
        XC_CHANNELLING,
        q0,
        p0,
        q_passive,
        p_passive
    );
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

    double bend_angle = length_eff / R;
    double p_final = px + bend_angle + tilt;

    if (is_vertical) {

        // Write back to y/yp for a vertical crystal.
        LocalParticle_set_y(part,  q_final + jaw_U);
        LocalParticle_set_yp(part, p_final);

        // x/xp are passive and remain unchanged in this approximation.
        LocalParticle_set_x(part,  x0);
        LocalParticle_set_xp(part, px0);
    }
    else {

        // Write back to x/xp for a horizontal crystal.
        LocalParticle_set_x(part,  q_final + jaw_U);
        LocalParticle_set_xp(part, p_final);

        // y/yp are passive and remain unchanged in this approximation.
        LocalParticle_set_y(part,  y0);
        LocalParticle_set_yp(part, py0);
    }
        
    
    LocalParticle_add_to_s(part, length_eff);
    LocalParticle_add_to_zeta(part, length_eff);
    }

#endif // TRACK_BENT_CHANNELLING_H



