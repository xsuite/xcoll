// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2026.                 #
// ######################################### #

#ifndef XCOLL_SWEEP_CAVITY_H
#define XCOLL_SWEEP_CAVITY_H

#include <xtrack/headers/track.h>
#include <xtrack/beam_elements/elements_src/track_rf.h>


GPUFUN
void SweepCavity_track_local_particle(SweepCavityData el, LocalParticle* part0){
    // Sweep data
    double const df_per_turn  = SweepCavityData_get_df_per_turn(el);
    CavityData parent = SweepCavityData_getp_cavity(el);
    int64_t const last_turn = SweepCavityData_get__last_turn(el);

    // Regular cavity data
    double const weight = 1;
    double length = CavityData_get_length(parent);
    double voltage = CavityData_get_voltage(parent);
    double frequency = CavityData_get_frequency(parent);
    double lag = CavityData_get_lag(parent);
    double harmonic = CavityData_get_harmonic(parent);
    // if (fabs(harmonic) < 1e-12){
    //     // Need to compute the frequency as correctly as possible to avoid issues with phase slip
    //     double const line_length = part0->line_length;
    //     double const beta0 = LocalParticle_get_beta0(part0);
    //     double const t_rev0 = line_length / (beta0 * C_LIGHT);
    //     double this_harmonic = roundf(frequency * t_rev0);
    //     frequency = this_harmonic / t_rev0;
    // }
    double transverse_voltage = 0.;
    double transverse_lag = 0;
    int64_t absolute_time = CavityData_get_absolute_time(parent);
    int64_t order = -1;
    GPUGLMEM const double* knl = NULL;
    GPUGLMEM const double* ksl = NULL;
    GPUGLMEM const double* pn = NULL;
    GPUGLMEM const double* ps = NULL;
    int64_t num_kicks = CavityData_get_num_kicks(parent);
    int8_t model = CavityData_get_model(parent);
    int8_t default_model = 6; // drift-kick-drift-expanded
    int8_t integrator = CavityData_get_integrator(parent);
    int8_t default_integrator = 3; // Uniform
    int64_t radiation_flag = 0;
    int64_t radiation_flag_parent = 0;
    double lag_taper = CavityData_get_lag_taper(parent);
    int64_t body_active = 1;

    double factor_knl_ksl = 1.0;
    uint8_t kill_energy_kick = LocalParticle_check_track_flag(
                    part0, XS_FLAG_KILL_CAVITY_KICK);

    // Backtracking
    double body_length;
    double factor_knl_ksl_body;

    #ifndef XTRACK_MULTIPOLE_NO_SYNRAD
        lag += lag_taper;
    #endif

    if (LocalParticle_check_track_flag(part0, XS_FLAG_BACKTRACK)) {
        body_length = -length;
        factor_knl_ksl_body = -factor_knl_ksl;
        voltage = -voltage;
        transverse_voltage = -transverse_voltage;
    } else {
        body_length = length;
        factor_knl_ksl_body = factor_knl_ksl;
    }

    if (integrator == 0){
        integrator = default_integrator;
    }
    if (model == 0){
        model = default_model;
    }
    if (model==-1){ // kick only
        integrator = 3; // uniform
        num_kicks = 1;
    }

    // Compute the number of kicks for auto mode
    if (num_kicks == 0) { // num_multipole_kicks = 0 means auto mode
        num_kicks = 1;
    }

    double k0_drift, k1_drift, h_drift, ks_drift;
    double k0_kick, k1_kick, h_kick;
    double k0_h_correction, k1_h_correction;
    int8_t kick_rot_frame;
    int8_t drift_model;
    configure_tracking_model(
        model,
        0, // k0
        0, // k1
        0, // h
        0, // ks
        &k0_drift,
        &k1_drift,
        &h_drift,
        &ks_drift,
        &k0_kick,
        &k1_kick,
        &h_kick,
        &k0_h_correction,
        &k1_h_correction,
        &kick_rot_frame,
        &drift_model
    );


    double const line_length = part0->line_length;
    double const beta0 = LocalParticle_get_beta0(part0);
    double const t_rev0 = line_length / (beta0 * C_LIGHT);
    START_PER_PARTICLE_BLOCK(part0, part);
        double df;
        int64_t at_turn = LocalParticle_get_at_turn(part);
        if (last_turn == -1){
            // Use df_per_turn mode
            df = df_per_turn * at_turn;
        } else {
            // Use df array mode (saturate at last value if beyond last_turn)
            if (at_turn > last_turn){
                at_turn = last_turn;
            }
            df = SweepCavityData_get_df(el, at_turn);
        }

        track_rf_body_single_particle(
            part,
            body_length * weight,
            voltage * weight,
            frequency + df,
            lag,
            harmonic,
            transverse_voltage * weight,
            transverse_lag,
            absolute_time,
            order,
            factor_knl_ksl_body * weight,
            knl,
            ksl,
            pn,
            ps,
            num_kicks,
            drift_model,
            integrator,
            kill_energy_kick
        );
    END_PER_PARTICLE_BLOCK;
}

#endif  // XCOLL_SWEEP_CAVITY_H