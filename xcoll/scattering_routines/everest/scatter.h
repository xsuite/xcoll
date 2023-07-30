// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_SCAT_H
#define XCOLL_EVEREST_SCAT_H
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// TODO:
//    Use drift function from xtrack Drift element (call its C function)
//    Use rotation function from xtrack XYRotation element (call its C function)

/*gpufun*/
void scatter(LocalParticle* part, double length, MaterialData material, RandomRutherfordData rng,
             struct ScatteringParameters scat, double aperture, double offset,
             double tilt_L, double tilt_R, double side, CollimatorImpactsData record, RecordIndex record_index){

    // Store initial coordinates for updating later
    double const rpp_in  = LocalParticle_get_rpp(part);
    double const rvv_in  = LocalParticle_get_rvv(part);
    double const e0      = LocalParticle_get_energy0(part) / 1e9; // Reference energy in GeV
    double const beta0   = LocalParticle_get_beta0(part);
    double const ptau_in = LocalParticle_get_ptau(part);
    double const x_in    = LocalParticle_get_x(part);
    double const px_in   = LocalParticle_get_px(part);
    double const y_in    = LocalParticle_get_y(part);
    double const py_in   = LocalParticle_get_py(part);
    double p0 = LocalParticle_get_p0c(part) / 1e9;

    // TODO: missing correction due to m/m0 (but also wrong in xpart...)
    double energy = p0*ptau_in + e0; // energy, not momentum, in GeV

    // Status flags
    int is_hit = 0;
    int is_abs = 0;

    double x = LocalParticle_get_x(part);

    // For one-sided collimators consider only positive X. For negative X jump to the next particle
    // TODO: need to drift!!
    if (side==0 || (side==1 && x>=0.) || (side==2 && x<=0.)) {

        double mirror = 1;
        double tiltangle = tilt_L;
        if (x < 0) {
            mirror = -1;
            tiltangle = -tilt_R;
            LocalParticle_scale_x(part, mirror);
            LocalParticle_scale_px(part, mirror);
        }

        // Shift with opening and offset
        LocalParticle_add_to_x(part, -aperture/2. - mirror*offset);

        // Include collimator tilt
        double rot_shift = YRotation_single_particle_rotate_only(part, 0., tiltangle);
        Drift_single_particle_4d(part, -rot_shift);
//         LocalParticle_add_to_px(part, -tiltangle/rpp_in);

        // Check if/where particle hits collimator
        double zlm; 
        x  = LocalParticle_get_x(part);
        double xp = LocalParticle_get_px(part)*rpp_in;
        if (x >= 0.) { // hit at front
            zlm = length;
            // val_part_impact = x;
            // store these in impact table, as interaction type 'hitting collimator' or something like that.
            // val_part_indiv  = xp;
        } else if (xp <= 0.) { // no hit
            zlm = 0.;
            Drift_single_particle_4d(part, length);
        } else { // hit from side
            double s = -x/xp;
            if (s < length) {
                zlm = length - s;
                Drift_single_particle_4d(part, s);
                // val_part_impact = 0.;
                // val_part_indiv  = xp;
            } else {
                zlm = 0.;
                Drift_single_particle_4d(part, length);
            }
        }

        // Now do the scattering part
        if (zlm > 0.) {
            is_hit = 1;

            double* jaw_result = jaw(part, material, rng, scat, energy, zlm, record, record_index);

            energy = jaw_result[0];
            if (jaw_result[1] == 1){
                is_abs = 1;
            }
            double s_out = jaw_result[2];
            free(jaw_result);

            if (is_abs != 1) {
            // Do the rest drift, if particle left collimator early
                Drift_single_particle_4d(part, zlm-s_out);
            }
        }

        // Transform back to particle coordinates with opening and offset
        // Include collimator tilt
        rot_shift = YRotation_single_particle_rotate_only(part, length, -tiltangle);
        Drift_single_particle_4d(part, length-rot_shift);

        // Transform back to particle coordinates with opening and offset
        LocalParticle_add_to_x(part, aperture/2. + mirror*offset);

        // Now mirror at the horizontal axis for negative X offset
        if (mirror < 0) {
            LocalParticle_scale_x(part, mirror);
            LocalParticle_scale_px(part, mirror);
        }
    } else {
        Drift_single_particle_4d(part, length);
    }

    // Cannot assign energy directly to LocalParticle as it would update dependent variables, but needs to be corrected first!

    // Update energy    ---------------------------------------------------
    // Only particles that hit the jaw and survived need to be updated
    if (is_hit>0 && is_abs==0){
        double ptau_out = (energy - e0) / (e0 * beta0);
        LocalParticle_update_ptau(part, ptau_out);
    }

    // Update 4D coordinates    -------------------------------------------
    // Absorbed particles get their coordinates set to the entrance of collimator
    if (is_abs>0){
        LocalParticle_set_x(part, x_in);
        LocalParticle_set_px(part, px_in);
        LocalParticle_set_y(part, y_in);
        LocalParticle_set_py(part, py_in);
    }

    // Update longitudinal coordinate zeta    -----------------------------
    // Absorbed particles keep coordinates at the entrance of collimator, others need correcting:
    // Non-hit particles are just drifting (zeta not yet drifted in K2, so do here)
    if (is_hit==0){
        LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, px_in*rpp_in, py_in*rpp_in, length) );
    }
    // Hit and survived particles need correcting:
    if (is_hit>0 && is_abs==0){
        double px  = LocalParticle_get_px(part);
        double py  = LocalParticle_get_py(part);
        double rvv = LocalParticle_get_rvv(part);
        double rpp = LocalParticle_get_rpp(part);
        // First we drift half the length with the old angles:
        LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, px_in*rpp_in, py_in*rpp_in, length/2) );
        // then half the length with the new angles:
        LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, px*rpp, py*rpp, length/2) );
    }

    // Update s    --------------------------------------------------------
    if (is_abs==0){
        LocalParticle_add_to_s(part, length);
    }

    // Update state    ----------------------------------------------------
    if (is_abs > 0){
        LocalParticle_set_state(part, XC_LOST_ON_EVEREST);
    }

    return;
}

#endif /* XCOLL_EVEREST_SCAT_H */
