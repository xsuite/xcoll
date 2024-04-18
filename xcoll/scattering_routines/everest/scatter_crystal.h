// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_SCAT_CRY_H
#define XCOLL_EVEREST_SCAT_CRY_H
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// TODO: we would like to pass EverestCrystalData as to avoid all variables, however, this part of code
//       is not yet generated when the EverestEngine code is generated, so the compiler does not find it.
//       Workaround: make a struct (with the definition in EverestEngine) that EverestCrystal uses...?

/*gpufun*/
void scatter_cry(EverestData restrict everest, LocalParticle* part, double length){

    // geometry values
    double aperture = everest->coll->aperture;
    double offset   = everest->coll->offset;
    double side     = everest->coll->side;
    double cry_tilt = everest->coll->tilt;
    double bend_r   = everest->coll->bend_r;
    double bend_ang = everest->coll->bend_ang;
    double xdim     = everest->coll->xdim;
    double miscut   = everest->coll->miscut;
    double const cry_cBend  = cos(bend_ang);
    double const cry_sBend  = sin(bend_ang);
    double const cry_cpTilt = cos(cry_tilt);
    double const cry_spTilt = sin(cry_tilt);

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
    double const p0      = LocalParticle_get_p0c(part) / 1e9;

    // TODO: missing correction due to m/m0 (but also wrong in xpart...)
    double energy = p0*ptau_in + e0; // energy in GeV

    // Status flags
    int is_hit = 0;

    double x = LocalParticle_get_x(part);

    // For one-sided collimators consider only positive X. For negative X jump to the next particle
    if (side==0 || (side==1 && x>=0.) || (side==2 && x<=0.)) {

        // Now mirror at the horizontal axis for negative X offset
        double mirror = 1;
        if (x < 0) {
            mirror = -1;
            LocalParticle_scale_x(part, mirror);
            LocalParticle_scale_px(part, mirror);
        }

        // Shift with opening and offset
        LocalParticle_add_to_x(part, -aperture/2.);

        // particle passing above the jaw are discarded => take new event
        // entering by the face, shorten the length (zlm) and keep track of
        // entrance longitudinal coordinate (keeps) for histograms

        // The definition is that the collimator jaw is at x>=0.

        // Move origin of x to inner front corner (transformation 4 in Figure 3.3 of thesis Valentina Previtali)
        double shift = 0;
        if (cry_tilt < 0) {
            shift = bend_r*(1 - cry_cpTilt);
            if (cry_tilt < -bend_ang) {
                shift = bend_r*(cry_cpTilt - cos(bend_ang - cry_tilt));
            }
            LocalParticle_add_to_x(part, -shift);
        } 

        // Rotate tilt (transformation 5 in Figure 3.3 of thesis Valentina Previtali)
        double s = YRotation_single_particle_rotate_only(part, 0., cry_tilt);

        // 3rd transformation: drift to the new coordinate s=0
        Drift_single_particle_4d(part, -s);

        // Check that particle hit the crystal
        double x = LocalParticle_get_x(part);
        double rpp_in  = LocalParticle_get_rpp(part);
        double xp = LocalParticle_get_xp(part);
        if (x >= 0. && x < xdim) {
            is_hit = 1;
            energy = interact(everest, part, energy, length);
            s = bend_r*cry_sBend;

        } else {
            double xp_tangent=0;
            if (x < 0) { // Crystal can be hit from below
                xp_tangent = sqrt((-(2.*x)*bend_r + pow(x,2.))/(pow(bend_r,2.)));
            } else {             // Crystal can be hit from above
                xp_tangent = asin((bend_r*(1. - cry_cBend) - x)/sqrt(((2.*bend_r)*(bend_r - x))*(1 - cry_cBend) + pow(x,2.)));
            }
            // If the hit is below, the angle must be greater or equal than the tangent,
            // or if the hit is above, the angle must be smaller or equal than the tangent
            if ((x < 0. && xp >= xp_tangent) || (x >= 0. && xp <= xp_tangent)) {

                // If it hits the crystal, calculate in which point and apply the transformation and drift to that point
                double a_eq  = 1 + pow(xp,2.);
                double b_eq  = (2.*xp)*(x - bend_r);
                double c_eq  = -(2.*x)*bend_r + pow(x,2.);
                double delta = pow(b_eq,2.) - 4*(a_eq*c_eq);
                double s_int = (-b_eq - sqrt(delta))/(2*a_eq);

                // MISCUT first step: P coordinates (center of curvature of crystalline planes)
                double s_P_tmp = (bend_r-xdim)*sin(-miscut);
                double x_P_tmp = xdim + (bend_r-xdim)*cos(-miscut);

                if (s_int < bend_r*cry_sBend) {
                    // Transform to a new reference system: shift and rotate
                    double tilt_int = s_int/bend_r;
                    double x_int  = xp*s_int + x;
                    LocalParticle_add_to_y(part, LocalParticle_get_py(part)*rpp_in*s_int);
                    LocalParticle_set_x(part, 0.);
                    LocalParticle_add_to_px(part, -tilt_int/rpp_in);

                    // MISCUT first step (bis): transform P in new reference system
                    // Translation
                    s_P_tmp = s_P_tmp - s_int;
                    x_P_tmp = x_P_tmp - x_int;
                    // Rotation
                    double s_P = s_P_tmp*cos(tilt_int) + x_P_tmp*sin(tilt_int);
                    double x_P = -s_P_tmp*sin(tilt_int) + x_P_tmp*cos(tilt_int);

                    is_hit = 1;
                    energy = interact(everest, part, energy, length-(tilt_int*bend_r));

                    s = bend_r*sin(bend_ang - tilt_int);

                    // un-rotate
                    s = YRotation_single_particle_rotate_only(part, s, -tilt_int);

                    // 2nd: shift back the 2 axis
                    LocalParticle_add_to_x(part, x_int);
                    s = s + s_int;
                } else {
                    // Drift
                    s = bend_r*sin(length/bend_r);
                    Drift_single_particle_4d(part, s);
                }
            } else {
                // Drift
                s = bend_r*sin(length/bend_r);
                Drift_single_particle_4d(part, s);
            }
        }

        // transform back from the crystal to the collimator reference system
        // 1st: un-rotate the coordinates
        s = YRotation_single_particle_rotate_only(part, length, -cry_tilt);
        Drift_single_particle_4d(part, length-s);

        // 2nd: shift back the reference frame
        if (cry_tilt < 0) {
            LocalParticle_add_to_px(part, shift);
        }

        // Transform back to particle coordinates with opening and offset
        LocalParticle_add_to_x(part, aperture/2.);

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
    if (is_hit>0 && LocalParticle_get_state(part) > 0){
        double ptau_out = (energy - e0) / (e0 * beta0);
        LocalParticle_update_ptau(part, ptau_out);
    }

    // Update 4D coordinates    -------------------------------------------
    // Absorbed particles get their coordinates set to the entrance of collimator
    if (LocalParticle_get_state(part) < 0){
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
    if (is_hit>0 && LocalParticle_get_state(part) > 0){
        double px  = LocalParticle_get_px(part);
        double py  = LocalParticle_get_py(part);
        double rvv = LocalParticle_get_rvv(part);
        double rpp = LocalParticle_get_rpp(part);
        // First we drift half the length with the old angles:
        LocalParticle_add_to_zeta(part, drift_zeta_single(rvv_in, px_in*rpp_in, py_in*rpp_in, length/2) );
        // then half the length with the new angles:
        LocalParticle_add_to_zeta(part, drift_zeta_single(rvv, px*rpp, py*rpp, length/2) );
    }

    return;
}

#endif /* XCOLL_EVEREST_SCAT_CRY_H */
