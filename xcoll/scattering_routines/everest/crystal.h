// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2023.                 #
// ######################################### #

#ifndef XCOLL_EVEREST_CRYSTAL_INTERACT_H
#define XCOLL_EVEREST_CRYSTAL_INTERACT_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/*gpufun*/
double interact(EverestData restrict everest, LocalParticle* part, double pc, double length) {

    double xdim     = everest->coll->xdim;
    double ydim     = everest->coll->ydim;
    double alayer   = everest->coll->amorphous_layer;

    double const x   = LocalParticle_get_x(part);
    double const xp  = LocalParticle_get_xp(part);
    double const y   = LocalParticle_get_y(part);

    double ymin = -ydim/2.;
    double ymax =  ydim/2.;

    // FIRST CASE: p don't interact with crystal
    if (y < ymin || y > ymax || x > xdim) {
        // TODO: is_hit is wrong
        Drift_single_particle_4d(part, length);
        return pc;

    } else if (x < alayer || y-ymin < alayer || ymax-y < alayer) {
    // SECOND CASE: p hits the amorphous layer
    // TODO: NOT SUPPORTED
        LocalParticle_kill_particle(part, XC_ERR_NOT_IMPLEMENTED);
        return pc;
//         double x0    = x;
//         double y0    = y;
//         double a_eq  = 1. + pow(xp,2.);
//         double b_eq  = (2.*x)*xp - (2.*xp)*bend_r;
//         double c_eq  = pow(x,2.) - (2.*x)*bend_r;
//         double delta = pow(b_eq,2.) - (4.*a_eq)*c_eq;
//         s = (-b_eq+sqrt(delta))/(2.*a_eq);
//         if (s >= length) {
//             s = length;
//         }
//         x   =  xp*s + x0;
//         double len_xs = sqrt(pow((x-x0),2.) + pow(s,2.));
//         double len_ys;
//         if (yp >= 0 && y + yp*s <= ymax) {
//             len_ys = yp*len_xs;
//         } else if (yp < 0 && y + yp*s >= ymin) {
//             len_ys = yp*len_xs;
//         } else {
//             s      = (ymax-y)/yp;
//             len_ys = sqrt(pow((ymax-y),2.) + pow(s,2.));
//             x   = x0 + xp*s;
//             len_xs = sqrt(pow((x-x0),2.) + pow(s,2.));
//         }
        
//         double am_len = sqrt(pow(len_xs,2.) + pow(len_ys,2.));
//         s     = s/2;
//         x  = x0 + xp*s;
//         y     = y0 + yp*s;

//         energy_loss = calcionloss(part, length, properties);

//         double* result_am = moveam(rng, part, am_len, energy_loss, pc, material);
//         pc = result_am[0];
//         iProc = result_am[1];
//         free(result_am);

//         x = x + xp*(length-s);
//         y = y + yp*(length-s);

//         result[0] = x;
//         result[1] = xp;
//         result[2] = y;
//         result[3] = yp;
//         result[4] = pc;
//         result[5] = iProc;
//         return result;

    } else if (x > xdim-alayer && x < xdim) {
    // TODO: NOT SUPPORTED
        LocalParticle_kill_particle(part, XC_ERR_NOT_IMPLEMENTED);
        return pc;
//         iProc = proc_AM;
        
//         energy_loss = calcionloss(part, length, properties);

//         double* result_am = moveam(rng, part, length, energy_loss, pc, material);
//         pc = result_am[0];
//         iProc = result_am[1];
//         free(result_am);

//         result[0] = x;
//         result[1] = xp;
//         result[2] = y;
//         result[3] = yp;
//         result[4] = pc;
//         result[5] = iProc;
//         return result;
    }

    CollimatorImpactsData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    CollimatorImpactsData_log(record, record_index, part, XC_ENTER_JAW);
    calculate_initial_angle(everest, part);
// printf("Start:  %f  %f  %f\n", xp, everest->t_I, everest->t_c);

    if (fabs(xp - everest->t_I) < everest->t_c) {
        pc = Channel(everest, part, pc, length);
    } else {
        pc = Amorphous(everest, part, pc, length);
    }
// printf("End\n\n");
    return pc;
}


/*gpufun*/
double* crystal(EverestData restrict everest, LocalParticle* part, double p, double length) {

    double* crystal_result = (double*)malloc(2 * sizeof(double));

    double cry_tilt = everest->coll->tilt;
    double bend_r   = everest->coll->bend_r;
    double bend_ang = everest->coll->bend_ang;
    double xdim     = everest->coll->xdim;
    double miscut   = everest->coll->miscut;
    double const cry_cBend  = cos(bend_ang);
    double const cry_sBend  = sin(bend_ang);
    double const cry_cpTilt = cos(cry_tilt);
    double const cry_spTilt = sin(cry_tilt);

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
    int is_hit = 0;
    double x = LocalParticle_get_x(part);
    double rpp_in  = LocalParticle_get_rpp(part);
    double xp = LocalParticle_get_px(part)*rpp_in;
    if (x >= 0. && x < xdim) {
        is_hit = 1;
        p = interact(everest, part, p, length);
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
                p = interact(everest, part, p, length-(tilt_int*bend_r));

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

    crystal_result[0] = is_hit;
    crystal_result[1] = p;
    return crystal_result;
}

#endif /* XCOLL_EVEREST_CRYSTAL_INTERACT_H */
