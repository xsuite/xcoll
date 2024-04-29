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

    InteractionRecordData record = everest->coll->record;
    RecordIndex record_index     = everest->coll->record_index;
    InteractionRecordData_log(record, record_index, part, XC_ENTER_JAW);
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

#endif /* XCOLL_EVEREST_CRYSTAL_INTERACT_H */
