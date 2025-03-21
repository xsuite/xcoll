// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_TRAJ_CIRCULAR_H
#define XCOLL_GEOM_TRAJ_CIRCULAR_H

#include <stdio.h>
#include <math.h>


// /*gpufun*/
// void CircularTrajectory_set_params(CircularTrajectory traj, double sR, double xR,
//                                         LocalParticle part){
//     CircularTrajectory_set_sR(traj, sR);
//     CircularTrajectory_set_xR(traj, xR);
//     double s0 = LocalParticle_get_s(part);
//     double x0 = LocalParticle_get_x(part);
//     double R = sqrt((s0-sR)*(s0-sR) + (x0-xR)*(x0-xR));
//     CircularTrajectory_set_R(traj, R);
//     CircularTrajectory_set_sin_tI(traj, (x0-xR)/R);
//     CircularTrajectory_set_cos_tI(traj, (s0-sR)/R);
//     CircularTrajectory_set_tan_tI(traj, (x0-xR)/(s0-sR));
// }


// TODO: maybe for this trajectory it is faster to use tI directly instead of sin_tI and cos_tI

/*gpufun*/
void CircularTrajectory_set_params(CircularTrajectory traj, double sR, double xR,
                                     double s0, double x0){
    CircularTrajectory_set_sR(traj, sR);
    CircularTrajectory_set_xR(traj, xR);
    double R = sqrt((s0-sR)*(s0-sR) + (x0-xR)*(x0-xR));
    CircularTrajectory_set_R(traj, R);
    CircularTrajectory_set_sin_tI(traj, (x0-xR)/R);
    CircularTrajectory_set_cos_tI(traj, (s0-sR)/R);
    CircularTrajectory_set_tan_tI(traj, (x0-xR)/(s0-sR));
}

/*gpufun*/
double CircularTrajectory_func_s(CircularTrajectory traj, double l){
    double R = CircularTrajectory_get_R(traj);
    double sR = CircularTrajectory_get_sR(traj);
    double sin_tI = CircularTrajectory_get_sin_tI(traj);
    double cos_tI = CircularTrajectory_get_cos_tI(traj);
    return sR + R*cos(l)*cos_tI - R*sin(l)*sin_tI; // s(𝜆) = sR + R cos(𝜆 + 𝜃I)
}

/*gpufun*/
double CircularTrajectory_func_x(CircularTrajectory traj, double l){
    double R = CircularTrajectory_get_R(traj);
    double xR = CircularTrajectory_get_xR(traj);
    double sin_tI = CircularTrajectory_get_sin_tI(traj);
    double cos_tI = CircularTrajectory_get_cos_tI(traj);
    return xR + R*sin(l)*cos_tI + R*cos(l)*sin_tI; // s(𝜆) = sR + R sin(𝜆 + 𝜃I)
}

/*gpufun*/
double CircularTrajectory_func_xp(CircularTrajectory traj, double l){
    double tan_tI = CircularTrajectory_get_tan_tI(traj);
    return (tan_tI + tan(l)) / (1. - tan_tI*tan(l)); // 𝜃(𝜆) = 𝜃I + 𝜆 + chan. effects
}

/*gpufun*/
double CircularTrajectory_deriv_s(CircularTrajectory traj, double l){
    double R = CircularTrajectory_get_R(traj);
    double sin_tI = CircularTrajectory_get_sin_tI(traj);
    double cos_tI = CircularTrajectory_get_cos_tI(traj);
    return -R*sin(l)*cos_tI - R*cos(l)*sin_tI; // s(𝜆) = sR + R cos(𝜆 + 𝜃I)
}

/*gpufun*/
double CircularTrajectory_deriv_x(CircularTrajectory traj, double l){
    double R = CircularTrajectory_get_R(traj);
    double sin_tI = CircularTrajectory_get_sin_tI(traj);
    double cos_tI = CircularTrajectory_get_cos_tI(traj);
    return R*cos(l)*cos_tI - R*sin(l)*sin_tI; // s(𝜆) = sR + R sin(𝜆 + 𝜃I)
}

/*gpufun*/
void CircularTrajectory_bounding_box(CircularTrajectory traj, double l1, double l2, BoundingBox* box){
    // Take note that these functions, bounding boxes, expect l1 and l2 in radians! 
    // The plot function will automatically change [0,1] to [0,2pi] for you.
    double s1 = CircularTrajectory_func_s(traj, l1);
    double s2 = CircularTrajectory_func_s(traj, l2);
    double x1 = CircularTrajectory_func_x(traj, l1);
    double x2 = CircularTrajectory_func_x(traj, l2);
    double sR = CircularTrajectory_get_sR(traj);
    double R  = CircularTrajectory_get_R(traj);
    double sin_t = CircularTrajectory_get_sin_tI(traj);
    double cos_t = CircularTrajectory_get_cos_tI(traj);
    double dx = x2 - x1;
    double ds = s2 - s1;
    double euclidean_length = sqrt(dx*dx + ds*ds);
    double sin_p, cos_p;

    if (sin_t < 0){   // if theta is larger than 180 degrees, theta = theta - 180
        sin_t = -sin_t;
        cos_t = -cos_t;
    }
    if (cos_t < 1){   // if theta is larger than 90 degrees, phi = theta + 90 
        sin_p = cos_t;
        cos_p = -sin_t;
    } else {          // if theta is between 0 and 90 degrees, phi = theta - 90
        sin_p = -cos_t;
        cos_p = sin_t;
    }
    if ((l2-l1) < M_PI){ // delta theta of box less than 180.
        box->l  = euclidean_length;   // length of the box
        box->w  = R - sqrt(R*R - box->l*box->l/4.);                  // width of the box, Sagitta
        // finding the first vertex. 
        if (x1 < x2){
            if (((l1 < M_PI && l2 > M_PI) && (cos_t < 0)) || ((-M_PI < l1 < 0.) && (-M_PI < l2 < 0.)) || ((l1 < 0 && l2 > 0) && (cos_t > 0))){                // l1 or l2 is lower vertex
                box->rC = sqrt( (s1+box->w*cos_p)*(s1+box->w*cos_p) +    // length of the position vector to the first vertex
                                (x1+box->w*sin_p)*(x1+box->w*sin_p) );
            } else {
                box->rC = sqrt(s1*s1 + x1*x1); // length of position vector to first vertex
            }
        } else {
            if (((l1 < M_PI && l2 > M_PI) && (cos_t < 0)) || ((-M_PI < l1 < 0.) && (-M_PI < l2 < 0.)) || ((l1 < 0 && l2 > 0) && (cos_t > 0))){                // l1 or l2 is lower vertex
                box->rC = sqrt( (s2+box->w*cos_p)*(s2+box->w*cos_p) +    // length of the position vector to the first vertex
                                (x2+box->w*sin_p)*(x2+box->w*sin_p) );
            } else {
                box->rC = sqrt(s2*s2 + x2*x2); // length of position vector to first vertex
            }
        }
    } else {
        box->l = 2*R;
        box->w = R + sqrt(R*R - box->l*box->l/4.); 
        double chord_side_p1_x, chord_side_p1_s, chord_side_p2_x, chord_side_p2_s;
        if (x1 < x2){
            chord_side_p1_x = x1 - (2*R - box->l)/2.*sin_t;
            chord_side_p1_s = s1 - (2*R - box->l)/2.*cos_t;
            chord_side_p2_x = x2 + (2*R - box->l)/2.*sin_t;
            chord_side_p2_s = s2 + (2*R - box->l)/2.*cos_t;
        } else {
            chord_side_p1_x = x2 - (2*R - box->l)/2.*sin_t;
            chord_side_p1_s = s2 - (2*R - box->l)/2.*cos_t;
            chord_side_p2_x = x1 + (2*R - box->l)/2.*sin_t;
            chord_side_p2_s = s1 + (2*R - box->l)/2.*cos_t;
        }
        double p3_x = chord_side_p1_x + box->w*sin_t;
        double p3_s = chord_side_p1_s + box->w*cos_t;
        double p4_x = chord_side_p2_x + box->w*sin_t;
        double p4_s = chord_side_p2_s + box->w*cos_t;
        // Compare with the other three points
        double min_x = chord_side_p1_x;
        double min_s = chord_side_p1_s;
        if (chord_side_p2_x < min_x) { 
            min_x = chord_side_p2_x; min_s = chord_side_p2_s;  
        }
        if (p3_x < min_x) {
            min_x = p3_x; min_s = p3_s;
        }
        if (p4_x < min_x) {
            min_x = p4_x; min_s = p4_s; 
        }

        box->rC = sqrt(min_s*min_s + min_x*min_x); // length of position vector to first vertex
    }
    box->sin_tC = x1 / box->rC;    // angle of position vector to first vertex
    box->cos_tC = s1 / box->rC;
    box->sin_tb = sin_t;           // orientation of the box (angle of length wrt horizontal)
    box->cos_tb = cos_t;
    box->proj_l = box->rC * (cos_t*box->cos_tC + sin_t*box->sin_tC); // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    box->proj_w = box->rC * (cos_t*box->sin_tC - sin_t*box->cos_tC); // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)

}


#endif /* XCOLL_GEOM_TRAJ_CIRCULAR_H */