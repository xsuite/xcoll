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
void CircularTrajectory_init_bounding_box(CircularTrajectory traj, BoundingBox box, double l1, double l2){
    double theta = atan2(CircularTrajectory_get_sin_tI(traj), CircularTrajectory_get_cos_tI(traj));
    double ll1 = theta + l1;
    double ll2 = ll1 + (l2-l1);
    double s1 = CircularTrajectory_func_s(traj, l1);
    double x1 = CircularTrajectory_func_x(traj, l1);
    printf("s1 = %f, x1 = %f\n", s1, x1);
    double s2 = CircularTrajectory_func_s(traj, l2);
    double x2 = CircularTrajectory_func_x(traj, l2);
    double sR = CircularTrajectory_get_sR(traj);
    double R  = CircularTrajectory_get_R(traj); 
    double dx = x2 - x1;
    double ds = s2 - s1;
    double chord_length = sqrt(dx*dx + ds*ds);
    double sin_chord = dx / chord_length;
    double cos_chord = ds / chord_length;
    double sin_t, cos_t;   // angle of the box wrt horizontal
    int8_t sign = 1;
    if ((cos_chord > 0.) && (sin_chord > 0.)){    // if 0 < chord angle < 90 deg, then chord angle = box angle
        sin_t = sin_chord;
        cos_t = cos_chord;
    } else {
        if (sin_chord < 0) {   // if theta is larger than 180 degrees, theta = theta - 180
            sin_chord = -sin_chord;
            cos_chord = -cos_chord;
        }
        if ((cos_chord > 0.) && (sin_chord > 0.)){    // if 0 < chord angle < 90 deg, then chord angle = box angle
            sin_t = sin_chord;
            cos_t = cos_chord;
        } else {
            sin_t = - cos_chord;   // box angle is 90 degrees less than chord angle
            cos_t = sin_chord;
            sign = -1;
        }
    }

    if ((ll2 - ll1) < M_PI){ // delta theta of box less than 180.
        BoundingBox_set_l(box, chord_length);   // length of the box
        BoundingBox_set_w(box, R - sqrt(R*R - chord_length*chord_length/4.));                  // width of the box, Sagitta
        double w = BoundingBox_get_w(box);
        // finding the first vertex. 
        if (x1 < x2){
            if (((l1 < M_PI && l2 > M_PI) && (cos_chord < 0)) || (( (-M_PI < l1) && (l1 < 0.) ) && ( (-M_PI < l2) && (l2 < 0.))) || ((l1 < 0 && l2 > 0) && (cos_chord > 0))){                // l1 or l2 is lower vertex
                BoundingBox_set_rC(box, sqrt( (s1+w*cos_t)*(s1+w*cos_t) +    // length of the position vector to the first vertex
                                              (x1+w*sin_t)*(x1+w*sin_t)) );
            } else {
                BoundingBox_set_rC(box, sqrt(s1*s1 + x1*x1)); // length of position vector to first vertex
            }
        } else {
            if (((l1 < M_PI && l2 > M_PI) && (cos_chord < 0)) || (( (-M_PI < l1) && (l1 < 0.) ) && ( (-M_PI < l2) && (l2 < 0.))) || ((l1 < 0 && l2 > 0) && (cos_chord > 0))){                // l1 or l2 is lower vertex
                BoundingBox_set_rC(box, sqrt( (s2+w*cos_t)*(s2+w*cos_t) +    // length of the position vector to the first vertex
                                              (x2+w*sin_t)*(x2+w*sin_t)) );
            } else {
                BoundingBox_set_rC(box, sqrt(s2*s2 + x2*x2)); // length of position vector to first vertex
            }
        }
    } else {
        BoundingBox_set_l(box, 2*R);                               // L is always on the side of the chord
        BoundingBox_set_w(box, R + sqrt(R*R - (chord_length*chord_length/4.))); // 2R - sagitta
        printf("2 r = %f\n", 2*R);
        printf("sagitta = %f\n", R - sqrt(R*R - chord_length*chord_length/4.));
        double chord_side_w1_x, chord_side_w1_s, chord_side_w2_x, chord_side_w2_s;
        double w3_x, w3_s, w4_x, w4_s;
        printf("w = %f\n", BoundingBox_get_w(box));
        printf("l = %f\n", BoundingBox_get_l(box));
        printf("\nchord length = %f\n", chord_length);
        printf("sin_chord = %f, cos_chord = %f\n", sin_chord, cos_chord);
        printf("x1 = %f, x2 = %f\n", x1, x2);
        // determining the (s,x) of the vertices on the chord side of box
        if (x1 < x2){
            chord_side_w1_x = x1 - (2*R - chord_length)/2.*sin_t;
            chord_side_w1_s = s1 - sign*(2*R - chord_length)/2.*cos_t;
            chord_side_w2_x = x2 + (2*R - chord_length)/2.*sin_t;
            chord_side_w2_s = s2 + sign*(2*R - chord_length)/2.*cos_t;
            w3_x = chord_side_w1_x - sign*BoundingBox_get_w(box)*sin_t; // sign because of the orientation of the box
            w3_s = chord_side_w1_s + BoundingBox_get_w(box)*cos_t;
            w4_x = chord_side_w2_x - sign*BoundingBox_get_w(box)*sin_t;
            w4_s = chord_side_w2_s + BoundingBox_get_w(box)*cos_t;
        } else {
            printf("x1 = %f, x2 = %f\n", x1, x2);
            printf("s1 = %f, s2 = %f\n", s1, s2);
            printf("chord_length = %f\n", chord_length);
            printf("2r - chord_length = %f\n", 2*R - chord_length);
            printf("2r - chord_length /2 = %f\n", (2*R - chord_length)/2.);
            printf("sint = %f, cost = %f\n", sin_t, cos_t);
            printf("L = %f\n", (2*R - chord_length)/2.*sin_t);
            printf("L = %f\n", (2*R - chord_length)/2.*cos_t);
            chord_side_w2_x = x2 - (2*R - chord_length)/2.*cos_t;
            chord_side_w2_s = s2 - sign*(2*R - chord_length)/2.*sin_t;
            chord_side_w1_x = x1 + (2*R - chord_length)/2.*cos_t;
            chord_side_w1_s = s1 + sign*(2*R - chord_length)/2.*sin_t;
            w3_x = chord_side_w1_x + sign*BoundingBox_get_w(box)*sin_t; 
            w3_s = chord_side_w1_s - BoundingBox_get_w(box)*cos_t;
            w4_x = chord_side_w2_x + sign*BoundingBox_get_w(box)*sin_t;
            w4_s = chord_side_w2_s - BoundingBox_get_w(box)*cos_t;
            printf("cos_ti = %f, sin_ti = %f\n", cos_t, sin_t);
            printf("chord_side_w1_x = %f, chord_side_w1_s = %f\n", chord_side_w1_x, chord_side_w1_s);
        }
        printf("sign = %d\n", sign);
        printf("w3_x = %f, w3_s = %f\n", w3_x, w3_s);
        printf("w4_x = %f, w4_s = %f\n", w4_x, w4_s);
        printf("chord_side_w2_x = %f, chord_side_w2_s = %f\n", chord_side_w2_x, chord_side_w2_s);
        // Compare with the other three points
        double min_x = chord_side_w1_x;
        double min_s = chord_side_w1_s;
        if (chord_side_w2_x < min_x) { 
            min_x = chord_side_w2_x; 
            min_s = chord_side_w2_s;  
        }
        if (w3_x < min_x) {
            min_x = w3_x; 
            min_s = w3_s;
        }
        if (w4_x < min_x) {
            min_x = w4_x; 
            min_s = w4_s; 
        }
        BoundingBox_set_rC(box, sqrt(min_s*min_s + min_x*min_x)); // length of position vector to first vertex
    }
    printf("BoundingBox_get_rC(box) = %f\n", BoundingBox_get_rC(box));
    double rC = BoundingBox_get_rC(box);
    BoundingBox_set_sin_tC(box, x1 / rC);   // angle of position vector to first vertex
    BoundingBox_set_cos_tC(box, s1 / rC);
    BoundingBox_set_sin_tb(box, sin_t);           // orientation of the box (angle of length wrt horizontal)
    BoundingBox_set_cos_tb(box, cos_t);
    BoundingBox_set_proj_l(box, rC * (cos_t*s1/rC + sin_t*x1/rC)); // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    BoundingBox_set_proj_w(box, rC * (cos_t*x1/rC - sin_t*s1/rC)); // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
}


#endif /* XCOLL_GEOM_TRAJ_CIRCULAR_H */