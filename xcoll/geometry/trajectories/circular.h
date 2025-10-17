// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

// #ifndef XCOLL_GEOM_TRAJ_CIRCULAR_H
// #define XCOLL_GEOM_TRAJ_CIRCULAR_H

// #include <stdio.h>
// #include <math.h>


// // /*gpufun*/
// // void CircularTrajectory_set_params(CircularTrajectory traj, double sR, double xR,
// //                                         LocalParticle part){
// //     CircularTrajectory_set_sR(traj, sR);
// //     CircularTrajectory_set_xR(traj, xR);
// //     double s0 = LocalParticle_get_s(part);
// //     double x0 = LocalParticle_get_x(part);
// //     double R = sqrt((s0-sR)*(s0-sR) + (x0-xR)*(x0-xR));
// //     CircularTrajectory_set_R(traj, R);
// //     CircularTrajectory_set_sin_tI(traj, (x0-xR)/R);
// //     CircularTrajectory_set_cos_tI(traj, (s0-sR)/R);
// //     CircularTrajectory_set_tan_tI(traj, (x0-xR)/(s0-sR));
// // }


// // TODO: maybe for this trajectory it is faster to use tI directly instead of sin_tI and cos_tI

// /*gpufun*/
// void CircularTrajectory_set_params(CircularTrajectory traj, double sR, double xR,
//                                      double s0, double x0){
//     CircularTrajectory_set_sR(traj, sR);
//     CircularTrajectory_set_xR(traj, xR);
//     double R = sqrt((s0-sR)*(s0-sR) + (x0-xR)*(x0-xR));
//     CircularTrajectory_set_R(traj, R);
//     CircularTrajectory_set_sin_tI(traj, (x0-xR)/R);
//     CircularTrajectory_set_cos_tI(traj, (s0-sR)/R);
//     CircularTrajectory_set_tan_tI(traj, (x0-xR)/(s0-sR));
// }

// /*gpufun*/
// double CircularTrajectory_func_s(CircularTrajectory traj, double l){
//     double R = CircularTrajectory_get_R(traj);
//     double sR = CircularTrajectory_get_sR(traj);
//     double sin_tI = CircularTrajectory_get_sin_tI(traj);
//     double cos_tI = CircularTrajectory_get_cos_tI(traj);
//     return sR + R*cos(l)*cos_tI - R*sin(l)*sin_tI; // s(𝜆) = sR + R cos(𝜆 + 𝜃I)
// }

// /*gpufun*/
// double CircularTrajectory_func_x(CircularTrajectory traj, double l){
//     double R = CircularTrajectory_get_R(traj);
//     double xR = CircularTrajectory_get_xR(traj);
//     double sin_tI = CircularTrajectory_get_sin_tI(traj);
//     double cos_tI = CircularTrajectory_get_cos_tI(traj);
//     return xR + R*sin(l)*cos_tI + R*cos(l)*sin_tI; // s(𝜆) = sR + R sin(𝜆 + 𝜃I)
// }

// /*gpufun*/
// double CircularTrajectory_func_xp(CircularTrajectory traj, double l){
//     double tan_tI = CircularTrajectory_get_tan_tI(traj);
//     return (tan_tI + tan(l)) / (1. - tan_tI*tan(l)); // 𝜃(𝜆) = 𝜃I + 𝜆 + chan. effects
// }

// /*gpufun*/
// double CircularTrajectory_deriv_s(CircularTrajectory traj, double l){
//     double R = CircularTrajectory_get_R(traj);
//     double sin_tI = CircularTrajectory_get_sin_tI(traj);
//     double cos_tI = CircularTrajectory_get_cos_tI(traj);
//     return -R*sin(l)*cos_tI - R*cos(l)*sin_tI; // s(𝜆) = sR + R cos(𝜆 + 𝜃I)
// }

// /*gpufun*/
// double CircularTrajectory_deriv_x(CircularTrajectory traj, double l){
//     double R = CircularTrajectory_get_R(traj);
//     double sin_tI = CircularTrajectory_get_sin_tI(traj);
//     double cos_tI = CircularTrajectory_get_cos_tI(traj);
//     return R*cos(l)*cos_tI - R*sin(l)*sin_tI; // s(𝜆) = sR + R sin(𝜆 + 𝜃I)
// }

// /*gpufun*/
// void CircularTrajectory_init_bounding_box(CircularTrajectory traj, BoundingBox box, double l1, double l2){
//     double theta = atan2(CircularTrajectory_get_sin_tI(traj), CircularTrajectory_get_cos_tI(traj));
//     double ll1 = theta + l1;
//     double ll2 = ll1 + (l2-l1);
//     if ((ll1 > M_PI) && (ll2 > M_PI)) {
//         ll1 = ll1 -2*M_PI;
//         ll2 = ll2 -2*M_PI;
//     }
//     double s1 = CircularTrajectory_func_s(traj, l1);
//     double x1 = CircularTrajectory_func_x(traj, l1);
//     double s2 = CircularTrajectory_func_s(traj, l2);
//     double x2 = CircularTrajectory_func_x(traj, l2);
//     double R  = CircularTrajectory_get_R(traj); 
//     double dx = x2 - x1;
//     double ds = s2 - s1;
//     double chord_length = sqrt(dx*dx + ds*ds);
//     double sin_chord = dx / chord_length;
//     double cos_chord = ds / chord_length;
//     double sin_t, cos_t;                                       // angle of the box wrt horizontal
//     double min_x, min_s;   
//     double sin_rot, cos_rot;                                   // angle of rotation
//     double l, w;
//     double rotate_box = -1.;
//     int8_t sign = 1;                                           // The orientation of the box can impact calculations
//     if (((cos_chord > 1e-10) && (sin_chord > 1e-10))){           // if 0 < chord angle < 90 deg, then chord angle = box angle
//         sin_t = sin_chord;
//         cos_t = cos_chord;
//         sin_rot = -cos_t;
//         cos_rot = sin_t;
//         sign = -1;
//     } else {
//         if (sin_chord < 1e-10) {                               // if theta is larger than 180 degrees, theta = theta - 180
//             sin_chord = -sin_chord;
//             cos_chord = -cos_chord;
//         }
//         if ( ((cos_chord > 1e-10) && (sin_chord > 1e-10))){       // if 0 < chord angle < 90 deg, then chord angle = box angle
//             sin_t = sin_chord;
//             cos_t = cos_chord;
//             sin_rot = -cos_t;
//             cos_rot = sin_t;
//             sign = -1;
//         } else if (((((cos_chord - 1.) < 1e-10) || ((cos_chord + 1.) > -1e-10)) && (sin_chord < 1e-10)) ||
//                      ((cos_chord < 1e-10) && (((sin_chord-1.) < 1e-10) || ((sin_chord+1) > -1e-10)))){
//             sin_t = sin_chord;
//             cos_t = cos_chord;
//             sin_rot = -cos_t;
//             cos_rot = sin_t;
//             sign = -1;
//         } else {
//             sin_t = - cos_chord;                               // box angle is 90 degrees less than chord angle
//             cos_t = sin_chord;
//             sin_rot = sin_t;    
//             cos_rot = cos_t;
//             sign = -1;
//         } 
//     }
//     if ((ll2 - ll1) < M_PI){                                   // delta theta of box less than 180.
//         w = R - sqrt(R*R - chord_length*chord_length/4.);      // sagitta
//         l = chord_length;
//         if (cos_chord > 0){
//             rotate_box = -1;
//         } else {
//             rotate_box = 1;
//         }
//         // finding the first vertex. 
//         if (x1 < x2){
//             if (((ll1 < M_PI && ll2 > M_PI) && (cos_chord < 0)) || (( (-M_PI < ll1) && (ll1 < 0.) ) && ( (-M_PI < ll2) && (ll2 < 0.))) || ((ll1 < 0 && ll2 > 0) && (cos_chord > 0))){                // l1 or l2 is lower vertex
//                 BoundingBox_set_rC(box,( sqrt( ((s1)+w*cos_rot)*((s1)+w*cos_rot) +
//                                                ((x1)+w*sin_rot)*((x1)+w*sin_rot)) ));
//                 double rC = BoundingBox_get_rC(box);
//                 BoundingBox_set_sin_tC(box,((x1)+w*sin_rot) / rC);
//                 BoundingBox_set_cos_tC(box,((s1)+w*cos_rot) / rC);
//             } else {
//                 double rC = (sqrt((s1)*(s1) + (x1)*(x1)));
//                 BoundingBox_set_rC(box, rC);
//                 BoundingBox_set_sin_tC(box, (x1) /rC);
//                 BoundingBox_set_cos_tC(box, (s1) /rC);
//             }
//         } else {
//             if (((ll1 < M_PI && ll2 > M_PI) && (cos_chord < 0)) || (( (-M_PI < ll1) && (ll1 < 0.) ) && ( (-M_PI < ll2) && (ll2 < 0.))) || ((ll1 < 0 && ll2 > 0) && (cos_chord > 0))){                // l1 or l2 is lower vertex
//                 BoundingBox_set_rC(box, (sqrt( ((s2)-w*cos_rot)*((s2)-w*cos_rot) +
//                                               ((x2)-w*sin_rot)*((x2)-w*sin_rot)) ));
//                 BoundingBox_set_sin_tC(box, (((x2)-w*sin_rot)) / BoundingBox_get_rC(box));
//                 BoundingBox_set_cos_tC(box, ((s2)-w*cos_rot) / BoundingBox_get_rC(box));

//             } else {
//                 BoundingBox_set_rC(box, sqrt((s2)*(s2) + (x2)*(x2)));
//                 BoundingBox_set_sin_tC(box, (x2) / BoundingBox_get_rC(box));
//                 BoundingBox_set_cos_tC(box, (s2) / BoundingBox_get_rC(box));
//             }
//         }
//     } else {
//         rotate_box = -1;
//         l = (2*R);                                            // L is always on the side of the chord in the code
//         w = (R + sqrt(R*R - (chord_length*chord_length/4.))); // 2R - sagitta
//         double chord_side_w1_x, chord_side_w1_s, chord_side_w2_x, chord_side_w2_s;
//         double w3_x, w3_s, w4_x, w4_s;
//         // determining the (s,x) of the vertices on the chord side of box
//         if ((x2-x1) > 1e-10){
//             rotate_box = 1;
//             chord_side_w1_x = x1 - (2*R - chord_length)/2.*sin_chord;
//             chord_side_w1_s = s1 - (2*R - chord_length)/2.*cos_chord;
//             chord_side_w2_x = x2 + (2*R - chord_length)/2.*sin_chord;
//             chord_side_w2_s = s2 + (2*R - chord_length)/2.*cos_chord;
//             w3_x = chord_side_w1_x - sign*w*sin_rot;
//             w3_s = chord_side_w1_s + w*cos_rot;
//             w4_x = chord_side_w2_x - sign*w*sin_rot;
//             w4_s = chord_side_w2_s + w*cos_rot;
//         } else {
//             rotate_box = -1;
//             chord_side_w2_x = x2 - (2*R - chord_length)/2.*sin_chord;
//             chord_side_w2_s = s2 - (2*R - chord_length)/2.*cos_chord;
//             chord_side_w1_x = x1 + (2*R - chord_length)/2.*sin_chord;
//             chord_side_w1_s = s1 + (2*R - chord_length)/2.*cos_chord;
//             w3_x = chord_side_w1_x + sign*w*sin_rot; 
//             w3_s = chord_side_w1_s - w*cos_rot;
//             w4_x = chord_side_w2_x + sign*w*sin_rot;
//             w4_s = chord_side_w2_s - w*cos_rot;
//         }
//         // Compare with the other three points to find the first vertex. Also taking into account if the box is horiztonal
//         min_x = chord_side_w1_x;
//         min_s = chord_side_w1_s;
//         if ((chord_side_w2_x < min_x) || (((chord_side_w2_x - min_x) < 1e-10) && (chord_side_w2_s < min_s))) {
//             if ((chord_side_w2_x < min_x) && (chord_side_w2_s > min_s)) {
//                 rotate_box = 1;
//             } else {
//                 rotate_box = -1;
//             }
//             min_x = chord_side_w2_x; 
//             min_s = chord_side_w2_s;
//         }
//         if (w3_x < min_x || (((w3_x - min_x) < 1e-10) && (w3_s < min_s))) {
//             min_x = w3_x; 
//             min_s = w3_s;   
//             rotate_box = -1;
//         }
//         if (w4_x < min_x || (((w4_x - min_x)<1e-10) && (w4_s <= min_s))) {
//             min_x = w4_x; 
//             min_s = w4_s;
//             rotate_box = 1;
//         }
//         BoundingBox_set_rC(box, sqrt((min_s)*(min_s) + (min_x)*(min_x)));
//         BoundingBox_set_sin_tC(box, min_x / BoundingBox_get_rC(box));
//         BoundingBox_set_cos_tC(box, min_s / BoundingBox_get_rC(box));
//     }
//     if (rotate_box == -1){
//         BoundingBox_set_l(box, l);   // length of the box
//         BoundingBox_set_w(box, w);   // width of the box
//     } else {
//         BoundingBox_set_l(box, w);   // length of the box
//         BoundingBox_set_w(box, l);   // width of the box
//     }
//     double rC = BoundingBox_get_rC(box);                           // length of the position vector to the first vertex
//     BoundingBox_set_sin_tb(box, sin_t);                            // orientation of the box (angle of length wrt horizontal)
//     BoundingBox_set_cos_tb(box, cos_t);
//     BoundingBox_set_proj_l(box, rC * (cos_t*s1/rC + sin_t*x1/rC)); // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
//     BoundingBox_set_proj_w(box, rC * (cos_t*x1/rC - sin_t*s1/rC)); // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
// }


// #endif /* XCOLL_GEOM_TRAJ_CIRCULAR_H */