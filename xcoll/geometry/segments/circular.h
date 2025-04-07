// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #
#include <math.h>
#ifndef XCOLL_GEOM__SEG_CIRCULAR_H
#define XCOLL_GEOM__SEG_CIRCULAR_H


/*gpufun*/
int8_t CircularSegment_bounded_below(CircularSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
int8_t CircularSegment_bounded_above(CircularSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
double CircularSegment_func_s(CircularSegment seg, double t){
    // t in [0,1] becomes tt in [theta1, theta2]
    double R = CircularSegment_get_R(seg);
    double sR = CircularSegment_get_sR(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    return sR + R*cos(theta1 + t*(theta2 - theta1));
}

/*gpufun*/
double CircularSegment_func_x(CircularSegment seg, double t){
    double R = CircularSegment_get_R(seg);
    double xR = CircularSegment_get_xR(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    return xR + R*sin(theta1 + t*(theta2 - theta1));
}

/*gpufun*/
double CircularSegment_deriv_s(CircularSegment seg, double t){
    double R = CircularSegment_get_R(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    return -R*sin(theta1 + t*(theta2 - theta1))*(theta2 - theta1);
}

/*gpufun*/
double CircularSegment_deriv_x(CircularSegment seg, double t){
    double R = CircularSegment_get_R(seg);
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    return R*cos(theta1 + t*(theta2 - theta1))*(theta2 - theta1);
}

/*gpufun*/
void CircularSegment_init_bounding_box(CircularSegment seg, BoundingBox box, double t1, double t2){
    double theta1 = CircularSegment_get__theta1(seg);
    double theta2 = CircularSegment_get__theta2(seg);
    double tt1    = theta1 + t1*(theta2 - theta1);
    double tt2    = theta1 + t2*(theta2 - theta1);
    if ((tt1 > M_PI) && (tt2 > M_PI)) {
        tt1 -= 2*M_PI;
        tt2 -= 2*M_PI;
    }
    double s1 = CircularSegment_func_s(seg, t1);
    double x1 = CircularSegment_func_x(seg, t1);
    double s2 = CircularSegment_func_s(seg, t2);
    double x2 = CircularSegment_func_x(seg, t2);
    double R  = CircularSegment_get_R(seg); 
    double dx = x2 - x1;
    double ds = s2 - s1;
    double chord_length = sqrt(dx*dx + ds*ds);
    double sin_chord = dx / chord_length;
    double cos_chord = ds / chord_length;
    double sin_t, cos_t;                                       // angle of the box wrt horizontal
    double min_x, min_s;   
    double sin_rot, cos_rot;
    double l, w;
    double rotate_box = -1.;
    int8_t sign = 1;                                           // The orientation of the box can impact calculations
    if ((cos_chord > 1e-10) && (sin_chord > 1e-10)){           // if 0 < chord angle < 90 deg, then chord angle = box angle
        sin_t = sin_chord;
        cos_t = cos_chord;
        sin_rot = -cos_t;
        cos_rot = sin_t;
        sign = -1;
    } else {
        if (sin_chord < 1e-10) {                               // if theta is larger than 180 degrees, theta = theta - 180
            sin_chord = -sin_chord;
            cos_chord = -cos_chord;
        }
        if ((cos_chord > 1e-10) && (sin_chord > 1e-10)){       // if 0 < chord angle < 90 deg, then chord angle = box angle
            sin_t = sin_chord;
            cos_t = cos_chord;
            sin_rot = -cos_t;
            cos_rot = sin_t;
            sign = -1;
        } else {
            sin_t = - cos_chord;                               // box angle is 90 degrees less than chord angle
            cos_t = sin_chord;
            sin_rot = sin_t;    
            cos_rot = cos_t;
            sign = -1;
        } 
    }
    if ((tt2 - tt1) < M_PI){                                   // delta theta of box less than 180.
        w = R - sqrt(R*R - chord_length*chord_length/4.);      // sagitta
        l = chord_length;
        if (cos_chord > 0){
            rotate_box = -1;
        } else {
            rotate_box = 1;
        }
        // finding the first vertex. 
        if (x1 < x2){
            if (((tt1 < M_PI && tt2 > M_PI) && (cos_chord < 0)) || (( (-M_PI < tt1) && (tt1 < 0.) ) && ( (-M_PI < tt2) && (tt2 < 0.))) || ((tt1 < 0 && tt2 > 0) && (cos_chord > 0))){                // t1 or t2 is lower vertex
                BoundingBox_set_rC(box,( sqrt( ((s1)+w*cos_rot)*((s1)+w*cos_rot) +
                                               ((x1)+w*sin_rot)*((x1)+w*sin_rot)) ));
                double rC = BoundingBox_get_rC(box);
                BoundingBox_set_sin_tC(box,((x1)+w*sin_rot) / rC);
                BoundingBox_set_cos_tC(box,((s1)+w*cos_rot) / rC);
            } else {
                double rC = (sqrt((s1)*(s1) + (x1)*(x1)));
                BoundingBox_set_rC(box, rC);
                BoundingBox_set_sin_tC(box, (x1) /rC);
                BoundingBox_set_cos_tC(box, (s1) /rC);
            }
        } else {
            if (((tt1 < M_PI && tt2 > M_PI) && (cos_chord < 0)) || (( (-M_PI < tt1) && (tt1 < 0.) ) && ( (-M_PI < tt2) && (tt2 < 0.))) || ((tt1 < 0 && tt2 > 0) && (cos_chord > 0))){                // t1 or t2 is lower vertex
                BoundingBox_set_rC(box, (sqrt( ((s2)-w*cos_rot)*((s2)-w*cos_rot) +
                                              ((x2)-w*sin_rot)*((x2)-w*sin_rot)) ));
                BoundingBox_set_sin_tC(box, (((x2)-w*sin_rot)) / BoundingBox_get_rC(box));
                BoundingBox_set_cos_tC(box, ((s2)-w*cos_rot) / BoundingBox_get_rC(box));

            } else {
                BoundingBox_set_rC(box, sqrt((s2)*(s2) + (x2)*(x2)));
                BoundingBox_set_sin_tC(box, (x2) / BoundingBox_get_rC(box));
                BoundingBox_set_cos_tC(box, (s2) / BoundingBox_get_rC(box));
            }
        }
    } else {
        rotate_box = -1;
        l = (2*R);                                            // L is always on the side of the chord in the code
        w = (R + sqrt(R*R - (chord_length*chord_length/4.))); // 2R - sagitta
        double chord_side_w1_x, chord_side_w1_s, chord_side_w2_x, chord_side_w2_s;
        double w3_x, w3_s, w4_x, w4_s;
        // determining the (s,x) of the vertices on the chord side of box
        if (x1 < x2){
            rotate_box = 1;
            chord_side_w1_x = x1 - (2*R - chord_length)/2.*sin_chord;
            chord_side_w1_s = s1 - (2*R - chord_length)/2.*cos_chord;
            chord_side_w2_x = x2 + (2*R - chord_length)/2.*sin_chord;
            chord_side_w2_s = s2 + (2*R - chord_length)/2.*cos_chord;
            w3_x = chord_side_w1_x - sign*w*sin_rot;
            w3_s = chord_side_w1_s + w*cos_rot;
            w4_x = chord_side_w2_x - sign*w*sin_rot;
            w4_s = chord_side_w2_s + w*cos_rot;
        } else {
            rotate_box = -1;
            chord_side_w2_x = x2 - (2*R - chord_length)/2.*sin_chord;
            chord_side_w2_s = s2 - (2*R - chord_length)/2.*cos_chord;
            chord_side_w1_x = x1 + (2*R - chord_length)/2.*sin_chord;
            chord_side_w1_s = s1 + (2*R - chord_length)/2.*cos_chord;
            w3_x = chord_side_w1_x + sign*w*sin_rot; 
            w3_s = chord_side_w1_s - w*cos_rot;
            w4_x = chord_side_w2_x + sign*w*sin_rot;
            w4_s = chord_side_w2_s - w*cos_rot;
        }
        // Compare with the other three points to find the first vertex. Also taking into account if the box is horiztonal
        min_x = chord_side_w1_x;
        min_s = chord_side_w1_s;
        if ((chord_side_w2_x < min_x) || (chord_side_w2_x == min_x && chord_side_w2_s < min_s)) {
            if ((chord_side_w2_x < min_x) && (chord_side_w2_s > min_s)) {
                rotate_box = 1;
            } else {
                rotate_box = -1;
            }
            min_x = chord_side_w2_x; 
            min_s = chord_side_w2_s;
        }
        if (w3_x < min_x || (w3_x == min_x && w3_s < min_s)) {
            min_x = w3_x; 
            min_s = w3_s;   
            rotate_box = -1;
        }
        if (w4_x < min_x || (((w4_x - min_x)<1e-10) && w4_s <= min_s)) {
            min_x = w4_x; 
            min_s = w4_s;
            rotate_box = 1;
        }
        BoundingBox_set_rC(box, sqrt((min_s)*(min_s) + (min_x)*(min_x)));
        BoundingBox_set_sin_tC(box, min_x / BoundingBox_get_rC(box));
        BoundingBox_set_cos_tC(box, min_s / BoundingBox_get_rC(box));
    }
    if (rotate_box == -1){
        BoundingBox_set_l(box, l);   // length of the box
        BoundingBox_set_w(box, w);   // width of the box
    } else {
        BoundingBox_set_l(box, w);   // length of the box
        BoundingBox_set_w(box, l);   // width of the box
    }
    double rC = BoundingBox_get_rC(box);                           // length of the position vector to the first vertex
    BoundingBox_set_sin_tb(box, sin_t);                            // orientation of the box (angle of length wrt horizontal)
    BoundingBox_set_cos_tb(box, cos_t);
    BoundingBox_set_proj_l(box, rC * (cos_t*s1/rC + sin_t*x1/rC)); // projection of position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    BoundingBox_set_proj_w(box, rC * (cos_t*x1/rC - sin_t*s1/rC)); // projection of position vector on width:  rC * (cos_t*sin_tC - sin_t*cos_tC)
}
// /*gpufun*/
// void CircularSegment_crossing_drift(CircularSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
//     // Get segment data
//     double R  = CircularSegment_get_R(seg);
//     double sC = CircularSegment_get_s(seg);
//     double xC = CircularSegment_get_x(seg);
//     double t1 = CircularSegment_get_t1(seg);
//     double t2 = CircularSegment_get_t2(seg);
//     // Move the angles to [-pi, pi]
//     int8_t reversed = 0, full_circle = 0;
//     if (fabs(fabs(t2 - t1) - 2*M_PI) < XC_EPSILON){
//         full_circle = 1;
//     }
//     while (t1 < -M_PI){
//         t1 += 2*M_PI;
//     }
//     while (t1 > M_PI){
//         t1 -= 2*M_PI;
//     }
//     while (t2 < -M_PI){
//         t2 += 2*M_PI;
//     }
//     while (t2 > M_PI){
//         t2 -= 2*M_PI;
//     }
//     if (t2 < t1){
//         reversed = 1;
//     }
//     // Calculate crossings
//     double a = 1 + xm*xm;
//     double bb = sC - xm*(x0 - xC - xm*s0); // This is -b/2 with b from the quadratic formula
//     double c = sC*sC + (x0 - xC - xm*s0)*(x0 - xC - xm*s0) - R*R;
//     double disc = bb*bb - a*c; // This is  2*discriminant**2
//     if (disc < 0){
//         // No crossing
//         return;
//     }
//     for (int8_t i = 0; i < 2; i++) {
//         double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
//         double new_s = (bb + sgnD*sqrt(fabs(disc)))/a;
//         double new_x = x0 + (new_s - s0)*xm;
//         double t = atan2(new_x - xC, new_s - sC);
//         if (full_circle){
//             // Full circle, so always hit
//             s[*n_hit] = new_s;
//             (*n_hit)++;
//         } else if (reversed){
//             // t2 < t1, so we are looking at the inverted region of angles
//             if (t1 <= t || t <= t2){
//                 s[*n_hit] = new_s;
//                 (*n_hit)++;
//             }
//         } else {
//             if (t1 <= t && t <= t2){
//                 s[*n_hit] = new_s;
//                 (*n_hit)++;
//             }
//         }
//     }
// }
// /*gpufun*/
// void CircularSegment_crossing_mcs(CircularSegment seg, int8_t* n_hit, double* s, double x, const double* Ax, const double Xo){
//     return grid_search_and_newton();
// }
#endif /* XCOLL_GEOM__SEG_CIRCULAR_H */