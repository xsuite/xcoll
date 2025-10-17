// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #
// #include <math.h>
// #ifndef XCOLL_GEOM__SEG_CIRCULAR_H
// #define XCOLL_GEOM__SEG_CIRCULAR_H


// /*gpufun*/
// int8_t CircularSegment_bounded_below(CircularSegment seg){
//     return 1;  // closed segment
// }

// /*gpufun*/
// int8_t CircularSegment_bounded_above(CircularSegment seg){
//     return 1;  // closed segment
// }

// /*gpufun*/
// double CircularSegment_func_s(CircularSegment seg, double t){
//     // t in [0,1] becomes tt in [theta1, theta2]
//     double R = CircularSegment_get_R(seg);
//     double sR = CircularSegment_get_sR(seg);
//     double theta1 = CircularSegment_get__theta1(seg);
//     double theta2 = CircularSegment_get__theta2(seg);
//     return sR + R*cos(theta1 + t*(theta2 - theta1));
// }

// /*gpufun*/
// double CircularSegment_func_x(CircularSegment seg, double t){
//     double R = CircularSegment_get_R(seg);
//     double xR = CircularSegment_get_xR(seg);
//     double theta1 = CircularSegment_get__theta1(seg);
//     double theta2 = CircularSegment_get__theta2(seg);
//     return xR + R*sin(theta1 + t*(theta2 - theta1));
// }

// /*gpufun*/
// double CircularSegment_deriv_s(CircularSegment seg, double t){
//     double R = CircularSegment_get_R(seg);
//     double theta1 = CircularSegment_get__theta1(seg);
//     double theta2 = CircularSegment_get__theta2(seg);
//     return -R*sin(theta1 + t*(theta2 - theta1))*(theta2 - theta1);
// }

// /*gpufun*/
// double CircularSegment_deriv_x(CircularSegment seg, double t){
//     double R = CircularSegment_get_R(seg);
//     double theta1 = CircularSegment_get__theta1(seg);
//     double theta2 = CircularSegment_get__theta2(seg);
//     return R*cos(theta1 + t*(theta2 - theta1))*(theta2 - theta1);
// }

// /*gpufun*/
// void CircularSegment_init_bounding_box(CircularSegment seg, BoundingBox box, double t1, double t2){
//     // The low-level angles always satisfy theta1 < theta2 and -pi <= theta1 <= pi (theta2 can be <= 3pi)
//     // interpolate between theta1 and theta2 to get ang. positions tt1 & tt2 corresponding to arc par. t1 and t1
//     double theta1 = CircularSegment_get__theta1(seg);
//     double theta2 = CircularSegment_get__theta2(seg);
//     double tt1    = theta1 + t1*(theta2 - theta1);
//     double tt2    = theta1 + t2*(theta2 - theta1);
//     // Getting (s,x) coord. for t1 and t2. Has to be t1 and t2 because of parametrization.
//     double s1 = CircularSegment_func_s(seg, t1);
//     double x1 = CircularSegment_func_x(seg, t1);
//     double s2 = CircularSegment_func_s(seg, t2);
//     double x2 = CircularSegment_func_x(seg, t2);
//     // Calculate the chord vector between two arc points and its angle wrt the horizontal.
//     double R  = CircularSegment_get_R(seg); 
//     double dx = x2 - x1;
//     double ds = s2 - s1;
//     double chord_length = sqrt(dx*dx + ds*ds);
//     if (chord_length < 1e-10){
//         // The two points are the same, so we have a 360 segment
//         double sR = CircularSegment_get_sR(seg);
//         double xR = CircularSegment_get_xR(seg);
//         double rC = sqrt((sR-R)*(sR-R) + (xR-R)*(xR-R)); // length of position vector to first vertex
//         BoundingBox_set_params(box, rC, (xR-R)/rC, (sR-R)/rC, 2*R, 2*R, 0., 1.);
//         return;
//     }
//     double sin_chord = dx / chord_length;
//     double cos_chord = ds / chord_length;
//     // Init. _t is the box angle, min and max is used to find lowest vertex,
//     double sin_t, cos_t;                                         // angle of the box wrt horizontal
//     double min_x, min_s;
//     double l, w;
//     double rC, sin_tC, cos_tC;
//     // rotate box is needed to make sure we sent the correct w and l to the box. When the box is tilted we tend to change what
//     // is the lowest vertex. In the box L is always the length from the first vertex to the second vertex.
//     // in this code L is always the length of the side with the chord. In which the chord is the length between the two points. 
//     // This is done to maintain some generality in the code. Without it we would have to change it IN the code, which would be messy. 
    
//     // Sign is needed due to the orientation of the box again. Bascially, when finding the lowest vertex etc we need to know if we are
//     // adding or subtracting stuff. This comes from sign. It is needed when having (tt2-tt1)>180 deg, because then we need to find all 4
//     // vertices and compare them to find the lowest one (slow, but i dont see how else we can do it since we dont know anything except (s1,x1)
//     // and (s2,x2), which in this case is not any of the vertices). 
//     double rotate_box = -1.;
//     int8_t sign = 1;                                             // The orientation of the box can impact calculations
//     // this next part is made for the angles. We want the angle to be [0, 180] for simplicity. We need to know if box angle = chord ang., 
//     // and we check what quadrant we're in (in which 1 and 3 gives equal, 2 and 4 does not). This is all done to find box angle
//     // and rot angle which we can only get from the chord angle. 
//     // We normalize to [0, pi] because the box angle is always from the lowest vertex wrt the horizontal. So, if the box angle is larger 
//     // than pi, then we also change the lowest vertex, and the box angle should adjust accordingly. So, it doesnt make sense to have it larger
//     // than pi. Therefore, we adjust after (s1,x1) and (s2,x2). Points --> chord angle --> box angle & rot angle.
    
//     // this if else check is to check if the chord angle is larger than 180 or not.
//     if (((cos_chord > 1e-10) && (sin_chord > 1e-10))){           // if 0 < chord angle < 90 deg, then chord angle = box angle
//         sin_t = sin_chord;
//         cos_t = cos_chord;
//         sign = -1;
//     } else {
//         // if larger than 180 degrees, then we need to adjust the angle.
//         if (sin_chord < 1e-10) {                               // if theta is larger than 180 degrees, theta = theta - 180
//             sin_chord = -sin_chord;
//             cos_chord = -cos_chord;
//         }
//         // then we need to recheck if 1) are we in the first quadrant, or
//         // 2) or are we at a border; 90, 180, etc..
//         // If neither of those then we know we're in the second quadrant, and we know that box != chord. 
//         if ( ((cos_chord > 1e-10) && (sin_chord > 1e-10))){       // if 0 < chord angle < 90 deg, then chord angle = box angle
//             sin_t = sin_chord;
//             cos_t = cos_chord;
//             sign = -1;
//         } else if ((((cos_chord - 1.) < 1e-10) || ((cos_chord + 1.) > -1e-10)) && ((sin_chord < 1e-10) && (sin_chord > -1e-10))){
//             sin_t = sin_chord;  // if cos_chord is 1 or -1, then box angle = chord angle
//             cos_t = cos_chord;
//             sign = -1;
//         } else {
//             sin_t = -cos_chord;                               // box angle is 90 degrees less than chord angle
//             cos_t = sin_chord;
//             sign = -1;
//         } 
//     }
//     // now that the angles are set, we need to figure out what the box looks like. First check is to check if 
//     // the box angle is less than 180 degrees. If its larger than 180 degrees, then we need to find the new width and length
//     if ((tt2 - tt1) <= M_PI){                                   // delta theta of box less than 180.
//         // this is the case where we have a box that is less than 180 degrees. Then the widt is the sagitta of the circle, and
//         // the length is the chord length.
//         w = R - sqrt(R*R - chord_length*chord_length/4.);       // sagitta
//         l = chord_length;
//         // this part is again to make sure that we switch w and l if needed. If cos_chord < 0, then the line between the lowest point.
//         // and the next point (counterclockwise) is actually what we define as W here. It is the shorter side of the box. But, instead of 
//         // mixing w and l in the code itself, we will just switch them at the end. Rotate box == 1 gives switch. 
//         // Try drawing it out. 
//         if ((cos_chord < 0. && cos_chord > -1.)){
//             rotate_box = 1;
//         } else {
//             rotate_box = -1;
//         }
//         // finding the first vertex. 
//         // this is the part where we find the first vertex. We need to check if the first vertex is the lower one or not. This
//         // depends on the angle and position of the box. We separate into cases of x1 < x2 and x1 > x2. This is done because
//         // the rC value depends on that; x1 lower --> x1 is first vertex (or closest to it), and opposite.
//         if (x1 < x2){
//             // this is where we check if (s1,x2) is the lowest vertex. Try drawing it out. Here you see that we use the rotation angle to
//             // find the lowest vertex from point (s1,x1). Remember the rotation angle depends on the box angle, and if box = chord angle, then
//             // we need to rotate 90 degrees to find the lowest vertex. If not, then it is equal to the box angle. It is generally always 
//             // 90 degress off the chord angle. That is why we are using the chord angle for this. 
//             if (((tt1 <= M_PI && tt2 > M_PI) && (cos_chord < 0)) || (( (-M_PI <= tt1) && (tt1 < 0.) ) 
//                 && ( (-M_PI <= tt2) && (tt2 < 0.))) || ((tt1 <= 0 && tt2 > 0) && (cos_chord > 0))){                // t1 or t2 is lower vertex
//                 rC = sqrt( ((s1)+w*sin_chord)*((s1)+w*sin_chord) +
//                            ((x1)+w*(-cos_chord))*((x1)+w*(-cos_chord)));
//                 sin_tC = ((x1)+w*(-cos_chord)) / rC;
//                 cos_tC = ((s1)+w*sin_chord) / rC;
//             } else {
//                 // in this case (s1,x1) is the lowest vertex. Easier case.
//                 rC = (sqrt((s1)*(s1) + (x1)*(x1)));
//                 sin_tC = (x1) /rC;
//                 cos_tC = (s1) /rC;
//             }
//         } else {
//             // this is the same as above, but for the other point.
//             if (((tt1 <= M_PI && tt2 >= M_PI) && (cos_chord < 0)) || (( (-M_PI <= tt1) && (tt1 <= 0.) ) 
//                 && ( (-M_PI <= tt2) && (tt2 <= 0.))) || ((tt1 <= 0 && tt2 >= 0) && (cos_chord > 0))){                // t1 or t2 is lower vertex
//                 rC = (sqrt(((s2) - w*sin_chord)*((s2) - w*sin_chord) +
//                             ((x2) - w*(-cos_chord))*((x2) - w*(-cos_chord))));
//                 sin_tC = (((x2) - w*(-cos_chord))) / rC;
//                 cos_tC = ((s2) - w*sin_chord) / rC;

//             } else {
//                 rC = sqrt((s2)*(s2) + (x2)*(x2));
//                 sin_tC = (x2) / rC;
//                 cos_tC = (s2) / rC;
//             }
//         }
//     } else {
//         // This is the case of when (tt2-tt1) > 180 degrees.
//         // the length of this box is always the diameter of the circle. This is because now the segment will always span > 180 deg, and
//         // therefore the maximum distance between two points is the diameter.
//         // The width of this circle is equal to the diameter - sagitta.
//         // The sagitta is the distance from the chord to the circle. 
//         l = (2*R);                                            // L is always on the side of the chord in the code
//         w = (R + sqrt(R*R - (chord_length*chord_length/4.))); // 2R - sagitta
//         // these variables are used to find the lowest vertex. In the other case (s1,x2) and (s2,x1) are two of the 4 vertices, 
//         // and we can "easier" check if either of the are the lowest of not. In this case it is a bit worse, seeing as we only have
//         // w, l and two points on the same side of the box (the chord side).
//         // The first variables are for the two points on the chord side of the box. The other two are for the other two points.
//         double chord_side_w1_x, chord_side_w1_s, chord_side_w2_x, chord_side_w2_s;
//         double w3_x, w3_s, w4_x, w4_s;
//         // determining the (s,x) of the vertices on the chord side of box
//         if ((x2-x1) > 1e-10){
//             // first we check if x2 > x1. Then, like before, we can assume that the box needs to be rotated. This time it is always the case, 
//             // and we dont need to check the angles. Might remove this tho seeing as we change rot_box later when we check the vertices.
//             rotate_box = 1;
//             // Here we find vertices on the chord side of the box. You take your point of the chord side (x1,s1) or (x2,s2), and
//             // you go the remaining distance of L to find the actual vertices. (Diameter - chord length)/2 is the distance 
//             // from (s1,x1) or (s2,x2) to the vertices
//             chord_side_w1_x = x1 - (2*R - chord_length)/2.*sin_chord;
//             chord_side_w1_s = s1 - (2*R - chord_length)/2.*cos_chord;
//             chord_side_w2_x = x2 + (2*R - chord_length)/2.*sin_chord;
//             chord_side_w2_s = s2 + (2*R - chord_length)/2.*cos_chord;
//             // here we find the other two vertices. This is where the orientation of the box matters, and we use the sign variable to
//             // determine if we need to add or subtract the distance from the chord side to the vertices. This comes from how the angles were
//             // defined in the beginning! Essentially, we need to rotate the box by 90 degrees (wrt chord angle) to find the other two vertices.
//             w3_x = chord_side_w1_x - sign*w*(-cos_chord);
//             w3_s = chord_side_w1_s + w*sin_chord;
//             w4_x = chord_side_w2_x - sign*w*(-cos_chord);
//             w4_s = chord_side_w2_s + w*sin_chord;
//         } else {
//             // this is the same as above, but for the other case.
//             rotate_box = -1;
//             chord_side_w2_x = x2 - (2*R - chord_length)/2.*sin_chord;
//             chord_side_w2_s = s2 - (2*R - chord_length)/2.*cos_chord;
//             chord_side_w1_x = x1 + (2*R - chord_length)/2.*sin_chord;
//             chord_side_w1_s = s1 + (2*R - chord_length)/2.*cos_chord;
//             w3_x = chord_side_w1_x + sign*w*(-cos_chord); 
//             w3_s = chord_side_w1_s - w*sin_chord;
//             w4_x = chord_side_w2_x + sign*w*(-cos_chord);
//             w4_s = chord_side_w2_s - w*sin_chord;
//         }
//         // Compare with the other three points to find the first vertex. Also taking into account if the box is horiztonal
//         // so, we start with one of the chord side vertices. We need to check if it is the lowest vertex or not.
//         // We do this by checking if the x values are lower than the other vertices.
//         min_x = chord_side_w1_x;
//         min_s = chord_side_w1_s;
//         // here several things are checked. 1) Is the other chord point lower? 
//         // 2) are they equal? If so the box is horizontal! And then 3) is p2 the leftmost point? 
//         if ((chord_side_w2_x < min_x) || (((chord_side_w2_x-min_x) < 1e-10) && (chord_side_w2_s < min_s))) {
//             if ((chord_side_w2_x < min_x) && (chord_side_w2_s > min_s)) {
//                 // if p2 is smaller than p1 wrt x, AND p2 has a larger s value, then we need to rotate the box. 
//                 // This is becaues the chord side will be on the left side of the lower point. And therefore, the line from
//                 // the lower point to the next point counterclockwise is actually the width of the box, but later on it should be the length
//                 // so we need to switch. 
//                 rotate_box = 1;
//             } else {
//                 rotate_box = -1;
//             }
//             min_x = chord_side_w2_x; 
//             min_s = chord_side_w2_s;
//         }
//         // we check if point three is lower (or horizontal and leftmost) than the current lowest point.
//         if (w3_x < min_x || (((w3_x-min_x)<1e-10) && (w3_s < min_s))) {
//             // if yes than point 3 is now the lowest point. Note that here we dont have to check that extra if regarding the 
//             // rotation, because w3 comes from point 1 (on the chord side), and therefore we know there will always be two sides of the box 
//             // between the chord side and w3. So, there will be a width side, w3 point, length, width, chord. All good.
//             min_x = w3_x; 
//             min_s = w3_s;   
//             rotate_box = -1;
//         }
//         if (w4_x < min_x || (((w4_x - min_x)<1e-10) && w4_s <= min_s)) {
//             // here we check the same for point 4, which comes from point 2 on the chord side. Here it is the opposite as w3, and there
//             // will always be a need to switch as there is only one side from w4 until chord side (counterclockwise), and since chord
//             // side is always the length side in this code, then this side in between is width. But if w4 is lowest, then this side has
//             // to be length. We need to switch. 
//             min_x = w4_x; 
//             min_s = w4_s;
//             rotate_box = 1;
//         }
//         rC = sqrt((min_s)*(min_s) + (min_x)*(min_x));
//         sin_tC = min_x / rC;
//         cos_tC = min_s / rC;
//     }
//     // this is the actual switching. 
//     if (rotate_box == 1){
//         BoundingBox_set_params(box, rC, sin_tC, cos_tC, w, l, sin_t, cos_t);
//     } else {
//         BoundingBox_set_params(box, rC, sin_tC, cos_tC, l, w, sin_t, cos_t);
//     }
// }
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