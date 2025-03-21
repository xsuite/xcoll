// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SEG_LINE_H
#define XCOLL_GEOM_SEG_LINE_H


/*gpufun*/
int8_t LineSegment_bounded_below(LineSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
int8_t LineSegment_bounded_above(LineSegment seg){
    return 1;  // closed segment
}

/*gpufun*/
double LineSegment_func_s(LineSegment seg, double t){
    double s1 = LineSegment_get_s1(seg);
    double s2 = LineSegment_get_s2(seg);
    return (1-t)*s1 + t*s2;
}

/*gpufun*/
double LineSegment_func_x(LineSegment seg, double t){
    double x1 = LineSegment_get_x1(seg);
    double x2 = LineSegment_get_x2(seg);
    return (1-t)*x1 + t*x2;
}

/*gpufun*/
double LineSegment_deriv_s(LineSegment seg, double t){
    UNUSED(t);
    return LineSegment_get_s2(seg) - LineSegment_get_s1(seg);
}

/*gpufun*/
double LineSegment_deriv_x(LineSegment seg, double t){
    UNUSED(t);
    return LineSegment_get_x2(seg) - LineSegment_get_x1(seg);
}

/*gpufun*/
void LineSegment_bounding_box(LineSegment seg, double t1, double t2, BoundingBox* box){
    //Calculate the bounding box of a line segment.
    //Theta is the angle of the line wrt the horizontal. 
    //Phi is the angle from s1 to the first vertex (in the frame of the box).
    double s1 = LineSegment_func_s(seg, t1);
    double s2 = LineSegment_func_s(seg, t2);
    double x1 = LineSegment_func_x(seg, t1);
    double x2 = LineSegment_func_x(seg, t2);
    double sin_t = (x2 - x1) / sqrt((x2 - x1)*(x2 - x1) + (s2 - s1)*(s2 - s1));
    double cos_t = (s2 - s1) / sqrt((x2 - x1)*(x2 - x1) + (s2 - s1)*(s2 - s1));
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
    box->l  = sqrt((s2 - s1)*(s2 - s1) + (x2 - x1)*(x2 - x1));   // length of the box
    box->w  = box->l/3.;                                         // width of the box 
    box->rC = sqrt( (s1+box->w/2.*cos_p)*(s1+box->w/2.*cos_p) + // length of the position vector to the first vertex
                    (x1+box->w/2.*sin_p)*(x1+box->w/2.*sin_p) );
    box->sin_tb = sin_t;  // orientation of the box (angle of length wrt horizontal)
    box->cos_tb = cos_t;
    box->sin_tC = x1 / box->rC;  // angle of the position vector to the first vertex
    box->cos_tC = s1 / box->rC;
    box->proj_l = box->rC * (box->cos_tb*box->cos_tC + box->sin_t*box->sin_tC); // projection of the position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    box->proj_w = box->rC * (box->cos_tb*box->sin_tC - box->sin_t*box->cos_tC); // projection of position vector on width: rC * (cos_t*sin_tC - sin_t*cos_tC)
}

// /*gpufun*/
// void LineSegment_crossing_drift(LineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
//     // Get segment data
//     double s1 = LineSegment_get_s1(seg);
//     double x1 = LineSegment_get_x1(seg);
//     double s2 = LineSegment_get_s2(seg);
//     double x2 = LineSegment_get_x2(seg);
//     double denom = (x2 - x1) - (s2 - s1)*xm;
//     if (fabs(denom) < XC_EPSILON){
//         // Trajectory is parallel to the segment
//         if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_EPSILON){
//             // Trajectory overlaps with the segment
//             // TODO: This is situational; we should return s1 if get_s_first and current_s if after current_s
//             //       For now we hit twice (because we go nor IN nor OUT)
//             s[*n_hit] = s1;
//             (*n_hit)++;
//             s[*n_hit] = s2;
//             (*n_hit)++;
//         } else {
//             // No crossing
//             return;
//         }
//     } else {
//         double t = (x0 - x1 - (s0 - s1)*xm) / denom;
//         if (t >= 0 && t <= 1){
//             s[*n_hit] = s1*(1-t) + s2*t;
//             (*n_hit)++;
//         }
//     }
// }


// /*gpufun*/ 
// void LineSegment_crossing_mcs(LineSegment seg, int8_t* n_hit, double* s, const double* Ax, const double Xo, void* params){
//     return grid_search_and_newton();
// }
#endif /* XCOLL_GEOM_SEG_LINE_H */