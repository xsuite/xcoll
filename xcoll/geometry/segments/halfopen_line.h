// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SEG_HALFOPENLINE_H
#define XCOLL_GEOM_SEG_HALFOPENLINE_H


/*gpufun*/
int8_t HalfOpenLineSegment_bounded_below(HalfOpenLineSegment seg){
    return 1;  // halfopen segment
}

/*gpufun*/
int8_t HalfOpenLineSegment_bounded_above(HalfOpenLineSegment seg){
    return 0;  // halfopen segment
}

/*gpufun*/
double HalfOpenLineSegment_func_s(HalfOpenLineSegment seg, double t){
    double s1 = HalfOpenLineSegment_get_s1(seg);
    double cos_t1 = HalfOpenLineSegment_get_cos_t1(seg);
    return s1 + t*cos_t1;
}

/*gpufun*/
double HalfOpenLineSegment_func_x(HalfOpenLineSegment seg, double t){
    double x1 = HalfOpenLineSegment_get_x1(seg);
    double sin_t1 = HalfOpenLineSegment_get_sin_t1(seg);
    return x1 + t*sin_t1;
}

/*gpufun*/
double HalfOpenLineSegment_deriv_s(HalfOpenLineSegment seg, double t){
    UNUSED(t);
    return HalfOpenLineSegment_get_cos_t1(seg);
}

/*gpufun*/
double HalfOpenLineSegment_deriv_x(HalfOpenLineSegment seg, double t){
    UNUSED(t);
    return HalfOpenLineSegment_get_sin_t1(seg);
}

/*gpufun*/
void HalfOpenLineSegment_init_bounding_box(HalfOpenLineSegment seg, BoundingBox box, double t1, double t2){
    double s1 = HalfOpenLineSegment_func_s(seg, t1);
    double s2 = HalfOpenLineSegment_func_s(seg, t2);
    double x1 = HalfOpenLineSegment_func_x(seg, t1);
    double x2 = HalfOpenLineSegment_func_x(seg, t2);
    double sin_t = HalfOpenLineSegment_get_sin_t1(seg);            // angle of the line wrt horizontal
    double cos_t = HalfOpenLineSegment_get_cos_t1(seg);
    double l = sqrt((s2 - s1)*(s2 - s1) + (x2 - x1)*(x2 - x1));    // length of the box
    double w = 0.0;                                                // width of the box 
    double rC = sqrt(s1*s1 + x1*x1);                               // length of the position vector to the first vertex
    double rC = BoundingBox_get_rC(box);
    double sin_tb = sin_t;                                               // orientation of the box (angle of length wrt horizontal)
    double cos_tb = cos_t;
    double sin_tC, cos_tC;                                   // angle of the position vector to the first vertex
    if (rC == 0.){
        sin_tC = 0.0;  // angle of the position vector to the first vertex
        cos_tC = 1.0;
    } else {
        sin_tC = x1 / rC;  // angle of the position vector to the first vertex
        cos_tC = s1 / rC;
    }
    BoundingBox_set_params(box, rC, sin_tC, cos_tC, l, w, sin_t, cos_t);
}





// /*gpufun*/
// void HalfOpenLineSegment_crossing_drift(HalfOpenLineSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
//     // Get segment data
//     double s1 = HalfOpenLineSegment_get_s(seg);
//     double x1 = HalfOpenLineSegment_get_x(seg);
//     double s2 = s1 + cos(HalfOpenLineSegment_get_t(seg));
//     double x2 = x1 + sin(HalfOpenLineSegment_get_t(seg));
//     double denom = (x2 - x1) - (s2 - s1)*xm;
//     if (fabs(denom) < XC_EPSILON){
//         // Trajectory is parallel to the segment
//         if (fabs((x0 - x1)/(s0 - s1) - xm) < XC_EPSILON){
//             // Trajectory overlaps with the segment
//             // TODO: This is situational; we should return s1 if get_s_first and current_s if after current_s
//             //       For now we hit twice (because we go nor IN nor OUT)
//             s[*n_hit] = s1;
//             (*n_hit)++;
//             s[*n_hit] = s1;
//             (*n_hit)++;
//         } else {
//             // No hit
//             return;
//         }
//     } else {
//         double t = (x0 - x1 - (s0 - s1)*xm) / denom;
//         if (t >= 0){  // We do not check for t<=1 as it is a half-open segment
//             s[*n_hit] = s1*(1-t) + s2*t;
//             (*n_hit)++;
//         }
//     }
// }

// /*gpufun*/
// void HalfOpenLineSegment_crossing_mcs(HalfOpenLineSegment seg, int8_t* n_hit, double* s, const double* Ax, const double Xo){
//    return grid_search_and_newton();
// }  
#endif /* XCOLL_GEOM_SEG_HALFOPENLINE_H */