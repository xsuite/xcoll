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
void LineSegment_update_box(LineSegment seg, double t1, double t2){
    // These ifs will be removed later when we know that the code works and never produces invalid t1, t2
    if (t1 >= t2){
        printf("t1 must be smaller than t2!\n");
        fflush(stdout);
        return;
    }
    if (t1 < 0 || t1 > 1){
        printf("t1 must be in [0, 1]!\n");
        fflush(stdout);
        return;
    }
    if (t2 < 0 || t2 > 1){
        printf("t2 must be in [0, 1]!\n");
        fflush(stdout);
        return;
    }
    // Calculate the bounding box of a line segment.
    // Theta is the angle of the line wrt the horizontal.
    // Phi is the angle from s1 to the first vertex (in the frame of the box).
    BoundingBox box = LineSegment_getp_box(seg);
    double s1 = LineSegment_func_s(seg, t1);
    double s2 = LineSegment_func_s(seg, t2);
    double x1 = LineSegment_func_x(seg, t1);
    double x2 = LineSegment_func_x(seg, t2);
    double sin_t = (x2 - x1) / sqrt((x2 - x1)*(x2 - x1) + (s2 - s1)*(s2 - s1));
    double cos_t = (s2 - s1) / sqrt((x2 - x1)*(x2 - x1) + (s2 - s1)*(s2 - s1));
    double l = sqrt((s2 - s1)*(s2 - s1) + (x2 - x1)*(x2 - x1));   // length of the box
    double w = 0.0;                                                // width of the box cannot be 0 for (0,0)
    double sin_tC, cos_tC;                                        // angle of the position vector to the first vertex
    double rC = sqrt(s1*s1 + x1*x1);                              // length of the position vector to the first vertex
    if (rC == 0.0){
        sin_tC = 0.0;                                             // angle of the position vector to the first vertex
        cos_tC = 1.0;
    } else {
        sin_tC = x1 / rC;                                         // angle of the position vector to the first vertex
        cos_tC = s1 / rC;
    }
    BoundingBox_set_params(box, rC, sin_tC, cos_tC, l, w, sin_t, cos_t);
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

/*gpufun*/
double LineSegment_prepare_newton(LineSegment seg, BoundingBox MCSbox, double tol){
    // Prepare initial guess for Newton-Raphson root finding
    double org_t1 = LineSegment_get__t1(seg);
    double org_t2 = LineSegment_get__t2(seg);
    while ((LineSegment_get__t2(seg) -  LineSegment_get__t1(seg)) > tol){
        double t1_old = LineSegment_get__t1(seg);
        double t2_old = LineSegment_get__t2(seg);
        double t_middle = 0.5 * (LineSegment_get__t2(seg) + LineSegment_get__t1(seg));

        // first half
        LineSegment_set__t2(seg, t_middle);
        double overlap_lower = BoundingBox_overlaps(MCSbox, 
                                                    LineSegment_getp_box(seg));
        // second half
        LineSegment_set__t2(seg, t2_old);
        LineSegment_set__t1(seg, t_middle);
        double overlap_upper = BoundingBox_overlaps(MCSbox, 
                                                    LineSegment_getp_box(seg));
        if (overlap_lower && !overlap_upper){
            LineSegment_set__t1(seg, t1_old);
            LineSegment_set__t2(seg, t_middle);
        }
    }
    double t = 0.5 * (LineSegment_get__t2(seg) + LineSegment_get__t1(seg));
    // Reset to original values
    LineSegment_set__t1(seg, org_t1);
    LineSegment_set__t2(seg, org_t2);
    return t;
}

// /*gpufun*/ 
// void LineSegment_crossing_mcs(LineSegment seg, int8_t* n_hit, double* s, const double* Ax, const double Xo, void* params){
//     return grid_search_and_newton();
// }
#endif /* XCOLL_GEOM_SEG_LINE_H */