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
    //UNUSED(t);
    return LineSegment_get_s2(seg) - LineSegment_get_s1(seg);
}

/*gpufun*/
double LineSegment_deriv_x(LineSegment seg, double t){
    //UNUSED(t);
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
    //BoundingBox box = LineSegment_get_box(seg);
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
    BoundingBox_set_params(LineSegment_getp_box(seg), rC, sin_tC, cos_tC, l, w, sin_t, cos_t);
}
#endif /* XCOLL_GEOM_SEG_LINE_H */