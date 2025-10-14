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
void LineSegment_update_box(LineSegment seg, BoundingBox box, double t1, double t2){
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
    double s1 = LineSegment_func_s(seg, t1);
    double s2 = LineSegment_func_s(seg, t2);
    double x1 = LineSegment_func_x(seg, t1);
    double x2 = LineSegment_func_x(seg, t2);
    box->sin_t = (x2 - x1) / sqrt((x2 - x1)*(x2 - x1) + (s2 - s1)*(s2 - s1));
    box->cos_t = (s2 - s1) / sqrt((x2 - x1)*(x2 - x1) + (s2 - s1)*(s2 - s1));
    box->l = sqrt((s2 - s1)*(s2 - s1) + (x2 - x1)*(x2 - x1));
    box->w = 0.0; // line segment has no width
    box->rC = sqrt(s1*s1 + x1*x1);
    if (box->rC == 0.0){
        box->sin_tC = 0.0;
        box->cos_tC = 1.0;
    } else {
        box->sin_tC = x1 / box->rC;
        box->cos_tC = s1 / box->rC;
    }
    BoundingBox_sync(box);
}


// Expose functions to Xobject test interface
// ------------------------------------------
void LineSegment_update_testbox(LineSegment seg, BoundingBoxTest box, double t1, double t2){
    BoundingBox_s box1;
    LineSegment_update_box(seg, &box1, t1, t2);
    BoundingBox_to_BoundingBoxTest(&box1, box);
}

#endif /* XCOLL_GEOM_SEG_LINE_H */