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
    //UNUSED(t);
    return HalfOpenLineSegment_get_cos_t1(seg);
}

/*gpufun*/
double HalfOpenLineSegment_deriv_x(HalfOpenLineSegment seg, double t){
    //UNUSED(t);
    return HalfOpenLineSegment_get_sin_t1(seg);
}

/*gpufun*/
void HalfOpenLineSegment_update_box(HalfOpenLineSegment seg, BoundingBox box, double t1, double t2){
    // These ifs will be removed later when we know that the code works and never produces invalid t1, t2
    if (t1 >= t2){
        printf("t1 must be smaller than t2!\n");
        fflush(stdout);
        return;
    }
    if (t1 < 0 || t1 > 10){
        printf("t1 must be in [0, 10]!\n");
        fflush(stdout);
        return;
    }
    if (t2 < 0 || t2 > 10){
        printf("t2 must be in [0, 10]!\n");
        fflush(stdout);
        return;
    }
    int8_t _SCALE_FACTOR = 10.;
    t2 = t2 /_SCALE_FACTOR;
    double s1 = HalfOpenLineSegment_func_s(seg, t1);
    double s2 = HalfOpenLineSegment_func_s(seg, t2);
    double x1 = HalfOpenLineSegment_func_x(seg, t1);
    double x2 = HalfOpenLineSegment_func_x(seg, t2);
    double sin_tb = HalfOpenLineSegment_get_sin_t1(seg);                // angle of the line wrt horizontal
    double cos_tb = HalfOpenLineSegment_get_cos_t1(seg);
    double l  = sqrt((s2 - s1)*(s2 - s1) + (x2 - x1)*(x2 - x1));     // length of the box
    double w  = 0.0;                                                 // width of the box 
    double rC = sqrt(s1*s1 + x1*x1);                                 // length of the position vector to the first vertex
    double sin_tC, cos_tC;
    if (rC == 0.){
        sin_tC = 0.0;  // angle of the position vector to the first vertex
        cos_tC = 1.0;
    } else {
        sin_tC = x1 / rC;  // angle of the position vector to the first vertex
        cos_tC = s1 / rC;
    }
    BoundingBox_set(box, rC, sin_tC, cos_tC, l, w, sin_tb, cos_tb);
}

#endif /* XCOLL_GEOM_SEG_HALFOPENLINE_H */