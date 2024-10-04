// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SEGMENTS_H
#define XCOLL_GEOM_SEGMENTS_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


// Parent type for all segments
// ----------------------------
typedef struct Segment_{
    int id;
} Segment_;
typedef Segment_* Segment;


// Line segment between two points (s1, x1) -- (s2, x2)
// -------------------------------------------------

#define XC_LINESEGMENT_ID 0
typedef struct LineSegment_ {
    int id;
    double point1_s;
    double point1_x;
    double point2_s;
    double point2_x;
} LineSegment_;
typedef LineSegment_* LineSegment;

/*gpufun*/
LineSegment create_line_segment(double point1_s, double point1_x, double point2_s, double point2_x){
    LineSegment seg = (LineSegment) malloc(sizeof(LineSegment_));
    seg->id = XC_LINESEGMENT_ID;
    seg->point1_s = point1_s;
    seg->point2_s = point2_s;
    seg->point1_x = point1_x;
    seg->point2_x = point2_x;
    return seg;
}


// Half-open line segment from a point (s, x) to +/-infinity along a slope
// -----------------------------------------------------------------------

#define XC_HALFOPENLINESEGMENT_ID 1
// For practical reasons, we store the inverse of the slope (inv_slope = -1/slope)
// of the segment. So inv_slope = 0 means the segment is vertical, inv_slope = 1
// implies a segment at an angle of 135 deg (tilt of +45 degrees), and a horizontal
// half-open segment is not allowed (and also not needed as it would go to infinity
// along the beam).
typedef struct HalfOpenLineSegment_ {
    int id;
    double point_s;
    double point_x;
    double inv_slope;
    int8_t sign;  // Does the segment go to +inf or -inf?
} HalfOpenLineSegment_;
typedef HalfOpenLineSegment_* HalfOpenLineSegment;

/*gpufun*/
HalfOpenLineSegment create_halfopen_line_segment(double point_s, double point_x, double inv_slope, int8_t sign){
    HalfOpenLineSegment seg = (HalfOpenLineSegment) malloc(sizeof(HalfOpenLineSegment_));
    seg->id = XC_HALFOPENLINESEGMENT_ID;
    seg->point_s = point_s;
    seg->point_x = point_x;
    seg->inv_slope = inv_slope;
    seg->sign = sign;
    return seg;
}


// Circular arc segment, defined by a centre and radius, and the starting/end angles
// ---------------------------------------------------------------------------------

#define XC_CIRCULARSEGMENT_ID 2
typedef struct CircularSegment_ {
    int id;
    double R;
    double centre_s;
    double centre_x;
    double point1_angle;
    double point2_angle;
} CircularSegment_;
typedef CircularSegment_* CircularSegment;

/*gpufun*/
CircularSegment create_circular_segment(double R, double centre_s, double centre_x, double point1_angle, double point2_angle){
    CircularSegment seg = (CircularSegment) malloc(sizeof(CircularSegment_));
    seg->id = XC_CIRCULARSEGMENT_ID;
    seg->R = R;
    seg->centre_s = centre_s;
    seg->centre_x = centre_x;
    seg->point1_angle = point1_angle;
    seg->point2_angle = point2_angle;
    return seg;
}


#endif /* XCOLL_GEOM_SEGMENTS_H */