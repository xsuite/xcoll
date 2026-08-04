// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SEGMENTS_H
#define XCOLL_GEOM_SEGMENTS_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#endif /* XO_CONTEXT_CPU */

#include "xobjects/headers/common.h"


// These functions compare a particle trajectory (straight line with slope part_tan going
// through the point [0, part_x]) with a given segment of specific type.
// The results are always stored in an array s, and n_hit keeps track of the number of hits.

#define XC_MAX_SEGMENTS 8
#define XC_MAX_CROSS_PER_SEGMENT 2  // Update if new segment type allows more crossings
#define XC_MAX_CROSSINGS (XC_MAX_SEGMENTS*XC_MAX_CROSS_PER_SEGMENT+1)

// Segment type tag, used for switch-based crossing dispatch.
// Device function pointers are unsupported on CUDA, so instead of a function-pointer
// vtable we store an int8_t tag as the first member of every segment struct and
// dispatch on it via a switch in find_crossing (only 3 subtypes -> cheap switch).
#define XC_SEGMENT_LINE     ((int8_t) 0)
#define XC_SEGMENT_HALFOPEN ((int8_t) 1)
#define XC_SEGMENT_CIRCULAR ((int8_t) 2)


// Segment value struct
// --------------------
// All three segment subtypes (line, half-open line, circular arc) share a single
// fixed-layout tagged-union struct so segments can be stored BY VALUE in a fixed-size
// array (no malloc/free, no device function pointers). The `type` tag selects which
// fields are meaningful; the field names of the three subtypes are disjoint, so we
// flatten them into one struct (a tagged union without per-member casts). A Segment is
// stored by value, so functions that scan an array of segments take a `Segment*`
// (pointer to the first element) exactly as before, but the array slots are values.
typedef struct Segment {
    int8_t type;
    // Line segment: two points (s1, x1), (s2, x2).
    double point1_s;
    double point1_x;
    double point2_s;
    double point2_x;
    // Half-open line segment: a point, a side (+-inf direction) and the perpendicular angle.
    double point_s;
    double point_x;
    int8_t side;
    double point_tan;
    // Circular arc segment: radius R, centre (Rs, Rx) and the two bounding angles.
    double R;
    double centre_s;
    double centre_x;
    double point1_angle;
    double point2_angle;
} Segment;


// Line segment
// ------------

// Line segment defined by two points (s1, x1), (s2, x2). The function will fill in the segment
// points in the trajectory equation; if the results have opposite sign, the two points lie on
// different sides of the trajectory and hence the segment is crossed.

GPUFUN void get_s_of_crossing_with_line_segment(int8_t* n_hit, double* s, double part_x, double part_tan, void* self) {
    // Get segment data
    Segment* seg = (Segment*) self;
    double s_p1 = seg->point1_s;
    double s_p2 = seg->point2_s;
    double x_p1 = seg->point1_x;
    double x_p2 = seg->point2_x;
    // Calculate crossings
    double trajectory_p1 = x_p1 - part_x - s_p1*part_tan;
    double trajectory_p2 = x_p2 - part_x - s_p2*part_tan;
    if (trajectory_p1*trajectory_p2 <= 0) {
        // It's a crossing
        if (fabs(s_p2 - s_p1) < 1.e-12) {
            s[*n_hit] = s_p1;
            (*n_hit)++;
        } else {
            double poly_tan = (x_p2 - x_p1)/(s_p2 - s_p1);
            if (fabs(poly_tan - part_tan) < 1.e-12) {
                s[*n_hit] = s_p1;
                (*n_hit)++;
            } else {
                s[*n_hit] = (part_x - x_p1 + s_p1*poly_tan)/(poly_tan - part_tan);
                (*n_hit)++;
            }
        }
    }
}

GPUFUN Segment create_line_segment(double point1_s, double point1_x, double point2_s, double point2_x) {
    Segment seg;
    seg.type = XC_SEGMENT_LINE;
    seg.point1_s = point1_s;
    seg.point2_s = point2_s;
    seg.point1_x = point1_x;
    seg.point2_x = point2_x;
    return seg;
}


// Half-open line segment
// ----------------------

// This function works as above, but considers a half-open segment. For this reason,
// the function needs to know whether this is a positive or negative jaw (representing
// to which infinite it points), and the overall tilt of the jaw.

// A half-open segment implies one of its points lies at +-inf.
// In practice we just add a polygon point at the wall overflow (at 1km for the x-coordinate).

GPUFUN void get_s_of_crossing_with_halfopen_line_segment(int8_t* n_hit, double* s, double part_x, double part_tan, void* self) {
    // Get segment data
    Segment* seg = (Segment*) self;
    double s_p1 = seg->point_s;
    double x_p1 = seg->point_x;
    double x_p2 = 1.e3*seg->side;
    double s_p2;
    if (fabs(seg->point_tan) < 1.e-12) {
        s_p2 = s_p1;
    } else {
        s_p2 = -(x_p2 - x_p1 - s_p1/seg->point_tan)*seg->point_tan;
    }
    // Calculate crossings
    double trajectory_p1 = x_p1 - part_x - s_p1*part_tan;
    double trajectory_p2 = x_p2 - part_x - s_p2*part_tan;
    if (trajectory_p1*trajectory_p2 <= 0) {
        // It's a crossing
        if (fabs(s_p2 - s_p1) < 1.e-12) {
            s[*n_hit] = s_p1;
            (*n_hit)++;
        } else {
            double poly_tan = (x_p2 - x_p1)/(s_p2 - s_p1);
            if (fabs(poly_tan - part_tan) < 1.e-12) {
                s[*n_hit] = s_p1;
                (*n_hit)++;
            } else {
                s[*n_hit] = (part_x - x_p1 + s_p1*poly_tan)/(poly_tan - part_tan);
                (*n_hit)++;
            }
        }
    }
}

GPUFUN Segment create_halfopen_line_segment(double point_s, double point_x, double point_tan, int8_t side) {
    Segment seg;
    seg.type = XC_SEGMENT_HALFOPEN;
    seg.point_s = point_s;
    seg.point_x = point_x;
    seg.side = side;
    seg.point_tan = point_tan;
    return seg;
}


// Circular arc segment
// --------------------

// This function finds the crossing points between a line defined by a point (s=0, x) and a tangent,
// and a circular arc segment defined by a radius R, a centre (Rs, Rx), and angles t1 and t2.
// The results are stored in an array s, and n_hit keeps track of the number of hits.

GPUFUN void get_s_of_crossing_with_circular_segment(int8_t* n_hit, double* s, double part_x, double part_tan, void* self) {
    // Get segment data
    Segment* seg = (Segment*) self;
    double R   = seg->R;
    double R_s = seg->centre_s;
    double R_x = seg->centre_x;
    double t1  = seg->point1_angle;
    double t2  = seg->point2_angle;
    // Calculate crossings
    int8_t reversed = 0;
    if (t2 < t1) {
        reversed = 1;
    }
    double a = 1 + part_tan*part_tan;
    double bb = R_s - part_tan*(part_x - R_x);
    double c = R_s*R_s + (part_x - R_x)*(part_x - R_x) - R*R;
    double disc = bb*bb - a*c;
    if (disc >= 0) {
        for (int8_t i = 0; i < 2; i++) {
            double sgnD = i*2-1; // negative and positive solutions
            double new_s = 1/a*(bb + sgnD*sqrt(bb*bb - a*c));
            double x = part_x + new_s*part_tan;
            double t = atan2(x - R_x, new_s - R_s);
            if (reversed) {
                // t2 < t1, so we are looking at the inverted region of angles
                if (t1 >= t || t >= t2) {
                    s[*n_hit] = new_s;
                    (*n_hit)++;
                }
            } else {
                if (t1 <= t && t <= t2) {
                    s[*n_hit] = new_s;
                    (*n_hit)++;
                }
            }
        }
    }
}

GPUFUN Segment create_circular_segment(double R, double centre_s, double centre_x, double point1_angle, double point2_angle) {
    Segment seg;
    seg.type = XC_SEGMENT_CIRCULAR;
    seg.R = R;
    seg.centre_s = centre_s;
    seg.centre_x = centre_x;
    seg.point1_angle = point1_angle;
    seg.point2_angle = point2_angle;
    return seg;
}


#endif /* XCOLL_GEOM_SEGMENTS_H */
