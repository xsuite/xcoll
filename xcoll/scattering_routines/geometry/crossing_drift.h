// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_CROSSING_DRIFT_H
#define XCOLL_GEOM_CROSSING_DRIFT_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


#define XC_DRIFT_MAX_CROSS_PER_SEGMENT 2  // Update if new segment type allows more crossings


// All segments have a function that calculates the s-coordinate of the crossing(s) with a
// particle trajectory in a drift (given by a line through (part_s, part_x) and a slope part_tan).
// The results are always stored in an array s, and n_hit keeps track of the number of hits over
// multiple segments.


// Line Segments
// -------------

/*gpufun*/
void _crossing_drift_line(int8_t* n_hit, double* s, double s1, double s2, double x1, double x2,
                          double part_s, double part_x, double part_tan){
    // We fill in the segment points in the particle trajectory equation; if the results have opposite sign,
    //the two points lie on different sides of the trajectory and hence the segment is crossed.
    double trajectory1 = x1 - part_x - (s1 - part_s)*part_tan;
    double trajectory2 = x2 - part_x - (s2 - part_s)*part_tan;
    if (trajectory1*trajectory2 <= 0){
        // It's a crossing
        if (fabs(s2 - s1) < 1.e-12){
            s[*n_hit] = s1;
            (*n_hit)++;
        } else {
            double poly_tan = (x2 - x1)/(s2 - s1);
            if (fabs(poly_tan - part_tan) < 1.e-12){
                // The only case where we have two hits is when the trajectory is parallel to the segment
                s[*n_hit] = s1;
                (*n_hit)++;
                s[*n_hit] = s2;
                (*n_hit)++;
            } else {
                // Normal crossing of two lines
                s[*n_hit] = (part_x - x1 + s1*poly_tan - part_s*part_tan)/(poly_tan - part_tan);
                (*n_hit)++;
            }
        }
    }
}

/*gpufun*/
void crossing_drift_line(void* self, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    LineSegment seg = (LineSegment) self;
    double s1 = seg->point1_s;
    double s2 = seg->point2_s;
    double x1 = seg->point1_x;
    double x2 = seg->point2_x;
    _crossing_drift_line(n_hit, s, s1, s2, x1, x2, part_s, part_x, part_tan);
}

/*gpufun*/
void crossing_drift_halfopenline(void* self, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    HalfOpenLineSegment seg = (HalfOpenLineSegment) self;
    double s1 = seg->point_s;
    double x1 = seg->point_x;
    // A half-open segment implies one of its points lies at +-inf.
    // In practice we just add a polygon point at the wall overflow (at 1000km for the x-coordinate).
    double x2 = 1.e6*seg->sign;
    double s2;
    if (fabs(seg->inv_slope) < 1.e-12){
        s2 = s1;
    } else {
        s2 = s1 - (x2 - x1)*seg->inv_slope;
    }
    _crossing_drift_line(n_hit, s, s1, s2, x1, x2, part_s, part_x, part_tan);
}


// Circular Segment
// ----------------

/*gpufun*/
void crossing_drift_circular(void* self, int8_t* n_hit, double* s, double part_s, double part_x, double part_tan){
    // Get segment data
    CircularSegment seg = (CircularSegment) self;
    double R   = seg->R;
    double sC = seg->centre_s;
    double xC = seg->centre_x;
    double t1  = seg->point1_angle;
    double t2  = seg->point2_angle;
    // Calculate crossings
    int8_t reversed = 0;
    if (t2 < t1){
        reversed = 1;
    }
    double a = 1 + part_tan*part_tan;
    double bb = sC - part_tan*(part_x - xC - part_tan*part_s); // This is -b/2 with b from the quadratic formula
    double c = sC*sC + (part_x - xC - part_tan*part_s)*(part_x - xC - part_tan*part_s) - R*R;
    double disc = bb*bb - a*c; // This is  2*discriminant**2
    if (disc >= 0){
        for (int8_t i = 0; i < 2; i++) {
            double sgnD = i*2-1; // negative and positive solutions
            double new_s = 1/a*(bb + sgnD*sqrt(disc));
            double new_x = part_x + (new_s - part_s)*part_tan;
            double t = atan2(new_x - xC, new_s - sC);
            if (reversed){
                // t2 < t1, so we are looking at the inverted region of angles
                if (t1 >= t || t >= t2){
                    s[*n_hit] = new_s;
                    (*n_hit)++;
                }
            } else {
                if (t1 <= t && t <= t2){
                    s[*n_hit] = new_s;
                    (*n_hit)++;
                }
            }
        }
    }
}


// Array of segments
// -----------------

/*gpufun*/
void crossing_drift(Segment* segments, int8_t n_segments, int8_t* n_hit, double* s, \
                    double part_s, double part_x, double part_tan){
    for (int8_t i=0; i<n_segments;i++) {
        int id = segments[i]->id;
        if (id == XC_LINESEGMENT_ID){
            crossing_drift_line(segments[i], n_hit, s, part_s, part_x, part_tan);
        } else if (id == XC_HALFOPENLINESEGMENT_ID){
            crossing_drift_halfopenline(segments[i], n_hit, s, part_s, part_x, part_tan);
        } else if (id == XC_CIRCULARSEGMENT_ID){
            crossing_drift_circular(segments[i], n_hit, s, part_s, part_x, part_tan);
        } // TODO: else throw fatal error
    }
    sort_array_of_double(s, (int64_t) *n_hit);
}

/*gpufun*/
void crossing_drift_vlimit(Segment* segments, int8_t n_segments, int8_t* n_hit, double* s, \
                           double part_s, double part_x, double part_tan_x, \
                           double part_y, double part_tan_y, \
                           double y_min, double y_max){
    if (fabs(part_tan_y) < 1.e-12){
        // Trajectory parallel to s axis
        if (part_y < y_min || part_y > y_max){
            // No crossing
            return;
        } else {
            // The particle is completely inside the vertical limits, so only check horizontal
            crossing_drift(segments, n_segments, n_hit, s, part_s, part_x, part_tan_x);
            return;
        }
    } else {
        crossing_drift(segments, n_segments, n_hit, s, part_s, part_x, part_tan_x);
        // restrict_s is the region [s0, s1] where the particle is inside the vertical limits
        double* restrict_s = (double*) malloc(2*sizeof(double));
        restrict_s[0] = (y_min - part_y)/part_tan_y + part_s;
        restrict_s[1] = (y_max - part_y)/part_tan_y + part_s;
        SWAP(restrict_s, 0, 1);   // To make sure these are sorted
        calculate_overlap_array_interval(s, n_hit, restrict_s);
        free(restrict_s);
    }
}

#endif /* XCOLL_GEOM_CROSSING_DRIFT_H */