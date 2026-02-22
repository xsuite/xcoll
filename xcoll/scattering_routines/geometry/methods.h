// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_METHODS_H
#define XCOLL_GEOM_METHODS_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdlib.h>
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU

#include <xobjects/headers/common.h>
#include <xcoll/scattering_routines/geometry/objects.h>


#define S_MAX 1.e21


GPUFUN
void find_crossing(int8_t* n_hit, double* s, double part_x, double part_tan, \
                       Segment* segments, int8_t n_segments){
    for (int8_t i=0; i<n_segments;i++) {
        segments[i]->crossing(n_hit, s, part_x, part_tan, segments[i]);
    }
    sort_array_of_double(s, (int64_t) *n_hit);
}


// IMPORTANT:
// The array and interval are assumed to be sorted!
// Furthermore, the array should have one extra slot allocated at the end, in case it needs to be expanded..
// This is always true for the arrays created by get_s, as we create them with 2*n_segments slots.
GPUFUN
void calculate_overlap_array_interval(double* arr, int8_t* length, double* interval){
    if (arr[0] > interval[1]){
        // No overlap
        *length = 0;
    }
    if ((*length)%2 == 1){
        // Special case: last interval of array is open until infinity,
        // so we add an extra point at the end to represent infinity.
        arr[*length] = S_MAX*0.1;
        (*length)++;
    } else if (arr[*length-1] < interval[0]){
        // No overlap
        *length = 0;
    }
    int8_t i_start = 0;
    // Find the start of overlap
    for (int8_t i=0; i<*length; i++){
        if (arr[i] >= interval[0]){
            if (i%2 == 0){
                // This is the first point of overlap
                i_start = i;
            } else {
                // The vertical restriction is the first point of overlap
                i_start = i-1;
                arr[i_start] = interval[0];
            }
            break;
        }
    }
    // Find the end of overlap
    int8_t i_stop = *length - 1;
    for (int8_t i=0; i<*length; i++){
        if (arr[i] >= interval[1]){
            if (i%2 == 0){
                // The previous point is the last point of overlap
                i_stop = i-1;
            } else {
                // The vertical restriction is the first point of overlap
                i_stop = i;
                arr[i_stop] = interval[1];
            }
            break;
        }
    }
    *length = i_stop - i_start + 1;
    if (i_start > 0){
        for (int8_t i=0; i<*length; i++){
            arr[i] = arr[i+i_start];
        }
    }
}


GPUFUN
void find_crossing_with_vlimit(int8_t* n_hit, double* s, double part_x, double part_tan_x, \
                               double part_y, double part_tan_y, Segment* segments, \
                               int8_t n_segments, double y_min, double y_max){
    if (fabs(part_tan_y) < 1.e-12){
        // Trajectory parallel to s axis
        if (part_y < y_min || part_y > y_max){
            // No crossing
            return;
        } else {
            // The particle is completely inside the vertical limits, so only check horizontal
            find_crossing(n_hit, s, part_x, part_tan_x, segments, n_segments);
            return;
        }
    } else {
        find_crossing(n_hit, s, part_x, part_tan_x, segments, n_segments);
        // restrict_s is the region [s0, s1] where the particle is inside the vertical limits
        double* restrict_s = (double*) malloc(2*sizeof(double));
        restrict_s[0] = (y_min - part_y)/part_tan_y;
        restrict_s[1] = (y_max - part_y)/part_tan_y;
        SWAP(restrict_s, 0, 1);   // To make sure these are sorted
        calculate_overlap_array_interval(s, n_hit, restrict_s);
        free(restrict_s);
    }
}

#endif /* XCOLL_GEOM_METHODS_H */
