// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_GET_S_H
#define XCOLL_GEOM_GET_S_H
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <headers/track.h>
#include <xcoll/scattering_routines/geometry/methods.h>

// IMPORTANT:
// These functions assume that the particle moves towards positive s!
// (hence no backscattering/backtracking is allowed)


GPUFUN
double get_s_of_first_crossing(double part_x, double part_tan, Segment* segments, \
                               int8_t n_segments){
    int8_t n_hit = 0;
    double* s = (double*) malloc(XC_MAX_CROSS_PER_SEGMENT*n_segments*sizeof(double));
    find_crossing(&n_hit, s, part_x, part_tan, segments, n_segments);
    if (n_hit==0){
        // No crossing
        free(s);
        return S_MAX;
    }
    double result = s[0];
    free(s);
    return result;
}

GPUFUN
double get_s_of_crossing_after_s(double part_x, double part_tan, Segment* segments, \
                                 int8_t n_segments, double current_s){
    int8_t n_hit = 0;
    double* s = (double*) malloc(XC_MAX_CROSS_PER_SEGMENT*n_segments*sizeof(double));
    find_crossing(&n_hit, s, part_x, part_tan, segments, n_segments);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= current_s){
            double result = s[i];
            free(s);
            return result;
        }
    }
    // No crossing
    free(s);
    return S_MAX;
}

GPUFUN
double get_s_of_first_crossing_with_vlimit(double part_x, double part_tan_x, \
                                double part_y, double part_tan_y, Segment* segments, \
                                int8_t n_segments, double y_min, double y_max){
    int8_t n_hit = 0;
    double* s = (double*) malloc(XC_MAX_CROSS_PER_SEGMENT*n_segments*sizeof(double));
    find_crossing_with_vlimit(&n_hit, s, part_x, part_tan_x, part_y, part_tan_y, \
                              segments, n_segments, y_min, y_max);
    if (n_hit==0){
        // No crossing
        free(s);
        return S_MAX;
    }
    double result = s[0];
    free(s);
    return result;
}

GPUFUN
double get_s_of_crossing_after_s_with_vlimit(double part_x, double part_tan_x, \
                                double part_y, double part_tan_y, Segment* segments, \
                                int8_t n_segments, double y_min, double y_max, double current_s){
    int8_t n_hit = 0;
    double* s = (double*) malloc(XC_MAX_CROSS_PER_SEGMENT*n_segments*sizeof(double));
    find_crossing_with_vlimit(&n_hit, s, part_x, part_tan_x, part_y, part_tan_y, \
                              segments, n_segments, y_min, y_max);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= current_s){
            double result = s[i];
            free(s);
            return result;
        }
    }
    // No crossing
    free(s);
    return S_MAX;
}


#endif /* XCOLL_GEOM_GET_S_H */