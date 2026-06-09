// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_GET_S_H
#define XCOLL_GEOM_GET_S_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU


// IMPORTANT:
// These functions assume that the particle moves towards positive s!
// (hence no backscattering/backtracking is allowed)


// Crossing-buffer allocation.
// find_crossing writes up to XC_MAX_CROSS_PER_SEGMENT*n_segments crossings, and
// calculate_overlap_array_interval (methods.h) may write one element past that
// (arr[*length]), so the buffer needs XC_MAX_CROSS_PER_SEGMENT*n_segments + 1 slots.
//  - On CPU any element/polygon may be tracked, including the geometry-test
//    generator's polygons with more than XC_MAX_SEGMENTS sides, so the buffer is
//    sized by the actual n_segments (heap).
//  - On GPU only EverestCrystal tracks, and it is bounded by XC_MAX_SEGMENTS, so a
//    fixed stack buffer avoids malloc (which is unavailable on the CUDA context).
#ifdef XO_CONTEXT_CPU
    #define XC_DECLARE_CROSSING_BUFFER(s, n_segments) \
        double* s = (double*) malloc((XC_MAX_CROSS_PER_SEGMENT*(n_segments) + 1)*sizeof(double))
    #define XC_FREE_CROSSING_BUFFER(s)  free(s)
#else
    #define XC_DECLARE_CROSSING_BUFFER(s, n_segments) \
        double s[XC_MAX_CROSS_PER_SEGMENT*XC_MAX_SEGMENTS + 1]
    #define XC_FREE_CROSSING_BUFFER(s)  ((void) 0)
#endif


/*gpufun*/
double get_s_of_first_crossing(double part_x, double part_tan, Segment* segments, \
                               int8_t n_segments){
    int8_t n_hit = 0;
    XC_DECLARE_CROSSING_BUFFER(s, n_segments);
    find_crossing(&n_hit, s, part_x, part_tan, segments, n_segments);
    if (n_hit==0){
        // No crossing
        XC_FREE_CROSSING_BUFFER(s);
        return S_MAX;
    }
    double result = s[0];
    XC_FREE_CROSSING_BUFFER(s);
    return result;
}

/*gpufun*/
double get_s_of_crossing_after_s(double part_x, double part_tan, Segment* segments, \
                                 int8_t n_segments, double current_s){
    int8_t n_hit = 0;
    XC_DECLARE_CROSSING_BUFFER(s, n_segments);
    find_crossing(&n_hit, s, part_x, part_tan, segments, n_segments);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= current_s){
            double result = s[i];
            XC_FREE_CROSSING_BUFFER(s);
            return result;
        }
    }
    // No crossing
    XC_FREE_CROSSING_BUFFER(s);
    return S_MAX;
}

/*gpufun*/
double get_s_of_first_crossing_with_vlimit(double part_x, double part_tan_x, \
                                double part_y, double part_tan_y, Segment* segments, \
                                int8_t n_segments, double y_min, double y_max){
    int8_t n_hit = 0;
    XC_DECLARE_CROSSING_BUFFER(s, n_segments);
    find_crossing_with_vlimit(&n_hit, s, part_x, part_tan_x, part_y, part_tan_y, \
                              segments, n_segments, y_min, y_max);
    if (n_hit==0){
        // No crossing
        XC_FREE_CROSSING_BUFFER(s);
        return S_MAX;
    }
    double result = s[0];
    XC_FREE_CROSSING_BUFFER(s);
    return result;
}

/*gpufun*/
double get_s_of_crossing_after_s_with_vlimit(double part_x, double part_tan_x, \
                                double part_y, double part_tan_y, Segment* segments, \
                                int8_t n_segments, double y_min, double y_max, double current_s){
    int8_t n_hit = 0;
    XC_DECLARE_CROSSING_BUFFER(s, n_segments);
    find_crossing_with_vlimit(&n_hit, s, part_x, part_tan_x, part_y, part_tan_y, \
                              segments, n_segments, y_min, y_max);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= current_s){
            double result = s[i];
            XC_FREE_CROSSING_BUFFER(s);
            return result;
        }
    }
    // No crossing
    XC_FREE_CROSSING_BUFFER(s);
    return S_MAX;
}

#undef XC_DECLARE_CROSSING_BUFFER
#undef XC_FREE_CROSSING_BUFFER

#endif /* XCOLL_GEOM_GET_S_H */
