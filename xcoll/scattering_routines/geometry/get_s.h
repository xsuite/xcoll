// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_GET_S_H
#define XCOLL_GEOM_GET_S_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdint.h>  // for int64_t etc
#endif  // XO_CONTEXT_CPU


// IMPORTANT:
// These functions assume that the particle moves towards positive s!
// (hence no backscattering/backtracking is allowed)


/*gpufun*/
double get_s_of_first_crossing(double part_x, double part_tan, Segment* segments, \
                               int8_t n_segments){
    int8_t n_hit = 0;
    // +1 spare slot: calculate_overlap_array_interval (methods.h) may write
    // arr[*length] one element past the logical end (see methods.h:43-44
    // contract). CPU tolerates the 1-element stack overrun; GPU traps it as
    // cudaErrorIllegalAddress.
    double s[XC_MAX_CROSS_PER_SEGMENT*XC_MAX_SEGMENTS + 1];
    find_crossing(&n_hit, s, part_x, part_tan, segments, n_segments);
    if (n_hit==0){
        // No crossing
        return S_MAX;
    }
    double result = s[0];
    return result;
}

/*gpufun*/
double get_s_of_crossing_after_s(double part_x, double part_tan, Segment* segments, \
                                 int8_t n_segments, double current_s){
    int8_t n_hit = 0;
    // +1 spare slot (see contract in methods.h:43-44 / first crossing above).
    double s[XC_MAX_CROSS_PER_SEGMENT*XC_MAX_SEGMENTS + 1];
    find_crossing(&n_hit, s, part_x, part_tan, segments, n_segments);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= current_s){
            double result = s[i];
            return result;
        }
    }
    // No crossing
    return S_MAX;
}

/*gpufun*/
double get_s_of_first_crossing_with_vlimit(double part_x, double part_tan_x, \
                                double part_y, double part_tan_y, Segment* segments, \
                                int8_t n_segments, double y_min, double y_max){
    int8_t n_hit = 0;
    // +1 spare slot (see contract in methods.h:43-44 / first crossing above).
    double s[XC_MAX_CROSS_PER_SEGMENT*XC_MAX_SEGMENTS + 1];
    find_crossing_with_vlimit(&n_hit, s, part_x, part_tan_x, part_y, part_tan_y, \
                              segments, n_segments, y_min, y_max);
    if (n_hit==0){
        // No crossing
        return S_MAX;
    }
    double result = s[0];
    return result;
}

/*gpufun*/
double get_s_of_crossing_after_s_with_vlimit(double part_x, double part_tan_x, \
                                double part_y, double part_tan_y, Segment* segments, \
                                int8_t n_segments, double y_min, double y_max, double current_s){
    int8_t n_hit = 0;
    // +1 spare slot (see contract in methods.h:43-44 / first crossing above).
    double s[XC_MAX_CROSS_PER_SEGMENT*XC_MAX_SEGMENTS + 1];
    find_crossing_with_vlimit(&n_hit, s, part_x, part_tan_x, part_y, part_tan_y, \
                              segments, n_segments, y_min, y_max);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= current_s){
            double result = s[i];
            return result;
        }
    }
    // No crossing
    return S_MAX;
}


#endif /* XCOLL_GEOM_GET_S_H */