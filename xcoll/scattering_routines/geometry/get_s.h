// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_GET_S_H
#define XCOLL_GEOM_GET_S_H


// IMPORTANT:
// These functions assume that the particle moves towards positive s!
// (hence no backscattering/backtracking is allowed)


// Find the s-coordinate of the first crossing of a drift with a set of segments
/*gpufun*/
double crossing_drift_first(Segment* segments, int8_t n_segments, \
                            double part_s, double part_x, double part_tan){
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift(segments, n_segments, &n_hit, s, part_s, part_x, part_tan);
    if (n_hit==0){
        // No crossing
        free(s);
        return XC_S_MAX;
    }
    double result = s[0];
    free(s);
    return result;
}

// Find the s-coordinate of the first crossing after a given s of a drift with a set of segments
/*gpufun*/
double crossing_drift_after_s(Segment* segments, int8_t n_segments, \
                              double part_s, double part_x, double part_tan, \
                              double current_s){
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift(segments, n_segments, &n_hit, s, part_s, part_x, part_tan);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= current_s){
            double result = s[i];
            free(s);
            return result;
        }
    }
    // No crossing
    free(s);
    return XC_S_MAX;
}

// Find the s-coordinate of the first crossing of a drift with a set of segments including a vertical restriction
/*gpufun*/
double crossing_drift_vlimit_first(Segment* segments, int8_t n_segments, \
                                   double part_s, double part_x, double part_tan_x, \
                                   double part_y, double part_tan_y, \
                                   double y_min, double y_max){
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift_vlimit(segments, n_segments, &n_hit, s, part_s, part_x, part_tan_x, \
                          part_y, part_tan_y, y_min, y_max);
    if (n_hit==0){
        // No crossing
        free(s);
        return XC_S_MAX;
    }
    double result = s[0];
    free(s);
    return result;
}

// Find the s-coordinate of the first crossing after a given s of a drift with a set of segments including a vertical restriction
/*gpufun*/
double crossing_drift_vlimit_after_s(Segment* segments, int8_t n_segments, \
                                     double part_s, double part_x, double part_tan_x, \
                                     double part_y, double part_tan_y, \
                                     double y_min, double y_max, double current_s){
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift_vlimit(segments, n_segments, &n_hit, s, part_s, part_x, part_tan_x, \
                          part_y, part_tan_y, y_min, y_max);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= current_s){
            double result = s[i];
            free(s);
            return result;
        }
    }
    // No crossing
    free(s);
    return XC_S_MAX;
}


#endif /* XCOLL_GEOM_GET_S_H */