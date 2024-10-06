// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_GET_S_H
#define XCOLL_GEOM_GET_S_H


// Find the s-coordinate of the first crossing of a drift with a set of segments
// Here, first means the smallest s-coordinate (needs to be adapted for back-scattering)
/*gpufun*/
double crossing_drift_first(Segment* segments, int8_t n_segments, \
                            double part_s, double part_x, double part_tan){
    // printf("crossing_drift_first\n");fflush(stdout);
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift(segments, n_segments, &n_hit, s, part_s, part_x, part_tan);
    if (n_hit==0){
        // No crossing
        free(s);
        // printf("crossing_drift_first done (no crossing)\n");fflush(stdout);
        return XC_S_MAX;
    }
    double result = s[0];
    free(s);
    // printf("crossing_drift_first done\n");fflush(stdout);
    return result;
}

// Find the s-coordinate of the first crossing after a given s of a drift with a set of segments
/*gpufun*/
double crossing_drift_after_s(Segment* segments, int8_t n_segments, \
                              double part_s, double part_x, double part_tan, \
                              double after_s){
    // printf("crossing_drift_after_s\n");fflush(stdout);
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift(segments, n_segments, &n_hit, s, part_s, part_x, part_tan);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= after_s){
            double result = s[i];
            free(s);
            // printf("crossing_drift_after_s done\n");fflush(stdout);
            return result;
        }
    }
    // No crossing
    free(s);
    // printf("crossing_drift_after_s done (no crossing)\n");fflush(stdout);
    return XC_S_MAX;
}

// Find the s-coordinate of the first crossing of a drift with a set of segments including a vertical restriction
/*gpufun*/
double crossing_drift_vlimit_first(Segment* segments, int8_t n_segments, \
                                   double part_s, double part_x, double part_tan_x, \
                                   double part_y, double part_tan_y, \
                                   double y_min, double y_max){
    // printf("crossing_drift_vlimit_first\n");fflush(stdout);
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift_vlimit(segments, n_segments, &n_hit, s, part_s, part_x, part_tan_x, \
                          part_y, part_tan_y, y_min, y_max);
    if (n_hit==0){
        // No crossing
        free(s);
        // printf("crossing_drift_vlimit_first done (no crossing)\n");fflush(stdout);
        return XC_S_MAX;
    }
    double result = s[0];
    free(s);
    // printf("crossing_drift_vlimit_first done\n");fflush(stdout);
    return result;
}

// Find the s-coordinate of the first crossing after a given s of a drift with a set of segments including a vertical restriction
/*gpufun*/
double crossing_drift_vlimit_after_s(Segment* segments, int8_t n_segments, \
                                     double part_s, double part_x, double part_tan_x, \
                                     double part_y, double part_tan_y, \
                                     double y_min, double y_max, double after_s){
    // printf("crossing_drift_vlimit_after_s\n");fflush(stdout);
    int8_t n_hit = 0;
    double* s = (double*) malloc(max_array_size_drift(segments, n_segments)*sizeof(double));
    crossing_drift_vlimit(segments, n_segments, &n_hit, s, part_s, part_x, part_tan_x, \
                          part_y, part_tan_y, y_min, y_max);
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= after_s){
            double result = s[i];
            free(s);
            // printf("crossing_drift_vlimit_after_s done\n");fflush(stdout);
            return result;
        }
    }
    // No crossing
    free(s);
    // printf("crossing_drift_vlimit_after_s done (no crossing)\n");fflush(stdout);
    return XC_S_MAX;
}


#endif /* XCOLL_GEOM_GET_S_H */