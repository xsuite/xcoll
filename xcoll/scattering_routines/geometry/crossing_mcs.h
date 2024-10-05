// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_CROSSING_MCS_H
#define XCOLL_GEOM_CROSSING_MCS_H


// Line Segments
// -------------

/*gpufun*/
void crossing_mcs_line(void* segment, int8_t* n_hit, double* s, __mcs_vars__);

/*gpufun*/
void crossing_mcs_halfopenline(void* segment, int8_t* n_hit, double* s, __mcs_vars__);


// Circular Segment
// ----------------

/*gpufun*/
void crossing_mcs_circular(void* segment, int8_t* n_hit, double* s, __mcs_vars__);


// Bézier Segment
// --------------

/*gpufun*/
void crossing_mcs_bezier(void* segment, int8_t* n_hit, double* s, __mcs_vars__);


// Array of segments
// -----------------

/*gpufun*/
void crossing_mcs(Segment* segments, int8_t n_segments, int8_t* n_hit, double* s, __mcs_vars__){}
    for (int8_t i=0; i<n_segments;i++) {
        int id = segments[i]->id;
        if (id == XC_LINESEGMENT_ID){
            crossing_mcs_line(segments[i], n_hit, s, __mcs_vars__);
        } else if (id == XC_HALFOPENLINESEGMENT_ID){
            crossing_mcs_halfopenline(segments[i], n_hit, s, __mcs_vars__);
        } else if (id == XC_CIRCULARSEGMENT_ID){
            crossing_mcs_circular(segments[i], n_hit, s, __mcs_vars__);
        } else if (id == XC_BEZIERSEGMENT_ID){
            crossing_mcs_bezier(segments[i], n_hit, s, __mcs_vars__);
        } // TODO: else throw fatal error
    }
    sort_array_of_double(s, (int64_t) *n_hit);
}

/*gpufun*/
void crossing_mcs_vlimit(Segment* segments, int8_t n_segments, int8_t* n_hit, double* s, \
                           __mcs_vars_V__, double y_min, double y_max){
    // Need implementation with MCS vertically
}

/*gpufun*/
int max_array_size_mcs(Segment* segments, int8_t n_segments){
    int size = 0;
    for (int8_t i=0; i<n_segments;i++) {
        int id = segments[i]->id;
        if (id == XC_LINESEGMENT_ID){
            size += 2;
        } else if (id == XC_HALFOPENLINESEGMENT_ID){
            size += 2;
        } else if (id == XC_CIRCULARSEGMENT_ID){
            size += 4;
        } else if (id == XC_BEZIERSEGMENT_ID){
            size += 5;
        } // TODO: else throw fatal error
    }
    return size;
}

#endif /* XCOLL_GEOM_CROSSING_MCS_H */
