// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_METHODS_H
#define XCOLL_GEOM_METHODS_H

#define XC_S_MAX 1.e21


// This function calculates the overlap between an array and a given interval.
// The array comes in pairs of points, e.g. in-out-in-out... or out-in-out-in...
// IMPORTANT:
// The array and interval are assumed to be sorted!
// Furthermore, the array should have one extra slot allocated at the end, in case it needs to be expanded..
// This is always true for the arrays created by get_s, as we create them with 2*n_segments slots.
/*gpufun*/
void calculate_overlap_array_interval(double* arr, int8_t* length, double* interval){
    if (arr[0] > interval[1]){
        // No overlap
        *length = 0;
    }
    if ((*length)%2 == 1){
        // Special case: last interval of array is open until infinity,
        // so we add an extra point at the end to represent infinity.
        arr[*length] = XC_S_MAX*0.1;
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


#endif /* XCOLL_GEOM_METHODS_H */