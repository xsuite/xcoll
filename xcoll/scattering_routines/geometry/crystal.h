// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_CRY_H
#define XCOLL_GEOM_CRY_H
#include <math.h>
#include <stdio.h>

// This function compares a particle trajectory to a crystal layout
double get_s_of_first_crossing_with_crystal(double part_x, double part_tan, double* poly_s, \
                                            double* poly_x, int8_t num_polys){
    double s = 1.e21;
    // double s_A = length/2*(1-cos_y);
    // double s_B = length/2*(1-cos_y);
    // double x_A =;
    // double x_B = jaw_U;
    // double trajectory_p0 = poly_x[0] - part_x - poly_s[0]*part_tan;
    // double trajectory_p1 = trajectory_p0;
    // double trajectory_p2, new_s;

    // for (int8_t i=0; i<num_polys-1; i++){
    //     trajectory_p2 = poly_x[i+1] - part_x - poly_s[i+1]*part_tan;
    //     if (trajectory_p1*trajectory_p2 <= 0){
    //         // It's a hit
    //         new_s = _get_s_of_crossing_with_segment(part_x, part_tan, poly_s[i], poly_s[i+1], \
    //                                                 poly_x[i], poly_x[i+1]);
    //         if (new_s < s){
    //             s = new_s;
    //         }
    //     }
    //     trajectory_p1 = trajectory_p2;
    // }
}


#endif /* XCOLL_GEOM_CRY_H */