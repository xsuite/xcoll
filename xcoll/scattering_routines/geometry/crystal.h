// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_CRY_H
#define XCOLL_GEOM_CRY_H
#include <math.h>
#include <stdio.h>
"""
Did not have time to finish at all, and it's sort of messy. Not sure it's useable. Good luck hehehe
"""


// This function compares a particle trajectory to a crystal layout
double get_s_of_first_crossing_with_crystal(double part_x, double part_tan, double theta_bend, double tilt, /
                                            double R, double width, double length, double jaw_LU, double jaw_LD){
    double s = 1.e21;
    double new_s;
    double delta;

    // if tilt is zero ~ ah can maybe remove if else 
    // probs not valid for +-R 
    if (tilt == 0){
        double AB_s[2] = {0.,0.};                                  // don't see how they could be L/2 
        double AB_x[2] = {jaw_LU + width, jaw_LU};                 // no R since we start at {0,0}
        double CD_s[2] = {length, length - width*sin(theta_bend)}  // cos pi/2 -theta = sin
        double CD_x[2] = {jaw_LD +, jaw_LD + } // do we have the height of the box (do we have width) xdim,ydim?
    } else { 
        double AB_s[2] = {length/2*(1-cos(tilt)) - sin(tilt)*width, length/2*(1-cos(tilt))};              
        double AB_x[2] = {jaw_LU + width*cos(tilt) , jaw_LU};                   
        double CD_s[2] = {length + R*sin(tilt), length + R*sin(tilt) - } // need the proj. of the width term  
        double CD_x[2] = {jaw_LD +, jaw_LD + } // do we have the height of the box (do we have width) xdim,ydim?
    }

    // A and B
    new_s = get_s_of_first_crossing_with_polygon(part_x, part_tan, AB_s, AB_x, 2);
    if (new_s < s){
        s = new_s
    }
    // BC
    double delta; //tbd, what is w, and delta R ? 
    if (delta >= 0) {
        new_s = get_new_s // calculate new s
        if (s >= AB_s[1] && s <= CD_s[0] && new_s < s){     // ensure we check the arc and not the circle
            s = new_s
        }
    }
    // C and D
    new_s = get_s_of_first_crossing_with_polygon(part_x, part_tan, CD_s, CD_x, 2);
    if (new_s < s){
        s = new_s
    }
    // DA
    double delta; // tbd
    if (delta >= 0) {
        new_s = get_new_s // calculate new s  
        if (s >= AB_s[0] && s <= CD_s[1] && new_s < s){
            s = new_s
        }
    }
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