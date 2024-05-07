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
    return s;
}


#endif /* XCOLL_GEOM_CRY_H */