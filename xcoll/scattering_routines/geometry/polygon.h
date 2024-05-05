// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_POLY_H
#define XCOLL_GEOM_POLY_H
#include <math.h>
#include <stdio.h>


double _get_s_of_crossing_with_segment(double part_x, double part_tan, double s_p1, double s_p2, \
                                       double x_p1, double x_p2){
    if (fabs(s_p2 - s_p1) < 1.e-12){
        return s_p1;
    } else {
        double poly_tan = (x_p2 - x_p1)/(s_p2 - s_p1);
        if (fabs(poly_tan - part_tan) < 1.e-12){
            return s_p1;
        } else {
            return (part_x - x_p1 + s_p1*poly_tan)/(poly_tan - part_tan);
        }
    }
}


// This function compares a particle trajectory (straight line with slope part_tan going through 
// [0, part_x]) with a polygon defined by a set of points (poly_s, poly_x). It will loop over the
// segments, and fill in the segment points in the trajectory equation. If they have opposite sign,
// they lie on different sides of the trajectory and hence the segment is crossed. The function
// returns the s position of the first intersection.
double get_s_of_first_crossing_with_polygon(double part_x, double part_tan, double* poly_s, \
                                            double* poly_x, int8_t num_polys, int8_t is_closed){
    double s = 1.e21;
    double trajectory_p0 = poly_x[0] - part_x - poly_s[0]*part_tan;
    double trajectory_p1 = trajectory_p0;
    double trajectory_p2, new_s;

    for (int8_t i=0; i<num_polys-1; i++){
        trajectory_p2 = poly_x[i+1] - part_x - poly_s[i+1]*part_tan;
        if (trajectory_p1*trajectory_p2 <= 0){
            // It's a hit
            new_s = _get_s_of_crossing_with_segment(part_x, part_tan, poly_s[i], poly_s[i+1], \
                                                    poly_x[i], poly_x[i+1]);
            if (new_s < s){
                s = new_s;
            }
        }
        trajectory_p1 = trajectory_p2;
    }
    // Last segment
    if (is_closed && trajectory_p1*trajectory_p0 <= 0){
        // It's a hit
        new_s = _get_s_of_crossing_with_segment(part_x, part_tan, poly_s[num_polys-1], poly_s[0], \
                                                poly_x[num_polys-1], poly_x[0]);
        if (new_s < s){
            s = new_s;
        }
    }
    return s;
}


// This function works as above, but considers a non-closed polygon, and adds half-open segments at
// the start and end (representing the half-infinite front and back faces of a jaw). For this
// reason, the function needs to know whether this is a positive or negative jaw (representing to
// which infinite it points), and the overall tilt of the jaw.
double get_s_of_first_crossing_with_open_polygon(double part_x, double part_tan, double* poly_s, \
                                                 double* poly_x, int8_t num_polys, double tan_tilt, \
                                                 int8_t side){
    double s = 1.e21;
    double new_s;
    double inf_poly_s[2];
    double inf_poly_x[2];

    // First segment is half-infinite; this implies one of the segment points lies at +-inf
    // In practice we just add a polygon point at the wall overflow (1km).
    // Important: we assume that this segment is perpendicular to the tilt (hence we use -1/tan_tilt)
    inf_poly_x[0] = 1.e3*side;
    if (fabs(tan_tilt) < 1.e-12){
        inf_poly_s[0] = poly_s[0];
    } else {
        inf_poly_s[0] = -(inf_poly_x[0] - poly_x[0] - poly_s[0]/tan_tilt)*tan_tilt;
    }
    inf_poly_s[1] = poly_s[0];
    inf_poly_x[1] = poly_x[0];
    new_s = get_s_of_first_crossing_with_polygon(part_x, part_tan, inf_poly_s, inf_poly_x, 2, 0);
    if (new_s < s){
        s = new_s;
    }

    // Middle segments
    new_s = get_s_of_first_crossing_with_polygon(part_x, part_tan, poly_s, poly_x, num_polys, 0);
    if (new_s < s){
        s = new_s;
    }

    // Last segment
    if (fabs(tan_tilt) < 1.e-12){
        inf_poly_s[0] = poly_s[num_polys-1];
    } else {
        inf_poly_s[0] = -(inf_poly_x[0] - poly_x[num_polys-1] - poly_s[num_polys-1]/tan_tilt)*tan_tilt;
    }
    inf_poly_s[1] = poly_s[num_polys-1];
    inf_poly_x[1] = poly_x[num_polys-1];
    new_s = get_s_of_first_crossing_with_polygon(part_x, part_tan, inf_poly_s, inf_poly_x, 2, 0);
    if (new_s < s){
        s = new_s;
    }
    return s;
}

#endif /* XCOLL_GEOM_POLY_H */
