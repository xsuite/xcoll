// copyright ############################### #
// This file is part of the Xcoll Package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_CRY_H
#define XCOLL_GEOM_CRY_H
#include <math.h>
#include <stdio.h>

#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

// This function compares a particle trajectory to a crystal layout. The four corners A, B, C, D are such that
// AB is the front face, BC the curve closest to the beam, CD the back face, and DA the curve furtest from the beam.
double get_s_of_first_crossing_with_crystal(double part_x, double part_tan, double R, double width, \
                                            double length, double jaw_U, double sin_tilt, double cos_tilt){
    double s = 1.e21;
    double new_s, x, t;
    double R_s, R_x;
    double pt_s[2];
    double pt_x[2];
    double sgnR = (R > 0) - (R < 0);
    double R_short  = sgnR*(fabs(R) - width);
    double sin_a = length/fabs(R);
    double cos_a = sqrt(1 - length*length/R/R);

    // First corner is what defines the crystal position
    double A_s = length/2*(1 - cos_tilt);
    double A_x = jaw_U;

    // Bending centre is defined w.r.t. this corner
    if (fabs(R) < 1.e-12){
        // straight crystal - not yet implemented 
        printf("Not yet implemented!"); //only_for_context cpu_serial
        fflush(stdout);                 //only_for_context cpu_serial
        return 0;

    } else if (R < 0){
        // This distinction is needed to keep the crystal at the same location when changing the bend direction
        double R_temp = R_short;
        R_short = R;
        R = R_temp;
    }
    R_s = A_s - R*sin_tilt;
    R_x = A_x + R*cos_tilt;

    // Three remaining corner points of crystal
    double B_s = R_s + R_short*sin_tilt;
    double C_s = R_s + fabs(R_short)*sin_a*cos_tilt + R_short*cos_a*sin_tilt;
    double D_s = R_s + fabs(R)*sin_a*cos_tilt + R*cos_a*sin_tilt;
    double B_x = R_x - R_short*cos_tilt;
    double C_x = R_x - cos_a*cos_tilt*R_short + sin_a*sin_tilt*fabs(R_short);
    double D_x = R_x - cos_a*cos_tilt*R + sin_a*sin_tilt*fabs(R);
    // printf("As: %f\n", A_s);
    // printf("Bs: %f\n", B_s);
    // printf("Cs: %f\n", C_s);
    // printf("Ds: %f\n", D_s);
    // printf("Ax: %f\n", A_x);
    // printf("Bx: %f\n", B_x);
    // printf("Cx: %f\n", C_x);
    // printf("Dx: %f\n", D_x);
    double A_t = atan2(A_x - R_x, A_s - R_s);
    double D_t = atan2(D_x - R_x, D_s - R_s);
    double t1 = MIN(A_t, D_t);
    double t2 = MAX(A_t, D_t);
    // printf("t1=%f  t2=%f\n", t1*180/3.1415, t2*180/3.1415);

    // Front segment
    pt_s[0] = A_s;
    pt_x[0] = A_x;
    pt_s[1] = B_s;
    pt_x[1] = B_x;
    new_s = get_s_of_first_crossing_with_polygon(part_x, part_tan, pt_s, pt_x, 2, 0);
    if (new_s < s){
        // printf("Hit front\n");
        s = new_s;
    }

    // Back segment
    pt_s[0] = C_s;
    pt_x[0] = C_x;
    pt_s[1] = D_s;
    pt_x[1] = D_x;
    new_s = get_s_of_first_crossing_with_polygon(part_x, part_tan, pt_s, pt_x, 2, 0);
    if (new_s < s){
        // printf("Hit back\n");
        s = new_s;
    }

    // Lower curve
    double a = 1 + part_tan*part_tan;
    double bb = R_s - part_tan*(part_x - R_x);
    double c = R_s*R_s + (part_x - R_x)*(part_x - R_x) - R*R;
    double disc = bb*bb - a*c;
    if (disc >= 0){
        new_s = 1/a*(bb + sqrt(bb*bb - a*c));  // positive solution
        x = part_x + new_s*part_tan;
        t = atan2(x - R_x, new_s - R_s);
        if (t1 <= t && t <= t2){
            // printf("Hit lower curve pos: s=%f  A_s=%f  D_s=%f  |  x=%f  A_x=%f  D_x=%f\n", new_s, A_s, D_s, x, A_x, D_x);
            if (new_s < s){
                s = new_s;
            }
        }
        new_s = 1/a*(bb - sqrt(bb*bb - a*c));  // negative solution
        x = part_x + new_s*part_tan;
        t = atan2(x - R_x, new_s - R_s);
        if (t1 <= t && t <= t2){
            // printf("Hit lower curve neg: s=%f  A_s=%f  D_s=%f  |  x=%f  A_x=%f  D_x=%f\n", new_s, A_s, D_s, x, A_x, D_x);
            if (new_s < s){
                s = new_s;
            }
        }
    }

    // Upper curve
    c = R_s*R_s + (part_x - R_x)*(part_x - R_x) - R_short*R_short;
    disc = bb*bb - a*c;
    if (disc >= 0){
        new_s = 1/a*(bb + sqrt(disc));  // positive solution
        x = part_x + new_s*part_tan;
        t = atan2(x - R_x, new_s - R_s);
        if (t1 <= t && t <= t2){
            // printf("Hit upper curve pos: s=%f  A_s=%f  D_s=%f  |  x=%f  A_x=%f  D_x=%f\n", new_s, A_s, D_s, x, A_x, D_x);
            if (new_s < s){
                s = new_s;
            }
        }
        new_s = 1/a*(bb - sqrt(disc));  // negative solution
        x = part_x + new_s*part_tan;
        t = atan2(x - R_x, new_s - R_s);
        if (t1 <= t && t <= t2){
            // printf("Hit upper curve neg: s=%f  A_s=%f  D_s=%f  |  x=%f  A_x=%f  D_x=%f\n", new_s, A_s, D_s, x, A_x, D_x);
            if (new_s < s){
                s = new_s;
            }
        }
    }
    // fflush(stdout);

    return s;
}


#endif /* XCOLL_GEOM_CRY_H */