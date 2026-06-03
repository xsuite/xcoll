// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_OBJECTS_H
#define XCOLL_GEOM_OBJECTS_H

#ifdef XO_CONTEXT_CPU
#include <math.h>
#include <stdio.h>
#include <stdint.h>  // for int64_t etc
#include <stdlib.h>  // for malloc and free
#endif  // XO_CONTEXT_CPU


// Assumption for all objects: the particle at -inf is outside the object (otherwise some comparisons might give wrong results)
// Segments are stored BY VALUE in a fixed-size array (Segment segments[N]); the create_*
// functions fill a caller-provided array (no malloc), and there is nothing to destroy.


// Collimator jaw  (CPU-only: collimators are not a GPU-track target in this PR,
// and they use heap-allocated variable-length segment arrays. The crystal path
// below uses the GPU-safe by-value fill.)
// --------------
#ifdef XO_CONTEXT_CPU

/*gpufun*/
Segment* create_jaw(double s_U, double x_U, double s_D, double x_D, double tilt_tan, int8_t side){
    Segment* segments = (Segment*) malloc(3*sizeof(Segment));
    segments[0] = create_halfopen_line_segment(s_U, x_U, tilt_tan, side);
    segments[1] = create_line_segment(s_U, x_U, s_D, x_D);
    segments[2] = create_halfopen_line_segment(s_D, x_D, tilt_tan, side);
    return segments;
}

/*gpufun*/
void destroy_jaw(Segment* segments){
    free(segments);
}


// Polygon
// -------

/*gpufun*/
Segment* create_polygon(double* s_poly, double* x_poly, int8_t num_polys){
    Segment* segments = (Segment*) malloc((unsigned int) num_polys*sizeof(Segment));
    for (int8_t i=0; i<num_polys-1; i++){
        segments[i] = create_line_segment(s_poly[i], x_poly[i], s_poly[i+1], x_poly[i+1]);
    }
    segments[num_polys-1] = create_line_segment(s_poly[num_polys-1], x_poly[num_polys-1], \
                                                s_poly[0], x_poly[0]);
    return segments;
}

/*gpufun*/
void destroy_polygon(Segment* segments, int8_t num_polys){
    free(segments);
}


// Open polygon
// ------------

/*gpufun*/
Segment* create_open_polygon(double* s_poly, double* x_poly, int8_t num_polys, double tilt_tan, int8_t side){
    Segment* segments = (Segment*) malloc((num_polys+1)*sizeof(Segment));
    segments[0] = create_halfopen_line_segment(s_poly[0], x_poly[0], tilt_tan, side);
    for (int8_t i=1; i<num_polys; i++){
        segments[i] = create_line_segment(s_poly[i-1], x_poly[i-1], s_poly[i], x_poly[i]);
    }
    segments[num_polys] = create_halfopen_line_segment(s_poly[num_polys-1], x_poly[num_polys-1], \
                                                       tilt_tan, side);
    return segments;
}

/*gpufun*/
void destroy_open_polygon(Segment* segments, int8_t num_polys){
    free(segments);
}

#endif  // XO_CONTEXT_CPU


// Crystal
// -------

// The four corners A, B, C, D are such that AB is the front face, BC the curve furthest from the beam,
// CD the back face, and DA the curve closest to the beam.
// Fills the caller-provided segments array (length XC_MAX_SEGMENTS == 4) by value; no
// allocation, so there is no matching destroy function. The straight-crystal case (R==0) is
// not implemented: a straight crystal is a per-element property, so instead of an (unsupported
// on GPU) fflush + bare return that leaves the segments uninitialised, the whole bunch is
// killed with the standard Xcoll error state -- mirroring EverestCrystal_init_geometry's own
// kill_all_particles(part0, XC_ERR_INVALID_XOFIELD) handling of the side==0 case.
/*gpufun*/
void create_crystal(LocalParticle* part0, Segment* segments, double R, double width, double length, double jaw_U, double tilt_sin, double tilt_cos){
    // First corner is what defines the crystal position
    double A_s = 0;
    double A_x = jaw_U;

    // Manipulate R in function of sign
    double sgnR = (R > 0) - (R < 0);
    double R_short  = sgnR*(fabs(R) - width);
    double sin_a = length/fabs(R);
    double cos_a = sqrt(1 - length*length/R/R);
    if (fabs(R) < 1.e-12){
        // straight crystal - not yet implemented
#ifdef XO_CONTEXT_CPU
        printf("Straight crystal not yet implemented!");
#endif  // XO_CONTEXT_CPU
        kill_all_particles(part0, XC_ERR_NOT_IMPLEMENTED);
        return;

    } else if (R < 0){
        // This distinction is needed to keep the crystal at the same location when changing the bend direction
        double R_temp = R_short;
        R_short = R;
        R = R_temp;
    }

    // Bending centre is defined w.r.t. A
    double R_s = A_s - R*tilt_sin;
    double R_x = A_x + R*tilt_cos;

    // Three remaining corner points of crystal
    double B_s = R_s + R_short*tilt_sin;
    double C_s = R_s + fabs(R_short)*sin_a*tilt_cos + R_short*cos_a*tilt_sin;
    double D_s = R_s + fabs(R)*sin_a*tilt_cos + R*cos_a*tilt_sin;
    double B_x = R_x - R_short*tilt_cos;
    double C_x = R_x - cos_a*tilt_cos*R_short + sin_a*tilt_sin*fabs(R_short);
    double D_x = R_x - cos_a*tilt_cos*R + sin_a*tilt_sin*fabs(R);
    double A_t = atan2(A_x - R_x, A_s - R_s);
    double D_t = atan2(D_x - R_x, D_s - R_s);
    double t1 = MIN(A_t, D_t);
    double t2 = MAX(A_t, D_t);

    // Fill segments (by value, into the caller-provided array)
    segments[0] = create_line_segment(A_s, A_x, B_s, B_x);
    segments[1] = create_circular_segment(R, R_s, R_x, t1, t2);
    segments[2] = create_line_segment(C_s, C_x, D_s, D_x);
    segments[3] = create_circular_segment(R_short, R_s, R_x, t1, t2);

    // printf("R: (%f, %f)   A: (%f, %f)   B: (%f, %f)   C: (%f, %f)   D: (%f, %f)   t1: %f   t2: %f\n", R_s,R_x,A_s,A_x,B_s,B_x,C_s,C_x,D_s,D_x,t1*180/3.141592653589793,t2*180/3.141592653589793); fflush(stdout);
}


#endif /* XCOLL_GEOM_OBJECTS_H */
