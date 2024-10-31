// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_CIRCULARSEG_H
#define XCOLL_COLL_GEOM_CIRCULARSEG_H

/*gpufun*/
void CircularSegment_crossing_drift(CircularSegment seg, int8_t* n_hit, double* s, double s0, double x0, double xm){
    // Get segment data
    double R  = CircularSegment_get_R(seg);
    double sC = CircularSegment_get_s(seg);
    double xC = CircularSegment_get_x(seg);
    double t1 = CircularSegment_get_t1(seg);
    double t2 = CircularSegment_get_t2(seg);
    // Move the angles to [-pi, pi]
    int8_t reversed = 0, full_circle = 0;
    if (fabs(fabs(t2 - t1) - 2*M_PI) < XC_EPSILON){
        full_circle = 1;
    }
    while (t1 < -M_PI){
        t1 += 2*M_PI;
    }
    while (t1 > M_PI){
        t1 -= 2*M_PI;
    }
    while (t2 < -M_PI){
        t2 += 2*M_PI;
    }
    while (t2 > M_PI){
        t2 -= 2*M_PI;
    }
    if (t2 < t1){
        reversed = 1;
    }
    // Calculate crossings
    double a = 1 + xm*xm;
    double bb = sC - xm*(x0 - xC - xm*s0); // This is -b/2 with b from the quadratic formula
    double c = sC*sC + (x0 - xC - xm*s0)*(x0 - xC - xm*s0) - R*R;
    double disc = bb*bb - a*c; // This is  2*discriminant**2
    if (disc < 0){
        // No crossing
        return;
    }
    for (int8_t i = 0; i < 2; i++) {
        double sgnD = i*2-1; // negative and positive solutions; if multiplicity 2, we add the same solution twice
        double new_s = (bb + sgnD*sqrt(fabs(disc)))/a;
        double new_x = x0 + (new_s - s0)*xm;
        double t = atan2(new_x - xC, new_s - sC);
        if (full_circle){
            // Full circle, so always hit
            s[*n_hit] = new_s;
            (*n_hit)++;
        } else if (reversed){
            // t2 < t1, so we are looking at the inverted region of angles
            if (t1 <= t || t <= t2){
                s[*n_hit] = new_s;
                (*n_hit)++;
            }
        } else {
            if (t1 <= t && t <= t2){
                s[*n_hit] = new_s;
                (*n_hit)++;
            }
        }
    }
}

#endif /* XCOLL_COLL_GEOM_CIRCULARSEG_H */