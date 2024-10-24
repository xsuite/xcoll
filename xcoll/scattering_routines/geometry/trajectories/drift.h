// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_COLL_GEOM_DRIFT_H
#define XCOLL_COLL_GEOM_DRIFT_H


/*gpufun*/
int8_t DriftTrajectory_vlimit(double* restrict_s, double s0, double y0, double ym, double ymin, double ymax){
    if (fabs(ym) < XC_EPSILON){
        // Trajectory parallel to s axis
        if (y0 < ymin || y0 > ymax){
            return 0;  // Completely outside - no crossing possible
        } else {
            return -1; // Completely inside - no vertical check needed
        }
    } else {
        restrict_s[0] = (ymin - y0)/ym + s0;
        restrict_s[1] = (ymax - y0)/ym + s0;
        SWAP(restrict_s, 0, 1);   // To make sure these are sorted
        return 1;  // Default behavior: check overlap with horizontal crossings
    }
}


// Curve length is int_s1^s2 sqrt(1 + (dx/ds)^2 + (dy/ds)^2) ds

/*gpufun*/
double DriftTrajectory_length(double s0, double x0, double xm, double y0, double ym, double s1, double s2){
    (void) s0;  // Avoid unused parameter warning
    (void) x0;  // Avoid unused parameter warning
    (void) y0;  // Avoid unused parameter warning
    return (s2-s1)*sqrt(1 + xm*xm + ym*ym);
}


#endif /* XCOLL_COLL_GEOM_DRIFT_H */