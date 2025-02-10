// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_TRAJ_DRIFT_H
#define XCOLL_GEOM_TRAJ_DRIFT_H

#include <stdio.h>
#include <math.h>


// /*gpufun*/
// void DriftTrajectory_set_params(DriftTrajectory traj, LocalParticle part){
//     DriftTrajectory_set_s0(traj, LocalParticle_get_s(part));
//     DriftTrajectory_set_x0(traj, LocalParticle_get_x(part));
//     double xp = LocalParticle_get_exact_xp(part);
//     DriftTrajectory_set_sin_t0(traj, xp / sqrt(1+xp*xp));
//     DriftTrajectory_set_cos_t0(traj, 1 / sqrt(1+xp*xp));
//     DriftTrajectory_set_tan_t0(traj, xp);
// }


/*gpufun*/
void DriftTrajectory_set_params(DriftTrajectory traj, double s0, double x0, double xp){
    DriftTrajectory_set_s0(traj, s0);
    DriftTrajectory_set_x0(traj, x0);
    DriftTrajectory_set_sin_t0(traj, xp / sqrt(1+xp*xp));
    DriftTrajectory_set_cos_t0(traj, 1 / sqrt(1+xp*xp));
    DriftTrajectory_set_tan_t0(traj, xp);
}

/*gpufun*/
double DriftTrajectory_func_s(DriftTrajectory traj, double l){
    double s0 = DriftTrajectory_get_s0(traj);
    double cos_t0 = DriftTrajectory_get_cos_t0(traj);
    return s0 + l*cos_t0;
}

/*gpufun*/
double DriftTrajectory_func_x(DriftTrajectory traj, double l){
    double x0 = DriftTrajectory_get_x0(traj);
    double sin_t0 = DriftTrajectory_get_sin_t0(traj);
    return x0 + l*sin_t0;
}

/*gpufun*/
double DriftTrajectory_func_xp(DriftTrajectory traj, double l){
    return DriftTrajectory_get_tan_t0(traj);
}

/*gpufun*/
double DriftTrajectory_deriv_s(DriftTrajectory traj, double l){
    return DriftTrajectory_get_cos_t0(traj);
}

/*gpufun*/
double DriftTrajectory_deriv_x(DriftTrajectory traj, double l){
    return DriftTrajectory_get_sin_t0(traj);
}


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


// The following functions do not need to be redefined for the other trajectories

/*gpufun*/
double Trajectory_get_first(int8_t n_hit, double* s){
    if (n_hit>0){
        return s[0];
    }
    return XC_S_MAX;
}

/*gpufun*/
double Trajectory_get_before_s(int8_t n_hit, double* s, double before_s){
    for (int8_t i=n_hit-1; i>=0; i--){
        if (s[i] <= before_s){
            return s[i];
        }
    }
    return XC_S_MAX;
}

/*gpufun*/
double Trajectory_get_after_s(int8_t n_hit, double* s, double after_s){
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= after_s){
            return s[i];
        }
    }
    return XC_S_MAX;
}

/*gpufun*/
double Trajectory_get_last(int8_t n_hit, double* s){
    if (n_hit>0){
        return s[n_hit-1];
    }
    return XC_S_MAX;
}

#endif /* XCOLL_GEOM_TRAJ_DRIFT_H */