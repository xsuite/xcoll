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
    UNUSED(l);
    return DriftTrajectory_get_tan_t0(traj);
}

/*gpufun*/
double DriftTrajectory_deriv_s(DriftTrajectory traj, double l){
    UNUSED(l);
    return DriftTrajectory_get_cos_t0(traj);
}

/*gpufun*/
double DriftTrajectory_deriv_x(DriftTrajectory traj, double l){
    UNUSED(l);
    return DriftTrajectory_get_sin_t0(traj);
}

/*gpufun*/
void DriftTrajectory_init_bounding_box(DriftTrajectory traj, BoundingBox box, double l1, double l2){
    double s1 = DriftTrajectory_get_s0(traj);
    double s2 = DriftTrajectory_func_s(traj, l2);
    double x1 = DriftTrajectory_get_x0(traj);
    double x2 = DriftTrajectory_func_x(traj, l2);
    double sin_t0 = DriftTrajectory_get_sin_t0(traj);
    double cos_t0 = DriftTrajectory_get_cos_t0(traj);
    double sin_p, cos_p;
    if (sin_t0 < 0){   // if theta is larger than 180 degrees, theta = theta - 180
        sin_t0 = -sin_t0;
        cos_t0 = -cos_t0;
    }
    if (cos_t0 < 1){   // if theta is larger than 90 degrees, phi = theta + 90 
        sin_p = cos_t0;
        cos_p = -sin_t0;
    } else {          // if theta is between 0 and 90 degrees, phi = theta - 90
        sin_p = -cos_t0;
        cos_p = sin_t0;
    }
    BoundingBox_set_l(box, sqrt((s2 - s1)*(s2 - s1) + (x2 - x1)*(x2 - x1)));   // length of the box
    BoundingBox_set_w(box, BoundingBox_get_l(box)/3.);       // width of the box 
    BoundingBox_set_rC(box, sqrt( (s1+BoundingBox_get_w(box)/2.*cos_p) * (s1+BoundingBox_get_w(box)/2.*cos_p) +  // length of the position vector to the first vertex
                                   (x1+BoundingBox_get_w(box)/2.*sin_p) * (x1+BoundingBox_get_w(box)/2.*sin_p) ));
    BoundingBox_set_sin_tb(box, sin_t0);  // orientation of the box (angle of length wrt horizontal)
    BoundingBox_set_cos_tb(box, cos_t0);
    BoundingBox_set_sin_tC(box, x1 / BoundingBox_get_rC(box));  // angle of the position vector to the first vertex
    BoundingBox_set_cos_tC(box, s1 / BoundingBox_get_rC(box));
    double sin_tC = BoundingBox_get_sin_tC(box);
    double cos_tC = BoundingBox_get_cos_tC(box);
    BoundingBox_set_proj_l(box, BoundingBox_get_rC(box) * (cos_t0*cos_tC + sin_t0*sin_tC)); // projection of the position vector on length: rC * (cos_t*cos_tC + sin_t*sin_tC)
    BoundingBox_set_proj_w(box, BoundingBox_get_rC(box) * (cos_t0*sin_tC - sin_t0*cos_tC)); // projection of position vector on width: rC * (cos_t*sin_tC - sin_t*cos_tC)
}



/*gpufun*/
int8_t DriftTrajectory_vlimit(double* restrict_s, double s0, double y0, double ym, double ymin, double ymax){
    if (fabs(ym) < XC_GEOM_EPSILON){
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
    return XC_GEOM_S_MAX;
}

/*gpufun*/
double Trajectory_get_before_s(int8_t n_hit, double* s, double before_s){
    for (int8_t i=n_hit-1; i>=0; i--){
        if (s[i] <= before_s){
            return s[i];
        }
    }
    return XC_GEOM_S_MAX;
}

/*gpufun*/
double Trajectory_get_after_s(int8_t n_hit, double* s, double after_s){
    for (int8_t i=0; i<n_hit; i++){
        if (s[i] >= after_s){
            return s[i];
        }
    }
    return XC_GEOM_S_MAX;
}

/*gpufun*/
double Trajectory_get_last(int8_t n_hit, double* s){
    if (n_hit>0){
        return s[n_hit-1];
    }
    return XC_GEOM_S_MAX;
}

#endif /* XCOLL_GEOM_TRAJ_DRIFT_H */