// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_SIMPSON_H
#define XCOLL_GEOM_SIMPSON_H
#include <stdio.h>
#include <math.h>
// Simpson's Rule function for numerical integration
// (b - a) / (3 * subintervals) * (f(l1) + 4f(subinterval1) + 2f(even subinterval) + 4f(odd subinterval)) + ... + f(l2))
/*gpufun*/
void simpson(FindRoot finder, LocalTrajectory traj, int32_t subintervals) {
    double l1 = 0;
    double l2 = FindRoot_get_solution_l(finder, 0);
    if (subintervals % 2 != 1) {
        subintervals++;
    }
    double step_size = (l2 - l1) / subintervals;
    double sum = sqrt(1 + LocalTrajectory_deriv_x(traj, l1)*LocalTrajectory_deriv_x(traj, l1)) + 
                 sqrt(1 + LocalTrajectory_deriv_x(traj, l2)*LocalTrajectory_deriv_x(traj, l2));  // f(l1) + f(l2)
    // Add subintervals to the sum
    for (int i = 1; i < subintervals; i++) {
        double l = l1 + i * step_size;
        // Even indices (except the endpoints) are multiplied by 2. Odd indices are multiplied by 4
        if (i % 2 == 0) {
            sum += 2.0 * sqrt(1 + LocalTrajectory_deriv_x(traj, l)*LocalTrajectory_deriv_x(traj, l));
        } else {
            sum += 4.0 * sqrt(1 + LocalTrajectory_deriv_x(traj, l)*LocalTrajectory_deriv_x(traj, l));
        }
    }
    sum *= step_size / 3.;
    if (sum < 0) {
        printf("Warning: Computed path length is negative. Setting to 1e21.\n");
        fflush(stdout);
        FindRoot_set_path_length(finder, 1e21);
        return;
    }
    FindRoot_set_path_length(finder, sum);
    return;
}

/*gpufun*/
void find_path_length_analytic(FindRoot finder, LocalTrajectory traj){
    // analytical drift does not save l2, so we need to compute it here
    double l2  = FindRoot_get_solution_l(finder, 0); // THIS DOES NOT EXIST FOR ANALYTICAL SOLUTIONS; EITHER FIX OR WE NEVER NEED THIS
    double s0  = LocalTrajectory_func_s(traj, 0);  // l1 = 0
    double x0  = LocalTrajectory_func_x(traj, 0);  
    double s1  = LocalTrajectory_func_s(traj, l2); // l2 = solution for l
    double x1  = LocalTrajectory_func_x(traj, l2);
    double path_length = sqrt((s1 - s0)*(s1 - s0) + (x1 - x0)*(x1 - x0));
    FindRoot_set_path_length(finder, path_length);
    return;
    // only for drift for now
}
/*gpufun*/
void FindRoot_find_path_length(FindRoot finder, LocalSegment seg, LocalTrajectory traj){
    if (LocalSegment_typeid(seg) == LocalSegment_BezierSegment_t){
        return simpson(finder, traj, XC_GEOM_SIMPSON_SUBINTERVALS);
    }
    // TODO: Check with Frederik. I assume we always use the first solution as the second solution might not
    // even be relevant (say out - in again, then we actually never go in again. IN-OUT -> we change to mcs before out)
    if (LocalTrajectory_typeid(traj) == LocalTrajectory_DriftTrajectory_t){
        return find_path_length_analytic(finder, traj);
    } else {
        return simpson(finder, traj, XC_GEOM_SIMPSON_SUBINTERVALS);
    }
}
#endif /* XCOLL_GEOM_SIMPSON_H */