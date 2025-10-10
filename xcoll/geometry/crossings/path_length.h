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
void simpson(FindRoot finder, LocalTrajectory traj, int subintervals) {
    double l1 = 0;
    double l2 = FindRoot_get_solution_l(finder, 0); ///// !!!!!!!!
    if (subintervals % 2 != 1) {
        subintervals++;                                         // Requires an even number of subintervals
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
    sum *= step_size / 3;
    FindRoot_set_path_length(finder, sum);
    return;
}

/*gpufun*/
void find_approximate_path_length(FindRoot finder, LocalTrajectory traj){
    
}

/*gpufun*/
void find_path_length_analytic(FindRoot finder, LocalTrajectory traj){
    
    // only for drift for now
}
/*gpufun*/
void FindRoot_find_path_length(FindRoot finder, LocalTrajectory traj){
    printf("Finding path length...\n");
    // TODO: we should find a find to classify if we are actually INSIDE or NOT. Drift outside dont need path length.
    switch (LocalTrajectory_typeid(traj)){
        case LocalTrajectory_DriftTrajectory_t:
            double s0  = DriftTrajectory_get_s0((DriftTrajectory) LocalTrajectory_member(traj));
            double x0  = DriftTrajectory_get_x0((DriftTrajectory) LocalTrajectory_member(traj));
            // TODO: Check with Frederik. I assume we always use the first solution as the second solution might not
            // even be relevant (say out - in again, then we actually never go in again. IN-OUT -> we change to mcs before out)
            double l2  = FindRoot_get_solution_l(finder, 0);
            find_path_length_analytic();
            break;
        default:
            printf("Using approximate crossing method.\n");
            return find_approximate_path_length(finder, traj);
        }
}
#endif /* XCOLL_GEOM_SIMPSON_H */