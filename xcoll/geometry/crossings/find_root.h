// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#include <math.h>
#include <stdio.h>

/*gpufun*/
void LocalCrossing_func(LocalSegment seg, LocalTrajectory traj, double TS[2], double t, double l){
    // Here we get the expr from segment and traj, and we connect them to create TSs and TSx
    // TSs = Ts(l) - Ss(t)
    // TSx = Tx(l) - Sx(t)
    TS[0] = LocalTrajectory_func_s(traj, l) - LocalSegment_func_s(seg, t);
    TS[1] = LocalTrajectory_func_x(traj, l) - LocalSegment_func_x(seg, t);
}
/*gpufun*/
int8_t LocalCrossing_inv_J(LocalSegment seg, LocalTrajectory traj, double J_inv[2][2], double t, double l){
    double J[2][2];
    // get derivatives
    J[0][0] = LocalTrajectory_deriv_s(traj, l); // get deriv. dsT/dl LocalTrajectory_deriv_s
    J[0][1] = -LocalSegment_deriv_s(seg, t);     // get deriv. dsS/dt LocalSegment_deriv_s
    J[1][0] = LocalTrajectory_deriv_x(traj, l); // get deriv. dxT/dl LocalTrajectory_deriv_x
    J[1][1] = -LocalSegment_deriv_x(seg, t);     // get deriv. dxS/dt LocalSegment_deriv_x
    double det = J[0][0] * J[1][1] - J[1][0]*J[0][1];
    if (fabs(det) < XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL){
        return 0; // Failure, Jacobian is singular
    }
    J_inv[0][0] =  J[1][1] / det;
    J_inv[0][1] = -J[0][1] / det;
    J_inv[1][0] = -J[1][0] / det;
    J_inv[1][1] =  J[0][0] / det;
    return 1; // Success
}

/*gpufun*/
void FindRoot_newton(FindRoot finder, LocalSegment seg, LocalTrajectory traj, double guess_t, double guess_l){
    // success = 0 means success, success = 1 means failure
    double J_inv[2][2];
    double TS[2];
    // double success[1] = {0};
    double corr0, corr1;
    double new_t, new_l;
    int8_t num_found = FindRoot_get_num_solutions(finder);
    double* solution_t = FindRoot_getp_solution_t(finder);
    double* solution_l = FindRoot_getp_solution_l(finder);
    int8_t* converged = FindRoot_getp_converged(finder);
    for (int i = 0; i < num_found; i++){
        guess_t = solution_t[i];
        guess_l = solution_l[i];
        for (int i = 0; i < XC_GEOM_ROOT_NEWTON_MAX_ITER; i++){
            LocalCrossing_func(seg, traj, TS, guess_t, guess_l);
            converged[i] = LocalCrossing_inv_J(seg, traj, J_inv, converged[i], guess_t, guess_l);
            if (converged[i] != 1){
                printf("Jacobian is singular. Stopping Newton's method.\n");
                solution_t[i] = 1.e21;
                solution_l[i] = 1.e21;
                continue;
            }
            // Residual convergence check
            if (fabs(TS[0]) < XC_GEOM_ROOT_NEWTON_EPSILON && fabs(TS[1]) < XC_GEOM_ROOT_NEWTON_EPSILON) {
                solution_t[i] = guess_t;
                solution_l[i] = guess_l;
                continue;
            }
            corr0 = J_inv[0][0]*TS[0] + J_inv[0][1]*TS[1]; // delta for l
            corr1 = J_inv[1][0]*TS[0] + J_inv[1][1]*TS[1]; // delta for t
            new_t = guess_t - corr1;
            new_l = guess_l - corr0;
            // Check for parameter update convergence
            if ((fabs(new_t -  guess_t) < XC_GEOM_ROOT_NEWTON_EPSILON) && (fabs(new_l - guess_l) < XC_GEOM_ROOT_NEWTON_EPSILON)){
                solution_t[i] = guess_t;
                solution_l[i] = guess_l;
                continue;
            }
            // Update the guesses for the next iteration
            guess_t = new_t;  // Keep *t updated
            guess_l = new_l;  // Keep *l updated
        }
        printf("Warning: Newton's method did not converge within the maximum number of iterations (%d).\n", XC_GEOM_ROOT_NEWTON_MAX_ITER);
        solution_t[i] = 1.e21;
        solution_l[i] = 1.e21;
        converged[i] = 0;
    }
}


XC_SLICING_NUM_STEPS = 4;
XC_SLICING_MAX_NEST_LEVEL = 8;


/*gpufun*/
void slice_before_newton(FindRoot finder, LocalSegment seg, LocalTrajectory traj,
                         double t1, double t2, double l1, double l2, int8_t nest_level){
    // Prepare initial guess for Newton-Raphson root finding
    box1 = LocalSegment_getp_box(seg);
    box2 = LocalTrajectory_getp_box(traj);
    t_step = (t2 - t1) / XC_SLICING_NUM_STEPS;
    l_step = (l2 - l1) / XC_SLICING_NUM_STEPS;
    double* solution_t = FindRoot_getp_solution_t(finder);
    double* solution_l = FindRoot_getp_solution_l(finder);
    double t, l;
    int16_t num;
    for (int i = 0; i < XC_SLICING_NUM_STEPS; i++){
        t = t1 + i * t_step;
        for (int j = 0; j < XC_SLICING_NUM_STEPS; j++){
            l = l1 + j * l_step;
            LocalSegment_update_box(seg, t, t + t_step);
            LocalTrajectory_update_box(traj, l, l + l_step);
            if (BoundingBox_overlaps(box1, box2)){
                if (nest_level >= XC_SLICING_MAX_NEST_LEVEL - 1){
                    // We reached the maximum nesting level, return the midpoint of the current t-interval
                    if (*num_found >= FindRoot_get_max_solutions(finder)){
                        printf("Warning: Maximum number of solutions (%d) reached in slice_before_newton. Some solutions may be missed.\n", _MAX_SOLUTIONS);
                        fflush(stdout);
                        return;
                    }
                    num = FindRoot_get_num_solutions(finder);
                    solution_t[num] = t + 0.5 * t_step;
                    solution_l[num] = l + 0.5 * l_step;
                    FindRoot_set_num_solutions(finder, num + 1);
                } else {
                    slice_before_newton(seg, traj, t, t + t_step, l, l + l_step, nest_level + 1);
                }
            }
        }
    }
    return;
}

/*gpufun*/
void find_crossing_approximate(FindRoot finder, LocalSegment seg, LocalTrajectory traj){
    slice_before_newton(finder, seg, traj, 0, 1, 0, 1, 0);
    int8_t num_found = FindRoot_get_num_solutions(finder);
    for (int i = 0; i < num_found; i++){
        guess_t = FindRoot_get_solution_t(finder, i);
        guess_l = FindRoot_get_solution_l(finder, i);
        find_root_newton(finder, traj, seg, guess_l, guess_t);
    }
}


void FindRoot_find_crossing(FindRoot finder, LocalTrajectory traj, LocalSegment seg){ //, double guess_l, double guess_t){
// this will act as the crossing main function - C magic
// First we check for the analytical solutions by checking what trajectory we have
    switch (LocalTrajectory_typeid(traj)){
        case LocalTrajectory_DriftTrajectory_t:
            switch (LocalSegment_typeid(seg)){
                // // may have to add #ifndef LOCALSEGMENT_SKIP_LINESEGMENT
                // case LocalSegment_LineSegment_t:
                //     return crossing_drift_line(traj, seg);
                //     break;
                // case LocalSegment_HalfOpenLineSegment_t:
                //     return crossing_drift_halfline(traj, seg);
                //     break;
                // // case LocalSegment_CircularSegment_t:
                // //     return crossing_drift_circular(traj, seg);
                // //     break;
                // case LocalSegment_BezierSegment_t:
                //     return crossing_drift_bezier(traj, seg);
                //     break;
                default:
                    // Custom segment
                    return find_crossing_approximate(finder, seg, traj);
            }
            break;
        case LocalTrajectory_CircularTrajectory_t:
            switch (LocalSegment_typeid(seg)){
                // case LocalSegment_LineSegment_t:
                //     return crossing_circ_line(traj, seg);
                //     break;
                // case LocalSegment_HalfOpenLineSegment_t:
                //     return crossing_circ_halfline(traj, seg);
                //     break;
                // // case LocalSegment_CircularSegment_t:
                // //     return crossing_circ_circular(traj, seg);
                // //     break;
                // case LocalSegment_BezierSegment_t:
                //     return crossing_circ_bezier(traj, seg);
                //     break;
                default:
                    // Custom segment
                    return find_crossing_approximate(finder, seg, traj);
            }
            break;
    default:
        return find_crossing_approximate(finder, seg, traj);
    }
} 
#endif /* XCOLL_GEOM_FIND_ROOT_H */