// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#include <math.h>
#include <stdio.h>

/*gpufun*/
void DriftTrajectory_analytical_solution_l(DriftTrajectory traj, double solution_l[3], int8_t num_solutions, double s, double x){
    if (DriftTrajectory_get_cos_t0(traj) - DriftTrajectory_get_cos_t0(traj) < XC_GEOM_EPSILON){
        solution_l[num_solutions] = (x - DriftTrajectory_get_x0(traj)) / DriftTrajectory_get_sin_t0(traj);
    } else {
        solution_l[num_solutions] = (s - DriftTrajectory_get_s0(traj)) / DriftTrajectory_get_cos_t0(traj);
    }
}

/*gpufun*/
void LocalCrossing_func(LocalSegment seg, LocalTrajectory traj, double TS[2], double t, double l){
    // Here we get the expr from segment and traj, and we connect them to create TSs and TSx
    // TSs = Ts(l) - Ss(t)
    // TSx = Tx(l) - Sx(t)
    //printf("Localcross: LocalSegment s: %f, x: %f\n", LocalSegment_func_s(seg, t), LocalSegment_func_x(seg, t));
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
void FindRoot_newton(FindRoot finder, LocalSegment seg, LocalTrajectory traj, double guess_t, double guess_l,
                     int8_t num){
    // success = 0 means success, success = 1 means failure
    double J_inv[2][2];
    double TS[2];
    double corr0, corr1;
    double new_t, new_l;
    int8_t converged;
    for (int i = 0; i < XC_GEOM_ROOT_NEWTON_MAX_ITER; i++){
        LocalCrossing_func(seg, traj, TS, guess_t, guess_l);
        //printf("Iteration %d: t = %f, l = %f, TSs = %f, TSx = %f\n", i, guess_t, guess_l, TS[0], TS[1]);
        converged = LocalCrossing_inv_J(seg, traj, J_inv, guess_t, guess_l);
        if (converged != 1){
            printf("Jacobian is singular. Stopping Newton's method.\n");
            FindRoot_set_solution_t(finder, num, 1.e21);
            FindRoot_set_solution_l(finder, num, 1.e21);
            FindRoot_set_converged(finder, num, converged);
            continue;
        }
        corr0 = J_inv[0][0]*TS[0] + J_inv[0][1]*TS[1]; // delta for l
        corr1 = J_inv[1][0]*TS[0] + J_inv[1][1]*TS[1]; // delta for t
        new_t = guess_t - corr1;
        new_l = guess_l - corr0;
        // For convergence study -----
        FindRoot_set_delta_t(finder, i, fabs(new_t - guess_t));
        FindRoot_set_delta_l(finder, i, fabs(new_l - guess_l));
        FindRoot_set_res_t(finder, i, fabs(TS[0]));
        FindRoot_set_res_l(finder, i, fabs(TS[1]));
        // ---------------------------
        // Residual convergence checks
        if (fabs(TS[0]) < XC_GEOM_ROOT_NEWTON_EPSILON && fabs(TS[1]) < XC_GEOM_ROOT_NEWTON_EPSILON) {
            FindRoot_set_solution_t(finder, num, guess_t);
            FindRoot_set_solution_l(finder, num, guess_l);
            FindRoot_set_converged(finder, num, converged);
            return;
        }
        // Check for parameter update convergence
        if ((fabs(new_t -  guess_t) < XC_GEOM_ROOT_NEWTON_EPSILON) && (fabs(new_l - guess_l) < XC_GEOM_ROOT_NEWTON_EPSILON)){
            FindRoot_set_solution_t(finder, num, guess_t);
            FindRoot_set_solution_l(finder, num, guess_l);
            FindRoot_set_converged(finder, num, converged);
            return;
        }
        // Update the guesses for the next iteration
        guess_t = new_t;  // Keep *t updated
        guess_l = new_l;  // Keep *l updated
    }
    printf("Warning: Newton's method did not converge within the maximum number of iterations (%d).\n", XC_GEOM_ROOT_NEWTON_MAX_ITER);
    FindRoot_set_solution_t(finder, num, 1.e21);
    FindRoot_set_solution_l(finder, num, 1.e21);
    FindRoot_set_num_solutions(finder, 0);
    FindRoot_set_converged(finder, num, converged);
}
// TODO: make this more clean pls
int8_t XC_SLICING_NUM_STEPS = 4;
int8_t XC_SLICING_MAX_NEST_LEVEL = 6;

/*gpufun*/
void slice_before_newton(FindRoot finder, LocalSegment seg, LocalTrajectory traj,
                         double t1, double t2, double l1, double l2, int8_t nest_level){
    // Prepare initial guesses for Newton-Raphson root finding
    BoundingBox_ _box_seg = {0};
    BoundingBox_ _box_traj = {0};
    BoundingBox box_seg = &_box_seg; // safe: BoundingBox_ has same layout as BoundingBox
    BoundingBox box_traj = &_box_traj;
    double t_step = (t2 - t1) / XC_SLICING_NUM_STEPS;
    double l_step = (l2 - l1) / XC_SLICING_NUM_STEPS;
    double t, l;
    int16_t num;
    for (int i = 0; i < XC_SLICING_NUM_STEPS; i++){
        t = t1 + i * t_step;
        LocalSegment_update_box(seg, box_seg, t, t + t_step);
        for (int j = 0; j < XC_SLICING_NUM_STEPS; j++){
            l = l1 + j * l_step;
            LocalTrajectory_update_box(traj, box_traj, l, l + l_step);
            if (BoundingBox_overlaps(box_seg, box_traj)){
                if (nest_level >= XC_SLICING_MAX_NEST_LEVEL - 1){
                    // We reached the maximum nesting level, return the midpoint of the current t-interval
                    if (FindRoot_get_num_solutions(finder) >= 100.){
                        printf("Warning: Maximum number of solutions (100) reached in slice_before_newton. Some solutions may be missed.\n");
                        fflush(stdout);
                        return;
                    }
                    num = FindRoot_get_num_solutions(finder);
                    fflush(stdout);
                    int8_t is_duplicate = 1;
                    for (int k = 0; k < num; k++) {
                        double dt = fabs(FindRoot_get_guess_t(finder, k) - (t + 0.5 * t_step));
                        double dl = fabs(FindRoot_get_guess_l(finder, k) - (l + 0.5 * l_step));
                        if (dt <= (0.5 * t_step) || dl <= (0.5 * l_step)) { // same region
                            is_duplicate = 0;
                            // we take the new guess instead
                            FindRoot_set_guess_t(finder, k, (t + 0.5 * t_step));
                            FindRoot_set_guess_l(finder, k, (l + 0.5 * l_step));
                            break;
                        }
                    }
                    if (is_duplicate == 0){
                        continue;
                    } else {
                        FindRoot_set_guess_t(finder, num, t + 0.5 * t_step);
                        FindRoot_set_guess_l(finder, num, l + 0.5 * l_step);
                        FindRoot_set_num_solutions(finder, num + 1);
                    }
                } else {
                    slice_before_newton(finder, seg, traj, t, t + t_step, l, l + l_step, nest_level + 1);
                }
            }
        }
    }
    return;
}

/*gpufun*/
void find_crossing_approximate(FindRoot finder, LocalSegment seg, LocalTrajectory traj){
    printf("Finding crossing approximately...\n");
    double t2;
    switch (LocalSegment_typeid(seg)){
        case LocalSegment_HalfOpenLineSegment_t:
            t2 = 10.0; // Half-open segment, we set an arbitrary large value
            break;
        default:
            t2 = 1.0;  // Line & Bezier segments
            break;
    }
    slice_before_newton(finder, seg, traj, 0, t2, 0, 10.0, 0);
    printf("Number of initial guesses found: %d\n", FindRoot_get_num_solutions(finder));
    for (int i = 0; i < FindRoot_get_num_solutions(finder); i++){ //Todo use something else than num solutions so that we can update in newton
        double guess_t = FindRoot_get_guess_t(finder, i);
        double guess_l = FindRoot_get_guess_l(finder, i);
        FindRoot_newton(finder, seg, traj, guess_t, guess_l, i);
    }
}
/*gpufun*/
void FindRoot_find_crossing(FindRoot finder, LocalSegment seg, LocalTrajectory traj){ // TODO: Temporary arrays
// this will act as the crossing main function - C magic
// First we check for the analytical solutions by checking what trajectory we have
    int8_t num_solutions = 0;
    double solution_t[3];
    double solution_l[3];
    switch (LocalTrajectory_typeid(traj)){
        case LocalTrajectory_DriftTrajectory_t:
            double xp  = DriftTrajectory_get_tan_t0((DriftTrajectory) LocalTrajectory_member(traj));
            double s0  = DriftTrajectory_get_s0((DriftTrajectory) LocalTrajectory_member(traj));
            double x0  = DriftTrajectory_get_x0((DriftTrajectory) LocalTrajectory_member(traj));
            switch (LocalSegment_typeid(seg)){
                case LocalSegment_LineSegment_t:
                    LineSegment_crossing_drift((LineSegment) LocalSegment_member(seg), solution_t, &num_solutions, s0, x0, xp);
                    for (int i = 0; i < num_solutions; i++){
                        DriftTrajectory_analytical_solution_l((DriftTrajectory) LocalTrajectory_member(traj), solution_l, i, 
                                                               LineSegment_func_s((LineSegment) LocalSegment_member(seg), solution_t[i]), 
                                                               LineSegment_func_x((LineSegment) LocalSegment_member(seg), solution_t[i]));
                        FindRoot_set_solution_t(finder, i, solution_t[i]); //Todo: temporary so we can transport to python
                        FindRoot_set_solution_l(finder, i, solution_l[i]);
                    }
                    return;
                    break;
                case LocalSegment_HalfOpenLineSegment_t:
                    HalfOpenLineSegment_crossing_drift((HalfOpenLineSegment) LocalSegment_member(seg), solution_t, &num_solutions, s0, x0, xp);
                    for (int i = 0; i < num_solutions; i++){
                        DriftTrajectory_analytical_solution_l((DriftTrajectory) LocalTrajectory_member(traj), solution_l, i, 
                                                               HalfOpenLineSegment_func_s((HalfOpenLineSegment) LocalSegment_member(seg), solution_t[i]), 
                                                               HalfOpenLineSegment_func_x((HalfOpenLineSegment) LocalSegment_member(seg), solution_t[i]));
                        FindRoot_set_solution_t(finder, i, solution_t[i]); //Todo: temporary so we can transport to python
                        FindRoot_set_solution_l(finder, i, solution_l[i]);
                    }   
                    return;
                    break;
                case LocalSegment_BezierSegment_t:
                    BezierSegment_crossing_drift((BezierSegment) LocalSegment_member(seg), solution_t, &num_solutions, s0, x0, xp);
                    for (int i = 0; i < num_solutions; i++){
                        DriftTrajectory_analytical_solution_l((DriftTrajectory) LocalTrajectory_member(traj), solution_l, i, 
                                                               BezierSegment_func_s((BezierSegment) LocalSegment_member(seg), solution_t[i]), 
                                                               BezierSegment_func_x((BezierSegment) LocalSegment_member(seg), solution_t[i]));
                        FindRoot_set_solution_t(finder, i, solution_t[i]); //Todo: temporary so we can transport to python
                        FindRoot_set_solution_l(finder, i, solution_l[i]);
                    }
                    return;
                    break;
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