// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2025.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#include <math.h>
#include <stdio.h>


// /*gpufun*/
// void grid_search_and_newton(LocalTrajectory traj, LocalSegment seg,
//                             double s_min, double s_max, double* roots, double max_crossings,
//                             int8_t* number_of_roots) {
//     /// Find the intervals where the function changes sign within the range [s_min, s_max]
//     //  in which later Newton's method can be applied to find the root(s) for each interval
//     double grid_step = (s_max - s_min) / XC_GRID_POINTS;
//     int interval_count = 0;

//     double prev_s   = s_min;
//     double prev_val = LocalTrajectory_func(traj, prev_s) - LocalSegment_func(seg, prev_s);

//     for (int i = 1; i <= XC_GRID_POINTS - 1; i++) {
//         if (interval_count >= max_crossings) break; // you cannot have more intervals than roots
//         double curr_s = s_min + i * grid_step;
//         double curr_val = LocalTrajectory_func(traj, curr_s) - LocalSegment_func(seg, curr_s);
//         if (prev_val * curr_val < 0) {
//             double initial_guess = 0.5 * (prev_s + curr_s);      // initial guess is midpoint
//             roots[interval_count] = newton(traj, seg, initial_guess);
//             interval_count++;
//         }
//         prev_s = curr_s;
//         prev_val = curr_val;
//     }
//     *number_of_roots = interval_count;
// }

// /*gpufun*/
// double newton(LocalTrajectory traj, LocalSegment seg, double guess_s) {
//     for (int i = 0; i < XC_NEWTON_MAX_ITER; i++) {
//         double f = LocalTrajectory_func(traj, guess_s) - LocalSegment_func(seg, guess_s);
//         double f_prime = LocalTrajectory_deriv(traj, guess_s) - LocalSegment_deriv(seg, guess_s);
//         if (fabs(f_prime) < XC_NEWTON_DERIVATIVE_TOL) return guess_s;

//         double guess_new = guess_s - f / f_prime;
//         if (fabs(guess_new - guess_s) < XC_NEWTON_EPSILON) return guess_new;
//         guess_s = guess_new;
//     }
//     return guess_s;
// }

// // --------------------------------------------------------------------------------------------
/*gpufun*/
void LocalCrossing_func(LocalTrajectory traj, LocalSegment seg, double TS[2], double l, double t){
    // Here we get the expr from segment and traj, and we connect them to create TSs and TSx
    // TSs = Ts(l) - Ss(t)
    // TSx = Tx(l) - Sx(t)
    TS[0] = LocalTrajectory_func_s(traj, l) - LocalSegment_func_s(seg, t);
    TS[1] = LocalTrajectory_func_x(traj, l) - LocalSegment_func_x(seg, t);
}
/*gpufun*/
void LocalCrossing_inv_J(LocalTrajectory traj, LocalSegment seg, double J_inv[2][2], double success[1], double l, double t){
    double J[2][2];
    // get derivatives
    J[0][0] = LocalTrajectory_deriv_s(traj, l); // get deriv. dsT/dl LocalTrajectory_deriv_s
    J[0][1] = -LocalSegment_deriv_s(seg, t);     // get deriv. dsS/dt LocalSegment_deriv_s
    J[1][0] = LocalTrajectory_deriv_x(traj, l); // get deriv. dxT/dl LocalTrajectory_deriv_x
    J[1][1] = -LocalSegment_deriv_x(seg, t);     // get deriv. dxS/dt LocalSegment_deriv_x
    double det = J[0][0] * J[1][1] - J[1][0]*J[0][1];
    if (fabs(det) < XC_GEOM_ROOT_NEWTON_DERIVATIVE_TOL){
        success[0] = 1;
        return; 
    }
    J_inv[0][0] =  J[1][1] / det;
    J_inv[0][1] = -J[0][1] / det;
    J_inv[1][0] = -J[1][0] / det;
    J_inv[1][1] =  J[0][0] / det;
    }

/*gpufun*/
void find_root_newton(find_root finder, LocalTrajectory traj, LocalSegment seg, double guess_l, double guess_t){
    // success = 0 means success, success = 1 means failure
    double J_inv[2][2];
    double TS[2];
    double success[1] = {0};
    double corr0, corr1;
    double new_t, new_l;
    for (int i = 0; i < XC_GEOM_ROOT_NEWTON_MAX_ITER; i++){
        LocalCrossing_func(traj, seg, TS, guess_l, guess_t);
        LocalCrossing_inv_J(traj, seg, J_inv, success, guess_l, guess_t);
        if (success[0]){
            printf("Jacobian is singular. Stopping Newton's method.\n");
            return;
        }
        // Residual convergence check
        if (fabs(TS[0]) < XC_GEOM_ROOT_NEWTON_EPSILON && fabs(TS[1]) < XC_GEOM_ROOT_NEWTON_EPSILON) {
            find_root_set_solution_l(finder, guess_l);
            find_root_set_solution_t(finder, guess_t);
            return;
        }
        corr0 = J_inv[0][0]*TS[0] + J_inv[0][1]*TS[1]; // delta for l
        corr1 = J_inv[1][0]*TS[0] + J_inv[1][1]*TS[1]; // delta for t
        new_t = guess_t - corr1;
        new_l = guess_l - corr0;
        // Check for parameter update convergence
        if ((fabs(new_t -  guess_t) < XC_GEOM_ROOT_NEWTON_EPSILON) && (fabs(new_l - guess_l) < XC_GEOM_ROOT_NEWTON_EPSILON)){
            find_root_set_solution_l(finder, guess_l);
            find_root_set_solution_t(finder, guess_t);
            return;
        }
        // Update the guesses for the next iteration
        guess_t = new_t;  // Keep *t updated
        guess_l = new_l;  // Keep *l updated
    }
    printf("Warning: Newton's method did not converge within the maximum number of iterations (%d).\n", XC_GEOM_ROOT_NEWTON_MAX_ITER);
    success[0] = 1;
}

// int8_t LocalCrossing_box_has_root(double TS_UL[2], double TS_UR[2], double TS_DL[2], double TS_DR[2]){

//     // Evaluate F at corners
//     evaluate_traj_seg(traj, seg, TS_DL, l1, t1);
//     evaluate_traj_seg(traj, seg, TS_DR, l2, t1);
//     evaluate_traj_seg(traj, seg, TS_UL, l1, t2);
//     evaluate_traj_seg(traj, seg, TS_UR, l2, t2);

//     // Check for sign changes for each component across the corners.
//     // (For simplicity, here we check if the min and max differ in sign.)
//     int8_t f1_has_sign_change = (min(TS_DL[0], TS_DR[0], TS_UL[0], TS_UR[0]) < 0 &&
//                                 max(TS_DL[0], TS_DR[0], TS_UL[0], TS_UR[0]) > 0);
//     int8_t f2_has_sign_change = (min(TS_DL[1], TS_DR[1], TS_UL[1], TS_UR[1]) < 0 &&
//                                 max(TS_DL[1], TS_DR[1], TS_UL[1], TS_UR[1]) > 0);
//     return (l1 < 0 && l2 > 0) || (l1 > 0 && l2 < 0) || (t1 < 0 && t2 > 0) || (t1 > 0 && t2 < 0);
// }


// void grid_search_and_newton(LocalTrajectory traj, LocalSegment seg, double* l, double* t, double l_min,
//                             double l_max, double t_min, double t_max, double* roots_l, double* roots_t,
//                             int max_crossings, int* number_of_roots){
//     int N_l = 100;
//     int N_t = 100;
//     double grid_step_l = (l_max - l_min) / N_l; // this is just for now. We need to get interval for t and l.
//     double grid_step_t = (t_max - t_min) / N_t;
//     int n_roots = 0;
//     double TS_prev, TS_curr;
//     double prev_t = t_min;
//     double prev_l = l_min;
//     double curr_t;
//     double curr_l;
//     int8_t no_crossing = 0;

//     LocalCrossing_func(traj, seg, TS_prev, prev_t, prev_l);

//     // Stage 1: Coarse grid search over (l, t)
//     for (int i=0; i < N_l-1; i++){
//         double l1 = l_min + i * grid_step_l;
//         double l2 = l_min + (i+1) * grid_step_l;
//         for (int j=0; j < N_t-1; j++){
//             double t1 = t_min + j * grid_step_t;
//             double t2 = t_min + (j+1) * grid_step_t;
//             F_G(traj, seg, FG_curr, curr_t, curr_l);
//             if (n_roots >= max_crossings) {
//                 return;  // all possible roots have been found
//             }
//             curr_t = t_min + t_step;
//             curr_l = l_min + l_step;
//             F_G(traj, seg, FG_curr, curr_t, curr_l);

//             if ((FG_prev[0] * FG_curr[0] < 0) || (FG_prev[1] * FG_curr[1] < 0)) {
//                 double initial_guess_t = 0.5*(curr_t - prev_t);
//                 double initial_guess_l = 0.5*(curr_l - prev_l);
//                 newton(traj, seg, &initial_guess_l, &initial_guess_t, &no_crossing);
//                 if (no_crossing){
//                     printf("No crossing found with Newton's method. \n");
//                 } else{
//                     roots_t[n_roots] = initial_guess_t;
//                     roots_l[n_roots] = initial_guess_l;
//                     n_roots++;
//                 }
//             }
//         }
//     }
// }

//  double prev_s   = s_min;
//     double prev_val = LocalTrajectory_func(traj, prev_s) - LocalSegment_func(seg, prev_s);

//     for (int i = 1; i <= XC_GRID_POINTS - 1; i++) {
//         if (interval_count >= max_crossings) break; // you cannot have more intervals than roots
//         double curr_s = s_min + i * grid_step;
//         double curr_val = LocalTrajectory_func(traj, curr_s) - LocalSegment_func(seg, curr_s);
//         if (prev_val * curr_val < 0) {
//             double initial_guess = 0.5 * (prev_s + curr_s);      // initial guess is midpoint
//             roots[interval_count] = newton(traj, seg, initial_guess);
//             interval_count++;
//         }
//         prev_s = curr_s;
//         prev_val = curr_val;


// int main() {
//     double s0 = 1.5;
//     double x0 = 0.0;
//     double theta = 20.0*M_PI/180.0;
//     double s1 = 1.0;
//     double x1 =-1.0;
//     double s2 =2.0;
//     double x2 =2.0;
//     double l = -1.0;
//     double t = -1.0;
//     Params p = {s0, x0, theta, s1, x1, s2, x2, l, t};
//     double guess_t = 0.5;
//     double guess_l = 0.5;
//     newton(&guess_t, &guess_l, &p);
//     return 0;
// }

#endif /* XCOLL_GEOM_FIND_ROOT_H */