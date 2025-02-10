// copyright ############################### #
// This file is part of the Xcoll package.   #
// Copyright (c) CERN, 2024.                 #
// ######################################### #

#ifndef XCOLL_GEOM_FIND_ROOT_H
#define XCOLL_GEOM_FIND_ROOT_H

#define XC_NEWTON_EPSILON 1.e-10
#define XC_NEWTON_MAX_ITER 100  // Maximum number of iterations
#define XC_NEWTON_DERIVATIVE_TOL 1e-10  // Threshold for small derivative
// #define XC_GRID_MAX_INTER 10     // Maximum number of intervals for grid search
#define XC_GRID_POINTS 1000     // Number of points to search in grid

#include <math.h>
#include <stdio.h>

/*gpufun*/
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

void F_G(LocalTrajectory traj, LocalSegment seg, double FG[2], double l, double t){
    // Here we get the expr from segment and traj, and we connect them to create F1 and F2
    // F1 = sS(t) - sT(l) 
    // F2 = xS(t) - xT(l)
    FG[0] = LocalSegment_func_s(seg, t) - LocalTrajectory_func_s(traj, l);
    FG[1] = LocalSegment_func_x(seg, t) - LocalTrajectory_func_x(traj, l);
}

void get_inv_J(LocalTrajectory traj, LocalSegment seg, double J_inv[2][2], int8_t* no_crossing, double l, double t){
    double J[2][2];
    // get derivatives
    J[0][0] = LocalSegment_deriv_s(seg, t);     // get deriv. dsS/dt LocalSegment_deriv_s 
    J[0][1] = -LocalTrajectory_deriv_s(traj,l); // get deriv. dsT/dl LocalTrajectory_deriv_s
    J[1][0] = LocalSegment_deriv_x(seg, t);     // get deriv. dxS/dt LocalSegment_deriv_x
    J[1][1] = -LocalTrajectory_deriv_x(traj,l); // get deriv. dxT/dl LocalTrajectory_deriv_x
    
    double det = J[0][0] * J[1][1] - J[1][0]*J[0][1];
    if (det < XC_NEWTON_DERIVATIVE_TOL){
        printf("There is no crossing. \n");
        *no_crossing = 1;
        return;
    }
    J_inv[0][0] =  J[1][1] / det;
    J_inv[0][1] = -J[0][1] / det;
    J_inv[1][0] = -J[1][0] / det;
    J_inv[1][1] =  J[0][0] / det;
}

void newton(LocalTrajectory traj, LocalSegment seg, double* guess_t, double* guess_l, int8_t* no_crossing) {
    double J_inv[2][2];
    double FG[2];
    

    for (int i = 0; i < XC_NEWTON_MAX_ITER; i++) {
        F_G(traj, seg, FG, *guess_l, *guess_t);
        get_inv_J(traj, seg, J_inv, no_crossing, *guess_l, *guess_t);
        if (no_crossing){
            return; 
        }
        double new_t = *guess_t - (J_inv[0][0]*FG[0] + J_inv[0][1]*FG[1]);
        double new_l = *guess_l - (J_inv[1][0]*FG[0] + J_inv[1][1]*FG[1]);

    // Check for convergence
        if ((fabs(new_t -  *guess_t) < XC_NEWTON_DERIVATIVE_TOL) && (fabs(new_l - *guess_l) < XC_NEWTON_DERIVATIVE_TOL)){
            return;
        }
        // Update the guesses for the next iteration
        *guess_t = new_t;  // Keep *t updated
        *guess_l = new_l;  // Keep *l updated
    }
}

void grid_search_and_newton(LocalTrajectory traj, LocalSegment seg, double* l, double* t, double t_min, 
                            double t_max, double l_min, double l_max, double* roots_t, double* roots_l, 
                            int max_crossings, int* number_of_roots) {
    double grid_step_t = (t_max - t_min) / 100;  // this is just for now. We need to get interval for t and l.
    double grid_step_l = (l_max - l_min) / 100;
    double FG_prev[2];
    double FG_curr[2]; 
    int interval_count = 0;
    double prev_t = t_min;  
    double prev_l = l_min;
    double curr_t;
    double curr_l;
    int8_t* no_crossing = 0;
    
    F_G(traj, seg, FG_prev, prev_t, prev_l);

    for (double t_step = t_min; t_step <= t_max; t_step += grid_step_t) {
        for (double l_step = l_min; l_step <= l_max; l_step += grid_step_l) {
            if (interval_count >= max_crossings) {
                return;                                // you cannot have more intervals than roots
            }
            curr_t = t_min + t_step;
            curr_l = l_min + l_step;
            F_G(traj, seg, FG_curr, curr_t, curr_l);
            
            if ((FG_prev[0] * FG_curr[0] < 0) || (FG_prev[1] * FG_curr[1] < 0)) {
                double initial_guess_t = 0.5*(curr_t - prev_t);
                double initial_guess_l = 0.5*(curr_l - prev_l);
                newton(traj, seg, &initial_guess_t, &initial_guess_l, no_crossing);
                if (no_crossing){
                    printf("No crossing found with Newton's method. \n");
                } else{
                    roots_t[interval_count] = initial_guess_t;
                    roots_l[interval_count] = initial_guess_l;
                    interval_count++;
                }
            }
        }
    }
}

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