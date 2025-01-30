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

// // Function to compute the value of the cubic equation at t
// /*gpufun*/
// double _bezier_cubic(double t, double a, double b, double c, double d) {
//     return a*t*t*t + b*t*t + c*t + d;
// }

// // Derivative of the cubic function
// /*gpufun*/
// double _bezier_cubic_prime(double t, double a, double b, double c) {
//     return 3*a*t*t + 2*b*t + c;
// }

// Newton's method to solve the cubic equation with safety checks
// TODO: This finds one root, but we want all. Need to find a smart combination of newton and bisection
// TODO: Need to find a way to catch oscillatory behaviour and start from another guess..
//       Example: f(x) = x3 − 2x + 2 will oscillate between 0 and 1
// TODO: Need a catch for undefined values (e.g. when derivative has a 1/x).
//       Similarily, the iteration can send itself beyond the valid domain
//       Example: f(x) = ln x  => x_(n+1) = xn(1- ln xn). If one starts at x >= e, it will get NaN


/*gpufun*/
double TrajectorySegmentFunc(int8_t traj_id, int8_t seg_id, double s, McsLineParams params) {
    switch (traj_id){
        case 0: 
            switch (seg_id) {
                case 0: return MultipleCoulomb_Line(s, params);
                case 1: return MultipleCoulomb_HalfOpenLine(s, params);
                case 2: return MultipleCoulomb_Circular(s, params);
                // case 3: return MultipleCoulomb_Bezier(s, params);
                default: 
                    printf("Segment ID not recognized\n");
                    fflush(stdout);
            }
        default: 
            printf("Trajectory ID not recognized\n");
            fflush(stdout);
    }
}

/*gpufun*/
double TrajectorySegmentDeriv(int8_t traj_id, int8_t seg_id, double s, McsLineParams params) {
    switch (traj_id){
        case 0: 
            switch (seg_id) {
                case 0: return MultipleCoulombDeriv_Line(s, params);
                case 1: return MultipleCoulombDeriv_HalfOpenLine(s, params);
                case 2: return MultipleCoulombDeriv_Circular(s, params);
                // case 3: return MultipleCoulomb_Bezier(s, params);
                default: 
                    printf("Segment ID not recognized\n");
                    fflush(stdout);
            }
        default: 
            printf("Trajectory ID not recognized\n");
            fflush(stdout);
    }
}

// TODO:  auto generate this function
/*gpufun*/
void grid_search_and_newton(int8_t traj_id, int8_t seg_id, double s_min, 
                            double s_max, double* roots, double max_crossings, McsLineParams params, int8_t* number_of_roots) {
    /// Find the intervals where the function changes sign within the range [s_min, s_max]
    //  in which later Newton's method can be applied to find the root(s) for each interval
    double grid_step = (s_max - s_min) / XC_GRID_POINTS;
    int interval_count = 0;

    double prev_s   = s_min;
    double prev_val = TrajectorySegmentFunc(traj_id, seg_id, prev_s, params);

    for (int i = 1; i <= XC_GRID_POINTS - 1; i++) {
        if (interval_count >= max_crossings) break; // you cannot have more intervals than roots
        double curr_s = s_min + i * grid_step;
        double curr_val = TrajectorySegmentFunc(traj_id, seg_id, curr_s, params);
        if (prev_val * curr_val < 0) {  
            double initial_guess = 0.5 * (prev_s + curr_s);      // initial guess is midpoint
            roots[interval_count] = newton(traj_id, seg_id, initial_guess, params);
            interval_count++;
        }
        prev_s = curr_s;
        prev_val = curr_val;
    }
    *number_of_roots = interval_count;
}

/*gpufun*/
double newton(int8_t traj_id, int8_t seg_id, double guess, McsLineParams params) {
    for (int i = 0; i < XC_NEWTON_MAX_ITER; i++) {
        double f = TrajectorySegmentFunc(traj_id, seg_id, guess, params);
        double f_prime = TrajectorySegmentDeriv(traj_id, seg_id, guess, params);
        if (fabs(f_prime) < XC_NEWTON_DERIVATIVE_TOL) return guess;

        double guess_new = guess - f / f_prime;
        if (fabs(guess_new - guess) < XC_NEWTON_EPSILON) return guess_new;
        guess = guess_new;
    }
    return guess;
}

// /*gpufun*/
// double newton_method(double a, double b, double c, double d, double initial_guess, int *status) {
//     double t = initial_guess;
//     double t_new;
//     int iter = 0;

//     while (fabs(_bezier_cubic(t, a, b, c, d)) > XC_NEWTON_EPSILON) {
//         double derivative = _bezier_cubic_prime(t, a, b, c);

//         // Check if the derivative is too small
//         if (fabs(derivative) < XC_NEWTON_DERIVATIVE_TOL) {
//             *status = -1;  // Indicate failure due to small derivative
//             return t;       // Return the best guess we have
//         }

//         t_new = t - _bezier_cubic(t, a, b, c, d) / derivative;

//         // Check if the change in t is very small
//         if (fabs(t_new - t) < XC_NEWTON_EPSILON) break;

//         t = t_new;
//         iter++;

//         // Check if we exceed maximum iterations
//         if (iter >= XC_NEWTON_MAX_ITER) {
//             *status = -2;  // Indicate failure due to max iterations
//             return t;
//         }
//     }

//     *status = 0;  // Indicate success
//     return t;
// }

// // Bisection method to find root in an interval [t1, t2]
// /*gpufun*/
// double bisection_method(double a, double b, double c, double d, double t1, double t2) {
//     double mid;

//     while ((t2 - t1) >= XC_NEWTON_EPSILON) {
//         mid = (t1 + t2) / 2;

//         // Check if the middle point is a root
//         if (fabs(_bezier_cubic(mid, a, b, c, d)) < XC_NEWTON_EPSILON) return mid;

//         // Decide the side to repeat the steps
//         if (_bezier_cubic(mid, a, b, c, d) * _bezier_cubic(t1, a, b, c, d) < 0)
//             t2 = mid;
//         else
//             t1 = mid;
//     }

//     return mid;
// }

#endif /* XCOLL_GEOM_FIND_ROOT_H */